#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$ROOT/.tmp_ambiguity"
LOG="$WORKDIR/allwmake.log"
LOCS="$WORKDIR/ambiguous_locations.txt"
FILES="$WORKDIR/ambiguous_files.txt"

mkdir -p "$WORKDIR"

usage() {
cat <<'EOF'
Usage:
  ./fix_ambiguous_typed.sh collect
  ./fix_ambiguous_typed.sh preview
  ./fix_ambiguous_typed.sh apply

Commands:
  collect  Run ./Allwmake with extra clang diagnostics and collect ambiguous tmp<Field> sites
  preview  Show dry-run diffs for collected sites
  apply    Apply the rewrite for collected sites

Notes:
  - Run this from a shell where OpenFOAM is already sourced.
  - Temporary artifacts are written to ./.tmp_ambiguity/.
  - Only compiler-reported lines are rewritten.
EOF
}

collect() {
    local diag_flags=(
        -ferror-limit=0
        -fno-color-diagnostics
        -fdiagnostics-absolute-paths
        -fdiagnostics-format=clang
    )

    : > "$LOG"
    : > "$LOCS"
    : > "$FILES"

    (
        cd "$ROOT"
        export WM_CXXFLAGS="${WM_CXXFLAGS-} ${diag_flags[*]}"
        export CXXFLAGS="${CXXFLAGS-} ${diag_flags[*]}"
        ./Allwmake -j1 2>&1 | tee "$LOG"
    ) || true

    awk '
        /error: conversion from / && /tmp</ && / is ambiguous/ {
            print $1
        }
    ' "$LOG" \
    | sed 's/:$//' \
    | awk -F: 'NF >= 2 { print $1 ":" $2 }' \
    | sort -u > "$LOCS"

    awk -F: 'NF >= 2 { print $1 }' "$LOCS" | sort -u > "$FILES"

    echo "log:   $LOG"
    echo "sites: $(wc -l < "$LOCS")"
    echo "files: $(wc -l < "$FILES")"

    if [[ -s "$LOCS" ]]; then
        echo
        echo "first few collected sites:"
        head -n 10 "$LOCS"
    fi
}

rewrite_file() {
    local src="$1"
    local target_lines="$2"
    local dst="$3"

    awk -v target_lines="$target_lines" '
        function ltrim(s) {
            sub(/^[[:space:]]+/, "", s)
            return s
        }

        function rtrim(s) {
            sub(/[[:space:]]+$/, "", s)
            return s
        }

        BEGIN {
            n = split(target_lines, lines, ",")
            for (i = 1; i <= n; ++i) {
                if (lines[i] != "") {
                    target[lines[i]] = 1
                }
            }
        }

        {
            original = $0

            if (!(FNR in target)) {
                print original
                next
            }

            code = original
            comment = ""

            if (match(code, /[[:space:]]*\/\/.*/)) {
                comment = substr(code, RSTART)
                code = substr(code, 1, RSTART - 1)
            }

            code = rtrim(code)

            if (substr(code, length(code), 1) != ";") {
                print original
                printf("%s:%d: skipped; expected a trailing semicolon\n", FILENAME, FNR) > "/dev/stderr"
                next
            }

            code = substr(code, 1, length(code) - 1)
            eq = index(code, "=")

            if (!eq) {
                print original
                printf("%s:%d: skipped; expected an assignment\n", FILENAME, FNR) > "/dev/stderr"
                next
            }

            lhs = rtrim(substr(code, 1, eq - 1))
            rhs = ltrim(substr(code, eq + 1))

            if (rhs == "") {
                print original
                printf("%s:%d: skipped; empty rhs\n", FILENAME, FNR) > "/dev/stderr"
                next
            }

            print lhs "((" rhs ")());" comment
        }
    ' "$src" > "$dst"
}

load_line_map() {
    declare -gA LINE_MAP=()

    while IFS=: read -r file line; do
        [[ -n "$file" && -n "$line" ]] || continue
        if [[ -n "${LINE_MAP[$file]:-}" ]]; then
            LINE_MAP["$file"]+=","
        fi
        LINE_MAP["$file"]+="$line"
    done < "$LOCS"
}

preview() {
    [[ -s "$LOCS" ]] || { echo "no collected sites; run '$0 collect' first"; exit 1; }

    load_line_map

    local shown=0
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue

        local tmp
        tmp="$(mktemp)"
        rewrite_file "$file" "${LINE_MAP[$file]}" "$tmp"

        if ! cmp -s "$file" "$tmp"; then
            echo "==== $file ===="
            diff -u "$file" "$tmp" || true
            shown=$((shown + 1))
        fi

        rm -f "$tmp"
    done < "$FILES"

    echo "previewed files: $shown"
}

apply() {
    [[ -s "$LOCS" ]] || { echo "no collected sites; run '$0 collect' first"; exit 1; }

    load_line_map

    local changed=0
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue

        local tmp
        tmp="$(mktemp)"
        rewrite_file "$file" "${LINE_MAP[$file]}" "$tmp"

        if ! cmp -s "$file" "$tmp"; then
            chmod --reference="$file" "$tmp"
            mv "$tmp" "$file"
            echo "patched $file"
            changed=$((changed + 1))
        else
            rm -f "$tmp"
        fi
    done < "$FILES"

    echo "patched files: $changed"
}

cmd="${1:-}"
case "$cmd" in
    collect) collect ;;
    preview) preview ;;
    apply) apply ;;
    *) usage; exit 1 ;;
esac
