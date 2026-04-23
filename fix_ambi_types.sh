#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${TMPDIR:-/tmp}/pmf_tmp_ambiguity"
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
collect  Run ./Allwmake with clang diagnostic flags and collect ambiguous tmp<Field> sites
preview  Show dry-run diffs for collected sites
apply    Apply the rewrite for collected sites

Assumptions:
- run this from the repo root
- OpenFOAM is already sourced
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

    TARGET_LINES="$target_lines" perl -0pe '
        my %target = map { $_ => 1 } grep { length } split /,/, $ENV{TARGET_LINES};
        my @lines = split /\n/, $_, -1;

        for my $i (0 .. $#lines) {
            my $ln = $i + 1;
            next unless $target{$ln};

            if (
                $lines[$i] =~ /^(\s*)(const\s+)?((?:Field<[^>]+>|[A-Za-z_][A-Za-z0-9_:<>]*Field))\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?);(\s*(?://.*)?)\s*$/
            ) {
                my ($indent, $const, $type, $name, $expr, $comment) =
                    ($1, $2 // q{}, $3, $4, $5, $6 // q{});
                $lines[$i] = "${indent}${const}${type} ${name}((${expr})());${comment}";
            }
            else {
                warn "$ARGV:$ln: skipped; expected single-line field initialization\n";
            }
        }

        $_ = join "\n", @lines;
    ' "$src" > "$dst"
}

preview() {
    [[ -s "$LOCS" ]] || { echo "no collected sites; run '$0 collect' first"; exit 1; }

    declare -A line_map=()
    while IFS=: read -r file line; do
        [[ -n "$file" && -n "$line" ]] || continue
        if [[ -n "${line_map[$file]:-}" ]]; then
            line_map["$file"]+=","
        fi
        line_map["$file"]+="$line"
    done < "$LOCS"

    local shown=0
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue
        local tmp
        tmp="$(mktemp)"
        rewrite_file "$file" "${line_map[$file]}" "$tmp"

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

    declare -A line_map=()
    while IFS=: read -r file line; do
        [[ -n "$file" && -n "$line" ]] || continue
        if [[ -n "${line_map[$file]:-}" ]]; then
            line_map["$file"]+=","
        fi
        line_map["$file"]+="$line"
    done < "$LOCS"

    local changed=0
    while IFS= read -r file; do
        [[ -n "$file" ]] || continue
        local tmp
        tmp="$(mktemp)"
        rewrite_file "$file" "${line_map[$file]}" "$tmp"

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
