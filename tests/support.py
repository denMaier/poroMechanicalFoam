from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def strip_c_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*", "", text)


def active_makefile_sources(makefile: Path) -> list[Path]:
    """Return active C/C++ source paths listed in an OpenFOAM Make/files file."""
    text = strip_c_comments(makefile.read_text())
    base = makefile.parents[1]
    sources: list[Path] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(("LIB", "EXE")):
            continue

        line = line.rstrip("\\").strip()
        if line.endswith((".C", ".cc", ".cpp", ".cxx")):
            sources.append(base / line)

    return sources


def header_typenames(paths: list[Path] | None = None) -> dict[str, set[Path]]:
    if paths is None:
        paths = list(ROOT.rglob("*.H"))

    typenames: dict[str, set[Path]] = {}
    for path in paths:
        if any(part in {"external", "oldLimitedHead", "templatePoroFluid"} for part in path.parts):
            continue

        text = path.read_text(errors="ignore")
        for match in re.finditer(r'TypeName\s*\(\s*"([^"]+)"\s*\)', text):
            typenames.setdefault(match.group(1), set()).add(path)

    return typenames


def dictionary_value(path: Path, key: str) -> str | None:
    text = strip_c_comments(path.read_text())
    match = re.search(rf"^\s*{re.escape(key)}\s+([^;]+);", text, flags=re.MULTILINE)
    return match.group(1).strip() if match else None
