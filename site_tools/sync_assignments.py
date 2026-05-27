"""Mirror /assignments into /course_notes/assignments before Quarto renders.

Assignment .qmd files live at the repository root so they are visible at the
top of the GitHub tree alongside data/, course_notes/, etc. Quarto, however,
cannot resolve links to files outside its project directory when publishing
to GitHub Pages (../assignments/ resolves above the site root). This script
runs as a pre-render hook to copy assignments into course_notes/assignments/
so Quarto sees them as a regular subdirectory. The copy is gitignored.
"""

from __future__ import annotations

import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "assignments"
DST = REPO_ROOT / "course_notes" / "assignments"


def main() -> None:
    if not SRC.is_dir():
        raise SystemExit(f"Expected assignments directory at {SRC}")
    if DST.exists():
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    print(f"Synced {SRC} -> {DST}")


if __name__ == "__main__":
    main()
