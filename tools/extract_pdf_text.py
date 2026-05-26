from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_pages(spec: str | None, total: int) -> list[int]:
    if not spec:
        return list(range(total))
    pages: set[int] = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_s, end_s = part.split("-", 1)
            start = max(1, int(start_s))
            end = min(total, int(end_s))
            pages.update(range(start - 1, end))
        else:
            page = int(part)
            if 1 <= page <= total:
                pages.add(page - 1)
    return sorted(pages)


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract PDF text into Markdown-ish private notes.")
    parser.add_argument("pdf", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--pages", help="1-based pages or ranges, e.g. 1-20,45,60-62")
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(".tools/pypdf").resolve()))
    from pypdf import PdfReader

    reader = PdfReader(str(args.pdf))
    selected = parse_pages(args.pages, len(reader.pages))
    args.out.parent.mkdir(parents=True, exist_ok=True)

    title = args.title or args.pdf.stem
    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# {title}\n\n")
        handle.write(f"Source PDF: `{args.pdf.as_posix()}`\n\n")
        handle.write(f"Extracted pages: {args.pages or 'all'}\n\n")
        handle.write("> Private extraction for indexing and instructor reference. Do not copy extracted prose into public course materials.\n\n")
        for page_idx in selected:
            handle.write(f"\n\n## Page {page_idx + 1}\n\n")
            text = reader.pages[page_idx].extract_text() or ""
            text = "\n".join(line.rstrip() for line in text.splitlines())
            handle.write(text.strip())
            handle.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

