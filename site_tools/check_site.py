"""Check the generated public site before deployment (standard library only)."""
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "_site"


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.path = path
        self.ids = set()
        self.links = []
        self.errors = []
        self.h1 = 0
        self.feed(path.read_text(encoding="utf-8"))

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            if attrs["id"] in self.ids:
                self.errors.append(f'Duplicate id: {attrs["id"]}')
            self.ids.add(attrs["id"])
        if tag == "h1":
            self.h1 += 1
        if tag == "img" and not attrs.get("alt"):
            self.errors.append("Image has no alternative text")
        for key in ("href", "src"):
            if key in attrs:
                self.links.append(attrs[key])


def main():
    pages = {p.resolve(): Page(p) for p in OUT.rglob("*.html")}
    errors = []
    assert pages, "Build the site first"
    for path, page in pages.items():
        name = path.relative_to(OUT)
        errors.extend(f"{name}: {error}" for error in page.errors)
        if page.h1 != 1:
            errors.append(f"{name}: expected one main heading, found {page.h1}")
        text = path.read_text(encoding="utf-8")
        if re.search(r"@(def|thm|prp|exm|exr|fig|eq|tbl)-[\w-]+", text):
            errors.append(f"{name}: unresolved cross-reference")
        if re.search(r"\{#(?:eq|fig|tbl)-", text):
            errors.append(f"{name}: unprocessed label")
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            if url.path.startswith("/"):
                errors.append(f"{name}: root-relative link breaks repository hosting: {link}")
                continue
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            if not target.is_relative_to(OUT):
                errors.append(f"{name}: link escapes public output: {link}")
            elif not target.exists():
                errors.append(f"{name}: missing target: {link}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                errors.append(f"{name}: missing fragment: {link}")
    index = json.loads((OUT / "search-index.json").read_text(encoding="utf-8"))
    indexed = {page["url"] for page in index}
    home = (OUT / "index.html").read_text(encoding="utf-8")
    for directory in (ROOT / "course_notes/lessons", ROOT / "assignments"):
        for source in directory.glob("*.md"):
            body = source.read_text(encoding="utf-8")
            route = directory.name + "/" + source.stem + ".html"
            if re.search(r"^published:\s*false\s*$", body, re.M):
                if route in indexed or f'href="{route}"' in home:
                    errors.append(f"Draft is linked or searchable: {route}")
                if "These notes are not available yet" not in (OUT / route).read_text(encoding="utf-8"):
                    errors.append(f"Draft does not show availability notice: {route}")
    for path in OUT.rglob("*"):
        if any(part in ("planning", "reference_docs", "source_digests", "syllabus_revised", ".git", ".tools") for part in path.relative_to(OUT).parts):
            errors.append(f"Private or non-site path in output: {path}")
    topics = json.loads((ROOT / "site/schedule.json").read_text(encoding="utf-8"))
    active_topics = [topic for topic in topics if topic["lessons"]]
    if home.count('class="topic-heading"') != len(active_topics):
        errors.append("Schedule does not have one heading per topic")
    if home.count('class="topic-lectures"') != len(active_topics):
        errors.append("Schedule does not group individual lectures under topics")
    lesson_numbers = [n for topic in topics for n in topic["lessons"]]
    if sorted(lesson_numbers) != list(range(1, 25)):
        errors.append("Expected each of the 24 lectures in exactly one topic")
    for topic in topics:
        if f'id="topic-{topic["id"]}"' not in home:
            errors.append("Missing schedule topic " + topic["title"])
    if home.count('Problem set') != len(active_topics):
        errors.append("Expected one problem-set entry per topic")
    for path in pages:
        text = path.read_text(encoding="utf-8")
        if any(label in text for label in ("View Markdown source", "Print this page", '<footer')):
            errors.append(f"Redundant page controls remain: {path.name}")
    for token in ("{{room}}", "{{meeting_days}}", "{{meeting_time}}"):
        if token in home:
            errors.append("Unexpanded course setting: " + token)
    if "Project workshop" in home or "Final project" in home:
        errors.append("Unsettled project activities in public schedule")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Checked {len(pages)} HTML pages, {len(index)} searchable pages, all local links and fragments, topic groupings, and draft exclusion.")


if __name__ == "__main__":
    main()
