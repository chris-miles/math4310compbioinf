"""Check the generated public site before deployment (standard library only)."""
from html.parser import HTMLParser
from datetime import date, timedelta
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
    schedule = json.loads((ROOT / "site/schedule.json").read_text(encoding="utf-8"))
    weeks = schedule["weeks"]
    if home.count('class="module-heading"') != len(schedule["modules"]):
        errors.append("Expected one heading per broad module")
    if home.count('class="week-row"') != sum(bool(w['meetings']) for w in weeks):
        errors.append("Expected one compact row per instructional week")
    lesson_numbers = [n for week in weeks for n in week["lessons"]]
    if sorted(lesson_numbers) != list(range(1, 24)):
        errors.append("Expected each lecture in exactly one weekly row")
    module_lessons = [n for module in schedule["modules"] for n in module["lessons"]]
    if sorted(module_lessons) != list(range(1, 24)):
        errors.append("Expected each lecture in exactly one module")
    for week in weeks:
        if f'id="week-{week["week"]}"' not in home:
            errors.append("Missing schedule week " + str(week["week"]))
    if len(weeks) != 16 or [week["week"] for week in weeks] != list(range(1, 17)):
        errors.append("Expected 16 dated calendar weeks through April 27")
    if weeks[0]["date_start"] != "2027-01-11" or weeks[-1]["date_end"] != "2027-04-27":
        errors.append("Schedule dates do not cover the official teaching term")
    if weeks[8]["date_start"] != "2027-03-06" or weeks[8]["date_end"] != "2027-03-14":
        errors.append("Spring break dates do not match the official calendar")
    if home.count('class="week-dates"') != 16:
        errors.append("A week is missing its displayed date range")
    if home.count('Problem set') != 10:
        errors.append("Expected 10 problem sets at topic boundaries")
    for path in pages:
        text = path.read_text(encoding="utf-8")
        if any(label in text for label in ("View Markdown source", "Print this page", '<footer')):
            errors.append(f"Redundant page controls remain: {path.name}")
    for token in ("{{room}}", "{{meeting_days}}", "{{meeting_time}}"):
        if token in home:
            errors.append("Unexpanded course setting: " + token)
    calendar = schedule["calendar"]
    current = date.fromisoformat(calendar["term_start"])
    end = date.fromisoformat(calendar["classes_end"])
    expected_dates = []
    while current <= end:
        day = current.isoformat()
        in_break = calendar["spring_break_start"] <= day <= calendar["spring_break_end"]
        if current.weekday() in (0, 2) and day not in calendar["holidays"] and not in_break:
            expected_dates.append(day)
        current += timedelta(days=1)
    meetings = [meeting for week in weeks for meeting in week["meetings"]]
    if sorted(m["date"] for m in meetings) != expected_dates:
        errors.append("Every available class date must have exactly one meeting")
    for week in weeks:
        if week["lessons"] != [m["lesson"] for m in week["meetings"] if "lesson" in m]:
            errors.append("Lesson list disagrees with meetings in week " + str(week["week"]))
        if any(not week["date_start"] <= m["date"] <= week["date_end"] for m in week["meetings"]):
            errors.append("Meeting lies outside its calendar week")
    events = [m for m in meetings if "event" in m]
    if len(events) != 4 or sum("Midterm" in m["event"] for m in events) != 2:
        errors.append("Expected two midterms and two presentation meetings")
    if [m["date"] for m in events if m["event"] == "Project presentations"] != ["2027-04-21", "2027-04-26"]:
        errors.append("Expected presentations on April 21 and April 26")
    if 'href="appendices/computing-math-reference.html">Reference</a>' not in home:
        errors.append("Missing reference page in main navigation")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Checked {len(pages)} HTML pages, {len(index)} searchable pages, all local links and fragments, module groupings and assignment placement, and draft exclusion.")


if __name__ == "__main__":
    main()
