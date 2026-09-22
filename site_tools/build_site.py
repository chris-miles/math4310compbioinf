"""Build the public course site from Markdown, Python examples, and HTML templates.

Only the explicit content directories below are published. Planning and reference
materials never enter the output. Run from any directory: python site_tools/build_site.py
"""
from __future__ import annotations

import ast
from collections import defaultdict
from datetime import date
from contextlib import redirect_stdout, redirect_stderr
import html
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "_site"
CONFIG = json.loads((ROOT / "site/config.json").read_text(encoding="utf-8"))
KINDS = {"def": "Definition", "thm": "Theorem", "lem": "Lemma",
         "cor": "Corollary", "prp": "Proposition", "exm": "Example",
         "exr": "Exercise", "fig": "Figure", "eq": "Equation", "tbl": "Table"}
REFERENCES = {}
CONTENT = {}


def pandoc(text: str, *args: str) -> str:
    result = subprocess.run(["pandoc", *args], input=text, text=True,
                            encoding="utf-8", capture_output=True, cwd=ROOT)
    if result.returncode or "[WARNING]" in result.stderr:
        raise RuntimeError(result.stderr)
    return result.stdout


def raw(text: str) -> dict:
    return {"t": "RawBlock", "c": ["html", text]}


def inlines(text: str) -> list:
    return [{"t": "Str", "c": text}]


def para(text: str) -> dict:
    return {"t": "Para", "c": [{"t": "Strong", "c": inlines(text)}]}


def walk(value, fn):
    if isinstance(value, list):
        return [walk(item, fn) for item in value]
    if isinstance(value, dict):
        value = {key: walk(item, fn) for key, item in value.items()}
        return fn(value)
    return value


def read_source(path: Path) -> tuple[dict, str]:
    text = path.read_text(encoding="utf-8")
    match = re.match(r"\A---\s*\n(.*?)\n---\s*\n", text, re.S)
    return (yaml.safe_load(match[1]) or {}, text[match.end():]) if match else ({}, text)


def execute_chunks(body: str, path: Path, route: str) -> str:
    namespace = {"__name__": "__main__"}
    count = 0
    previous = Path.cwd()
    sys_path = sys.path[:]
    os.chdir(path.parent)
    plt.close("all")
    defaults = plt.rcParams.copy()

    def chunk(match):
        nonlocal count
        count += 1
        lines = match[1].splitlines()
        options = yaml.safe_load("\n".join(line[2:].strip() for line in lines if line.startswith("#|"))) or {}
        code = "\n".join(line for line in lines if not line.startswith("#|")).strip()
        label = options.get("label", f"chunk-{count}")
        if not re.fullmatch(r"[a-zA-Z0-9_-]+", label):
            raise ValueError(f"Invalid chunk label: {label}")
        plt.rcParams["figure.figsize"] = [options.get("fig-width", CONFIG["figure"]["width"]), options.get("fig-height", CONFIG["figure"]["height"])]
        plt.rcParams["figure.dpi"] = CONFIG["figure"]["dpi"]
        output = io.StringIO()
        if options.get("eval", True):
            tree = ast.parse(code, filename=str(path))
            final = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
            with redirect_stdout(output), redirect_stderr(output):
                exec(compile(tree, str(path), "exec"), namespace)
                if final:
                    value = eval(compile(ast.Expression(final.value), str(path), "eval"), namespace)
                    if value is not None:
                        print(repr(value))
        parts = []
        if options.get("echo", True):
            block = "\n".join(["~~~python", code, "~~~"])
            if options.get("code-fold"):
                summary = html.escape(options.get("code-summary", "Show code"))
                block = f"<details class='code-fold'><summary>{summary}</summary>\n\n{block}\n\n</details>"
            parts.append(block)
        if output.getvalue().strip():
            parts.append("~~~{.output}\n" + output.getvalue().strip() + "\n~~~")
        for i, number in enumerate(plt.get_fignums()):
            figure = plt.figure(number)
            filename = f"{Path(route).stem}-{label}-{i}.svg"
            figure.savefig(OUT / "assets/figures" / filename, format="svg")
            plt.close(figure)
            alt = options.get("fig-alt")
            caption = options.get("fig-cap")
            if not alt or not caption:
                raise ValueError(f"Figure needs fig-alt and fig-cap: {path}: {label}")
            prefix = "../" * len(Path(route).parent.parts)
            alt = html.escape(str(alt), quote=True)
            image = f'<img src="{prefix}assets/figures/{filename}" alt="{alt}" loading="lazy">'
            parts.append(f"::: {{#{label} .course-figure}}\n\n{image}\n\n{caption}\n\n:::")
        return "\n\n".join(parts) + "\n"

    try:
        body = re.sub(r"^\x60{3}\{\.python \.execute\}\s*\n(.*?)^\x60{3}\s*$", chunk, body, flags=re.M | re.S)
    except Exception as error:
        raise RuntimeError(f"Python example failed in {path.relative_to(ROOT)}: {error}") from error
    finally:
        os.chdir(previous)
        sys.path[:] = sys_path
        plt.rcParams.update(defaults)
        plt.close("all")
    # Give labeled equations a real HTML anchor before Pandoc parses them.
    return re.sub(r"(\$\$[ \t]*\n(?:(?!\$\$).)*?\n\$\$)[ \t]*\{#(eq-[\w-]+)\}",
                  r"\n\n::: {#\2 .equation}\n\n\1\n\n:::\n", body, flags=re.S)


def prepare(path: Path, route: str, number: str = "") -> dict:
    meta, body = read_source(path)
    if path.parent.name in ("lessons", "assignments"):
        meta.setdefault("published", False)
    for key in ("room", "canvas", "email", "instructor", "term", "institution", "section", "class_number", "credits", "meeting_days", "meeting_time", "format"):
        body = body.replace("{{" + key + "}}", html.escape(CONFIG[key]))
    if meta.get("published") is False:
        body = "These notes are not available yet. Return to the [weekly schedule](" + ("../" * len(Path(route).parent.parts)) + "index.html#course-schedule)."
    for marker, value in CONTENT.items():
        body = body.replace("<!-- " + marker + " -->", value)
    body = execute_chunks(body, path, route)
    citation_header = ""
    if meta.get("nocite"):
        citation_header = "---\n" + yaml.safe_dump({"nocite": meta["nocite"]}) + "---\n"
    doc = json.loads(pandoc(citation_header + body, "-f", "markdown", "-t", "json"))
    counters = defaultdict(int)

    def labels(node):
        if node.get("t") == "Table":
            caption = node["c"][1]
            match = re.search(r"\{#(tbl-[\w-]+)\}", plain(caption))
            if match:
                ident = match[1]
                counters["tbl"] += 1
                label = f'Table {number + "." if number else ""}{counters["tbl"]}'
                if ident in REFERENCES:
                    raise ValueError(f"Duplicate reference: {ident}")
                REFERENCES[ident] = (route, label)
                node["c"][0][0] = ident
                def clean_caption(value):
                    if value.get("t") == "Str":
                        value["c"] = value["c"].replace("{#" + ident + "}", "")
                    return value
                node["c"][1] = walk(caption, clean_caption)
                if node["c"][1][1]:
                    node["c"][1][1][0]["c"] = inlines(label + ". ") + node["c"][1][1][0]["c"]
        if node.get("t") == "Div":
            attr, blocks = node["c"]
            ident, classes, pairs = attr
            kind = ident.split("-")[0]
            if kind in KINDS:
                counters[kind] += 1
                label = f'{KINDS[kind]} {number + "." if number else ""}{counters[kind]}'
                if ident in REFERENCES:
                    raise ValueError(f"Duplicate reference: {ident}")
                REFERENCES[ident] = (route, label)
                classes.append("numbered-block")
                classes.append(kind)
                heading = inlines(label)
                if blocks and blocks[0]["t"] == "Header":
                    heading += inlines(". ") + blocks.pop(0)["c"][2]
                blocks.insert(0, {"t": "Para", "c": [{"t": "Strong", "c": heading}]})
            options = dict(pairs)
            if "proof" in classes:
                blocks.insert(0, para("Proof."))
            if any(c.startswith("callout-") for c in classes):
                title = options.get("title", "Note")
                if options.get("collapse") == "true":
                    blocks.insert(0, raw(f"<details><summary>{html.escape(title)}</summary>"))
                    blocks.append(raw("</details>"))
                else:
                    blocks.insert(0, para(title))
        if node.get("t") == "Header":
            ident = node["c"][1][0]
            if ident.startswith("sec-"):
                REFERENCES[ident] = (route, "Section " + plain(node["c"][2]))
        return node
    doc = walk(doc, labels)
    title = meta.get("title")
    if not title and doc["blocks"] and doc["blocks"][0]["t"] == "Header":
        title = plain(doc["blocks"].pop(0)["c"][2])
    return {"path": path, "route": route, "meta": meta, "doc": doc,
            "title": title or path.stem, "number": number, "body": body}


def plain(value) -> str:
    if isinstance(value, list):
        return "".join(plain(v) for v in value)
    if isinstance(value, dict):
        if value.get("t") in ("Str", "Code", "Math"):
            c = value["c"]
            return c if isinstance(c, str) else c[-1]
        if value.get("t") in ("Space", "SoftBreak", "LineBreak"):
            return " "
        return plain(value.get("c", []))
    return ""


def relative(target: str, route: str) -> str:
    return os.path.relpath(target, str(Path(route).parent)).replace(os.sep, "/")


def render(page: dict, lectures: list[dict]) -> None:
    route = page["route"]
    prefix = "../" * len(Path(route).parent.parts)

    def links(node):
        if node.get("t") == "Cite":
            citations = node["c"][0]
            if any(c["citationId"] in REFERENCES for c in citations):
                if len(citations) != 1:
                    raise ValueError("Use separate links for grouped cross-references")
                ident = citations[0]["citationId"]
                target, label = REFERENCES[ident]
                return {"t": "Link", "c": [["", [], []], inlines(label), [relative(target, route) + "#" + ident, ""]]}
            for citation in citations:
                if citation["citationId"].split("-")[0] in KINDS:
                    raise ValueError(f'Unresolved reference {citation["citationId"]} in {route}')
        if node.get("t") == "Link":
            url = node["c"][2][0]
            if not re.match(r"[a-z]+:", url):
                node["c"][2][0] = re.sub(r"\.(?:qmd|md)(?=#|$)", ".html", url)
        return node
    doc = walk(page["doc"], links)
    nav = [("Schedule", "index.html#course-schedule"),
           ("Course info", "course-info.html"), ("Search", "search.html")]
    section = "index.html#course-schedule" if route == "index.html" or route.startswith(("lessons/", "assignments/")) else route
    navhtml = "".join(f'<a href="{prefix}{url}"' + (' aria-current="page"' if section == url else '') + f'>{label}</a>' for label, url in nav)
    navhtml += f'<a class="canvas-link" href="{html.escape(CONFIG["canvas"])}">Canvas ↗</a>'
    breadcrumb = ""
    pagination = ""
    if page["number"] and page["meta"].get("published", True):
        breadcrumb = f'<a href="{prefix}index.html#course-schedule">Schedule</a> / Lesson {page["number"]}'
        i = lectures.index(page)
        for neighbor, label in [(i - 1, "← Previous"), (i + 1, "Next →")]:
            if 0 <= neighbor < len(lectures):
                other = lectures[neighbor]
                pagination += f'<a href="{relative(other["route"], route)}"><span>{label}</span>{html.escape(other["title"])}</a>'
    meta = {"title": page["title"], "subtitle": page["meta"].get("subtitle", ""),
            "root": prefix, "repository": CONFIG["repository"], "nav": navhtml,
            "pageclass": "home" if route == "index.html" else "reading",
            "eyebrow": "University of Utah · Spring 2027" if route == "index.html" else "MATH 4310 · Computational Bioinformatics",
            "breadcrumb": breadcrumb, "pagination": pagination,
            "source": CONFIG["repository"] + "/blob/main/" + page["path"].relative_to(ROOT).as_posix()}
    doc["meta"].update({key: {"t": "MetaString", "c": str(value)} for key, value in meta.items()})
    for key in ("nav", "breadcrumb", "pagination"):
        doc["meta"][key] = {"t": "MetaBlocks", "c": [raw(meta[key])]}
    args = ["-f", "json", "-t", "html5", "--standalone", "--template", str(ROOT / "site/template.html"),
            "--mathjax=https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js", "--highlight-style=tango",
            "--citeproc", "--bibliography", str(ROOT / "course_notes/references.bib")]
    if page["meta"].get("toc", True):
        args += ["--toc", "--toc-depth=2"]
    target = OUT / route
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(pandoc(json.dumps(doc), *args), encoding="utf-8")


def catalogs() -> dict[str, str]:
    lessons = {}
    assignments = {}
    for path in sorted((ROOT / "course_notes/lessons").glob("*.md")):
        meta, _ = read_source(path)
        lessons[int(re.search(r"lesson(\d+)", path.stem)[1])] = (path, meta)
    for path in sorted((ROOT / "assignments").glob("*.md")):
        assignments[path.stem] = (path, read_source(path)[0])

    def unavailable(label):
        return '<span class="unreleased">' + html.escape(label) + '<span class="sr-only"> (not yet available)</span></span>'

    def entry(item, directory, label=None):
        path, meta = item
        label = label or meta["title"]
        if meta.get("published") is True:
            return f'<a href="{directory}/{path.stem}.html">{html.escape(label)}</a>'
        return unavailable(label)

    def problem_set(group):
        assignment = group["assignment"]
        if not assignment:
            return ""
        item = assignments.get(assignment["id"])
        label = entry(item, "assignments", assignment["label"]) if item else unavailable(assignment["label"])
        due = date.fromisoformat(assignment["due"])
        return label + f'<span class="assignment-due">Due {due:%b} {due.day}</span>'

    def week_label(week):
        start = date.fromisoformat(week["date_start"])
        end = date.fromisoformat(week["date_end"])
        dates = f'{start:%b} {start.day}–'
        dates += str(end.day) if start.month == end.month else f'{end:%b} {end.day}'
        return f'<span class="week-number">{week["week"]:02}</span><span class="week-dates">{dates}</span>'

    schedule_data = json.loads((ROOT / "site/schedule.json").read_text(encoding="utf-8"))
    rows = []
    groups = []
    assignment_items = []
    current_module = None
    for week in schedule_data["weeks"]:
        if not week["lessons"]:
            rows.append(f'<tr class="schedule-break" id="week-{week["week"]}"><th scope="row">{week_label(week)}</th><td colspan="2">{html.escape(week["topic"])}</td></tr>')
            continue
        module = next(m for m in schedule_data["modules"] if week["lessons"][0] in m["lessons"])
        if module != current_module:
            rows.append(f'<tr class="module-heading" id="module-{module["id"]}"><th colspan="3">{html.escape(module["title"])}</th></tr>')
            current_module = module
        meeting_items = []
        for meeting in week["meetings"]:
            day = date.fromisoformat(meeting["date"])
            date_label = f'{day:%a %b} {day.day}'
            if "lesson" in meeting:
                n = meeting["lesson"]
                content = f'<span class="lesson-number">{n:02}</span>' + entry(lessons[n], "lessons")
            else:
                content = '<span class="schedule-event">' + html.escape(meeting["event"]) + '</span>'
            meeting_items.append('<div class="schedule-item"><span class="meeting-date">' + date_label + '</span>' + content + '</div>')
        lecture_list = "".join(meeting_items)
        if week.get("note"):
            lecture_list += '<div class="calendar-note">' + html.escape(week["note"]) + '</div>'
        rows.append(f'<tr class="week-row" id="week-{week["week"]}"><th scope="row">{week_label(week)}</th><td>{lecture_list}</td><td>{problem_set(week)}</td></tr>')
        if week["assignment"]:
            assignment_items.append('<li><h2>Week ' + str(week["week"]) + ' · ' + html.escape(week["topic"]) + '</h2>' + problem_set(week) + '</li>')
    for module in schedule_data["modules"]:
        groups.append("## " + module["title"] + "\n\n" + "\n".join(f'<div class="catalog-item"><span class="lesson-number">{n:02}</span>{entry(lessons[n], "lessons")}</div>' for n in module["lessons"]))
    schedule = '<div class="table-scroll"><table class="course-table"><thead><tr><th scope="col">Week / dates</th><th scope="col">Lectures</th><th scope="col">Assignment</th></tr></thead><tbody>' + "".join(rows) + '</tbody></table></div>'
    assignment_list = '<ul class="resource-list">' + "".join(assignment_items) + '</ul>'
    files = subprocess.check_output(["git", "ls-files", "-z", "data"], cwd=ROOT).decode().split("\0")
    datasets = [f for f in files if f and Path(f).suffix.lower() in (".csv", ".tsv", ".fasta", ".fa", ".json", ".txt")]
    data_list = '<ul>' + "".join(f'<li><a href="{html.escape(f)}">{html.escape(Path(f).name)}</a></li>' for f in datasets) + '</ul>' if datasets else ''
    return {"WEEKLY_SCHEDULE": schedule, "LECTURE_CATALOG": "\n\n".join(groups), "ASSIGNMENT_CATALOG": assignment_list, "DATA_CATALOG": data_list, "DATA_STATUS": "Download the datasets below. Small examples are also included directly in the lecture code." if datasets else "There are no downloadable datasets yet. Small examples are included directly in the lecture code."}


def main() -> None:
    if OUT.resolve() != (ROOT / "_site").resolve():
        raise ValueError("Unexpected output directory")
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "assets/figures").mkdir(parents=True)
    for filename in ("site.css", "site.js"):
        shutil.copy2(ROOT / "site" / filename, OUT / "assets" / filename)
    shutil.copy2(ROOT / "course_notes/img/UU Logo-CTR_RGB.svg", OUT / "assets/utah.svg")
    if (ROOT / "site/assets").exists():
        shutil.copytree(ROOT / "site/assets", OUT / "assets", dirs_exist_ok=True)
    # Only public data already tracked by Git is copied; local files stay local.
    data_files = subprocess.check_output(["git", "ls-files", "-z", "data"], cwd=ROOT).decode().split("\0")
    for filename in filter(None, data_files):
        source = ROOT / filename
        if source.is_symlink():
            raise ValueError(f"Data symlinks are not published: {filename}")
        target = OUT / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    CONTENT.update(catalogs())
    lectures = [prepare(p, "lessons/" + p.stem + ".html", str(i)) for i, p in enumerate(sorted((ROOT / "course_notes/lessons").glob("*.md")), 1)]
    pages = list(lectures)
    for directory, output in [("course_notes/appendices", "appendices/"), ("assignments", "assignments/")]:
        pages += [prepare(p, output + p.stem + ".html") for p in sorted((ROOT / directory).glob("*.md")) if p.name != "README.md"]
    pages.append(prepare(ROOT / "course_notes/using-notes.md", "using-notes.html"))
    for p in sorted((ROOT / "site/pages").glob("*.md")):
        pages.append(prepare(p, p.stem + ".html"))
    for page in pages:
        render(page, [p for p in lectures if p["meta"].get("published", True)])
        print("Built", page["route"], flush=True)
    search = [{"title": p["title"], "url": p["route"], "text": re.sub(r"\s+", " ", plain(p["doc"]["blocks"]))} for p in pages if p["meta"].get("published", True) and p["meta"].get("search", True)]
    (OUT / "search-index.json").write_text(json.dumps(search, ensure_ascii=False), encoding="utf-8")
    (OUT / ".nojekyll").touch()
    print(f"Built {len(pages)} pages in {OUT}")


if __name__ == "__main__":
    main()
