# Editing the course website

## Where to edit

| Content | Source |
|---|---|
| Homepage and course information | site/pages/*.md |
| Topics, week ranges, and problem-set links | site/schedule.json |
| Course links and default figure dimensions | site/config.json |
| Lecture content | course_notes/lessons/*.md |
| Assignment content | assignments/*.md |
| Layout and visual presentation | site/template.html, site/site.css |
| Browser interactions | site/site.js |

The prose, code, figure, and schematic style guides in course_notes/style/ still govern the notes. Their historical references to Quarto describe the source conventions retained below; the website is now built with Pandoc and Python. The palette, plot typography, helper functions, and per-chunk figure dimensions are unchanged.

## Release a lecture or assignment

Each lecture and assignment has YAML front matter:

~~~yaml
---
title: "Forward–Backward and Posterior Decoding"
subtitle: "Lesson 13 · Week 7"
published: false
---
~~~

Change **published: false** to **published: true** to make the title clickable in the schedule and catalog. Released pages enter search and lecture navigation. Draft text is excluded from the generated website: old URLs show a short availability notice. **The repository is public, so draft sources are still public.** Keep instructor-only material outside tracked files.

Only set this flag after reviewing the page. Existing lessons 1–12 are initially linked; the explicit scaffolds in lessons 13–24 and the assignment drafts are not linked.

The schedule has broad modules and compact weekly rows. Edit site/schedule.json: the modules list groups lecture numbers by subject, and the weeks list assigns lectures to weeks and places 11 problem sets at topic boundaries. Set assignment to null for a week without a new problem set. Each row also carries date_start and date_end; holiday notes are optional. The calendar source and term boundaries are recorded in the same file. The short final week has a class-meeting entry without a new problem set. An assignment ID matches its filename without .md; entries without a source file remain unlinked. The homepage, lecture index, and assignment catalog share this schedule.

## Markdown and mathematics

Use ordinary headings, links, lists, tables, and dollar-delimited LaTeX. Citations remain [@needleman1970], with entries in course_notes/references.bib.

Numbered blocks and references retain the existing syntax:

~~~md
::: {#def-example}
## A named definition
The definition goes here.
:::

See @def-example.
~~~

Definitions, theorems, propositions, examples, exercises, figures, and labeled equations are numbered per lecture. Keep IDs unique across the course. Unresolved references fail the build. Standard Markdown footnotes, raw HTML, fenced divs, and collapsible callouts also work.

## Python examples and figures

A code block marked with {.python .execute} runs during the build. A plain python code fence is displayed without execution. Chunks on a page share a namespace, run in order with that page's directory as the working directory, and display their final expression and printed output. Use the existing lecture files as examples.

Keep using coursefigs.use_style(), figstyle.mplstyle, and the existing #| options: label, fig-cap, fig-alt, fig-width, fig-height, code-fold, code-summary, echo, and eval. Figure captions and alt text are required. Default dimensions are 6 × 4 inches at 150 dpi. Generated figures are SVGs and are not committed. The builder executes trusted course code; never build unreviewed code with access to credentials.

## HTML and interactive content

Raw HTML can be embedded directly in Markdown. For a reusable interactive activity, add a JavaScript file under site/assets/, then reference it from the page with a relative URL. The builder copies this public asset directory into _site/assets/. Use semantic controls, keyboard support, descriptive labels, and a useful static explanation. No framework is required.

For example, a native disclosure works without JavaScript:

~~~html
<details>
<summary>Reveal the answer</summary>
<p>The Hamming distance is 1.</p>
</details>
~~~

## Data

Place only public-safe data in data/ and add it to Git. The build copies tracked files only; local files are not published by accident. CSV, TSV, FASTA, JSON, and text files appear automatically on the data page. Link a dataset from a lecture with ../data/filename.csv.

## Build and publish

From the repository root:

~~~sh
python -m pip install -r site_tools/requirements.txt
python site_tools/build_site.py
python site_tools/check_site.py
python -m http.server 8000 --directory _site
~~~

Install Pandoc separately. The build replaces only the ignored _site/ output directory. The checker validates internal links and fragments, draft exclusion, figure alt text, and the publication boundary. GitHub Actions uses the same build and checker before deploying to GitHub Pages.

Commit only intended public changes and push to main to deploy. The shared HTML template uses relative URLs, including when the site is hosted under the repository's GitHub Pages path. Existing lecture .html URLs remain valid.

## Legacy Quarto files

The original .qmd files and Quarto configuration are retained for reference. They are not used by the new build and do not update automatically. Make future content edits in the .md files; edit the course homepage in site/pages/index.md.
