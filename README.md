# MATH 4310 Computational Bioinformatics

Course materials for **MATH 4310 Computational Bioinformatics** at the University of Utah, taught by Chris Miles, Spring 2027.

**Course website:** <https://chris-miles.github.io/math4310compbioinf/>

## Course materials

- [Lecture sources](course_notes/lessons/) — editable Markdown with math and Python examples.
- [Reference pages](course_notes/appendices/) — computing conventions, math refreshers, and notation.
- [Assignments](assignments/) — problem sets and project workshops; unreleased drafts are marked in front matter.
- [Data](data/) — small public datasets for lessons and assignments.

The website groups lectures and assignments by week. Unreleased material appears as a muted title until it is ready.

## Build and preview

Install [Pandoc](https://pandoc.org/installing.html) and Python 3.11 or newer, then run from the repository root:

~~~sh
python -m pip install -r site_tools/requirements.txt
python site_tools/build_site.py
python site_tools/check_site.py
python -m http.server 8000 --directory _site
~~~

Open <http://localhost:8000>. The output is plain HTML, CSS, and JavaScript in the ignored _site/ directory. Quarto is not required. GitHub Actions builds, checks, and deploys the site on pushes to main.

See [the authoring guide](site/AUTHORING.md) for editing pages, releasing lectures, adding datasets, and embedding interactive content. The existing [prose, code, and figure conventions](course_notes/style/) remain in force.
