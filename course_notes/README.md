# MATH 4310 Computational Bioinformatics Notes

This folder is the working area for developing the course notes for MATH 4310.

Structure:

- `_quarto.yml` — Quarto book configuration.
- `index.qmd` — book landing page with a rough week-by-week topic timeline.
- `using-notes.qmd` — student-facing guide to reading lessons and working with code chunks.
- `lessons/` — numbered student-facing Quarto lessons, grouped by conceptual module in the book navigation.
- `assignments/` — student-facing assignment and project-workshop shells.
- `appendices/` — computing conventions and notation appendices.
- `data/` — small public-safe datasets for examples and assignments.
- `img/` — image assets, currently including the University of Utah logo used in the sidebar.
- `style/` — **style guides (read first)**: prose, figures, code, and the shared matplotlib style file. Any human or LLM drafting chapter content must read these before writing.
- `planning/` — instructor-facing course architecture, daily schedule, homework plan, project plan, and reference map.

The immediate goal is to replace lesson scaffolds with short student-facing readings, starting with the foundations and pairwise-alignment modules and the first two problem sets.

## Build Status

The Quarto book is HTML-only. The active HTML parameters live in `_quarto.yml`:

- `output-dir: _book`
- sidebar logo: `img/UU Logo-CTR_RGB.svg`
- author: Chris Miles
- footer email: `chris.miles@utah.edu`
- `number-sections: true`
- `code-copy: true`
- `code-overflow: wrap`
- default figure size: 6 by 4 inches at 150 dpi

The public site is built by GitHub Actions from the repository root and publishes `course_notes/_book` to GitHub Pages.
