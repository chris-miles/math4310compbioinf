# MATH 4310 Computational Bioinformatics

Course notes and planning materials for MATH 4310 Computational Bioinformatics.

The student-facing Quarto book lives in `course_notes/`. It is built as an HTML-only site; PDF output is intentionally not configured.

Public repository:

- https://github.com/chris-miles/math4310compbioinf

Public site:

- https://chris-miles.github.io/math4310compbioinf/

## Local Preview

```powershell
quarto render .\course_notes
```

Open `course_notes/_book/index.html`.

## Public Site

GitHub Pages is configured to deploy from GitHub Actions. The workflow in `.github/workflows/publish.yml`:

1. installs Quarto,
2. installs the Python packages in `course_notes/requirements.txt`,
3. renders `course_notes/`, and
4. publishes `course_notes/_book`.

Generated Quarto output, local caches, private reference materials, PDFs, DOCX files, and large archives are ignored by git. Small public-safe datasets can go in `data/`.

## Current Book Status

The book scaffold is organized as numbered lessons inside conceptual modules, with assignments and appendices separated from the main lesson sequence. The landing page now gives a rough week-by-week topic timeline, and `using-notes.qmd` gives students a small code example for practicing with executable chunks.
