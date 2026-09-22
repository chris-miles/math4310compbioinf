# MATH 4310 — Code Style

Lecture examples use Python and are **intended to be copied and adapted by students**. Assignments may be completed in R, Julia, MATLAB, or another language unless a problem has a specific restriction. Python packages and plotting helpers describe the notes' implementation, not student submission requirements. Write small, readable algorithms for undergraduates. Chris will rarely read submissions line by line; the assessed evidence is reasoning, validation, actual output, and interpretation.

## Stack
- Python ≥ 3.11.
- Core: `numpy`, `scipy`, `pandas`, `matplotlib`.
- Domain as needed: `scikit-learn`, `biopython`, `networkx`, `umap-learn`.
- Never use a library where a 20-line implementation reveals the algorithm. Use the library version only after the from-scratch version exists in the chapter.

## Style
- PEP 8, four-space indent, 88-char soft limit.
- `snake_case` functions and variables, `PascalCase` classes, `ALL_CAPS` module constants.
- Single-letter names only where they match the math ($i, j, k$ indices; $n, m$ lengths; $x, y$ sequences when consistent with chapter notation).
- Type hints on every function signature; none inside function bodies.
- One function per algorithmic idea. No clever one-liners.
- Favor explicit loops, small helper functions, and named intermediate variables over compact comprehensions or chained one-liners. Students should be able to step through copied code by hand and modify a single input without untangling a dense expression.

## Comments
Comments explain **why**, the **biological meaning**, or the **link to the math** — never restate the code. Aim for one comment per logical block, not per line. Every non-trivial function gets a one-line docstring; longer docstrings only when arguments are non-obvious.

## Pseudocode vs Python
Introduce every algorithm with pseudocode first in a `.algorithm` div, then give the Python implementation in a code chunk. Pseudocode uses mathematical notation; Python uses our naming conventions.

## Figures from Code
Apply `figstyle.mplstyle` once at lesson top (see `STYLE_FIGURES.md`). **Never** set `figsize`, `dpi`, or other figure-size rcParams in chunk code — The site builder controls those via `site/config.json` and chunk options.

## Executable Markdown chunks
- Active lesson sources are `.md` files. Use a fenced block marked `{.python .execute}` for code the website runs, and a plain `python` fence for nonexecuting examples such as deliberate assignment bugs. The builder retains `#|` options; `.qmd` files are legacy.
- Chunk labels: `code-{chapter}-{shortname}` for code, `fig-{chapter}-{shortname}` when the chunk produces a figure.
- Figure-producing chunks set `fig-cap`, `fig-alt`, and `label`.
- `echo: true` by default — students see the code.
- Imports at chapter top in a single `setup` chunk; later chunks assume those imports.

## Reproducibility
- Set `np.random.seed(0)` (or a chapter-specific seed) at the top of any chunk using randomness.
- Tiny illustrative inputs inline; larger inputs loaded from `data/`. Never download at render time.

## Snippets for Assignments
Mark starter chunks students will modify:
```python
# --- starter for Problem Set N, Problem M ---
def viterbi(...):
    raise NotImplementedError
```
Put incomplete starter functions and intentional bugs in nonexecuting code fences. Every executable block must run during a clean site build. Show expected and actual values for checks; avoid assertions that merely repeat the implementation.
