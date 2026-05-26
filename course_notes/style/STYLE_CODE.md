# MATH 4310 — Code Style

Python only. Code is in service of understanding the algorithm and is **explicitly intended to be copied by students** into assignments. Write for a strong undergraduate, not for production.

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

## Comments
Comments explain **why**, the **biological meaning**, or the **link to the math** — never restate the code. Aim for one comment per logical block, not per line. Every non-trivial function gets a one-line docstring; longer docstrings only when arguments are non-obvious.

## Pseudocode vs Python
Introduce every algorithm with pseudocode first in a `.algorithm` div, then give the Python implementation in a code chunk. Pseudocode uses mathematical notation; Python uses our naming conventions.

## Figures from Code
Apply `figstyle.mplstyle` once at chapter top (see `STYLE_FIGURES.md`). **Never** set `figsize`, `dpi`, or other figure-size rcParams in chunk code — Quarto controls those via `_quarto.yml` and chunk options.

## Chunks (Quarto)
- Chunk labels: `code-{chapter}-{shortname}` for code, `fig-{chapter}-{shortname}` when the chunk produces a figure.
- Figure-producing chunks set `fig-cap`, `fig-alt`, and `label`.
- `echo: true` by default — students see the code.
- Imports at chapter top in a single `setup` chunk; later chunks assume those imports.

## Reproducibility
- Set `np.random.seed(0)` (or a chapter-specific seed) at the top of any chunk using randomness.
- Tiny illustrative inputs inline; larger inputs loaded from `course_notes/data/`. Never download at render time.

## Snippets for Assignments
Mark starter chunks students will modify:
```python
# --- starter for Problem Set N, Problem M ---
def viterbi(...):
    raise NotImplementedError
```
These must run as-is (raising at the right line is acceptable).
