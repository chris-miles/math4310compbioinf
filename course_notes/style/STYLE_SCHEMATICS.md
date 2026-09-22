# MATH 4310 — Schematics Style

Most figures in these notes are not data plots. They are dynamic-programming tables, HMM trellises and state diagrams, k-mer and neighborhood graphs, trees, printed alignments, and phase planes. `STYLE_FIGURES.md` covers data plots. This file covers everything else, and it is strict: every schematic is drawn in Python through `style/coursefigs.py`, from the same instance the prose uses, and nothing is drawn by hand or in an external editor.

The reason is consistency at scale. Twenty-four lessons will draw something like sixty schematics. If each one is composed ad hoc, fonts, arrow styles, and the meaning of red drift within a week, and students lose the visual vocabulary that lets them read a Viterbi trellis as the alignment table it is.

## The one rule

Call a `coursefigs` helper. If no helper does what the lesson needs, add one to `coursefigs.py` (small, typed, docstring, same color roles) rather than composing matplotlib calls inside the lesson. A lesson chunk should read as "here is the instance, draw it," not as drawing code.

```python
import sys
sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

That goes in the lesson's setup chunk once. `cf.use_style()` applies `figstyle.mplstyle`, so data plots and schematics share fonts and line weights.

## Color roles

Schematics use roles, not a palette. A role means the same thing on every page.

| Role | Hex | Meaning |
|---|---|---|
| `INK` | `#000000` | Text, node outlines, anything that is "the object" |
| `RULE` | `#707271` | Structure: cell borders, ordinary edges, back pointers, branches |
| `LIGHT` | `#E2E6E6` | Background fill for "the cell being computed" and for stream lines behind a phase plane |
| `PAPER` | `#FFFFFF` | Node and cell fill |
| `EMPHASIS` | `#BE0000` | The thing the figure is about: the optimal path, the selected edge, the fixed point |
| `SECOND` | `#0072B2` | A genuinely different second thing (the $x$-nullcline; the other of two models) |
| `THIRD` | `#E69F00` | A third thing, rarely |

Rules that follow from this:

- Red means "look here." It never decorates. A figure with two red things is making two points and should be two figures.
- Grey structure is recessive. If a grey element carries meaning, that is a sign it wants to be red or wants a label.
- `SECOND` and `THIRD` are validated against each other and against `EMPHASIS` for colorblind separation; do not introduce other hues in a schematic. Data plots use the categorical series in `figstyle.mplstyle`, which is a different job.
- Text is always ink or paper (on a red fill). Text is never colored to match a series.

## Typography

One font for everything, set by `figstyle.mplstyle` (Arial, then Helvetica, then DejaVu Sans). Sequence letters and k-mers are monospace (`DejaVu Sans Mono`) so that columns line up and a `-` gap reads as a gap. Math goes through mathtext, so `$F(i,j)$` renders the same in a figure as in prose. Sizes: labels 10 pt, values inside cells 8–9 pt, never smaller.

## Conventions by object

These are fixed so that the same picture means the same thing in every lesson. The helper enforces most of them; the rest is on the author.

**Dynamic-programming tables** (`dp_table`). Rows are $x$, columns are $y$, and the boundary row and column for the empty prefix are always shown with a `$-$` label. Cell $(i,j)$ in the figure is $F(i,j)$ in the prose, with the same indices. The optimal path is red-filled cells. The cell currently being computed is `LIGHT` with red arrows coming in from its predecessors. Show at most one of those two things at a time unless the point is their relationship.

**Trellises** (`trellis`). Rows are states, columns are observations, entries are circles. Same conventions as tables (red path, grey back pointers) on purpose: students should see Viterbi as the table with a different rule in each cell. Show scores in log space when the lesson works in log space.

**Graphs** (`digraph`). Nodes are white circles with ink outlines and monospace labels. Edges are grey with arrowheads; parallel edges are drawn as separate arcs because multiplicity is the point in a de Bruijn graph. Highlight one path or one node in red. Neighborhood graphs pass `directed=False`. Fix `pos` by hand for small graphs so the layout matches how the instructor draws it on the board; use the returned positions for a second panel of the same graph.

**HMM state diagrams** (`state_diagram`). States in a row, self-loops above, emission tables below in monospace, transition probabilities on the arcs. The same B/I toy model is used from Lesson 11 through Lesson 13 and its diagram is identical on all three pages.

**Trees** (`tree`). Root at the left, leaves at the right in the order given, branch lengths as horizontal extent so that an ultrametric tree visibly lines its leaves up. Highlight a leaf or clade in red when the lesson is about it. Internal node names, if any, are small and grey.

**Alignments** (`alignment`). Two monospace rows with a match line between them: `|` for a match, nothing for a gap column, a red dot for a mismatch. Row labels are the sequence names from the prose.

**Phase planes** (`phase_plane`). Stream lines in `LIGHT` so they recede, the $\dot{x}=0$ nullcline in `SECOND`, the $\dot{y}=0$ nullcline in `THIRD`, fixed points in red (filled stable, open unstable), example trajectories in ink. Axes are labeled with the state variables in the lesson's notation.

## Sizing

Schematics are the exception to the "never set figure size" rule in `STYLE_FIGURES.md`: their aspect ratio is dictated by the object, so set `fig-width` and `fig-height` per chunk. A DP table for strings of length 4 is roughly square; a state diagram is wide and short. Keep width at or below 6 inches. Never set `figsize` in Python.

## Captions and alt text

Same rule as data plots: the caption says what the figure shows and then, in one sentence, what to take from it. `fig-alt` describes the structure ("a five-by-four grid with a red path along the diagonal and one shaded cell"), not the takeaway. Label the chunk `fig-{lesson-slug}-{shortname}`.

## Same instance, everywhere

A lesson's figure draws the instance the prose is working through, with the same strings, scores, and indices. Do not invent a second instance for the picture. Across lessons, reuse instances deliberately: the `ACGT`/`AGT` table from Lesson 3 is traced back in Lesson 4; the B/I HMM is defined in Lesson 11, decoded in Lesson 12, and given posteriors in Lesson 13. The example banks in `planning/source_digests/` list which instances are shared.

## Checklist before a figure ships

- Drawn through a `coursefigs` helper; no bare matplotlib in the lesson chunk.
- Exactly one red idea.
- Same instance as the surrounding prose, same index convention.
- Monospace for sequence letters; mathtext for symbols.
- `label`, `fig-cap`, `fig-alt` set; `fig-width`/`fig-height` set for schematics.
- Rendered and looked at. The helper cannot check for a label sitting on an arrow.
