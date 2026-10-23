# MATH 4310 — Figure Style

Scope: data plots (lines, points, bars, heatmaps, scree plots). Schematics (DP tables, trellises, graphs, trees, alignments, phase planes) follow `STYLE_SCHEMATICS.md` and are drawn through `style/coursefigs.py`; the two share `figstyle.mplstyle`.

Target aesthetic: **clean scientific** — minimal chrome, single sans-serif, colorblind-safe palette, every mark earns its place. All figures use matplotlib with the shared style file `figstyle.mplstyle`.

## Setup
```python
import matplotlib.pyplot as plt
plt.style.use("course_notes/style/figstyle.mplstyle")
```

## Sizes — Do Not Set in Python
**Never** call `plt.figure(figsize=...)`, pass `figsize=` to `plt.subplots`, or override `figure.figsize` in rcParams. The website builder controls canvas size via chunk-level `fig-width` / `fig-height`, with project defaults in `site/config.json`:

```json
"figure": {"width": 6, "height": 4, "dpi": 150}
```

Override per-chunk only for genuinely wide figures (multi-panel comparisons) and for schematics, whose aspect ratio follows the object (see `STYLE_SCHEMATICS.md`). HTML responsive scaling is handled by the site CSS  regardless of these values — they exist to set the rendered text-to-plot ratio, not to control display width.

## Palette
Categorical: a six-color Okabe–Ito subset in fixed order (blue, orange, green, pink, vermillion, sky), validated for colorblind separation. The order lives in `figstyle.mplstyle`; do not redeclare hex values in chapter code, and never cycle past six series (facet or aggregate instead). The first color is blue, so a single-series plot is blue, not black. For emphasis inside a data plot use `coursefigs.EMPHASIS` (Utah red) and nothing else.

**Assignment rule:** when a chapter uses categorical color, fix the mapping on first use and reuse the same color for the same category across every figure in that chapter.

Sequential colormap: `viridis`. Diverging: `RdBu_r` (centered at 0). Never use `jet`, `rainbow`, or `hsv`.

## Type
- Family: sans-serif (Arial / Helvetica / DejaVu Sans fallback), set in `figstyle.mplstyle`.
- Axis labels 10 pt, tick labels 9 pt, panel labels (A, B, C) 11 pt bold top-left outside axes.
- Math rendered via mathtext (no LaTeX dependency).

## Axes
- Spines: left and bottom only.
- Ticks: outward, on left and bottom only.
- No gridlines unless conveying a real reference (e.g. probability $= 0.5$).
- Let data dictate limits unless a comparison demands fixed limits.

## Marks
- Lines: 1.5 pt default.
- Points: filled, 5 pt; encode categories by **marker shape *and* color** (accessibility).
- Bars: `edgecolor='none'`, `alpha=0.9`.
- Heatmaps: square cells when conceptually symmetric (DP tables, score matrices); always include a colorbar.

## Captions
Set `fig-cap` in the chunk header — that is the caption mechanism. Every caption has two parts: (1) what the figure shows technically, (2) one short sentence on the biological or algorithmic takeaway. Also set `fig-alt` (accessibility) and `label: fig-{chapter}-{shortname}` for cross-references.

## When to Make a Figure
Every figure must answer a specific question. Skip ornamentation. A DP-table heatmap is worth it; a generic "schematic of the pipeline" usually is not. Schematics that *are* worth it (HMM state diagrams, de Bruijn graphs) go through matplotlib + networkx, not external editors, so they regenerate from source. The most valuable figure shows why something is true: a path through the table, evidence accumulating along a sequence, points projected onto an axis. A table of computed values shows the result, not the idea. Never add a picture for decoration.

## Borrowed Figures
A published figure can show real data or an instrument better than a simulation: a sequencing run, a genome-browser track, a single-cell atlas. Use one only when its license permits reuse (CC BY, CC BY-SA, CC0, or public domain such as NIH/NHGRI images); textbook figures are not reusable. Save it under `site/assets/borrowed/`, record its source URL and license in `site/assets/borrowed/README.md`, and place it as the builder places generated figures: a `::: {#fig-... .course-figure}` div holding an `<img>` tag with `alt` text, then the caption paragraph. The caption ends with the credit and license: "Reproduced from @key (CC BY 4.0)." When the license is unclear, redraw the idea from public data with `coursefigs` and credit it "after @key."

## File Conventions
- Figures generated inline in executable Markdown chunks; no committed PNG/PDF unless hand-drawn or third-party.
- Chunk labels: `fig-{chapter}-{shortname}`, kebab-case (e.g. `fig-align-nw-table`).
