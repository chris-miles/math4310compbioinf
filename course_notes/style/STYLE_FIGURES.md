# MATH 4310 — Figure Style

Target aesthetic: **clean scientific** — minimal chrome, single sans-serif, colorblind-safe palette, every mark earns its place. All figures use matplotlib with the shared style file `figstyle.mplstyle`.

## Setup
```python
import matplotlib.pyplot as plt
plt.style.use("course_notes/style/figstyle.mplstyle")
```

## Sizes — Do Not Set in Python
**Never** call `plt.figure(figsize=...)`, pass `figsize=` to `plt.subplots`, or override `figure.figsize` in rcParams. Quarto controls canvas size via chunk-level `fig-width` / `fig-height`, with project defaults in `_quarto.yml`:

```yaml
format:
  html:
    fig-width: 6
    fig-height: 4
    fig-dpi: 150
  pdf:
    fig-width: 5.5
    fig-height: 3.5
    fig-dpi: 300
```

Override per-chunk only for genuinely wide figures (multi-panel comparisons). HTML responsive scaling is handled by Quarto/Bootstrap CSS (`.img-fluid { max-width: 100%; height: auto; }`) regardless of these values — they exist to set the rendered text-to-plot ratio, not to control display width.

## Palette
Categorical: **Okabe–Ito** (colorblind-safe, 8 colors). The order is canonical and lives in `figstyle.mplstyle` — do not redeclare hex values in chapter code.

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
Every figure must answer a specific question. Skip ornamentation. A DP-table heatmap is worth it; a generic "schematic of the pipeline" usually is not. Schematics that *are* worth it (HMM state diagrams, de Bruijn graphs) go through matplotlib + networkx, not external editors, so they regenerate from source.

## File Conventions
- Figures generated inline in Quarto chunks; no committed PNG/PDF unless hand-drawn or third-party.
- Chunk labels: `fig-{chapter}-{shortname}`, kebab-case (e.g. `fig-align-nw-table`).
