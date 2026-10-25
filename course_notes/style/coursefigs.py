"""Shared drawing helpers for MATH 4310 course-note figures.

Every schematic in the notes (dynamic-programming tables, HMM trellises and
state diagrams, de Bruijn and neighborhood graphs, trees, alignments, phase
planes) is drawn through these functions so that fonts, colors, line weights,
and arrow styles are identical across lessons. The functions are deliberately
small; copy and modify one if a lesson needs something they do not do.

Usage at the top of a lesson:

    import sys
    sys.path.insert(0, "../style")
    import coursefigs as cf
    cf.use_style()

Never set figure sizes here or in lesson code; Quarto controls them.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

STYLE_PATH = Path(__file__).with_name("figstyle.mplstyle")

# --- Color roles -------------------------------------------------------------
# A schematic has at most three colored roles on top of ink and grey structure.
# EMPHASIS is the thing the figure is about (the optimal path, the current cell).
# SECOND and THIRD are for a genuinely different second and third thing.
# The categorical SERIES list is for data plots with several series; it is the
# same order as the prop_cycle in figstyle.mplstyle.
INK = "#000000"
RULE = "#707271"
LIGHT = "#E2E6E6"
PAPER = "#FFFFFF"
EMPHASIS = "#BE0000"
SECOND = "#0072B2"
THIRD = "#E69F00"
SERIES = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9"]

MONO = "DejaVu Sans Mono"


def use_style() -> None:
    """Apply the shared matplotlib style. Call once per lesson."""
    plt.style.use(STYLE_PATH)


def sequence_tracks(
    ax: plt.Axes,
    rows: dict[str, str],
    *,
    highlighted_row: str | None = None,
    highlighted_positions: Sequence[int] = (),
    show_positions: bool = True,
) -> None:
    """Align observed symbols and proposed labels at the same sequence positions."""
    _bare(ax)
    lengths = {len(sequence) for sequence in rows.values()}
    if len(lengths) != 1 or not rows:
        raise ValueError("Sequence tracks need nonempty rows of equal length.")
    length = lengths.pop()
    highlighted_positions = set(highlighted_positions)
    for row, (label, sequence) in enumerate(rows.items()):
        y = -row
        ax.text(-0.7, y, label, ha="right", va="center", fontsize=10)
        for i, symbol in enumerate(sequence):
            selected = label == highlighted_row and i in highlighted_positions
            if selected:
                ax.add_patch(mpatches.Rectangle(
                    (i - 0.36, y - 0.32), 0.72, 0.64,
                    facecolor=EMPHASIS, edgecolor="none", zorder=1,
                ))
            ax.text(i, y, symbol, ha="center", va="center", fontsize=11,
                    family=MONO, color=PAPER if selected else INK, zorder=2)
    if show_positions:
        for i in range(length):
            ax.text(i, 0.75, str(i + 1), ha="center", va="center", fontsize=9)
        ax.text(-0.7, 0.75, "Position", ha="right", va="center", fontsize=9)
    ax.set_xlim(-3.5, max(1, length - 0.4))
    ax.set_ylim(-len(rows) + 0.4, 1.2 if show_positions else 0.6)


def read_assignment(
    ax: plt.Axes,
    references: dict[str, str],
    reads: Sequence[str],
    assignments: Sequence[str],
) -> None:
    """Show exact toy cDNA-read matches and the resulting sample-by-gene count row."""
    _bare(ax)
    if len(reads) != len(assignments) or not references:
        raise ValueError("Supply a gene assignment for every read.")
    grouped = {gene: [] for gene in references}
    for read, gene in zip(reads, assignments):
        if gene not in references or read not in references[gene]:
            raise ValueError("Every assigned read must match its reference segment.")
        grouped[gene].append(read)
    pitch = 0.24
    width = max(len(sequence) for sequence in references.values()) * pitch
    spacing = width + 0.9
    max_reads = max(1, max(len(group) for group in grouped.values()))
    count_y = -max_reads - 1.1
    ax.text(-0.55, 0, "Reference\ncDNA segments", ha="right",
            va="center", fontsize=10)
    ax.text(-0.55, -(max_reads + 1) / 2, "Read records\nfrom sample s1",
            ha="right", va="center", fontsize=10)
    ax.text(-0.55, count_y, "Count row s1", ha="right", va="center", fontsize=10)
    for column, (gene, sequence) in enumerate(references.items()):
        start = column * spacing
        center = start + (len(sequence) - 1) * pitch / 2
        ax.text(center, 0.9, gene, ha="center", va="center", fontsize=10)
        for i, base in enumerate(sequence):
            ax.text(start + pitch * i, 0, base, ha="center", va="center",
                    family=MONO, fontsize=10)
        for row, read in enumerate(grouped[gene]):
            offset = sequence.index(read)
            for i, base in enumerate(read):
                ax.text(start + pitch * (offset + i), -row - 1.1, base,
                        ha="center", va="center", family=MONO, fontsize=10)
        if not grouped[gene]:
            ax.text(center, -1.1, "no assigned reads", ha="center",
                    va="center", fontsize=9)
        ax.annotate("", xy=(center, count_y + 0.45),
                    xytext=(center, -max_reads - 0.35),
                    arrowprops=dict(arrowstyle="->", color=RULE, lw=1))
        ax.add_patch(mpatches.Rectangle(
            (center - 0.42, count_y - 0.32), 0.84, 0.64,
            facecolor=LIGHT, edgecolor=RULE, linewidth=0.8,
        ))
        ax.text(center, count_y, str(len(grouped[gene])), ha="center",
                va="center", fontsize=11, color=INK)
    ax.set_xlim(-2.5, (len(references) - 1) * spacing + width + 0.25)
    ax.set_ylim(count_y - 0.7, 1.4)


def pca_projection(
    ax: plt.Axes,
    points: np.ndarray,
    direction: Sequence[float],
    labels: Sequence[str] = (),
) -> None:
    """Draw centered two-gene samples and their perpendicular line projections."""
    points = np.asarray(points, dtype=float)
    direction = np.asarray(direction, dtype=float)
    direction = direction / np.linalg.norm(direction)
    scores = points @ direction
    projected = scores[:, None] * direction
    extent = max(1.0, float(np.max(np.abs(points)))) * 1.4
    line = np.array([-extent, extent])[:, None] * direction
    ax.plot(line[:, 0], line[:, 1], color=RULE, linewidth=1)
    for i, (point, projection) in enumerate(zip(points, projected)):
        ax.plot([point[0], projection[0]], [point[1], projection[1]],
                color=RULE, linewidth=1)
        if labels:
            ax.annotate(labels[i], point, xytext=(4, 5),
                        textcoords="offset points", fontsize=9)
    ax.scatter(points[:, 0], points[:, 1], color=INK, s=25,
               label="sample", zorder=3)
    ax.scatter(projected[:, 0], projected[:, 1], color=EMPHASIS,
               marker="x", s=30, linewidths=1.5, label="projection", zorder=4)
    ax.plot(0, 0, marker="+", color=RULE, markersize=7, markeredgewidth=1)
    ax.set(xlim=(-extent, extent), ylim=(-extent, extent),
           xlabel="Centered gene 1", ylabel="Centered gene 2")
    ax.set_aspect("equal", adjustable="box")


def _bare(ax: plt.Axes) -> None:
    """Remove axes chrome for schematics that are not data plots."""
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_aspect("equal")


# --- Dynamic-programming tables ---------------------------------------------

def dp_table(
    ax: plt.Axes,
    rows: Sequence[str],
    cols: Sequence[str],
    values: Sequence[Sequence[float | int | str | None]],
    path: Iterable[tuple[int, int]] = (),
    arrows: Iterable[tuple[tuple[int, int], tuple[int, int]]] = (),
    highlight: Iterable[tuple[int, int]] = (),
    fmt: str = "{}",
    corner: str = "",
) -> None:
    """Draw a (len(rows)+1) x (len(cols)+1) table with row and column labels.

    Cell (i, j) is row i, column j with i, j starting at 0 for the boundary
    row and column, matching the recurrence indices in the lessons. `path`
    cells are filled in the emphasis color; `arrows` are drawn from the first
    cell to the second (predecessor to successor); `highlight` cells get a
    light fill for "the cell we are computing now".
    """
    _bare(ax)
    n_rows = len(rows) + 1
    n_cols = len(cols) + 1
    values = np.asarray(values, dtype=object)
    path = set(path)
    highlight = set(highlight)

    for i in range(n_rows):
        for j in range(n_cols):
            x, y = j, -i
            if (i, j) in path:
                fill, text = EMPHASIS, PAPER
            elif (i, j) in highlight:
                fill, text = LIGHT, INK
            else:
                fill, text = PAPER, INK
            ax.add_patch(mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                            facecolor=fill, edgecolor=RULE,
                                            linewidth=0.8))
            value = values[i][j]
            if value is not None:
                ax.text(x, y, fmt.format(value), ha="center", va="center",
                        color=text, fontsize=9)

    # labels: the boundary row/column stand for the empty prefix
    ax.text(-1, 0, corner or "$-$", ha="center", va="center", fontsize=10)
    for i, label in enumerate(rows, start=1):
        ax.text(-1, -i, label, ha="center", va="center", fontsize=10, family=MONO)
    ax.text(0, 1, "$-$", ha="center", va="center", fontsize=10)
    for j, label in enumerate(cols, start=1):
        ax.text(j, 1, label, ha="center", va="center", fontsize=10, family=MONO)

    for (i0, j0), (i1, j1) in arrows:
        ax.annotate("", xy=(j1, -i1), xytext=(j0, -i0),
                    arrowprops=dict(arrowstyle="->", color=EMPHASIS,
                                    lw=1.5, shrinkA=9, shrinkB=9))

    ax.set_xlim(-1.6, n_cols - 0.4)
    ax.set_ylim(-n_rows + 0.4, 1.6)


def trellis(
    ax: plt.Axes,
    states: Sequence[str],
    observations: Sequence[str],
    values: Sequence[Sequence[float | str | None]],
    backpointers: Iterable[tuple[tuple[int, int], tuple[int, int]]] = (),
    path: Iterable[tuple[int, int]] = (),
    fmt: str = "{:.3g}",
) -> None:
    """Draw an HMM trellis: one row per state, one column per observation.

    `values[k][t]` is the entry for state k at position t (0-based). Back
    pointers and path use the same (k, t) indices. Same look as `dp_table`
    so that students see Viterbi as the alignment table with a different
    recurrence.
    """
    _bare(ax)
    values = np.asarray(values, dtype=object)
    path = set(path)
    n_states, n_obs = len(states), len(observations)

    nodes = {}
    for k in range(n_states):
        for t in range(n_obs):
            x, y = t, -k
            on_path = (k, t) in path
            node = mpatches.Circle((x, y), 0.34,
                                  facecolor=EMPHASIS if on_path else PAPER,
                                  edgecolor=EMPHASIS if on_path else RULE,
                                  linewidth=1.0, zorder=3)
            nodes[k, t] = node
            ax.add_patch(node)
            value = values[k][t]
            if value is not None:
                ax.text(x, y, fmt.format(value) if not isinstance(value, str) else value,
                        ha="center", va="center", fontsize=8,
                        color=PAPER if on_path else INK)

    for k, label in enumerate(states):
        ax.text(-1, -k, label, ha="center", va="center", fontsize=10)
    for t, label in enumerate(observations):
        ax.text(t, 1, label, ha="center", va="center", fontsize=10, family=MONO)

    for (k0, t0), (k1, t1) in backpointers:
        ax.annotate("", xy=(t1, -k1), xytext=(t0, -k0),
                    arrowprops=dict(arrowstyle="->", color=RULE, lw=1.0,
                                    patchA=nodes[k0, t0], patchB=nodes[k1, t1],
                                    shrinkA=2, shrinkB=2), zorder=2)

    ax.set_xlim(-1.6, n_obs - 0.4)
    ax.set_ylim(-n_states + 0.4, 1.6)


# --- Graphs ------------------------------------------------------------------

def digraph(
    ax: plt.Axes,
    edges: Sequence[tuple[str, str]],
    pos: dict[str, tuple[float, float]] | None = None,
    edge_labels: dict[tuple[str, str], str] | None = None,
    highlight_edges: Iterable[tuple[str, str] | tuple[str, str, int]] = (),
    highlight_nodes: Iterable[str] = (),
    node_size: int = 900,
    directed: bool = True,
) -> dict[str, tuple[float, float]]:
    """Draw a graph with circular white nodes and grey edges.

    Directed by default, with arrowheads; pass `directed=False` for a kNN or
    similarity graph. Repeated edges in `edges` are drawn as parallel arcs.
    Highlighted edges and nodes use the emphasis color. A highlight (u, v)
    selects all parallel copies; (u, v, key) selects one copy, numbered from
    zero in input order for that pair. Returns node positions for reuse.
    """
    _bare(ax)
    graph = nx.MultiDiGraph() if directed else nx.MultiGraph()
    graph.add_edges_from(edges)
    if pos is None:
        pos = nx.kamada_kawai_layout(graph) if len(graph) > 2 else nx.circular_layout(graph)
    highlight_edges = list(highlight_edges)
    highlight_nodes = set(highlight_nodes)

    node_colors = [EMPHASIS if n in highlight_nodes else PAPER for n in graph.nodes]
    edge_colors = [EMPHASIS if n in highlight_nodes else INK for n in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=node_size,
                           node_color=node_colors, edgecolors=edge_colors,
                           linewidths=1.0)
    for n in graph.nodes:
        ax.text(*pos[n], n, ha="center", va="center", fontsize=9, family=MONO,
                color=PAPER if n in highlight_nodes else INK, zorder=5)

    # parallel edges get increasing curvature; opposite-direction pairs curve too
    seen: dict[tuple[str, str], int] = {}
    for u, v, key in graph.edges(keys=True):
        count = seen.get((u, v), 0)
        seen[(u, v)] = count + 1
        reverse = directed and graph.has_edge(v, u)
        rad = 0.18 * (count + (1 if reverse else 0))
        if u == v:
            rad = 0.0
        on = (u, v) in highlight_edges or (u, v, key) in highlight_edges
        if not directed:
            on = on or (v, u) in highlight_edges or (v, u, key) in highlight_edges
        color = EMPHASIS if on else RULE
        width = 1.8 if on else 1.2
        nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=[(u, v, key)],
                               edge_color=color, width=width, arrows=True,
                               arrowstyle="-|>" if directed else "-",
                               arrowsize=12, node_size=node_size,
                               connectionstyle=f"arc3,rad={rad}")
    if edge_labels:
        nx.draw_networkx_edge_labels(graph, pos, ax=ax, edge_labels=edge_labels,
                                     font_size=8, font_color=INK,
                                     bbox=dict(facecolor=PAPER, edgecolor="none", pad=0.5))
    ax.margins(0.15)
    return pos


def state_diagram(
    ax: plt.Axes,
    states: Sequence[str],
    transitions: dict[tuple[str, str], float],
    emissions: dict[str, dict[str, float]] | None = None,
    pos: dict[str, tuple[float, float]] | None = None,
) -> None:
    """Draw an HMM state diagram with labeled transition arcs and self-loops.

    `transitions[(a, b)]` is the probability of moving from state a to b.
    `emissions[a]` is a dict symbol -> probability, printed under state a.
    """
    _bare(ax)
    n = len(states)
    if pos is None:
        pos = {s: (2.2 * i, 0.0) for i, s in enumerate(states)}
    nodes = {}
    for s in states:
        node = mpatches.Circle(pos[s], 0.5, facecolor=PAPER,
                                edgecolor=INK, linewidth=1.0, zorder=3)
        nodes[s] = node
        ax.add_patch(node)
        ax.text(*pos[s], s, ha="center", va="center", fontsize=10, zorder=4)
        if emissions and s in emissions:
            lines = [f"{sym}: {p:g}" for sym, p in emissions[s].items()]
            ax.text(pos[s][0], pos[s][1] - 0.75, "\n".join(lines),
                    ha="center", va="top", fontsize=8, family=MONO, color=INK)
    for (a, b), p in transitions.items():
        (xa, ya), (xb, yb) = pos[a], pos[b]
        if a == b:
            # self-loop drawn as an arc above the state
            loop = mpatches.FancyArrowPatch((xa - 0.25, ya + 0.43), (xa + 0.25, ya + 0.43),
                                            connectionstyle="arc3,rad=-1.6",
                                            arrowstyle="-|>", mutation_scale=10,
                                            color=RULE, lw=1.2)
            ax.add_patch(loop)
            ax.text(xa, ya + 1.25, f"{p:g}", ha="center", va="bottom", fontsize=8)
        else:
            reverse = (b, a) in transitions
            rad = 0.25 if reverse else 0.0
            arrow = mpatches.FancyArrowPatch((xa, ya), (xb, yb),
                                             connectionstyle=f"arc3,rad={rad}",
                                             arrowstyle="-|>", mutation_scale=10,
                                             color=RULE, lw=1.2,
                                             patchA=nodes[a], patchB=nodes[b],
                                             shrinkA=2, shrinkB=2, zorder=2)
            ax.add_patch(arrow)
            mx, my = (xa + xb) / 2, (ya + yb) / 2
            # offset the label to the side the arc bends toward
            side = -1 if xb > xa else 1
            ax.text(mx, my + 0.38 * side * (1 if reverse else 0) + (0.15 if not reverse else 0),
                    f"{p:g}", ha="center", va="bottom" if side > 0 else "top", fontsize=8)
    ax.set_xlim(-1.2, 2.2 * (n - 1) + 1.2)
    ax.set_ylim(-1.9, 1.9)


# --- Trees -------------------------------------------------------------------

def tree(
    ax: plt.Axes,
    root: dict,
    highlight: Iterable[str] = (),
    show_lengths: bool = False,
    highlight_edges: Iterable[tuple[str, str]] = (),
) -> None:
    """Draw a rooted tree with branch lengths as horizontal extent.

    `root` is a nested dict: {"name": str, "length": float, "children": [...]}.
    Leaves are drawn at the right, in the order they appear. Names of nodes
    in `highlight` are colored with the emphasis color. Directed parent-child
    name pairs in `highlight_edges` emphasize the corresponding branches.
    """
    _bare(ax)
    highlight = set(highlight)
    highlight_edges = set(highlight_edges)
    leaf_y: dict[int, float] = {}
    counter = [0]

    def place(node: dict, x0: float, parent: str = "") -> tuple[float, float]:
        x1 = x0 + node.get("length", 0.0)
        children = node.get("children", [])
        name = node.get("name", "")
        branch_color = EMPHASIS if (parent, name) in highlight_edges else RULE
        if not children:
            y = -counter[0]
            counter[0] += 1
            ax.plot([x0, x1], [y, y], color=branch_color, lw=1.2,
                    solid_capstyle="butt")
            name = node.get("name", "")
            ax.text(x1 + 0.03, y, name, ha="left", va="center", fontsize=10,
                    color=EMPHASIS if name in highlight else INK)
            return x1, y
        ys = []
        for child in children:
            _, cy = place(child, x1, name)
            ys.append(cy)
        y = float(np.mean(ys))
        ax.plot([x1, x1], [min(ys), max(ys)], color=RULE, lw=1.2)
        ax.plot([x0, x1], [y, y], color=branch_color, lw=1.2)
        name = node.get("name", "")
        if name:
            ax.text(x1 - 0.03, y + 0.15, name, ha="right", va="bottom", fontsize=8,
                    color=EMPHASIS if name in highlight else RULE)
        if show_lengths and node.get("length"):
            ax.text((x0 + x1) / 2, y + 0.12, f"{node['length']:g}", ha="center",
                    va="bottom", fontsize=8, color=RULE)
        return x1, y

    place(root, 0.0)
    ax.set_aspect("auto")
    ax.margins(x=0.25, y=0.2)


# --- Alignments --------------------------------------------------------------

def alignment(
    ax: plt.Axes,
    top: str,
    bottom: str,
    labels: tuple[str, str] = ("x", "y"),
    mark_mismatches: bool = True,
) -> None:
    """Print two aligned strings in monospace with a match line between them.

    Matches get a vertical bar, gaps a blank, mismatches a dot (in the
    emphasis color when `mark_mismatches`).
    """
    _bare(ax)
    assert len(top) == len(bottom), "aligned strings must have equal length"
    for i, (a, b) in enumerate(zip(top, bottom)):
        ax.text(i, 1, a, ha="center", va="center", fontsize=11, family=MONO)
        ax.text(i, -1, b, ha="center", va="center", fontsize=11, family=MONO)
        if a == b:
            mark, color = "|", RULE
        elif "-" in (a, b):
            mark, color = "", RULE
        else:
            mark, color = "·", EMPHASIS if mark_mismatches else RULE
        ax.text(i, 0, mark, ha="center", va="center", fontsize=11,
                family=MONO, color=color)
    ax.text(-1.2, 1, f"${labels[0]}$", ha="center", va="center", fontsize=10)
    ax.text(-1.2, -1, f"${labels[1]}$", ha="center", va="center", fontsize=10)
    ax.set_xlim(-2, len(top))
    ax.set_ylim(-2, 2)


def sequence_windows(
    ax: plt.Axes,
    sequence: str,
    *,
    word: str,
    label: str = "x",
) -> None:
    """Print a sequence and underline every occurrence of a word, including overlaps."""
    if not word:
        raise ValueError("The highlighted word must not be empty.")
    _bare(ax)
    for i, base in enumerate(sequence):
        ax.text(i, 0, base, ha="center", va="center", fontsize=11,
                family=MONO, color=INK)
    ax.text(-1, 0, label, ha="right", va="center", fontsize=10, color=INK)

    # Overlapping words need separate lanes so each occurrence stays visible.
    lane_ends = []
    for start in range(len(sequence) - len(word) + 1):
        if sequence[start:start + len(word)] != word:
            continue
        lane = 0
        while lane < len(lane_ends) and lane_ends[lane] >= start:
            lane += 1
        end = start + len(word) - 1
        if lane == len(lane_ends):
            lane_ends.append(end)
        else:
            lane_ends[lane] = end
        height = -0.4 - 0.3 * lane
        ax.plot([start - 0.25, start - 0.25, end + 0.25, end + 0.25],
                [height + 0.1, height, height, height + 0.1],
                color=EMPHASIS, linewidth=1.5)
    ax.set_xlim(-2, max(1, len(sequence)))
    ax.set_ylim(-0.8 - 0.3 * max(0, len(lane_ends) - 1), 0.6)


def replicate_sampling(
    ax: plt.Axes,
    biological_values: Sequence[float],
    technical_values: Sequence[float],
    culture_labels: Sequence[str],
) -> None:
    """Contrast assays of three independent cultures with three assays of one."""
    if len(biological_values) != 3 or len(technical_values) != 3:
        raise ValueError("This comparison shows three measurements per panel.")
    if len(culture_labels) != 3:
        raise ValueError("Supply a name for each independent culture.")
    _bare(ax)
    ax.set_xlim(0, 10.3)
    ax.set_ylim(-0.8, 4.8)
    ax.text(2.3, 4.5, "Separate cultures", ha="center", fontsize=10)
    ax.text(7.6, 4.5, "Repeated assays", ha="center", fontsize=10)

    def dish(x: float, y: float, label: str, emphasize: bool = False) -> None:
        """Draw the physical source of a measurement as a small culture dish."""
        ax.add_patch(mpatches.Circle((x, y), 0.38, facecolor=PAPER,
                                    edgecolor=EMPHASIS if emphasize else RULE))
        for dx, dy in ((-0.15, 0.08), (0.12, 0.12), (0.02, -0.13)):
            ax.add_patch(mpatches.Ellipse((x + dx, y + dy), 0.16, 0.07,
                                         angle=25, facecolor=RULE))
        ax.text(x, y - 0.6, label, ha="center", fontsize=9)

    def measurement(x: float, y: float, value: float) -> None:
        """Draw one assay result, without treating its box as a biological unit."""
        ax.add_patch(mpatches.Rectangle((x - 0.7, y - 0.28), 1.4, 0.56,
                                        facecolor=LIGHT, edgecolor=RULE))
        ax.text(x, y, f"{value:g}", ha="center", va="center", fontsize=10)

    ax.text(3.5, 4.0, "Measured value", ha="center", fontsize=9)
    ax.text(9.0, 4.0, "Measured value", ha="center", fontsize=9)
    dish(6.3, 2.0, culture_labels[0], emphasize=True)
    for i, y in enumerate((3.5, 2.0, 0.5)):
        dish(1.1, y, culture_labels[i])
        measurement(3.5, y, biological_values[i])
        measurement(9.0, y, technical_values[i])
        ax.annotate("", xy=(2.7, y), xytext=(1.55, y),
                    arrowprops={"arrowstyle": "->", "color": RULE, "lw": 1})
        ax.annotate("", xy=(8.2, y), xytext=(6.75, 2.0),
                    arrowprops={"arrowstyle": "->", "color": RULE, "lw": 1})
    ax.text(2.3, -0.6, "3 independent cultures", ha="center", fontsize=10)
    ax.text(7.6, -0.6, "1 independent culture", ha="center", fontsize=10)


def cell_barcode_records(
    ax: plt.Axes,
    counts: Sequence[Sequence[int]],
    *,
    cells: Sequence[str],
    genes: Sequence[str],
    barcodes: Sequence[str],
) -> None:
    """Connect two cells' captured RNA, mixed tagged records, and count rows.

    This small conceptual diagram omits capture chemistry. Each record is one
    distinct captured molecule after duplicate handling, not one raw read.
    """
    counts = np.asarray(counts)
    if counts.shape != (2, 2) or len(cells) != 2 or len(genes) != 2:
        raise ValueError("Use two cells and two genes for this small schematic.")
    if len(barcodes) != 2 or len(set(barcodes)) != 2:
        raise ValueError("Each cell needs a distinct barcode.")
    if np.any(counts < 0) or np.any(counts != counts.astype(int)):
        raise ValueError("Molecule counts must be nonnegative integers.")
    if np.any(counts.sum(axis=1) > 6) or counts.sum() == 0:
        raise ValueError("Draw at most six captured molecules per cell.")
    _bare(ax)
    ax.set_xlim(-0.2, 12.1)
    ax.set_ylim(-1.1, 5.1)
    ax.text(1.1, 4.7, "RNA in separate droplets", ha="center", fontsize=10)
    ax.text(5.9, 4.7, "Mixed molecule records", ha="center", fontsize=10)
    ax.text(10.2, 4.7, "Counts by cell", ha="center", fontsize=10)

    records = []
    for i, center_y in enumerate((3.2, 0.7)):
        ax.add_patch(mpatches.Circle((1.1, center_y), 0.8,
                                    edgecolor=RULE, facecolor=PAPER))
        molecule_genes = []
        for j in range(2):
            molecule_genes.extend([j] * int(counts[i, j]))
        cell_records = []
        for number, gene in enumerate(molecule_genes):
            x = 0.58 + 0.55 * (number % 2)
            y = center_y + 0.32 - 0.31 * (number // 2)
            strand_x = np.linspace(x, x + 0.43, 25)
            strand_y = y + 0.045 * np.sin(np.linspace(0, 4 * np.pi, 25))
            ax.plot(strand_x, strand_y, color=(SECOND, THIRD)[gene], lw=1.5,
                    linestyle=("-", "--")[gene])
            cell_records.append((barcodes[i], genes[gene], f"m{number + 1}"))
        records.append(cell_records)
        ax.text(1.1, center_y - 1.0, f"{cells[i]}: tag {barcodes[i]}",
                ha="center", fontsize=9)
        ax.annotate("", xy=(3.7, 2.1), xytext=(2.0, center_y),
                    arrowprops={"arrowstyle": "->", "color": RULE, "lw": 1})

    mixed = []
    for number in range(max(map(len, records))):
        for row in records:
            if number < len(row):
                mixed.append(row[number])
    ax.text(4.2, 4.2, "cell tag", ha="center", fontsize=8)
    ax.text(5.8, 4.2, "gene", ha="center", fontsize=8)
    ax.text(7.2, 4.2, "molecule", ha="center", fontsize=8)
    for y, (barcode, gene, molecule) in zip(np.linspace(3.8, 0.0, len(mixed)), mixed):
        ax.add_patch(mpatches.Rectangle((3.8, y - 0.14), 0.8, 0.28,
                                        facecolor=PAPER, edgecolor=EMPHASIS))
        ax.text(4.2, y, barcode, ha="center", va="center", family=MONO,
                fontsize=9)
        ax.text(5.8, y, gene, ha="center", va="center", fontsize=9)
        ax.text(7.2, y, molecule, ha="center", va="center", family=MONO,
                fontsize=9)
    ax.annotate("", xy=(8.5, 2.1), xytext=(7.8, 2.1),
                arrowprops={"arrowstyle": "->", "color": RULE, "lw": 1})

    for j, gene in enumerate(genes):
        ax.text(10.0 + j * 1.3, 3.6, gene, ha="center", fontsize=9)
    for i, y in enumerate((2.9, 1.2)):
        ax.text(8.6, y + 0.45, f"{cells[i]} ({barcodes[i]})", fontsize=9)
        for j in range(2):
            x = 10.0 + j * 1.3
            ax.add_patch(mpatches.Rectangle((x - 0.5, y - 0.3), 1, 0.6,
                                            facecolor=LIGHT, edgecolor=RULE))
            ax.text(x, y, str(int(counts[i, j])), ha="center", va="center",
                    fontsize=10)
    for j, gene in enumerate(genes):
        ax.plot([0.0, 0.45], [-0.85 + j * 0.27] * 2,
                color=(SECOND, THIRD)[j], linestyle=("-", "--")[j], lw=1.5)
        ax.text(0.55, -0.85 + j * 0.27, gene, va="center", fontsize=8)


# --- Phase planes ------------------------------------------------------------

def phase_plane(
    ax: plt.Axes,
    f: Callable[[np.ndarray, np.ndarray], np.ndarray],
    g: Callable[[np.ndarray, np.ndarray], np.ndarray],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    fixed_points: Iterable[tuple[float, float, bool]] = (),
    trajectories: Iterable[np.ndarray] = (),
    labels: tuple[str, str] = ("$x$", "$y$"),
    density: float = 0.8,
) -> None:
    """Draw a 2-D vector field with nullclines, fixed points, and trajectories.

    `f` and `g` are the right-hand sides dx/dt and dy/dt, vectorized. Fixed
    points are (x, y, stable): stable ones are filled, unstable ones open.
    Trajectories are arrays of shape (n, 2).
    """
    xs = np.linspace(*xlim, 200)
    ys = np.linspace(*ylim, 200)
    X, Y = np.meshgrid(xs, ys)
    U, V = f(X, Y), g(X, Y)
    ax.streamplot(X, Y, U, V, color=LIGHT, density=density, linewidth=0.8,
                  arrowsize=0.8)
    ax.contour(X, Y, U, levels=[0], colors=[SECOND], linewidths=1.5)
    ax.contour(X, Y, V, levels=[0], colors=[THIRD], linewidths=1.5)
    for traj in trajectories:
        traj = np.asarray(traj)
        ax.plot(traj[:, 0], traj[:, 1], color=INK, lw=1.2)
        ax.plot(traj[0, 0], traj[0, 1], marker="o", color=INK, ms=3)
    for x, y, stable in fixed_points:
        ax.plot(x, y, marker="o", ms=7, color=EMPHASIS,
                markerfacecolor=EMPHASIS if stable else PAPER, markeredgewidth=1.5)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])


def read_layout(
    ax: plt.Axes,
    reads: Sequence[str],
    *,
    starts: Sequence[int] | None = None,
    reference: str | None = None,
    labels: Sequence[str] | None = None,
) -> None:
    """Draw read records alone or at proposed starts beneath a candidate sequence."""
    if not reads:
        raise ValueError("A read layout needs at least one read.")
    if labels is None:
        labels = [f"r{i + 1}" for i in range(len(reads))]
    if len(labels) != len(reads):
        raise ValueError("Supply one label per read.")
    if starts is None:
        if reference is not None:
            raise ValueError("A candidate sequence requires proposed read starts.")
        starts = [0] * len(reads)
    if len(starts) != len(reads) or any(start < 0 for start in starts):
        raise ValueError("Supply one nonnegative start per read.")
    if reference is not None:
        for read, start in zip(reads, starts):
            if reference[start:start + len(read)] != read:
                raise ValueError("Each proposed read must match the candidate sequence.")

    _bare(ax)
    rows = []
    if reference is not None:
        rows.append(("candidate", reference, 0))
    rows.extend(zip(labels, reads, starts))
    for row, (label, sequence, start) in enumerate(rows):
        ax.text(-0.8, -row, label, ha="right", va="center",
                fontsize=9, color=INK)
        for offset, base in enumerate(sequence):
            ax.text(start + offset, -row, base, ha="center", va="center",
                    fontsize=11, family=MONO, color=INK)
    extent = max(start + len(sequence) for _, sequence, start in rows)
    ax.set_xlim(-3.4, extent - 0.3)
    ax.set_ylim(-len(rows) + 0.4, 0.7)


def independent_pair_area(
    ax: plt.Axes,
    probability: float,
    *,
    base: str = "A",
) -> None:
    """Show a matching pair as the intersection of two independent base draws."""
    if not 0 < probability < 1:
        raise ValueError("The base probability must lie strictly between zero and one.")
    if len(base) != 1:
        raise ValueError("Use a one-letter base label.")
    _bare(ax)
    ax.add_patch(mpatches.Rectangle((0, 0), 1, 1, facecolor=PAPER,
                                   edgecolor=RULE, linewidth=0.8))
    ax.add_patch(mpatches.Rectangle((0, 0), probability, probability,
                                   facecolor=EMPHASIS, edgecolor=RULE,
                                   linewidth=0.8))
    ax.plot([probability, probability], [0, 1], color=RULE, linewidth=0.8)
    ax.plot([0, 1], [probability, probability], color=RULE, linewidth=0.8)
    for midpoint, label in [(probability / 2, base),
                            ((1 + probability) / 2, "other")]:
        ax.text(midpoint, -0.06, label, ha="center", va="top", fontsize=9)
        ax.text(-0.06, midpoint, label, ha="right", va="center", fontsize=9)
    ax.text(0.5, -0.2, "First base", ha="center", va="top", fontsize=10)
    ax.text(-0.35, 0.5, "Second base", ha="center", va="center",
            rotation=90, fontsize=10)
    ax.text(probability / 2, probability / 2, base + "/" + base,
            ha="center", va="center", color=PAPER, fontsize=9)
    ax.set_xlim(-0.45, 1.05)
    ax.set_ylim(-0.3, 1.08)


def read_pileup(
    ax: plt.Axes,
    reference: str,
    reads: Sequence[str],
    starts: Sequence[int],
    *,
    highlight_column: int | None = None,
    reference_label: str = "reference",
) -> None:
    """Stack placed reads under a reference, greying agreeing bases.

    A read base that disagrees with the reference is drawn in ink; at
    ``highlight_column`` (0-based reference position) disagreeing bases get a
    red fill, marking the difference the figure is about.
    """
    if len(reads) != len(starts):
        raise ValueError("Supply one start per read.")
    for read, start in zip(reads, starts):
        if start < 0 or start + len(read) > len(reference):
            raise ValueError("Each read must lie within the reference.")
    _bare(ax)
    ax.text(-0.8, 0, reference_label, ha="right", va="center", fontsize=9)
    for i, base in enumerate(reference):
        ax.text(i, 0, base, ha="center", va="center", fontsize=11,
                family=MONO, color=INK, fontweight="bold")
    ax.plot([-0.5, len(reference) - 0.5], [-0.55, -0.55], color=RULE,
            linewidth=0.8)
    for row, (read, start) in enumerate(zip(reads, starts), start=1):
        y = -row - 0.3
        for offset, base in enumerate(read):
            i = start + offset
            differs = base != reference[i]
            selected = differs and i == highlight_column
            if selected:
                ax.add_patch(mpatches.Rectangle(
                    (i - 0.36, y - 0.32), 0.72, 0.64,
                    facecolor=EMPHASIS, edgecolor="none", zorder=1,
                ))
            color = PAPER if selected else (INK if differs else RULE)
            ax.text(i, y, base, ha="center", va="center", fontsize=11,
                    family=MONO, color=color, zorder=2,
                    fontweight="bold" if differs else "normal")
    ax.set_xlim(-3.4, len(reference) - 0.4)
    ax.set_ylim(-len(reads) - 0.8, 0.6)


# --- One-variable rate balance and phase line -------------------------------

def rate_balance(
    ax: plt.Axes,
    production: Callable[[np.ndarray], np.ndarray],
    removal: Callable[[np.ndarray], np.ndarray],
    xlim: tuple[float, float],
    *,
    fixed_points: Iterable[tuple[float, bool]] = (),
    labels: tuple[str, str] = ("Production", "Removal"),
    xlabel: str = "$x$",
    ylabel: str = "Rate",
) -> None:
    """Plot production and removal rates against one state, above a phase line.

    Production is drawn in SECOND and removal in THIRD. Fixed points are
    (x, stable) at crossings of the two curves: red filled if stable, open if
    unstable. Below the curves, a grey phase line carries one ink arrow per
    interval between fixed points, pointing the way the state moves (to the
    right wherever production exceeds removal).
    """
    xs = np.linspace(*xlim, 400)
    # A constant rate function may return a scalar, so broadcast to the grid.
    produced = np.broadcast_to(production(xs), xs.shape)
    removed = np.broadcast_to(removal(xs), xs.shape)
    top = max(produced.max(), removed.max())
    line_y = -0.12 * top
    ax.plot(xs, produced, color=SECOND, label=labels[0])
    ax.plot(xs, removed, color=THIRD, label=labels[1])
    ax.plot(xlim, [line_y, line_y], color=RULE, lw=0.8)
    points = sorted(fixed_points)
    bounds = [xlim[0]] + [x for x, _ in points] + [xlim[1]]
    for left, right in zip(bounds[:-1], bounds[1:]):
        middle = (left + right) / 2
        net = float(np.asarray(production(middle)) - np.asarray(removal(middle)))
        if net == 0:
            continue
        half = 0.2 * (right - left)
        start, end = (middle - half, middle + half) if net > 0 else (
            middle + half, middle - half)
        ax.annotate("", xy=(end, line_y), xytext=(start, line_y),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2,
                                    mutation_scale=12))
    for x, stable in points:
        level = float(np.asarray(production(x)))
        face = EMPHASIS if stable else PAPER
        ax.plot([x, x], [line_y, level], color=RULE, lw=0.8, ls=":")
        for y in (level, line_y):
            ax.plot(x, y, marker="o", ms=7, color=EMPHASIS,
                    markerfacecolor=face, markeredgewidth=1.5, zorder=3)
    ax.set_xlim(*xlim)
    ax.set_ylim(1.9 * line_y, 1.1 * top)
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= 1.1 * top])
    ax.spines["left"].set_bounds(0, 1.1 * top)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper left")


# --- Clustering and neighbor graphs on measured coordinates -------------------

CLUSTER_MARKERS = ["o", "^", "s", "D", "v", "P"]


def cluster_step(
    ax: plt.Axes,
    points: Sequence[Sequence[float]],
    labels: Sequence[int],
    centers: Sequence[Sequence[float]],
    *,
    names: Sequence[str] = (),
    previous_centers: Sequence[Sequence[float]] | None = None,
) -> None:
    """Draw one K-means state: points by cluster, centers as crosses.

    Cluster g uses SERIES[g] and the g-th marker shape. When previous_centers
    is given, red arrows show how each center moved in the update step.
    Axes stay visible because positions are measured coordinates.
    """
    points = np.asarray(points, dtype=float)
    centers = np.asarray(centers, dtype=float)
    labels = np.asarray(labels)
    for g in range(len(centers)):
        members = points[labels == g]
        ax.scatter(members[:, 0], members[:, 1], s=36, color=SERIES[g],
                   marker=CLUSTER_MARKERS[g], zorder=3)
    for name, point in zip(names, points):
        ax.annotate(name, point, xytext=(5, 4), textcoords="offset points",
                    fontsize=9, family=MONO)
    if previous_centers is not None:
        for old, new in zip(np.asarray(previous_centers, dtype=float), centers):
            if np.allclose(old, new):
                continue
            ax.annotate("", xy=new, xytext=old, zorder=4,
                        arrowprops={"arrowstyle": "-|>", "color": EMPHASIS,
                                    "lw": 1.5, "shrinkA": 0, "shrinkB": 5})
    for g, center in enumerate(centers):
        ax.scatter(*center, s=90, marker="X", color=SERIES[g],
                   edgecolors=INK, linewidths=1.0, zorder=5)
    ax.margins(0.08)
    ax.set_aspect("equal", adjustable="box")


def neighbor_graph(
    ax: plt.Axes,
    points: Sequence[Sequence[float]],
    edges: Iterable[tuple[int, int]],
    *,
    names: Sequence[str] = (),
    highlight: Iterable[int] = (),
) -> None:
    """Draw an undirected neighbor graph with nodes at measured coordinates.

    Edges are grey; highlighted nodes, and edges joining two of them, are red.
    With names, nodes are labeled white circles; without, small dots.
    """
    points = np.asarray(points, dtype=float)
    highlight = set(highlight)
    for i, j in edges:
        on = i in highlight and j in highlight
        ax.plot(points[[i, j], 0], points[[i, j], 1],
                color=EMPHASIS if on else RULE, lw=1.6 if on else 0.8,
                zorder=1)
    colors = [EMPHASIS if i in highlight else PAPER for i in range(len(points))]
    edgecolors = [EMPHASIS if i in highlight else INK for i in range(len(points))]
    if names:
        ax.scatter(points[:, 0], points[:, 1], s=260, c=colors,
                   edgecolors=edgecolors, linewidths=1.0, zorder=3)
        for i, name in enumerate(names):
            ax.text(*points[i], name, ha="center", va="center", fontsize=9,
                    family=MONO, color=PAPER if i in highlight else INK,
                    zorder=4)
    else:
        ax.scatter(points[:, 0], points[:, 1], s=14, c=colors,
                   edgecolors=edgecolors, linewidths=0.6, zorder=3)
    ax.margins(0.12 if names else 0.05)
    ax.set_aspect("equal", adjustable="box")


def multiple_alignment(ax: plt.Axes, rows: dict[str, str]) -> None:
    """Print a multiple alignment in monospace columns, shading columns where rows differ.

    Row labels are the sample names from the prose; column numbers run from 1.
    """
    _bare(ax)
    lengths = {len(sequence) for sequence in rows.values()}
    if len(lengths) != 1 or not rows:
        raise ValueError("A multiple alignment needs nonempty rows of equal length.")
    length = lengths.pop()
    sequences = list(rows.values())
    for i in range(length):
        if len({sequence[i] for sequence in sequences}) > 1:
            ax.add_patch(mpatches.Rectangle(
                (i - 0.45, -len(rows) + 0.5), 0.9, len(rows),
                facecolor=LIGHT, edgecolor="none", zorder=0,
            ))
        ax.text(i, 0.85, str(i + 1), ha="center", va="center", fontsize=9,
                color=RULE)
    for row, (label, sequence) in enumerate(rows.items()):
        ax.text(-0.9, -row, label, ha="right", va="center", fontsize=10)
        for i, symbol in enumerate(sequence):
            ax.text(i, -row, symbol, ha="center", va="center", fontsize=11,
                    family=MONO, color=INK, zorder=2)
    ax.set_xlim(-3, length - 0.4)
    ax.set_ylim(-len(rows) + 0.3, 1.3)


# --- Several alignment paths converging on one table cell -------------------

def path_bundle(
    ax: plt.Axes,
    rows: Sequence[str],
    cols: Sequence[str],
    paths: Sequence[Sequence[tuple[int, int]]],
    *,
    keep: int = 0,
    continuation: Sequence[tuple[int, int]] = (),
    labels: Sequence[str | None] | None = None,
) -> None:
    """Draw several alignment paths through the cells of an alignment grid.

    Cells are light dots at the `dp_table` positions, with the same (i, j)
    indices and row and column labels. Path `keep` (the best one, which the
    table retains) is red; the others are grey. A shared `continuation` from
    the common end cell is drawn dashed in ink. Optional `labels`, one per
    path, become a monospace legend in path order; a None label leaves that
    path out of the legend, so several alike paths can share one entry.
    """
    if not 0 <= keep < len(paths):
        raise ValueError("keep must index one of the paths.")
    if labels is not None and len(labels) != len(paths):
        raise ValueError("Supply one label (or None) per path.")
    _bare(ax)
    for i in range(len(rows) + 1):
        for j in range(len(cols) + 1):
            ax.plot(j, -i, marker="o", ms=7, color=LIGHT,
                    markeredgecolor=RULE, markeredgewidth=0.6, zorder=1)
    ax.text(-1, 0, "$-$", ha="center", va="center", fontsize=10)
    for i, label in enumerate(rows, start=1):
        ax.text(-1, -i, label, ha="center", va="center", fontsize=10, family=MONO)
    ax.text(0, 1, "$-$", ha="center", va="center", fontsize=10)
    for j, label in enumerate(cols, start=1):
        ax.text(j, 1, label, ha="center", va="center", fontsize=10, family=MONO)
    # draw the kept path last so it sits on top of shared grey segments
    order = [k for k in range(len(paths)) if k != keep] + [keep]
    lines = {}
    for k in order:
        path = paths[k]
        best = k == keep
        (lines[k],) = ax.plot([j for _, j in path], [-i for i, _ in path],
                              color=EMPHASIS if best else RULE,
                              lw=2.5 if best else 1.5, marker="o", ms=4,
                              zorder=4 if best else 3)
    if continuation:
        ax.plot([j for _, j in continuation], [-i for i, _ in continuation],
                color=INK, lw=1.5, ls="--", marker="o", ms=4, zorder=2)
    if labels is not None:
        # list the paths in the order given, not the drawing order
        shown = [k for k in range(len(paths)) if labels[k] is not None]
        ax.legend([lines[k] for k in shown], [labels[k] for k in shown],
                  loc="center left", bbox_to_anchor=(1.0, 0.5),
                  prop={"family": MONO, "size": 9}, handlelength=1.5)
    ax.set_xlim(-1.6, len(cols) + 0.4)
    ax.set_ylim(-len(rows) - 0.4, 1.6)


# --- Experimental setting: an induced fluorescent reporter -----------------

def reporter_setting(ax: plt.Axes) -> None:
    """Sketch an induced reporter gene in one cell beside time-lapse images.

    Left: an inducer reaches the promoter, the reporter gene (gfp) is copied
    to mRNA and translated into fluorescent protein (red dots, concentration
    x) at rate beta, and protein is removed at rate alpha x. Production and
    removal arrows use SECOND and THIRD, the colors of `rate_balance`.
    Right: one lineage imaged 0, 30, and 60 minutes after induction with a
    30-minute division; the followed cell has an ink outline and its sister
    a grey one. Dot density stands for concentration, so it rises from frame
    to frame although a division halves the molecules in each cell.
    """
    _bare(ax)

    def arrow(start: tuple[float, float], end: tuple[float, float],
              color: str, lw: float = 1.5) -> None:
        """Draw one straight arrow with a filled head."""
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                    mutation_scale=11, shrinkA=0, shrinkB=0))

    def cell(x0: float, y0: float, width: float, height: float, n_dots: int,
             edge: str = RULE, lw: float = 1.2) -> None:
        """Draw a rod-shaped bacterium holding n_dots reporter molecules."""
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0, y0), width, height,
            boxstyle=f"round,pad=0,rounding_size={height / 2}",
            facecolor=PAPER, edgecolor=edge, linewidth=lw))
        pad = 0.3
        step = (width - 2 * pad) / max(n_dots, 1)
        for i in range(n_dots):
            ax.add_patch(mpatches.Circle(
                (x0 + pad + (i + 0.5) * step, y0 + height / 2 + 0.14 * (-1) ** i),
                0.09, facecolor=EMPHASIS, edgecolor="none"))

    def inducer(x: float, y: float) -> None:
        """Draw one inducer molecule as a small ink diamond."""
        ax.add_patch(mpatches.RegularPolygon((x, y), 4, radius=0.13,
                                             facecolor=INK, edgecolor="none"))

    # Left: what happens inside one cell after induction.
    ax.text(4.5, 5.75, "Inside one cell", ha="center", fontsize=10)
    cell(0.4, 0.9, 8.2, 3.9, 0)
    for x, y in ((0.35, 5.2), (0.9, 5.4), (1.45, 5.15)):
        inducer(x, y)
    ax.text(1.8, 5.25, "inducer", ha="left", va="center", fontsize=9)
    arrow((1.2, 4.95), (2.55, 2.8), RULE, lw=0.8)
    ax.plot([1.7, 6.2], [2.1, 2.1], color=INK, lw=1.5, zorder=1,
            solid_capstyle="butt")
    ax.add_patch(mpatches.Rectangle((2.2, 1.85), 1.0, 0.5, facecolor=LIGHT,
                                    edgecolor=RULE, linewidth=0.8, zorder=2))
    ax.text(2.7, 1.65, "promoter", ha="center", va="top", fontsize=9)
    ax.add_patch(mpatches.Rectangle((3.3, 1.85), 2.2, 0.5, facecolor=PAPER,
                                    edgecolor=INK, linewidth=1.0, zorder=2))
    ax.text(4.4, 2.1, "gfp", ha="center", va="center", fontsize=9,
            style="italic", zorder=3)
    inducer(2.7, 2.55)
    wave = np.linspace(0, 6 * np.pi, 60)
    ax.plot(np.linspace(3.6, 5.2, 60), 3.35 + 0.08 * np.sin(wave),
            color=INK, lw=1.2)
    ax.text(4.4, 3.6, "mRNA", ha="center", va="bottom", fontsize=9)
    arrow((4.4, 2.45), (4.4, 3.15), SECOND)
    for dx, dy in ((-0.35, 0.3), (0.05, 0.45), (0.4, 0.2), (-0.25, -0.15),
                   (0.15, -0.05), (0.5, -0.3), (-0.05, -0.45)):
        ax.add_patch(mpatches.Circle((6.7 + dx, 3.35 + dy), 0.11,
                                     facecolor=EMPHASIS, edgecolor="none"))
    ax.text(6.7, 3.98, "GFP, $x$", ha="center", va="bottom", fontsize=9)
    arrow((5.35, 3.35), (6.1, 3.35), SECOND)
    ax.text(5.72, 3.1, r"$\beta$", ha="center", va="top", fontsize=10)
    arrow((6.8, 2.75), (6.8, 1.95), THIRD)
    ax.text(7.0, 2.35, r"$\alpha x$", ha="left", va="center", fontsize=10)
    ax.text(6.8, 1.8, "removed", ha="center", va="top", fontsize=9)
    ax.text(4.5, 0.15, r"$\dot x = \beta - \alpha x$", ha="center",
            va="center", fontsize=10)

    # Right: the same lineage photographed as reporter accumulates.
    ax.text(14.1, 5.75, "Images of one lineage", ha="center", fontsize=10)
    unit, height, gap, axis_y = 1.4, 0.85, 0.12, 1.6
    frames = [(11.0, "0", [(unit, 0, INK)]),
              (14.0, "30", [(unit, 3, INK), (unit, 3, RULE)]),
              (17.2, "60", [(2 * unit, 9, INK)])]
    for center, label, cells in frames:
        x = center - (sum(w for w, _, _ in cells) + gap * (len(cells) - 1)) / 2
        for width, n_dots, edge in cells:
            cell(x, 3.1, width, height, n_dots, edge=edge,
                 lw=1.2 if edge == INK else 0.9)
            x += width + gap
        ax.plot([center, center], [axis_y, axis_y - 0.15], color=RULE, lw=0.8)
        ax.text(center, axis_y - 0.3, label, ha="center", va="top", fontsize=9)
    arrow((10.0, axis_y), (19.2, axis_y), RULE, lw=0.8)
    ax.text(14.1, 0.15, "minutes after induction", ha="center", va="center",
            fontsize=9)
    ax.set_xlim(0.0, 19.3)
    ax.set_ylim(-0.2, 6.1)


# --- Opening cartoons: where the data come from ------------------------------

def sequencing_run(
    ax: plt.Axes,
    reference: str,
    reads: Sequence[str],
    *,
    query: int = 0,
    offset: int = 3,
) -> None:
    """Cartoon of a short-read experiment, from sample DNA to unplaced reads.

    Left to right: copies of one genome region, the same copies broken into
    fragments, a sequencer holding a flow cell, and up to four read strings
    with no positions. Below, read `query` (red, also red in the list) hovers
    over the reference at `offset` with arrows both ways: its position is
    the unknown.
    """
    if not 1 <= len(reads) <= 4 or not 0 <= query < len(reads):
        raise ValueError("Show one to four reads and pick one of them as the query.")
    if offset < 0 or offset + len(reads[query]) > len(reference):
        raise ValueError("The hovering read must lie over the reference.")
    _bare(ax)
    title_y, rows = 4.45, [3.6, 3.15, 2.7, 2.25]
    for x, title in [(1.1, "DNA sample"), (4.3, "Fragments"),
                     (7.55, "Sequencer"), (10.75, "Reads")]:
        ax.text(x, title_y, title, ha="center", va="center", fontsize=10)
    for x0 in (2.35, 5.55, 8.95):
        ax.annotate("", xy=(x0 + 0.6, 2.92), xytext=(x0, 2.92),
                    arrowprops=dict(arrowstyle="-|>", color=RULE, lw=1.2,
                                    mutation_scale=12))
    # copies of one region, then the same copies broken at different places
    cuts = [(0.75, 1.45), (1.05,), (0.45, 1.25), (0.9, 1.7)]
    jitter = [(0, 0.05, -0.04), (0.04, -0.05), (-0.03, 0.04, 0), (0.05, -0.03, 0.04)]
    for y, row_cuts, row_jitter in zip(rows, cuts, jitter):
        ax.plot([0.1, 2.1], [y, y], color=INK, lw=2.2, solid_capstyle="round")
        edges = [0.1, *row_cuts, 2.1]
        for a, b, dy in zip(edges[:-1], edges[1:], row_jitter):
            ax.plot([3.26 + a, 3.14 + b], [y + dy] * 2, color=INK, lw=2.2,
                    solid_capstyle="round")
    # the instrument: a box holding a flow cell with fragments attached
    ax.add_patch(mpatches.FancyBboxPatch(
        (6.4, 2.0), 2.3, 1.85, boxstyle="round,pad=0,rounding_size=0.18",
        facecolor=PAPER, edgecolor=INK, linewidth=1.2))
    ax.add_patch(mpatches.FancyBboxPatch(
        (6.65, 2.3), 1.8, 0.85, boxstyle="round,pad=0,rounding_size=0.08",
        facecolor=LIGHT, edgecolor=RULE, linewidth=0.8))
    for k, y in enumerate((2.5, 2.72, 2.94)):
        for x in (6.85 + 0.15 * k, 7.55 - 0.1 * k):
            ax.plot([x, x + 0.45], [y, y], color=INK, lw=1.6,
                    solid_capstyle="round")
    ax.text(7.55, 3.45, "flow cell", ha="center", va="center", fontsize=9)
    for k, (y, read) in enumerate(zip(rows, reads)):
        ax.text(10.75, y, read, ha="center", va="center", fontsize=10,
                family=MONO, color=EMPHASIS if k == query else INK)
    # the question: where on the reference does the read belong?
    pitch, ref_y, read_y = 0.27, 0.4, 1.0
    x0 = 12.0 - pitch * (len(reference) - 1)
    ax.text(x0 - 0.3, ref_y, "reference", ha="right", va="center", fontsize=10)
    for i, base in enumerate(reference):
        ax.text(x0 + pitch * i, ref_y, base, ha="center", va="center",
                fontsize=10, family=MONO, fontweight="bold")
    read = reads[query]
    for i, base in enumerate(read):
        ax.text(x0 + pitch * (offset + i), read_y, base, ha="center",
                va="center", fontsize=10, family=MONO, color=EMPHASIS)
    for end, step in [(x0 + pitch * offset - 0.25, -0.55),
                      (x0 + pitch * (offset + len(read) - 1) + 0.25, 0.55)]:
        ax.annotate("", xy=(end + step, read_y), xytext=(end, read_y),
                    arrowprops=dict(arrowstyle="-|>", color=EMPHASIS, lw=1.2,
                                    mutation_scale=10))
    ax.text(0.0, 0.7, "Where does each\nread belong?", ha="left",
            va="center", fontsize=10)
    ax.set_xlim(-0.1, 12.2)
    ax.set_ylim(0.05, 4.75)


def shared_domain(
    ax: plt.Axes,
    lengths: tuple[int, int],
    starts: tuple[int, int],
    width: int,
    *,
    labels: tuple[str, str] = ("protein 1", "protein 2"),
    domain: str = "homeodomain",
) -> None:
    """Two protein bars sharing one red domain, paired end to end and locally.

    Left: a grey band joins the full lengths, as a global alignment must, so
    each domain lies opposite unrelated residues. Right: the band joins only
    the two domains. Lengths, starts, and width are in residues (0-based
    starts); bars are drawn to one scale.
    """
    for n, s in zip(lengths, starts):
        if s < 0 or s + width > n:
            raise ValueError("Each domain must lie inside its protein.")
    _bare(ax)
    ax.set_aspect("auto")
    scale, h, ys = 1.0 / max(lengths), 0.24, (1.0, 0.0)
    for x_shift, title, local in [(0.0, "Global alignment", False),
                                  (1.25, "Local alignment", True)]:
        ax.text(x_shift + 0.5, 1.55, title, ha="center", va="center",
                fontsize=10)
        for y, n, s in zip(ys, lengths, starts):
            ax.add_patch(mpatches.FancyBboxPatch(
                (x_shift, y - h / 2), n * scale, h,
                boxstyle="round,pad=0,rounding_size=0.04",
                facecolor=PAPER, edgecolor=RULE, linewidth=1.0, zorder=3))
            ax.add_patch(mpatches.Rectangle(
                (x_shift + s * scale, y - h / 2), width * scale, h,
                facecolor=EMPHASIS, edgecolor="none", zorder=4))
        if local:
            a0, b0 = starts[0] * scale, starts[1] * scale
            a1, b1 = a0 + width * scale, b0 + width * scale
        else:
            a0 = b0 = 0.0
            a1, b1 = lengths[0] * scale, lengths[1] * scale
        top, bottom = ys[0] - h / 2, ys[1] + h / 2
        ax.add_patch(mpatches.Polygon(
            [(x_shift + a0, top), (x_shift + a1, top),
             (x_shift + b1, bottom), (x_shift + b0, bottom)],
            closed=True, facecolor=LIGHT, edgecolor="none", zorder=1))
    for y, label in zip(ys, labels):
        ax.text(-0.04, y, label, ha="right", va="center", fontsize=10)
    ax.text((starts[0] + width / 2) * scale, 1.27, domain, ha="center",
            va="center", fontsize=9)
    ax.set_xlim(-0.3, 2.27)
    ax.set_ylim(-0.25, 1.75)


def torn_copies(
    ax: plt.Axes,
    cuts: Sequence[Sequence[float]],
    *,
    bridge: tuple[int, int],
    length: float = 10.0,
) -> None:
    """Copies of one genome, each broken at different places.

    Copy r is a bar of `length` broken at the positions in `cuts[r]`. Piece
    `bridge = (r, k)`, the k-th piece of copy r, is red; a dotted line marks
    the first break of copy 1 (row 0), which that piece must span. The piece
    joins the two sides of that break, as overlapping reads do in assembly.
    """
    pieces = {}
    for row, row_cuts in enumerate(cuts):
        edges = [0.0, *row_cuts, length]
        for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            pieces[row, k] = (a, b)
    if bridge not in pieces or bridge[0] == 0:
        raise ValueError("bridge must name a piece of a copy other than copy 1.")
    a, b = pieces[bridge]
    if not a < cuts[0][0] < b:
        raise ValueError("The bridging piece must span the first break of copy 1.")
    _bare(ax)
    ax.set_aspect("auto")
    gap, pitch = 0.15, 0.9
    for (row, k), (a, b) in pieces.items():
        on = (row, k) == bridge
        ax.plot([a + gap, b - gap], [-row * pitch] * 2, lw=3,
                color=EMPHASIS if on else INK, solid_capstyle="round", zorder=3)
    for row in range(len(cuts)):
        ax.text(-0.3, -row * pitch, f"copy {row + 1}", ha="right",
                va="center", fontsize=10)
    ax.plot([cuts[0][0]] * 2, [0.3, -bridge[0] * pitch - 0.3],
            color=RULE, lw=0.8, ls=":", zorder=1)
    ax.set_xlim(-2.0, length + 0.2)
    ax.set_ylim(-(len(cuts) - 1) * pitch - 0.4, 0.4)


# --- Orientation figures: from a gene to the data that read it ---------------

# The standard genetic code (NCBI translation table 1). Codons are ordered with
# T, C, A, G at each position, so TTT is index 0, TTC is 1, ..., GGG is 63;
# "*" marks a stop codon.
STANDARD_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def translate(dna: str) -> str:
    """Translate a coding DNA string codon by codon with the standard code."""
    if len(dna) % 3 or set(dna) - set("ACGT"):
        raise ValueError("Translate an A/C/G/T string whose length is a multiple of 3.")
    index = {base: k for k, base in enumerate("TCAG")}
    return "".join(STANDARD_CODE[16 * index[c[0]] + 4 * index[c[1]] + index[c[2]]]
                   for c in (dna[i:i + 3] for i in range(0, len(dna), 3)))


def gene_to_protein(
    ax: plt.Axes,
    exons: Sequence[str],
    introns: Sequence[str],
    *,
    codon: int = 0,
    data: Sequence[str] = (),
) -> None:
    """A two-strand gene, its spliced mRNA, and its protein, one row each.

    The coding strand is `exons` and `introns` interleaved (exon, intron,
    exon, ...), written 5' to 3', with the complementary strand below it.
    Grey bands carry each exon down to the mRNA, where U replaces T and the
    introns are gone; brackets group the mRNA into codons, and each codon's
    amino acid (standard code) sits below its bracket. Codon number `codon`
    (0-based) is red-filled in all three rows, so the reader can follow it.
    `data`, if given, holds one short label per row (DNA, mRNA, protein),
    printed at the right, such as the measurement that reads that molecule.
    """
    if not exons or len(introns) != len(exons) - 1:
        raise ValueError("Give one intron between each pair of exons.")
    coding = "".join(exons)
    protein = translate(coding)
    if not 0 <= codon < len(protein):
        raise ValueError("codon must index one codon of the spliced mRNA.")
    if data and len(data) != 3:
        raise ValueError("Give one data label per row, or none.")
    _bare(ax)
    pair = dict(zip("ACGT", "TGCA"))
    y_top, y_bottom, y_rna, y_aa = 0.0, -1.15, -4.8, -7.5
    half = 0.42

    # Lay out the gene; columns[i] is the gene column of mRNA letter i.
    gene, columns, spans = "", [], []
    for k, exon in enumerate(exons):
        spans.append((len(gene), len(gene) + len(exon), "exon"))
        columns += range(len(gene), len(gene) + len(exon))
        gene += exon
        if k < len(introns):
            spans.append((len(gene), len(gene) + len(introns[k]), "intron"))
            gene += introns[k]
    rna_x0 = (len(gene) - len(coding)) / 2
    red_rna = set(range(3 * codon, 3 * codon + 3))
    red_gene = {columns[i] for i in red_rna}

    def letter(x: float, y: float, text: str, red: bool) -> None:
        """One monospace letter, paper on a red fill when highlighted."""
        if red:
            ax.add_patch(mpatches.Rectangle((x - half, y - 0.5), 2 * half, 1.0,
                                            facecolor=EMPHASIS, edgecolor="none",
                                            zorder=2))
        ax.text(x, y, text, ha="center", va="center", fontsize=10, family=MONO,
                color=PAPER if red else INK, zorder=3)

    # DNA: coding strand, its partner strand, and exon/intron brackets above.
    for i, base in enumerate(gene):
        letter(i, y_top, base, i in red_gene)
        letter(i, y_bottom, pair[base], False)
    for y, left_end, right_end in ((y_top, "5′", "3′"), (y_bottom, "3′", "5′")):
        ax.text(-0.9, y, left_end, ha="right", va="center", fontsize=9)
        ax.text(len(gene) - 0.1, y, right_end, ha="left", va="center", fontsize=9)
    exon_number = 0
    for a, b, kind in spans:
        ax.plot([a - 0.35, a - 0.35, b - 0.65, b - 0.65],
                [0.75, 0.95, 0.95, 0.75], color=RULE, lw=0.8)
        if kind == "exon":
            exon_number += 1
        name = f"exon {exon_number}" if kind == "exon" else "intron"
        ax.text((a + b - 1) / 2, 1.15, name, ha="center", va="bottom", fontsize=9)

    # Splicing: one grey band per exon, from the gene to its place in the mRNA.
    rna_start = 0
    for a, b, kind in spans:
        if kind != "exon":
            continue
        c, d = rna_x0 + rna_start, rna_x0 + rna_start + (b - a)
        ax.add_patch(mpatches.Polygon(
            [(a - 0.45, y_bottom - 0.6), (b - 0.55, y_bottom - 0.6),
             (d - 0.55, y_rna + 0.6), (c - 0.45, y_rna + 0.6)],
            closed=True, facecolor=LIGHT, edgecolor="none", zorder=1))
        rna_start += b - a

    # mRNA letters, codon brackets, and amino acids.
    for i, base in enumerate(coding):
        letter(rna_x0 + i, y_rna, "U" if base == "T" else base, i in red_rna)
    for k, residue in enumerate(protein):
        a, b = rna_x0 + 3 * k - 0.35, rna_x0 + 3 * k + 2.35
        ax.plot([a, a, b, b], [y_rna - 0.65, y_rna - 0.85, y_rna - 0.85,
                               y_rna - 0.65], color=RULE, lw=0.8)
        x = rna_x0 + 3 * k + 1
        if residue == "*":
            ax.text(x, y_aa, "stop", ha="center", va="center", fontsize=9)
        else:
            letter(x, y_aa, residue, k == codon)

    # Row names at the left, with the step between rows in italics.
    left = -2.6
    for y, name in (((y_top + y_bottom) / 2, "DNA"), (y_rna, "mRNA"),
                    (y_aa, "protein")):
        ax.text(left, y, name, ha="right", va="center", fontsize=10)
    for y, step in (((y_bottom + y_rna) / 2, "transcription\nand splicing"),
                    ((y_rna + y_aa) / 2, "translation")):
        ax.text(left, y, step, ha="right", va="center", fontsize=9,
                style="italic")
    right = len(gene) + 1.6
    for y, label in zip(((y_top + y_bottom) / 2, y_rna, y_aa), data):
        ax.text(right, y, label, ha="left", va="center", fontsize=9)
    ax.set_xlim(left - 7.2, right + (9.5 if data else 0))
    ax.set_ylim(y_aa - 0.8, 2.0)
