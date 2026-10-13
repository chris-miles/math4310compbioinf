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
