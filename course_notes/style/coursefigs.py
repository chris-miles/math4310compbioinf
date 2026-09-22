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
    highlight_edges: Iterable[tuple[str, str]] = (),
    highlight_nodes: Iterable[str] = (),
    node_size: int = 900,
    directed: bool = True,
) -> dict[str, tuple[float, float]]:
    """Draw a graph with circular white nodes and grey edges.

    Directed by default, with arrowheads; pass `directed=False` for a kNN or
    similarity graph. Repeated edges in `edges` are drawn as parallel arcs.
    Highlighted edges and nodes use the emphasis color. Returns the node
    positions so a second panel can reuse the same layout.
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
        on = (u, v) in highlight_edges or (not directed and (v, u) in highlight_edges)
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
) -> None:
    """Draw a rooted tree with branch lengths as horizontal extent.

    `root` is a nested dict: {"name": str, "length": float, "children": [...]}.
    Leaves are drawn at the right, in the order they appear. Names of nodes
    in `highlight` are colored with the emphasis color.
    """
    _bare(ax)
    highlight = set(highlight)
    leaf_y: dict[int, float] = {}
    counter = [0]

    def place(node: dict, x0: float) -> tuple[float, float]:
        x1 = x0 + node.get("length", 0.0)
        children = node.get("children", [])
        if not children:
            y = -counter[0]
            counter[0] += 1
            ax.plot([x0, x1], [y, y], color=RULE, lw=1.2, solid_capstyle="butt")
            name = node.get("name", "")
            ax.text(x1 + 0.03, y, name, ha="left", va="center", fontsize=10,
                    color=EMPHASIS if name in highlight else INK)
            return x1, y
        ys = []
        for child in children:
            _, cy = place(child, x1)
            ys.append(cy)
        y = float(np.mean(ys))
        ax.plot([x1, x1], [min(ys), max(ys)], color=RULE, lw=1.2)
        ax.plot([x0, x1], [y, y], color=RULE, lw=1.2)
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
