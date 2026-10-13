---
published: false
title: "Distances, Clustering, and Neighborhood Graphs"
subtitle: "Lesson 19 · 2027-04-05 · Week 13"
nocite: |
  @holmeshuber2019
---

## Core question

Which cells have similar measurements, and when does that similarity support a biological group?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Finding groups without condition labels

Lessons 17 and 18 began with known control and treatment labels. Suppose instead that we have a mixture of cells and want to find groups with similar properties. *Flow cytometry* measures fluorescent signals from individual cells passing through an instrument. Fluorescent tags can target selected proteins, called markers, so a record contains a cell identifier and an intensity for each measured marker. Cells with similar combinations of marker measurements are candidates for a shared cell type or state [@holmeshuber2019]. Those identities still need biological evidence; the instrument supplies measurements, not type labels.

In a toy two-marker experiment, suppose the recorded values are:

| cell | marker 1 intensity | marker 2 intensity |
|---|---:|---:|
| A | 0 | 0 |
| B | 2 | 0 |
| C | 0 | 10 |
| D | 3 | 1 |
| E | 1 | 9 |

On these raw coordinates, B is nearer to A than C is: their Euclidean distances from A are $2$ and $10$. If we divide marker 2 by $10$ before comparing, C becomes nearer, with distance $1$. Neither calculation changed a cell. We changed how strongly each measured protein contributes to similarity. A clustering result therefore depends on the biological comparison encoded by the scale and distance, as well as on the grouping algorithm.

## From measured cells to distances and neighbors

Let $X\in\mathbb R^{n\times p}$ contain a fixed, stated representation of the measurements, with cells as rows and marker features as columns. Its row $x_i$ is the feature vector of cell $i$. A distance function supplies $d(x_i,x_j)$ for each pair. Clustering assigns group labels using these comparisons; a neighborhood graph instead connects each cell to selected nearby cells. We will construct both kinds of output and ask how changing the comparison changes their biological interpretation.

For a $k$-nearest-neighbor graph, first select the $k$ nearest other cells to each cell. This selection is directed: cell $i$ can select $j$ without $j$ selecting $i$. Here an undirected edge joins $i$ and $j$ when either selects the other, the *union* convention. Stating this rule matters because a mutual-selection convention produces a different graph from the same distances.

The two graphs below use $k=2$ and the same five cells. Their drawing positions are fixed for comparison, so position in the diagram does not represent a measured distance. The highlighted connections show how changing the relative weight of marker 2 changes A's neighbors and therefore the local evidence for grouping A with other cells.

::: {.algorithm}
**Algorithm: construct a union k-nearest-neighbor graph.**

For each cell, calculate its distance to every other cell, sort them from nearest to farthest, and select the first $k$. Break ties by the stated cell order. Include one undirected edge for every selected pair, retaining an edge if either cell selected the other.
:::

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

```{.python .execute}
#| label: fig-neighbors-marker-scale
#| code-fold: true
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "Union two-nearest-neighbor graphs from the same five cells before and after rescaling marker 2; red edges connect A to its neighbors. Marker weighting changes the relationships supplied to a clustering method."
#| fig-alt: "Two graphs with the same five fixed node positions: A connects to B and D in the raw graph, and to C and E after marker 2 is divided by ten."
def union_neighbors(X: np.ndarray, names: list[str], k: int) -> list[tuple[str, str]]:
    """Connect a pair if either cell selects the other among its k nearest."""
    edges = set()
    for i in range(len(X)):
        distances = np.linalg.norm(X - X[i], axis=1)
        distances[i] = np.inf
        order = np.argsort(distances, kind="stable")
        for j in order[:k]:
            pair = (min(i, int(j)), max(i, int(j)))
            edges.add(pair)
    return [(names[i], names[j]) for i, j in sorted(edges)]


names = ["A", "B", "C", "D", "E"]
raw_markers = np.array([[0, 0], [2, 0], [0, 10], [3, 1], [1, 9]], dtype=float)
scaled_markers = raw_markers / np.array([1, 10])
positions = {"A": (0, 0), "B": (2, 0), "C": (0, 2),
             "D": (3, 1), "E": (1, 2)}
fig, axes = plt.subplots(1, 2)
for ax, markers, title in [
    (axes[0], raw_markers, "A  Raw marker intensities"),
    (axes[1], scaled_markers, "B  Marker 2 divided by 10"),
]:
    edges = union_neighbors(markers, names, k=2)
    a_edges = [edge for edge in edges if "A" in edge]
    cf.digraph(ax, edges, pos=positions, directed=False,
               highlight_edges=a_edges, node_size=500)
    ax.set_title(title)
plt.show()
```

## The model

- A distance is a modeling choice: Euclidean, correlation, cosine each assume something different.

## The technique

- k-means and hierarchical clustering briefly; kNN graph construction; the curse of dimensionality by a small computation.

## Worked example

- Use the five cells above to construct a union kNN graph for $k = 2$. Compare raw Euclidean distance with Euclidean distance after dividing marker 2 by $10$; use alphabetical cell order to break any distance ties.

## Why it works, and what it costs

- kNN graph cost; why distances concentrate as $p$ grows.

## Where the model breaks

- Clusters found by every method on data with no clusters.

## Exercises

- Build a kNN graph by hand.
- Give a dataset where Euclidean and correlation distance disagree about the nearest neighbor.

## Schedule

PS8 due. PS9 out.

::: {#refs}
**References**
:::
