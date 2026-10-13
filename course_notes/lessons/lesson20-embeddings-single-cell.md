---
published: false
title: "Single-Cell Data and Embeddings"
subtitle: "Lesson 20 · 2027-04-07 · Week 13"
nocite: |
  @holmeshuber2019, @macosko2015
---

## Core question

Which cells show similar RNA expression, and how faithfully can a two-dimensional picture display those relationships?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Seeing the cells inside a tissue

A bulk RNA measurement combines signals from many cells. If a tissue contains several kinds of cells, its average expression can hide differences among them. Single-cell RNA sequencing aims to measure RNA separately for individual cells, so we can ask which cells share expression patterns and which genes distinguish them. In droplet methods such as Drop-seq, cells are separated into small droplets and RNA receives sequence tags called *barcodes*, which allow pooled sequencing reads to be assigned back to their cell of origin [@macosko2015].

The resulting table contains a cell identifier, a gene identifier, and a count of captured RNA molecules assigned to that pair after processing. Arrange these counts as $C\in\mathbb Z_{\geq0}^{n\times p}$, with $n$ cells as rows and $p$ genes as columns. Many entries are zero. A zero can mean that a gene had no RNA in that cell, or that its RNA was not captured and detected. We therefore compare a stated normalized representation of the counts, with sampling depth and technical variation in mind. A visual group of cells is a hypothesis about shared expression, not yet a cell-type assignment.

The toy capture-and-tag picture below follows two cells into two count rows. Suppose cell A contributes three captured molecules from gene 1 and one from gene 2, while cell B contributes one and three. Their average row is $(2,2)$, which hides their opposite expression patterns. A cell barcode preserves which cell supplied a molecule; gene assignment supplies the column. A separate molecule tag distinguishes captured molecules, allowing repeated sequencing copies of the same molecule to be collapsed before counting [@macosko2015]. The labels m1 to m4 stand for distinct molecule tags within each cell; the displayed records represent unique captured molecules rather than raw sequencing reads.

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
#| label: fig-single-cell-barcode-counts
#| code-fold: true
#| fig-width: 6
#| fig-height: 3.4
#| fig-cap: "A toy capture-and-tag schematic connects two cells' RNA molecules to barcode-labeled records and count rows. Retaining cell identity reveals expression differences hidden by the average row (2, 2)."
#| fig-alt: "Two cells contribute four molecules each, with opposite three-to-one gene patterns; tags AA and CC identify their records, which become count rows (3, 1) and (1, 3)."
captured_counts = np.array([[3, 1], [1, 3]])
fig, ax = plt.subplots()
cf.cell_barcode_records(
    ax, captured_counts, cells=["Cell A", "Cell B"],
    genes=["gene 1", "gene 2"], barcodes=["AA", "CC"],
)
plt.show()
```

## From expression neighborhoods to a picture

As in Lesson 19, start from feature vectors $x_i\in\mathbb R^d$ in a chosen expression representation, often a selected set of PCA scores. Distances between these vectors define local neighbors. The new task is to place each cell at a point $z_i\in\mathbb R^2$ so that a useful part of those relationships is visible. The collection $Z\in\mathbb R^{n\times2}$ is an *embedding*: new display coordinates for the same cells. Unlike PCA, a nonlinear embedding need not assign its axes to fixed linear combinations of genes.

Two dimensions cannot generally preserve all the distances among high-dimensional expression vectors. Methods such as t-SNE emphasize local similarities, so nearby cells in the input should tend to remain neighbors in the picture. That emphasis can make local structure visible while changing apparent gaps, group sizes, and distances between groups [@holmeshuber2019]. We will compare displays of the same cells under different settings and check their neighborhoods against the input representation. Evidence for a biological population must return to the gene measurements and independent metadata, rather than rest on the outline of a plotted island.

To isolate the layout question, consider a separate six-cell toy neighbor graph: cells A1, A2, A3 select one another, and B1, B2, B3 select one another. The two constructed layouts below keep every within-group connection and its length unchanged, yet have different gaps between groups. These are explanatory drawings, not outputs of a fitted t-SNE model. They show why a local neighbor graph alone cannot tell us how far apart the groups should appear, or whether their separation represents a biological barrier.

```{.python .execute}
#| label: fig-single-cell-layout-gaps
#| code-fold: true
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "Two constructed layouts of one six-cell neighbor graph, drawn on the same coordinate scale. Identical local relationships do not determine the displayed gap between groups."
#| fig-alt: "Two rows of graphs containing the same pair of three-node triangles; the triangles are close in the upper row and much farther apart in the lower row, with their shapes and edges unchanged."
neighbor_edges = [("A1", "A2"), ("A2", "A3"), ("A1", "A3"),
                  ("B1", "B2"), ("B2", "B3"), ("B1", "B3")]
fig, axes = plt.subplots(2, 1)
for ax, gap, title in [(axes[0], 3, "A  Smaller display gap"),
                       (axes[1], 8, "B  Larger display gap")]:
    positions = {"A1": (0, 0), "A2": (0.5, 1), "A3": (1, 0),
                 "B1": (gap, 0), "B2": (gap + 0.5, 1), "B3": (gap + 1, 0)}
    cf.digraph(ax, neighbor_edges, pos=positions, directed=False, node_size=450)
    ax.set(xlim=(-0.5, 9.5), ylim=(-0.5, 1.7), title=title)
plt.show()
```

## The model

- Local-similarity objectives emphasize neighborhood preservation. Distances between separated plotted groups need not measure their expression difference on a common scale.

## The technique

- t-SNE and UMAP at the level of objectives; the standard workflow matrix $\to$ PCA $\to$ kNN $\to$ clusters $\to$ labels.

## Worked example

- The same small dataset embedded under two hyperparameter settings.

## Why it works, and what it costs

- What is optimized, and why the optimum is not unique.

## Where the model breaks

- Cluster size and spacing are artifacts.
- Cell-type labels are claims layered on clusters, not data.

## Exercises

- Explain which steps of the workflow use methods from this course.
- Critique a provided embedding figure.

::: {#refs}
**References**
:::
