---
published: false
title: "Single-Cell Data and Embeddings"
subtitle: "Lesson 20 · 2027-04-07 · Week 13"
nocite: |
  @holmeshuber2019, @macosko2015, @vandermaaten2008, @kobak2019, @wattenberg2016, @mcinnes2018
---

## Core question

Which cells show similar RNA expression, and how faithfully can a two-dimensional picture display those relationships?

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## Seeing the cells inside a tissue

A mouse retina contains rods, cones, bipolar cells, and several other neuron and support-cell classes, each with a different role in turning light into a nerve signal. A bulk RNA measurement grinds the tissue together and reports one average profile, so it cannot say which genes are used by which cells. Single-cell RNA sequencing measures RNA separately in each of thousands of cells. @fig-single-cell-retina-atlas shows the kind of picture such an experiment produces: 44,808 retina cells measured by the droplet method Drop-seq [@macosko2015], each drawn as a dot, positioned by a method called t-SNE, and colored by the cell class the original study assigned.

::: {#fig-single-cell-retina-atlas .course-figure}

<img src="../assets/borrowed/kobak-berens-2019-retina-tsne.png" alt="A scatter plot of tens of thousands of small dots forming colored islands. A very large light-blue island is labeled Rods, a smaller light-blue island at the left is labeled Cones, dark-blue islands at the lower right are labeled Bipolar cells, green islands Amacrine cells, and smaller islands are labeled Ganglion cells, Horizontal cells, MG, and single letters." loading="lazy">

A t-SNE embedding of 44,808 mouse retina cells from a Drop-seq experiment, colored by the cell classes assigned in the original study; single letters mark smaller groups named there. Cells of one class form compact islands, while the axes, island areas, and gaps have no units. Reproduced, cropped to panel a, from @kobak2019 (CC BY 4.0); data from @macosko2015.

:::

Reading such a picture requires knowing how the experiment turned cells into numbers and how the method turned those numbers into positions. Lesson 19 compared cells through distances and neighbor graphs; t-SNE uses the same neighbor relationships to place every cell in the plane, and some features of the result are evidence about the cells while others are not.

### From droplets to a count matrix

In Drop-seq, a microfluidic device encloses single cells with single beads in nanoliter droplets. Each bead carries many DNA capture strands. All strands on one bead share a *cell barcode* sequence that differs between beads, and each strand also carries its own random *unique molecular identifier* (UMI). Inside a droplet the cell bursts, its messenger RNA sticks to the bead's strands, and copying each captured RNA into DNA attaches the barcode and one UMI to it. The droplets are then broken, and everything is amplified and sequenced together [@macosko2015].

Each read therefore records a gene (by matching to gene sequences, as in Lesson 15), a cell (by its barcode), and a captured molecule (by its UMI). Amplification makes many reads from one molecule, all sharing barcode, gene, and UMI, so counting distinct UMIs rather than reads counts molecules. @fig-single-cell-barcode-counts follows two toy cells through this bookkeeping. Cell A contributes three molecules from gene 1 and one from gene 2, cell B the reverse, and their average row $(2,2)$ would hide both patterns.

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

The result is a count matrix $C\in\mathbb Z_{\geq0}^{n\times p}$, where $C_{ij}$ is the number of distinct UMIs for gene $j$ in cell $i$. Cells take the rows that samples occupied in Lessons 15 and 16, and an experiment can have thousands to millions of them.

### Why most entries are zero

A cell yields relatively few captured molecules, spread over thousands of genes. Recall Lesson 16's model $\mathbb E[C_{ij}\mid L_i]=L_i\rho_{ij}$, with $L_i=\sum_jC_{ij}$ the cell's total and $\rho_{ij}$ the expected fraction from gene $j$. Suppose the cell holds $N$ copies of gene $j$'s RNA and each is captured independently with small probability $q$. Then no copy is captured with probability
$$
(1-q)^N=e^{N\ln(1-q)}\approx e^{-Nq}=e^{-m},
$$
where $m=Nq$ is the expected count and we used $\ln(1-q)\approx-q$. A gene expected to supply one molecule in a thousand, in a cell with $L_i=300$, has $m=0.3$ and reads zero with probability $e^{-0.3}\approx0.74$. A zero, then, usually says that the expected count was small, not that the gene was off. Such sampling zeros are often called *dropouts*.

To have data with a known answer, we simulate an experiment. Six cell types share a baseline profile over $p=300$ genes. Types 1 to 4 each raise their own block of 12 marker genes eightfold; types 5a and 5b share a fifth block and differ only in four genes each. The types have 150, 80, 50, 40, 20, and 20 cells, every cell has the same amount of random variation around its type's profile, and each captures about 300 molecules by the independent-capture model.

```{.python .execute}
#| label: code-single-cell-simulate
#| code-fold: true
#| code-summary: "Simulate UMI counts for six cell types"
TYPE_NAMES = ["1", "2", "3", "4", "5a", "5b"]
TYPE_MARKERS = ["o", "s", "^", "D", "v", "P"]


def simulate_counts(seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Simulated UMI counts (cells x genes) and the true type of each cell."""
    rng = np.random.default_rng(seed)
    n_genes, block, fold = 300, 12, np.log(8)
    log_profile = np.tile(rng.normal(0, 1.5, n_genes), (6, 1))
    for t in range(5):  # types 1-4, then the program shared by 5a and 5b
        log_profile[t, block * t:block * (t + 1)] += fold
    log_profile[5, 4 * block:5 * block] += fold
    log_profile[4, 60:64] += fold  # 5a and 5b differ in four genes each
    log_profile[5, 64:68] += fold
    counts, labels = [], []
    for t, size in enumerate([150, 80, 50, 40, 20, 20]):
        for _ in range(size):
            rho = np.exp(log_profile[t] + rng.normal(0, 0.25, n_genes))
            depth = rng.lognormal(np.log(300), 0.4)
            # independent capture of many copies gives Poisson counts
            counts.append(rng.poisson(depth * rho / rho.sum()))
            labels.append(t)
    return np.array(counts), np.array(labels)


counts, cell_type = simulate_counts()
```

```{.python .execute}
#| label: fig-single-cell-count-heatmap
#| code-fold: true
#| fig-width: 6
#| fig-height: 3.2
#| fig-cap: "Simulated UMI counts for five cells of each type over the first 72 genes, which contain every marker block; white entries are zeros. Even marker genes are often recorded as zero, so a cell's type shows only as a pattern across many genes."
#| fig-alt: "A heat map with 30 rows grouped by cell type and 72 gene columns. Many entries are white; colored counts concentrate in diagonal blocks where each type's marker genes lie, with white gaps inside the blocks."
shown_cells = np.concatenate([np.flatnonzero(cell_type == t)[:5] for t in range(6)])
shown = np.ma.masked_equal(counts[shown_cells, :72], 0)
fig, ax = plt.subplots()
image = ax.imshow(shown, aspect="auto", interpolation="nearest",
                  norm=matplotlib.colors.LogNorm(vmin=1))
ax.set_yticks(np.arange(2, 30, 5), [f"type {name}" for name in TYPE_NAMES])
ax.set_xlabel("gene")
fig.colorbar(image, ax=ax, label="UMI count (log scale)")
plt.show()
```

In @fig-single-cell-count-heatmap the marker blocks form a diagonal pattern, but each cell shows only scattered counts within its block, and types 5a and 5b look nearly identical. Real data over twenty thousand genes is sparser still. The task is to recover groups of similar cells from many such noisy rows and to draw them faithfully.

## From counts to a cell-by-cell comparison

As in Lesson 16, we first put cells on a common scale:
$$
X_{ij}=\log_2\!\Bigl(1+10^4\,\frac{C_{ij}}{L_i}\Bigr).
$$
The factor $10^4$ replaces the per-million factor used for bulk samples because a cell has only hundreds to thousands of molecules. We then center the columns to form $Z$ and keep the first $r=20$ PCA components. Cell $i$ is represented by its score vector $t_i=z_iV_r\in\mathbb R^{r}$, the $r$-direction version of Lesson 15's score $t_i=z_iv$, and distances $d_{ij}=\|t_i-t_j\|$ define its neighbors. (Lesson 19 wrote $x_i$ for whatever representation was chosen.)

```{.python .execute}
#| label: code-single-cell-pca
totals = counts.sum(axis=1)
X = np.log2(1 + 1e4 * counts / totals[:, None])
Z = X - X.mean(axis=0)
_, singular_values, Vt = np.linalg.svd(Z, full_matrices=False)
T = Z @ Vt[:20].T  # row i is the score vector t_i
explained = singular_values**2 / np.sum(singular_values**2)
print("fraction of zero counts:", round(np.mean(counts == 0), 3))
print("explained by PC1-PC5:", np.round(explained[:5], 3))
```

About two-thirds of the counts are zero, and no component dominates: PC1 and PC2 explain under 5% each, with no sharp drop afterward. Five expression programs raise disjoint gene blocks, so the type means differ along about five directions of comparable size. A PCA plot shows two of them, and types separated only along the other three land on top of each other, as in panel A of @fig-single-cell-embeddings below. Most of the remaining variance is sampling noise spread over all 300 genes.

The 20 scores still recognize the types: on average 92% of each cell's ten nearest neighbors in $\mathbb R^{20}$ share its type (computed later). The information is in the neighbor relationships, and a linear projection to the plane discards it. We therefore seek points $y_1,\ldots,y_n\in\mathbb R^2$, an *embedding* $Y\in\mathbb R^{n\times2}$, that keep each cell's near neighbors close while letting large distances distort.

## Neighbor probabilities

The first step turns each cell's distances into a statement about who its neighbors are. Instead of the yes-or-no choice of a $k$-nearest-neighbor graph, t-SNE gives every other cell a probability of being chosen that falls quickly with distance [@vandermaaten2008].

::: {#def-neighbor-probabilities}
## Neighbor probabilities and perplexity
For score vectors $t_1,\ldots,t_n$ with distances $d_{ij}$ and a bandwidth $\sigma_i>0$ for each cell, the *conditional neighbor probabilities* are
$$
p_{j\mid i}=\frac{\exp\bigl(-d_{ij}^2/2\sigma_i^2\bigr)}{\sum_{k\ne i}\exp\bigl(-d_{ik}^2/2\sigma_i^2\bigr)}\quad(j\ne i),
\qquad p_{i\mid i}=0.
$$
The *perplexity* of row $i$ is $\operatorname{Perp}_i=\exp\bigl(-\sum_{j}p_{j\mid i}\ln p_{j\mid i}\bigr)$, with $0\ln0=0$. The *joint neighbor probabilities* are $p_{ij}=(p_{j\mid i}+p_{i\mid j})/(2n)$.
:::

Read $p_{j\mid i}$ as two operations: weight every other cell by a Gaussian function of its distance from cell $i$, then divide by the total so the row sums to one. Perplexity converts a row into an effective number of neighbors. A row with probability $1/k$ on each of $k$ cells has $-\sum_jp_{j\mid i}\ln p_{j\mid i}=\ln k$, so its perplexity is $k$; a row concentrated on one cell has perplexity near one.

t-SNE fixes a target perplexity, typically between 5 and 50, and chooses each $\sigma_i$ so that every row reaches it. A cell in a dense region gets a small $\sigma_i$ and a cell in a sparse region a large one. Perplexity thus plays the role of $k$ in a neighbor graph, and the per-cell bandwidth deliberately divides out differences in local density. The joint probability combines both directions of choice, as Lesson 19's union rule kept an edge if either cell chose the other, and the divisor $2n$ makes all $p_{ij}$ sum to one.

::: {#exm-tsne-three-cells}
## Three cells on a line
Place cells at $t=0,1,3$ and, to compute by hand, fix $\sigma_i=1$ instead of calibrating. From cell 1 the squared distances are 1 and 9, with weights $e^{-1/2}=0.607$ and $e^{-9/2}=0.011$; dividing by their sum $0.618$ gives $p_{2\mid1}=0.982$ and $p_{3\mid1}=0.018$. The other rows follow in the same way:

| from cell | to cell 1 | to cell 2 | to cell 3 |
|---|---:|---:|---:|
| 1 | 0 | 0.982 | 0.018 |
| 2 | 0.818 | 0 | 0.182 |
| 3 | 0.076 | 0.924 | 0 |

Symmetrizing gives $p_{12}=(0.982+0.818)/6=0.300$, $p_{23}=(0.182+0.924)/6=0.184$, and $p_{13}=0.016$. Each unordered pair appears twice in the sum over $i\ne j$, and $2(0.300+0.184+0.016)=1$.
:::

Cell 3's nearest cell is two units away, yet it gives that cell 92% of its probability: neighbor probabilities are relative to each cell's surroundings. The code below computes a row and its perplexity, then calibrates $\sigma_i$ by bisection on $\ln\sigma_i$. Perplexity increases with $\sigma_i$, so a trial bandwidth with too many effective neighbors is shrunk and one with too few is enlarged.

```{.python .execute}
#| label: code-single-cell-neighbor-probabilities
def neighbor_row(sq_dists: np.ndarray, sigma: float) -> np.ndarray:
    """Conditional probabilities p_{j|i} from cell i's squared distances to others."""
    # subtracting the minimum rescales every weight equally and avoids underflow
    weights = np.exp(-(sq_dists - sq_dists.min()) / (2 * sigma**2))
    return weights / weights.sum()


def perplexity(p: np.ndarray) -> float:
    """Effective number of neighbors: exp of the entropy of a probability row."""
    p = p[p > 0]
    return float(np.exp(-np.sum(p * np.log(p))))


def calibrate_row(sq_dists: np.ndarray, target: float) -> np.ndarray:
    """Bisect on log(sigma) until the row's perplexity matches the target."""
    low, high = -10.0, 10.0
    for _ in range(50):
        row = neighbor_row(sq_dists, np.exp((low + high) / 2))
        if perplexity(row) > target:
            high = (low + high) / 2
        else:
            low = (low + high) / 2
    return row


def squared_distances(A: np.ndarray) -> np.ndarray:
    """All pairwise squared Euclidean distances between the rows of A."""
    sq = np.sum(A**2, axis=1)
    return np.maximum(sq[:, None] + sq[None, :] - 2 * A @ A.T, 0)


def joint_probabilities(T: np.ndarray, target: float) -> np.ndarray:
    """Symmetric neighbor probabilities p_ij, summing to one over all pairs."""
    n = len(T)
    D2 = squared_distances(T)
    conditional = np.zeros((n, n))
    for i in range(n):
        others = np.arange(n) != i
        conditional[i, others] = calibrate_row(D2[i, others], target)
    return (conditional + conditional.T) / (2 * n)


row_from_cell_1 = neighbor_row(np.array([1.0, 9.0]), sigma=1.0)
print(np.round(row_from_cell_1, 3), round(perplexity(row_from_cell_1), 3))
```

The row agrees with the hand calculation, and its perplexity of 1.09 says that cell 1 effectively has one neighbor.

## Placing cells in the plane

To choose the display points $y_i$, t-SNE describes the plane with similar probabilities but a different kernel:
$$
q_{ij}=\frac{w_{ij}}{\sum_{k\ne l}w_{kl}},
\qquad w_{ij}=\frac{1}{1+\|y_i-y_j\|^2}\quad(i\ne j),
\qquad q_{ii}=0.
$$
First compute a similarity $w_{ij}$ for every pair of display points, then normalize over all ordered pairs. This kernel decays like $\|y_i-y_j\|^{-2}$, far more slowly than a Gaussian, and the heavy tail addresses a shortage of room. In 20 dimensions a cell can have dozens of cells at moderate, nearly equal distances, which the plane cannot reproduce. With a slowly decaying kernel, moderately distant pairs can be pushed well apart at little cost, which opens gaps between groups [@vandermaaten2008].

::: {#def-tsne-cost}
## t-SNE embedding
**Input:** joint neighbor probabilities $p_{ij}$ for $n$ cells.

**Output:** a configuration $Y=(y_1,\ldots,y_n)$ in $\mathbb R^2$, a *t-SNE embedding*, that locally minimizes
$$
\mathcal C(Y)=\sum_{i\ne j,\;p_{ij}>0}p_{ij}\ln\frac{p_{ij}}{q_{ij}}.
$$
:::

$\mathcal C$, called the Kullback-Leibler divergence, is zero when every $q_{ij}=p_{ij}$. Each term is weighted by $p_{ij}$. Input neighbors have large $p_{ij}$, and placing them far apart makes $q_{ij}$ small and the term large. When $p_{ij}\approx0$ the term is negligible wherever the two cells sit. Crowding unrelated cells together is penalized only indirectly: because the $q_{ij}$ sum to one, probability spent on unrelated pairs is taken from neighbor pairs.

Continuing @exm-tsne-three-cells, try $y=(0,1,3)$, a copy of the input. The kernel values for pairs 12, 13, and 23 are $1/2$, $1/10$, and $1/5$, which total $1.6$ over ordered pairs, so $q_{12}=0.3125$, $q_{13}=0.0625$, and $q_{23}=0.125$. The pair terms are $0.300\ln(0.300/0.3125)=-0.0122$, $0.016\ln(0.016/0.0625)=-0.0218$, and $0.184\ln(0.184/0.125)=0.0712$; counting each pair twice gives $\mathcal C\approx0.074$. The code compares three other candidates, using the rounded probabilities of the table.

```{.python .execute}
#| label: code-single-cell-cost
def low_dim_affinities(Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Kernel values w_ij and normalized probabilities q_ij for an embedding."""
    W = 1 / (1 + squared_distances(Y))
    np.fill_diagonal(W, 0)
    return W, W / W.sum()


def tsne_cost(P: np.ndarray, Y: np.ndarray) -> float:
    """C(Y): sum of p_ij log(p_ij / q_ij) over pairs with p_ij > 0."""
    _, Q = low_dim_affinities(Y)
    keep = P > 0
    return float(np.sum(P[keep] * np.log(P[keep] / Q[keep])))


P_line = np.array([[0, 0.300, 0.016], [0.300, 0, 0.184], [0.016, 0.184, 0]])
for y in [(0, 1, 3), (0, 2, 6), (0, 1, 1.5), (0, 3, 1)]:
    print(y, round(tsne_cost(P_line, np.array(y, dtype=float)[:, None]), 3))
```

Doubling every distance raises the cost only from $0.074$ to $0.087$. Pulling cell 3 in to $1.5$, closer to cell 1 than the input suggests, triples it. Swapping cells 2 and 3, so that cell 2 sits far from both of its neighbors, multiplies it by thirteen. The objective cares a great deal about who is next to whom and much less about how far apart things are.

### How the cost moves the points

To minimize $\mathcal C$ we follow its gradient [@vandermaaten2008]:
$$
\frac{\partial\mathcal C}{\partial y_i}=4\sum_{j\ne i}\bigl(p_{ij}-q_{ij}\bigr)\,w_{ij}\,(y_i-y_j).
$$ {#eq-tsne-gradient}
Gradient descent moves $y_i$ against this vector, so each cell $j$ exerts a force along the line joining $y_i$ and $y_j$. The $p_{ij}$ part pulls $y_i$ toward input neighbors; the $q_{ij}$ part pushes it away from every cell it currently sits near. For two cells in well-separated groups, $p_{ij}\approx0$ and only the push remains, with size proportional to $w_{ij}^2\|y_i-y_j\|$, which falls like $\|y_i-y_j\|^{-3}$. Once two islands have separated, the forces between them are negligible, so their final distance reflects where the optimization left them rather than how different their expression is [@kobak2019].

Rotating, reflecting, or shifting every point by $y_i'=Ry_i+b$, with $R$ an orthogonal $2\times2$ matrix, leaves every distance $\|y_i-y_j\|$ unchanged, hence every $w_{ij}$, every $q_{ij}$, and the cost. A minimizer is therefore never unique, and the axes of a t-SNE plot mean nothing, unlike PCA axes, which are fixed gene combinations. The cost is also not convex, so different random starts can reach different local minima with islands arranged differently.

::: {.algorithm}
**Algorithm: t-SNE by gradient descent.**
Compute the $p_{ij}$ at the target perplexity and start from small random positions with a fixed seed. Repeatedly compute $w_{ij}$, $q_{ij}$, and @eq-tsne-gradient, then move each $y_i$ a step against the gradient plus a fraction of its previous move (momentum). For the first 250 iterations multiply every $p_{ij}$ by 12 (early exaggeration), so attraction dominates while groups form.
:::

In the code, $\sum_jF_{ij}(y_i-y_j)$ with $F_{ij}=(p_{ij}-q_{ij})w_{ij}$ is computed for all cells at once as $\bigl(\sum_jF_{ij}\bigr)y_i-(FY)_i$.

```{.python .execute}
#| label: code-single-cell-tsne
#| code-fold: true
#| code-summary: "t-SNE by gradient descent (implementation)"
def tsne(T: np.ndarray, target: float, seed: int = 0,
         n_iter: int = 1000, step: float = 200.0) -> np.ndarray:
    """Two-dimensional t-SNE embedding of the rows of T by gradient descent."""
    rng = np.random.default_rng(seed)
    P = joint_probabilities(T, target)
    Y = rng.normal(0, 1e-4, size=(len(T), 2))
    velocity = np.zeros_like(Y)
    for iteration in range(n_iter):
        boost, momentum = (12.0, 0.5) if iteration < 250 else (1.0, 0.8)
        W, Q = low_dim_affinities(Y)
        F = (boost * P - Q) * W
        gradient = 4 * (F.sum(axis=1)[:, None] * Y - F @ Y)
        velocity = momentum * velocity - step * gradient
        Y = Y + velocity
    return Y


embeddings = {target: tsne(T, target, seed=0) for target in (5, 30, 100)}
```

## t-SNE on the simulated cells

@fig-single-cell-embeddings compares the first two PCA scores with t-SNE embeddings computed from the same 20 scores at three perplexities. At perplexity 30 (panel C), types 1 to 4 each form an island, and 5a and 5b share a fifth, as their construction suggests: two subtypes differing in eight genes are each other's closest relatives.

```{.python .execute}
#| label: fig-single-cell-embeddings
#| code-fold: true
#| fig-width: 6
#| fig-height: 5.2
#| fig-cap: "The simulated cells drawn by their first two PCA scores (A) and by t-SNE from 20 PCA scores at perplexity 5, 30, and 100 (B to D); color and marker show the true type. The linear projection overlaps several types, while every embedding separates the same five groups in a different arrangement."
#| fig-alt: "Four scatter plots of 360 points in six colors and marker shapes. In panel A, type 1 separates at positive PC1 while the other types overlap at negative PC1. Panels B, C, and D each show five colored islands, with 5a and 5b sharing one island; a few points sit between islands at perplexity 5, and the islands are placed differently in each panel."
fig, axes = plt.subplots(2, 2)
panels = [(T[:, :2], "A  PCA")] + [
    (Y, f"{label}  t-SNE, perplexity {target}")
    for label, (target, Y) in zip("BCD", embeddings.items())]
for ax, (coords, title) in zip(axes.flat, panels):
    for t, name in enumerate(TYPE_NAMES):
        chosen = cell_type == t
        ax.scatter(coords[chosen, 0], coords[chosen, 1], s=6, color=f"C{t}",
                   marker=TYPE_MARKERS[t], label=f"type {name}")
    ax.set_title(title, loc="left", fontweight="bold")
    if title.startswith("A"):
        ax.set(xlabel=f"PC1 ({100 * explained[0]:.1f}%)",
               ylabel=f"PC2 ({100 * explained[1]:.1f}%)")
    else:
        ax.set(xticks=[], yticks=[])
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside lower center", ncol=6,
           handletextpad=0.2, columnspacing=0.8)
plt.show()
```

Because we simulated the data, we can test the picture against the truth. The next chunk asks two questions of each representation: what fraction of each cell's ten nearest neighbors share its type, and how far apart the centers of types 1 to 4 lie.

```{.python .execute}
#| label: code-single-cell-checks
def same_type_fraction(A: np.ndarray, k: int = 10) -> float:
    """Average fraction of each cell's k nearest neighbors that share its type."""
    D2 = squared_distances(A)
    np.fill_diagonal(D2, np.inf)
    neighbors = np.argsort(D2, axis=1, kind="stable")[:, :k]
    return float(np.mean(cell_type[neighbors] == cell_type[:, None]))


def center_gaps(A: np.ndarray) -> dict[str, float]:
    """Distances between the centers of types 1-4, keyed by the pair of types."""
    centers = [A[cell_type == t].mean(axis=0) for t in range(4)]
    gaps = {}
    for a in range(4):
        for b in range(a + 1, 4):
            gaps[f"{a + 1}-{b + 1}"] = round(float(np.linalg.norm(centers[a] - centers[b])), 1)
    return gaps


for name, A in [("20 PCs", T), ("PC1-PC2", T[:, :2]), ("t-SNE 30", embeddings[30])]:
    print(f"{name:9s} same-type {same_type_fraction(A):.2f}; gaps {center_gaps(A)}")
```

Neighbors in the embedding share a type as often as neighbors in the 20 scores, and far more often than neighbors in the PCA plane, so membership in an island is supported by the input neighborhoods. The center gaps are not. In the input the four types are nearly equidistant, with center distances between 17.6 and 20.5, and types 3 and 4 are the farthest pair. In the embedding, types 3 and 4 are the closest pair, and types 1 and 2 are two and a half times as far apart. The spacing of islands is a by-product of the optimization, as the gradient predicted.

## Reading an embedding

Perplexity changes the picture without any change in the data. All three t-SNE panels of @fig-single-cell-embeddings recover the same five groups, but the shapes and arrangement differ. At perplexity 5 each cell keeps only a few neighbors, and a few cells end up stranded between islands; in noisier real data, small perplexities can split one population into fragments, and even structureless random data can form convincing clumps [@wattenberg2016]. No setting is the true one. The five groups are worth reporting because they appear at every perplexity we tried; their arrangement, which changes with perplexity, is not.

Island area is unreliable too. Every simulated type has the same variation around its profile, yet type 1, with 150 cells, covers a far larger area than type 5b, with 20. Because each $\sigma_i$ was tuned to give every cell the same number of effective neighbors, differences in local density were removed in the first step, and area tracks the number of cells more than their variability [@wattenberg2016]. The large rod island in @fig-single-cell-retina-atlas reflects at least partly how many rods were captured.

::: {.callout-warning title="Gaps and areas in an embedding are not measurements"}
A t-SNE or UMAP plot can support the claim that cells in one island were each other's neighbors in the input. It does not show how different two islands are, how variable a group is, or that an empty region contains no intermediate cells.
:::

UMAP, the other widely used embedding, builds the $k$-nearest-neighbor graph of Lesson 19 explicitly, with $k$ playing the role of perplexity, and lays it out by pulling neighbors together and pushing sampled non-neighbors apart [@mcinnes2018]. The warning above applies to it unchanged.

In a typical analysis, clusters are computed on the neighbor graph built from PCA scores, as in Lesson 19, and the embedding is drawn alongside for display, not used to define them. Naming a cluster as a cell type is a further biological claim, made from marker genes; the colors in @fig-single-cell-retina-atlas are such claims [@macosko2015].

## Exercises

::: {#exr-tsne-three-cells}
Three cells lie on a line at $t=0,2,3$. With $\sigma_i=1$ for every cell, compute all $p_{j\mid i}$ and $p_{ij}$. Which pair has the largest joint probability, and is it the closest pair?
:::

::: {#exr-tsne-perplexity}
Compute the perplexity of the conditional row $(1/2,1/4,1/4)$ and explain why it lies between 2 and 3.
:::

::: {#exr-tsne-gradient}
Derive @eq-tsne-gradient. Write $S=\sum_{k\ne l}w_{kl}$ and $\ln q_{kl}=\ln w_{kl}-\ln S$, use $\sum_{k\ne l}p_{kl}=1$, and show that the repulsion comes entirely from $\ln S$.
:::

::: {#exr-tsne-neighbor-check}
In real data the types are unknown. Write a label-free check: for each cell, the fraction of its ten nearest neighbors in the 20 PCA scores that remain among its ten nearest in an embedding. Compare the first two PCA scores with the perplexity-30 embedding, then rerun t-SNE with `seed=1` and compare the center gaps.
:::

::: {#exr-tsne-atlas-critique}
For each statement about @fig-single-cell-retina-atlas, decide whether the figure supports it, and name the additional evidence needed if not. (a) Cones are more similar in expression to rods than bipolar cells are. (b) Bipolar cells comprise several distinct groups of similar cells. (c) Rods vary more in expression than amacrine cells. (d) The cells labeled amacrine are amacrine cells.
:::

::: {#refs}
**References**
:::
