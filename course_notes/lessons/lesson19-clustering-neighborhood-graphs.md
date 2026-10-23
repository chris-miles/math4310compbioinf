---
published: false
title: "Distances, Clustering, and Neighborhood Graphs"
subtitle: "Lesson 19 · 2027-04-05 · Week 13"
nocite: |
  @holmeshuber2019
---

## Core question

Which cells have similar measurements, and when does a group of similar cells support a biological population?

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## Finding cell populations without labels

Lessons 17 and 18 compared conditions whose labels we knew in advance. Many experiments supply no labels at all. A blood sample mixes helper T cells, cytotoxic T cells, B cells, and other types, and nothing in the tube says which cell is which. We want to find groups of cells with similar measurements, then ask whether each group is a known cell type, an unrecognized one, or nothing in particular.

*Flow cytometry* measures cells one at a time as a fluid stream carries them past a laser. The cells are first stained with antibodies that bind particular surface proteins, called *markers*; each antibody carries a fluorescent dye, and detectors record how brightly each cell glows in each dye's color. One record is a cell identifier followed by one intensity per marker, usually after an $\operatorname{asinh}$ transformation that behaves like a logarithm for large values. A typical experiment measures about ten markers on many thousands of cells [@holmeshuber2019, sec. 5.5].

Helper T cells carry the marker CD4, and cytotoxic T cells carry CD8. Nobody has labeled the simulated cells in @fig-clustering-cytometry-cloud, yet the eye finds three clouds. Immunologists traditionally *gate* such data by drawing a polygon around each cloud. With forty markers, or the thousands of genes per cell measured by single-cell RNA sequencing, there is no single plot on which to draw, and we need a computational definition of a group of similar cells.

```{.python .execute}
#| label: fig-clustering-cytometry-cloud
#| code-fold: true
#| fig-width: 4.5
#| fig-height: 4
#| fig-cap: "Simulated CD4 and CD8 intensities for 750 unlabeled cells; each point is one cell. The instrument supplies only coordinates, and the grey annotations are what an analyst would infer from the three clouds."
#| fig-alt: "Scatter plot of CD4 against CD8 with three dense clouds: high CD4 only, high CD8 only, and low values of both."
rng = np.random.default_rng(19)
cloud_centers = np.array([[4.0, 0.8], [0.9, 3.8], [0.8, 0.7]])
clouds = []
for center, size in zip(cloud_centers, [320, 180, 250]):
    clouds.append(center + 0.45 * rng.standard_normal((size, 2)))
cytometry = np.vstack(clouds)

fig, ax = plt.subplots()
ax.scatter(cytometry[:, 0], cytometry[:, 1], s=5, alpha=0.5, edgecolors="none")
for (x, y), text in zip(cloud_centers, ["CD4+ CD8−", "CD4− CD8+", "CD4− CD8−"]):
    ax.text(x + 0.6, y + 1.0, text, fontsize=9, color=cf.RULE)
ax.set(xlabel="CD4 intensity (asinh scale)", ylabel="CD8 intensity (asinh scale)")
ax.set_aspect("equal")
plt.show()
```

The input is a matrix $X\in\mathbb R^{n\times p}$ whose row $x_i$ holds the $p$ transformed measurements of cell $i$, oriented like the expression matrices of Lessons 15 and 16. The rows could equally be bulk samples; clustering tumors by expression has split tumors that looked alike under the microscope into molecularly distinct groups [@holmeshuber2019, sec. 5.2]. The output is a label $c_i\in\{1,\ldots,K\}$ for each cell, or a graph recording which cells are neighbors. Both depend on two choices the instrument does not make for us: a *distance* between rows and a *criterion* for a good group. A cluster is therefore a hypothesis about a population, to be checked against biology such as a marker known to define the type.

## Cells as points: what "near" means

### Euclidean distance and the scale of each marker

Each row is a point in $\mathbb R^p$, and the first comparison to try is the Euclidean distance between cells $i$ and $i'$,
$$
d(x_i,x_{i'})=\|x_i-x_{i'}\|=\Bigl(\sum_{j=1}^p (X_{ij}-X_{i'j})^2\Bigr)^{1/2}.
$$
Every feature adds its squared difference, so a feature whose numbers span a wide range outweighs one with a narrow range, even when the range reflects the instrument rather than the biology. Suppose six cells were measured on two markers, and marker 2's dye is brighter, so its intensities come out about ten times larger (@tbl-clustering-six-cells).

| cell | marker 1 | marker 2 | marker 2 ÷ 10 |
|---|---:|---:|---:|
| A | 0 | 10 | 1 |
| B | 1 | 0 | 0 |
| C | 2 | 20 | 2 |
| D | 5 | 20 | 2 |
| E | 6 | 30 | 3 |
| F | 7 | 10 | 1 |

: Six toy cells measured on two markers; the last column rescales marker 2. {#tbl-clustering-six-cells}

Which cell is nearest to A? On the recorded values, $d(A,F)=7$ and $d(A,B)=\sqrt{1^2+10^2}\approx10.05$, so F is nearest. After dividing marker 2 by ten, $d(A,B)=\sqrt{1^2+1^2}\approx1.41$ while $d(A,F)=7$ is unchanged, so B is nearest. No cell changed; only the weight of each marker did. A weighted Euclidean distance makes that choice explicit:
$$
d_w(x_i,x_{i'})^2=\sum_{j=1}^p w_j\,(X_{ij}-X_{i'j})^2,\qquad w_j>0.
$$
The raw comparison uses $w=(1,1)$ and the rescaled one $w=(1,1/100)$. Standardizing, as in Lesson 16, sets $w_j=1/\sigma_j^2$. That is a convention rather than a neutral default: when a marker's spread comes from the difference between two populations, standardizing shrinks the very difference we hope to find. From here on, the six cells use the rescaled coordinates.

### Correlation distance compares profile shape

Weights decide how much each feature counts; a second choice decides which aspect of a profile is compared. Suppose three cells have RNA counts $u=(1,2,3)$, $v=(2,4,6)$, and $w=(3,2,1)$ for three genes. Cell $v$ has the composition of $u$ with twice as many captured molecules, the depth effect of Lesson 16, while $w$ reverses the pattern. Euclidean distance gives $d(u,v)=\sqrt{14}\approx3.74$ and $d(u,w)=\sqrt8\approx2.83$, so it calls $w$ the nearer cell.

Correlation compares shapes. Center each profile on its own mean across genes (not the gene-by-gene centering of PCA) and scale it to unit length:
$$
\tilde x=\frac{x-\bar x\,\mathbf 1}{\|x-\bar x\,\mathbf 1\|},
\qquad r(x,y)=\tilde x\cdot\tilde y,
\qquad d_r(x,y)=1-r(x,y),
$$
where $\bar x$ is the mean of the entries of $x$, $\mathbf 1$ is the all-ones vector, and $r$ is the Pearson correlation. Centering gives $(-1,0,1)$ for $u$ and $(-2,0,2)$ for $v$, which normalize to the same unit vector, so $d_r(u,v)=0$. Centering $w$ gives $(1,0,-1)$, the opposite direction, so $d_r(u,w)=2$, the largest possible value. Because the normalized profiles have unit length, $\|\tilde x-\tilde y\|^2=2-2\,\tilde x\cdot\tilde y=2(1-r)$: correlation distance is half the squared Euclidean distance between centered, normalized profiles, so it ranks pairs exactly as that Euclidean distance does.

Which answer is right depends on why $v$ is larger. Deeper sampling is technical, and correlation, like depth normalization, rightly removes it. A larger cell making more of every RNA is biology, and correlation discards it. Correlation is popular for comparing genes, which can rise and fall together from different baseline levels; for cells, Euclidean distance on normalized, log-transformed data is a frequent choice.

::: {.callout-note title="Key idea"}
Each distance in this lesson is a Euclidean distance, or a function of one, after some transformation of the rows: rescaling features, normalizing depth, or centering and normalizing profiles. Choosing a distance is choosing which measured differences count as biological.
:::

## K-means: groups around centers

Gating draws a region around each cloud. A computational version picks $K$ centers, gives each cell to its nearest center, and judges the result by how close cells sit to their own center. We write $K$ for the number of clusters and reserve lowercase $k$ for neighbors.

::: {#def-kmeans-objective}
## Within-cluster sum of squares
For labels $c_i\in\{1,\ldots,K\}$ and centers $m_1,\ldots,m_K\in\mathbb R^p$,
$$
W(c,m)=\sum_{i=1}^n\|x_i-m_{c_i}\|^2 .
$$ {#eq-clustering-within}
The *K-means problem* minimizes $W$ over labels and centers. Cluster $\ell$ is $C_\ell=\{i:c_i=\ell\}$, with $n_\ell$ members.
:::

Small $W$ means cells sit close to their own centers. $W$ cannot choose $K$ by itself, since one cluster per cell gives $W=0$. Trying all $K^n$ labelings is hopeless, but $W$ is easy to minimize over labels with the centers held fixed, or over centers with the labels held fixed.

*Best labels for fixed centers.* Cell $i$ appears in one term of @eq-clustering-within, which depends only on its own label, so the best label is the nearest center.

*Best centers for fixed labels.* The sum splits into one sum per cluster. For any set $C$ of cells with mean $\bar x_C$ and any point $m$,
$$
\sum_{i\in C}\|x_i-m\|^2=\sum_{i\in C}\|x_i-\bar x_C\|^2+|C|\,\|\bar x_C-m\|^2 .
$$ {#eq-clustering-mean-split}
To see this, write $x_i-m=(x_i-\bar x_C)+(\bar x_C-m)$ and expand the square. The cross term is $2(\bar x_C-m)\cdot\sum_{i\in C}(x_i-\bar x_C)$, which vanishes because deviations from a mean sum to zero, as the centered columns of Lesson 15 did. Only the last term depends on $m$, and it is zero exactly at $m=\bar x_C$, so the best center is the cluster mean.

::: {.algorithm}
**Algorithm: Lloyd's algorithm for K-means.**

Choose $K$ starting centers, for example $K$ distinct cells. Repeat: (1) assign each cell to its nearest center, keeping its current center if that one ties for nearest; (2) move each center to the mean of its cells. Stop when an assignment changes no label. If a cluster empties, its center stays where it was.
:::

::: {#thm-lloyd-descent}
## Lloyd's algorithm descends and stops
Neither step of Lloyd's algorithm increases $W$. The algorithm stops after finitely many steps, at labels and centers that neither step can improve.
:::

::: {.proof}
The assignment step gives each cell the smallest term available to it. By @eq-clustering-mean-split, moving center $\ell$ to its cluster mean lowers that cluster's sum by $n_\ell\|\bar x_\ell-m_\ell\|^2\geq0$. After each update the centers are the means of the current partition, so $W$ is then determined by the partition alone. A label changes only when another center is strictly closer, so a change lowers $W$ strictly, and the following update cannot raise it. Each new partition therefore has a smaller value than every earlier one, none recurs, and there are at most $K^n$ partitions.
:::

The guarantee is about descent, not about the best answer: the algorithm stops at the first partition that neither step improves, which need not minimize $W$ (@exr-clustering-local-optimum). Software therefore runs it from several starts and keeps the smallest $W$. Each iteration computes $nK$ distances in $p$ dimensions, $\mathcal O(nKp)$ operations.

::: {#exm-kmeans-six-cells}
## Two iterations on the six cells
Use the rescaled cells A $(0,1)$, B $(1,0)$, C $(2,2)$, D $(5,2)$, E $(6,3)$, F $(7,1)$ with $K=2$, and start the centers at B and C. Both lie in the left group, a poor start that the algorithm must repair.

*Assign.* For A, $\|A-B\|^2=1+1=2$ and $\|A-C\|^2=4+1=5$, so A joins B. For D, $\|D-B\|^2=16+4=20$ and $\|D-C\|^2=9$, so D joins C, as do E and F. With clusters $\{A,B\}$ and $\{C,D,E,F\}$, $W=2+0+0+9+17+26=54$.

*Update.* The means are $m_1=(0.5,0.5)$ and $m_2=\bigl(\tfrac{2+5+6+7}{4},\tfrac{2+2+3+1}{4}\bigr)=(5,2)$, and $W$ falls to 17. Center 2 has left C for the middle of D, E, and F, and C now contributes 9 of the total.

*Assign again.* For C, $\|C-m_1\|^2=1.5^2+1.5^2=4.5<9=\|C-m_2\|^2$, so C switches and $W$ falls to 12.5.

*Update again.* The centers move to $(1,1)$ and $(6,2)$, and $W=8$. The next assignment changes nothing, so the algorithm stops with $\{A,B,C\}$ and $\{D,E,F\}$, the cells with low and high marker 1.
:::

```{.python .execute}
#| label: code-clustering-kmeans
def assign(X: np.ndarray, centers: np.ndarray,
           labels: np.ndarray | None = None) -> np.ndarray:
    """Label each cell with its nearest center, keeping its label on a tie."""
    new_labels = np.zeros(len(X), dtype=int)
    for i in range(len(X)):
        squared = np.sum((centers - X[i]) ** 2, axis=1)
        best = int(np.argmin(squared))
        if labels is not None and squared[labels[i]] == squared[best]:
            best = int(labels[i])
        new_labels[i] = best
    return new_labels


def update(X: np.ndarray, labels: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Move each center to the mean of its cells; an empty cluster keeps its center."""
    new_centers = centers.astype(float).copy()
    for g in range(len(centers)):
        members = X[labels == g]
        if len(members) > 0:
            new_centers[g] = members.mean(axis=0)
    return new_centers


def kmeans(X: np.ndarray, centers: np.ndarray) -> tuple[np.ndarray, np.ndarray, list]:
    """Lloyd's algorithm; also returns W after every step."""
    labels = assign(X, centers)
    history = [np.sum((X - centers[labels]) ** 2)]
    while True:
        centers = update(X, labels, centers)
        history.append(np.sum((X - centers[labels]) ** 2))
        new_labels = assign(X, centers, labels)
        if np.array_equal(new_labels, labels):
            return labels, centers, history
        labels = new_labels
        history.append(np.sum((X - centers[labels]) ** 2))


names = ["A", "B", "C", "D", "E", "F"]
X = np.array([[0, 1], [1, 0], [2, 2], [5, 2], [6, 3], [7, 1]], dtype=float)
labels, centers, history = kmeans(X, X[[1, 2]])
print("W after each step:", [float(W) for W in history])
print("labels:", labels, " centers:", centers.tolist())
```

The values $54, 17, 12.5, 8$ match the hand calculation. In the second panel of @fig-clustering-kmeans-steps, the orange center travels most of the way to D, E, and F, which leaves C closer to the blue center.

```{.python .execute}
#| label: fig-clustering-kmeans-steps
#| code-fold: true
#| fig-width: 6
#| fig-height: 3.7
#| fig-cap: "Lloyd's algorithm on the six rescaled cells from centers at B and C. Crosses are centers, red arrows show the update step moving them, and color and shape give each cell's cluster. W falls at every step, and C changes cluster once the orange center has moved toward D, E, and F."
#| fig-alt: "Four panels of six labeled points. First: A, B blue and C to F orange, centers at B and C. Second: red arrows move the centers to the means. Third: C turns blue. Fourth: centers move to (1,1) and (6,2)."
step_labels = [assign(X, X[[1, 2]])]
step_centers = [X[[1, 2]], update(X, step_labels[0], X[[1, 2]])]
step_labels.append(assign(X, step_centers[1], step_labels[0]))
step_centers.append(update(X, step_labels[1], step_centers[-1]))
# (labels, centers, previous centers) for each panel, in the order of the steps
states = [(0, 0, None), (0, 1, 0), (1, 1, None), (1, 2, 1)]
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
for ax, (a, c, p) in zip(axes.flat, states):
    previous = None if p is None else step_centers[p]
    cf.cluster_step(ax, X, step_labels[a], step_centers[c], names=names,
                    previous_centers=previous)
    W = np.sum((X - step_centers[c][step_labels[a]]) ** 2)
    ax.set_title(f"W = {W:g}", loc="left")
for ax in axes[1]:
    ax.set_xlabel("marker 1")
for ax in axes[:, 0]:
    ax.set_ylabel("marker 2 ÷ 10")
plt.show()
```

The last panel also shows what kind of group K-means can report. A cell $x$ goes to $m_1$ rather than $m_2$ when $\|x-m_1\|^2<\|x-m_2\|^2$; expanding, the $\|x\|^2$ terms cancel and leave
$$
2\,x\cdot(m_2-m_1)<\|m_2\|^2-\|m_1\|^2 .
$$
The condition is linear in $x$, so the boundary between two clusters is the line (in $p$ dimensions, the hyperplane) perpendicular to $m_2-m_1$ through the midpoint of the centers. Each cluster's region is an intersection of half-planes, a convex polygon around its center. K-means therefore works best when the true groups are compact, round, and of similar size [@holmeshuber2019, sec. 5.4.1].

## Hierarchical clustering revisits UPGMA

K-means needs $K$ in advance. UPGMA in Lesson 13 instead merged the two closest clusters again and again, averaging distances over all cross pairs. Run on cells with any distance, the same algorithm is *average-linkage hierarchical clustering*, and without a molecular clock its merge heights are average distances, not ages. On the six rescaled cells, the pairs A,B and D,E merge at $\sqrt2$, C and F join them at $\sqrt5$, and the final merge comes at about $5.23$; cutting between $\sqrt5$ and $5.23$ returns the K-means clusters. Every dataset yields a tree, so the cut height is our choice, not a discovery [@holmeshuber2019, sec. 5.6].

## Neighborhood graphs

K-means and average linkage summarize each group by averages over all its cells. Cells changing gradually between states may form a curved chain rather than a compact cloud, so one group center can describe them poorly. Keeping each cell's nearest neighbors lets a group follow that chain. A neighborhood graph keeps only that.

::: {#def-knn-graph}
## k-nearest-neighbor graph
For cell $i$, let $N_k(i)$ be the $k$ other cells nearest to $x_i$, breaking ties by row order. The *union $k$-nearest-neighbor graph* has one vertex per cell and an undirected edge $\{i,i'\}$ whenever $i'\in N_k(i)$ or $i\in N_k(i')$.
:::

Cell $i$ can select $i'$ without being selected in return. The union convention keeps an edge if either cell selects the other; the *mutual* convention requires both and gives a sparser graph, so an analysis must say which it uses. Under the union convention each cell has at least $k$ edges, and the graph has at most $nk$.

::: {.algorithm}
**Algorithm: construct a union k-nearest-neighbor graph.**

For each cell, compute its distances to all other cells, sort them, and select the nearest $k$, breaking ties by row order. Add an undirected edge for each selected pair, once even if both ends select it.
:::

::: {#exm-knn-six-cells}
## Neighbors of the six cells for $k=2$
On the rescaled coordinates, A's distances to B, C, D, E, F are $1.41$, $2.24$, $5.10$, $6.32$, and $7.00$, so $N_2(A)=\{B,C\}$. B selects A and C; C, at distance $\sqrt5$ from both, selects A and B; and D, E, F select one another. The graph is two separate triangles.

On the raw intensities, A's nearest cells are F at distance $7$ and B at $10.05$. B selects A and F, C selects D and A, D selects C and E, E selects D and C, and F selects A and D. The eight edges $AB, AC, AF, BF, CD, CE, DE, DF$ form triangles $\{A,B,F\}$ and $\{C,D,E\}$ joined by $AC$ and $DF$: this graph groups cells by marker 2.
:::

```{.python .execute}
#| label: code-clustering-knn
def knn_edges(X: np.ndarray, k: int) -> list[tuple[int, int]]:
    """Union k-nearest-neighbor graph: join i and j if either selects the other."""
    edges = set()
    for i in range(len(X)):
        distances = np.linalg.norm(X - X[i], axis=1)
        distances[i] = np.inf  # a cell is not its own neighbor
        for j in np.argsort(distances, kind="stable")[:k]:  # ties by row order
            edges.add((min(i, int(j)), max(i, int(j))))
    return sorted(edges)


raw_edges = knn_edges(X * np.array([1, 10]), k=2)  # undo the rescaling
rescaled_edges = knn_edges(X, k=2)
print("raw:     ", [names[i] + names[j] for i, j in raw_edges])
print("rescaled:", [names[i] + names[j] for i, j in rescaled_edges])
```

@fig-clustering-scaled-neighbors draws both graphs at the rescaled positions, so only the distance used to choose neighbors differs. The raw graph's edges run across the gap in marker 1, joining cells at similar heights.

```{.python .execute}
#| label: fig-clustering-scaled-neighbors
#| code-fold: true
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "Union two-nearest-neighbor graphs of the six cells, both drawn at the rescaled coordinates, with neighbors chosen from raw intensities (left) or after dividing marker 2 by ten (right). Red marks the highest-modularity community containing A, which changes with the marker scale."
#| fig-alt: "Two panels with the same six labeled nodes. Left: edges cross between the left and right groups, and A, B, F are red. Right: two separate triangles, with A, B, C red."
fig, axes = plt.subplots(1, 2, sharey=True)
panels = [(axes[0], raw_edges, [0, 1, 5], "Neighbors from raw intensities"),
          (axes[1], rescaled_edges, [0, 1, 2], "Neighbors after marker 2 ÷ 10")]
for ax, edges, community, title in panels:
    cf.neighbor_graph(ax, X, edges, names=names, highlight=community)
    ax.set(title=title, xlabel="marker 1")
axes[0].set_ylabel("marker 2 ÷ 10")
plt.show()
```

### Communities in a graph

In a graph, a group of cells is a *community*: vertices with many edges among themselves and few to the rest. The rescaled graph's two triangles share no edge, so its communities are plain. In the raw graph every split leaves some edges crossing, and we need a number to compare candidate partitions.

::: {#def-modularity}
## Modularity
Let a graph have $m$ edges, and let $\deg(i)$ count the edges at vertex $i$. For a partition into communities $\ell=1,\ldots,L$, let $e_\ell$ be the number of edges with both ends in community $\ell$ and $D_\ell=\sum_{i\in\ell}\deg(i)$. The *modularity* of the partition is
$$
Q=\sum_{\ell=1}^L\left[\frac{e_\ell}{m}-\left(\frac{D_\ell}{2m}\right)^2\right].
$$
:::

The first term is the fraction of edges inside community $\ell$. The second is the fraction expected there if edges were rewired at random while every cell kept its number of edges: a fraction $D_\ell/2m$ of the $2m$ edge ends belong to community $\ell$, and an edge lies inside when both its ends do, with probability about $(D_\ell/2m)^2$. Modularity measures how many more edges lie inside the communities than random wiring would put there.

In the raw graph, $m=8$. Grouping by marker 2 as $\{A,B,F\}$ and $\{C,D,E\}$ gives two triangles, each with $e_\ell=3$ and $D_\ell=8$, so
$$
Q=2\left[\frac38-\left(\frac{8}{16}\right)^2\right]=\frac14 .
$$
Grouping by marker 1 as $\{A,B,C\}$ and $\{D,E,F\}$ keeps only $AB, AC$ and $DE, DF$ inside, so $Q=2(2/8-1/4)=0$, no better than random wiring. In the rescaled graph the marker 1 grouping has $Q=2(3/6-1/4)=1/2$. The marker 2 grouping is the best two-group partition of the raw graph and the marker 1 grouping is the best of the rescaled graph, so the distance chosen before the graph was built decided which grouping won. For thousands of cells, single-cell software searches for high-modularity partitions with greedy algorithms such as Louvain and Leiden, which move cells between communities while $Q$ increases.

## Why cells are clustered on graphs

Single-cell analyses usually cluster a neighbor graph rather than run K-means directly. Three properties of cell data explain the choice.

*Shape.* A K-means cluster is a convex region around a center, while a graph community grows along chains of neighbors. In @fig-clustering-shape, simulated cells spread along a curve, as cells partway through a gradual change of state might be, bend around a compact population. K-means with $K=2$ must draw a straight boundary through them.

```{.python .execute}
#| label: fig-clustering-shape
#| code-fold: true
#| fig-width: 6
#| fig-height: 2.9
#| fig-cap: "Simulated cells along a curve around a compact population. Left: the best of ten K-means runs with K = 2; its straight boundary mixes the populations. Right: the union ten-nearest-neighbor graph, with the highest-modularity community containing the compact cells in red. Local neighbors follow a shape that center-based clusters cannot."
#| fig-alt: "Two panels of the same 90 points, an arc of 60 above a round cloud of 30. Left: blue and orange groups split by a diagonal through both. Right: grey edges join nearby points, and the round cloud is red."
rng = np.random.default_rng(3)
angles = rng.uniform(np.pi / 2 - 1.2, np.pi / 2 + 1.2, 60)
curve = 3 * np.column_stack([np.cos(angles), np.sin(angles)])
curve += 0.3 * rng.standard_normal((60, 2))
compact = np.array([0, 1.05]) + 0.45 * rng.standard_normal((30, 2))
cells = np.vstack([curve, compact])
is_compact = np.arange(90) >= 60

best_W = np.inf
for seed in range(10):  # several random starts; keep the smallest W
    starts = np.random.default_rng(seed).choice(90, 2, replace=False)
    run_labels, run_centers, run_history = kmeans(cells, cells[starts])
    if run_history[-1] < best_W:
        best_W, km_labels, km_centers = run_history[-1], run_labels, run_centers
shape_edges = knn_edges(cells, k=10)
communities = nx.community.louvain_communities(nx.Graph(shape_edges), seed=0)
print("K-means (curve, compact) per cluster:",
      [(int(np.sum(~is_compact & (km_labels == g))),
        int(np.sum(is_compact & (km_labels == g)))) for g in range(2)])
print("graph (curve, compact) per community:",
      [(int(np.sum(~is_compact[list(c)])), int(np.sum(is_compact[list(c)])))
       for c in communities])

compact_community = max(communities, key=lambda c: np.sum(is_compact[list(c)]))
fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
cf.cluster_step(axes[0], cells, km_labels, km_centers)
cf.neighbor_graph(axes[1], cells, shape_edges, highlight=compact_community)
axes[0].set(title="K-means, K = 2", xlabel="feature 1", ylabel="feature 2")
axes[1].set(title="10-nearest-neighbor graph", xlabel="feature 1")
plt.show()
```

K-means puts 27 of the 30 compact cells in a cluster with 27 curve cells, while every graph community lies within one population. The curve itself is split into four communities: modularity has no notion of a continuous process, so a gradual change of state still gets cut into pieces.

*Scale.* For $n=10^5$ cells, average linkage needs all $n(n-1)/2\approx5\times10^9$ pairwise distances, while a $k$-nearest-neighbor graph stores at most $nk$ edges, $10^6$ for $k=10$. Large analyses find neighbors approximately, without comparing every pair, and run community search on the sparse graph.

*Dimension.* With thousands of noisy genes, distances between cells tend to look alike, so noise can decide which cell is nearest. Single-cell workflows therefore reduce the matrix with PCA first, keeping 10 to 50 component scores in which coordinated variation across many genes is concentrated, and build the neighbor graph on those scores. That graph and its communities are where Lesson 20 begins, when it asks how to draw the cells in two dimensions.

:::: {.callout-note collapse="true" title="Going deeper: why distances concentrate"}
A cell measured on thousands of genes accumulates noise in every coordinate. Suppose two cells differ only by independent noise, each coordinate a standard normal value. Each coordinate difference then has variance 2, and its square has mean 2 and variance 8, so the squared distance over $p$ coordinates has mean $2p$ and standard deviation $\sqrt{8p}$. The relative spread
$$
\frac{\sqrt{8p}}{2p}=\sqrt{\frac2p}
$$
is $1$ for $p=2$ but $0.01$ for $p=20{,}000$: all distances look alike, and noise largely decides which cell is "nearest."
::::

## What a cluster supports

Every method in this lesson returns groups whether or not the cells form groups. K-means applied to cells spread evenly over a square still reports $K$ clusters, and a neighbor graph of such cells still has a highest-modularity partition [@holmeshuber2019, sec. 5.7]. The curve in @fig-clustering-shape shows the opposite risk, one continuous process cut into apparent types.

A cluster becomes evidence when it survives reasonable changes to the analysis and agrees with something the clustering did not use. Rerunning with other starts, other $k$, other distances, or subsets of cells tests stability; a marker known from biology, such as CD4 for helper T cells, or an independent measurement tests meaning. A label like "helper T cells" is a claim layered on the representation, the distance, and the criterion. It is not part of the data.

## Exercises

::: {#exr-clustering-knn-components}
Five cells have PCA scores A $(0,0)$, B $(2,0)$, C $(0,3)$, D $(7,7)$, E $(9,8)$. Build the union $k$-nearest-neighbor graph for $k=1$ and $k=2$, and count the connected components of each. What decided whether these cells form one group or two?
:::

::: {#exr-clustering-correlation}
Construct three cells measured on three genes for which Euclidean and correlation distance disagree about which cell is nearest the first. Describe a biological situation in which each distance gives the better answer.
:::

::: {#exr-clustering-local-optimum}
Four cells lie at $(0,0)$, $(0,2)$, $(6,0)$, and $(6,2)$, with $K=2$. Show that the left/right and bottom/top partitions, each with centers at its means, are both stopping points of Lloyd's algorithm. Compute $W$ for each, and give starting centers that lead to the worse one.
:::

::: {#exr-clustering-choose-k}
Let $W_K^*$ be the smallest $W$ achievable with $K$ clusters. Show that $W_{K+1}^*\leq W_K^*$ for $K<n$ and that $W_n^*=0$. Why can minimizing $W$ not choose the number of cell populations?
:::

::: {#exr-clustering-single-linkage}
Single linkage measures the distance between clusters by the smallest distance between their members. Show that the single-linkage clusters formed by merges at heights at most $r$ are the connected components of the graph joining two cells whenever their distance is at most $r$.
:::

::: {#refs}
**References**
:::
