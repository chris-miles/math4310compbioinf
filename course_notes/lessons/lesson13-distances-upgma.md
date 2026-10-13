---
published: true
title: "Evolutionary Distances and UPGMA"
subtitle: "Lesson 13 · 2027-03-03 · Week 8"
nocite: |
  @compeau2015, @durbin1998
---

## Core question

Which samples share a more recent ancestor?

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

## From sampled sequences to shared ancestry

Suppose we have DNA from several organisms or viral samples and want to know which are more closely related. Their common ancestors are usually unavailable to sequence. We instead compare the descendants: sequencing gives a DNA letter string from each sample, and substitutions accumulated since their ancestors can leave different letters at corresponding positions. A *phylogeny*, or evolutionary tree, describes the inferred pattern of shared ancestry. Its branches can also represent the amount of sequence change along each lineage.

During the 2003 SARS investigation, researchers compared coronavirus sequences sampled from humans and other animals to investigate their evolutionary relationships [@compeau2015]. Sequence relatedness can help assess an origin hypothesis. A tree of sampled sequences alone does not identify a direct transmission event or prove which animal infected a particular person. Here we develop the distance-based inference using synthetic DNA examples.

The sequences must be *homologous*: descended from a corresponding sequence in a common ancestor. Our input is a supplied *multiple alignment*, with one sequence per row and corresponding positions in the same columns. For example:

```text
A  ACGTACGT
B  ACGTTCGT
C  ATGTTCGA
```

A and B differ at one of eight columns; A and C differ at three. Alignment matters because comparing unrelated positions would turn errors in correspondence into apparent evolutionary changes. We will use only columns containing a known DNA base in every row. This common set of columns keeps all pairwise comparisons on the same footing; choosing or computing the alignment is a separate task.

Unlike the HMM problems, we now infer relationships among samples rather than hidden labels along one sequence. We first summarize their differences in a matrix, then ask what evolutionary tree could explain that matrix.

## Observed differences and evolutionary distance

For $L$ retained columns, define the mismatch fraction
$$
f_{ij}=\frac{1}{L}\sum_{\ell=1}^L
\mathbf{1}\{x^{(i)}_\ell\ne x^{(j)}_\ell\}.
$$
The indicator contributes one at a mismatch and zero otherwise. In the alignment above, $f_{AB}=1/8$, $f_{AC}=3/8$, and $f_{BC}=2/8$. The matrix is symmetric, has a zero diagonal, and is nonnegative. Counting mismatches on a common set of columns also gives the triangle inequality $f_{ij}\leq f_{ik}+f_{kj}$: if $i$ and $j$ differ at a column, at least one must differ from $k$ there. Distinct samples can nevertheless have identical retained sequences and distance zero.

Mismatch fractions describe the observed letters, not the full history of substitutions. Along a lineage, A could change to G and later back to A, leaving no mismatch. A change to G followed by T leaves one mismatch after two substitutions. This loss of historical information becomes substantial for divergent sequences [@durbin1998].

An *evolutionary distance* estimates the expected number of substitutions per site along the path between two samples. Unlike a mismatch fraction, it can exceed one. We will supply these distances to the tree algorithm; they may come from a fitted substitution model rather than raw mismatch counts.

### A correction for hidden substitutions
The Jukes–Cantor model assumes equal base frequencies and equal rates for changes between each pair of different bases. Its predicted mismatch fraction after total path distance $d$ is
$$
f=\frac34\left(1-e^{-4d/3}\right).
$$
Taking this model relation as given, rearranging and taking the natural logarithm gives
$$
d=-\frac34\ln\left(1-\frac43f\right),\qquad 0\leq f<\frac34.
$$
For $f=0.1$, this gives $d\approx0.1073$, slightly more than the observed mismatch fraction. As the lineages diverge, their bases approach independent uniform draws, so the probability of a mismatch approaches $3/4$. The observed difference then *saturates*: large additional evolutionary distances produce small changes in $f$. The correction diverges as $f$ approaches $3/4$ and is not a finite real estimate at or above that value. Near saturation, sampling error in $f$ becomes large uncertainty in $d$ [@durbin1998].

A corrected matrix estimated from finite sequences need not fit any tree exactly: correcting for hidden substitutions still leaves sampling noise and model assumptions.

```{.python .execute}
#| label: fig-upgma-saturation
#| fig-cap: "Expected mismatch fraction under the Jukes-Cantor model. Additional substitutions become harder to detect as the curve approaches three quarters, so observed mismatches increasingly underestimate evolutionary distance."
#| fig-alt: "A rising concave curve of mismatch fraction against substitutions per site approaches the horizontal level 0.75; a dashed no-hidden-changes line rises diagonally."
distance = np.linspace(0, 3, 301)
mismatch = 0.75 * (1 - np.exp(-4 * distance / 3))
fig, ax = plt.subplots()
ax.plot(distance, mismatch, label="Expected mismatch fraction")
ax.plot([0, 1], [0, 1], "--", color=cf.RULE, label="Every change remains visible")
ax.axhline(0.75, color=cf.RULE, linewidth=0.8, linestyle=":")
ax.set(xlabel="Evolutionary distance (substitutions per site)",
       ylabel="Mismatch fraction", ylim=(0, 1))
ax.legend()
plt.show()
```

## What the tree must explain

We need to distinguish a branching pattern from the distances it predicts. A rooted tree has a distinguished ancestral node, the *root*. Its *leaves* or *tips* are the sampled sequences; internal nodes represent inferred common ancestors. Each edge $e$ has a nonnegative branch length $\ell_e$. For the unique path $P_T(i,j)$ between two tips in tree $T$, define
$$
d_T(i,j)=\sum_{e\in P_T(i,j)}\ell_e.
$$
These are *additive tree distances*. Given an $n\times n$ matrix $D=(d_{ij})$, our exact reconstruction question is to find a tree with $d_T(i,j)=d_{ij}$ for every pair. Arbitrary matrices need not admit such a tree; we will establish a guarantee under a stronger assumption.

Pairwise distances alone generally do not locate the root. UPGMA adds a *molecular-clock assumption*: all lineages accumulate substitutions at the same rate, and the tips are sampled at the same time. Under this model every root-to-tip path has the same length [@durbin1998]. A node's height $h_u$ is its distance down to any descendant tip. Thus tips have height zero and an edge from parent $u$ to child $v$ has length $\ell_{uv}=h_u-h_v$.

If $u$ is the most recent common ancestor of tips $i,j$, their paths to $u$ have equal length. Therefore
$$
d_{ij}=2h_u.
$$
This equation explains both the factor of two in UPGMA and its biological interpretation: a smaller distance means a more recent common ancestor *under the clock assumption*. Heights are measured in substitutions per site; converting them to calendar ages requires a calibrated substitution rate.

::: {#def-ultrametric}
## Ultrametric distances
A symmetric, nonnegative matrix with zero diagonal is **ultrametric** if
$$
d_{ij}\leq\max(d_{ik},d_{jk})
\quad\text{for every }i,j,k.
$$
Equivalently, among the three distances for any three tips, the largest two are equal. All three may be equal.
:::

A clock tree has this property. For three tips, either all share the same most recent common ancestor, or two meet below the ancestor they share with the third. In the latter case, the two distances to the third tip are equal and larger. This gives an observable constraint on the matrix, stronger than the ordinary triangle inequality. Real estimated distances rarely satisfy it exactly; UPGMA still returns a clock tree, but exact recovery then has no guarantee.

## UPGMA: merge the closest clusters

Start with one cluster per tip. After a merge, a cluster represents an inferred ancestral group. How should we measure its separation from another group? UPGMA uses the average over all original cross-cluster tip pairs:
$$
d(C,D)=\frac{1}{|C||D|}\sum_{i\in C}\sum_{j\in D}d_{ij}.
$$
Every tip pair has equal weight. The method's name is *unweighted pair group method with arithmetic mean*; equal weighting of tips requires weighting cluster updates by their sizes.

When $C$ and $D$ merge, partitioning the sum into pairs originating in $C$ and pairs originating in $D$ gives
$$
d(C\cup D,E)
=\frac{|C|d(C,E)+|D|d(D,E)}{|C|+|D|}.
$$
This avoids revisiting every original pair. A simple average of the two old cluster distances is correct only when their sizes are equal, or when those distances happen to coincide.

::: {.algorithm}
**Algorithm: UPGMA.**
Initialize singleton clusters, each with size one and height zero.

1. Choose two distinct active clusters $C,D$ with smallest distance $\delta$.
2. Create their parent at height $\delta/2$. Its branches to their roots have lengths $\delta/2-h_C$ and $\delta/2-h_D$.
3. Replace $C,D$ by $C\cup D$, whose size is $|C|+|D|$. Update its distances to the other clusters by the size-weighted formula.
4. Repeat until one cluster remains; its node is the root. Break ties by a stated rule.
:::

The merge distances never decrease: every new distance is a weighted average of distances at least as large as the current minimum, and retained distances also exceed or equal it. Parent heights therefore never lie below child heights, so the algorithm produces nonnegative branches even when its input is not ultrametric.

## Four samples, three merges

Suppose a distance estimator supplies the following synthetic matrix. For convenient hand arithmetic, entries are in units of $0.01$ substitutions per site; an entry 4 means a distance of 0.04.

| | A | B | C | D |
|---|---:|---:|---:|---:|
| A | 0 | 4 | 8 | 16 |
| B | 4 | 0 | 8 | 16 |
| C | 8 | 8 | 0 | 16 |
| D | 16 | 16 | 16 | 0 |

The closest pair is A,B at distance 4. Their parent U has height 2, and both terminal branches have length 2. The new distances are
$$
d(AB,C)=\frac{8+8}{2}=8,
\qquad d(AB,D)=\frac{16+16}{2}=16.
$$
The next merge is AB with C at distance 8. Its parent V has height 4; the branch V→U has length $4-2=2$, while V→C has length 4. Now the clusters have unequal sizes:
$$
d(ABC,D)=\frac{2\,d(AB,D)+1\,d(C,D)}{3}
=\frac{2(16)+16}{3}=16.
$$
The final root R has height 8. Its branch to V has length 4 and its branch to D length 8. Every root-to-tip distance is 8. For instance, the path A→U→V→C has length $2+2+4=8$, reproducing $d_{AC}$.

Under the clock model, A and B share the most recent ancestor, C joins their lineage earlier, and D separates earlier still. A is not inferred to be the ancestor of B: both are observed tips descending from U. The heights here are 0.02, 0.04, and 0.08 substitutions per site, not dates.

The equal distances to D conceal the numerical effect of the weights. If instead $d(AB,D)=16$ and $d(C,D)=20$, the update would be $(2\cdot16+20)/3=17\frac13$, not 18. The two-member cluster contributes twice as many cross-cluster pairs.

## Computing and drawing the result

The implementation stores a size and height for each active cluster and a dictionary of distances between pairs of cluster identifiers. The returned nested tree records each node's name, children, and branch length from its parent. It also returns the merge history so we can compare execution with the hand calculation.

```{.python .execute}
#| label: code-upgma
def upgma(labels: list[str], matrix: np.ndarray) -> tuple[dict, list[tuple]]:
    """Return an average-linkage clock tree and its merge history."""
    matrix = np.asarray(matrix, dtype=float)
    n = len(labels)
    if n == 0 or len(set(labels)) != n or matrix.shape != (n, n):
        raise ValueError("Use distinct labels and a matching nonempty square matrix.")
    if (not np.all(np.isfinite(matrix)) or np.any(matrix < 0)
            or not np.allclose(matrix, matrix.T)
            or not np.allclose(np.diag(matrix), 0)):
        raise ValueError("Distances must be finite, symmetric, nonnegative, and zero-diagonal.")

    active = list(range(n))
    sizes = dict.fromkeys(active, 1)
    heights = dict.fromkeys(active, 0.0)
    members = {i: (labels[i],) for i in active}
    nodes = {i: {"name": labels[i], "length": 0.0, "children": []} for i in active}
    distances = {}
    for i in range(n):
        for j in range(i + 1, n):
            distances[i, j] = matrix[i, j]
    history = []
    next_id = n

    while len(active) > 1:
        # Sorted identifiers make the first minimizing pair the tie rule.
        pairs = []
        for position, i in enumerate(active):
            for j in active[position + 1:]:
                pairs.append((i, j))
        i, j = min(pairs, key=lambda pair: distances[pair])
        delta = distances[i, j]
        height = delta / 2
        history.append((members[i], members[j], float(delta), float(height)))

        nodes[i]["length"] = height - heights[i]
        nodes[j]["length"] = height - heights[j]
        nodes[next_id] = {"name": "", "length": 0.0,
                          "children": [nodes[i], nodes[j]]}
        sizes[next_id] = sizes[i] + sizes[j]
        heights[next_id] = height
        members[next_id] = members[i] + members[j]
        remaining = [k for k in active if k not in (i, j)]
        for k in remaining:
            distance_i = distances[tuple(sorted((i, k)))]
            distance_j = distances[tuple(sorted((j, k)))]
            distances[k, next_id] = (
                sizes[i] * distance_i + sizes[j] * distance_j
            ) / sizes[next_id]
        active = remaining + [next_id]
        next_id += 1

    return nodes[active[0]], history


labels = ["A", "B", "C", "D"]
distances = np.array([
    [0, 4, 8, 16],
    [4, 0, 8, 16],
    [8, 8, 0, 16],
    [16, 16, 16, 0],
])
clock_tree, merge_history = upgma(labels, distances)
assert [step[2] for step in merge_history] == [4, 8, 16]
assert [step[3] for step in merge_history] == [2, 4, 8]
merge_history
```

```{.python .execute}
#| label: fig-upgma-clock-tree
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "UPGMA tree for the four-sample matrix. Parentheses give node heights in units of 0.01 substitutions per site; all tips have height zero, and internal branches carry their lengths. The aligned tips show the clock constraint imposed on the inferred ancestry."
#| fig-alt: "Root at the left branches to D and a cluster containing C and the sister tips A and B; all four tips have the same horizontal coordinate."
fig, ax = plt.subplots()
# Name the internal nodes of this worked instance, leaving the algorithm general.
clock_tree["name"] = "R (8)"
clock_tree["children"][1]["name"] = "V (4)"
clock_tree["children"][1]["children"][1]["name"] = "U (2)"
cf.tree(ax, clock_tree, show_lengths=True)
plt.show()
```

The figure's horizontal distances encode branch lengths, not merely the order of merges. Rotating the children of an internal node changes the drawing's tip order but leaves the inferred ancestry and all path distances unchanged.

## Why exact ultrametric input is recovered

Under a clock, the closest tips meet at the youngest ancestral height because their distance is twice that height. Merging them therefore recovers the most recent split. Once this group is represented by a cluster, every member has the same distance to a tip outside it: all those paths meet at the same older ancestor. Averaging preserves that distance, so the same reasoning applies at the next merge. The proof makes this repeated step precise.

::: {#thm-upgma-recovery}
## Exact distance recovery
Given an ultrametric distance matrix, UPGMA constructs a rooted tree whose tip-to-tip path distances equal every input distance. Equal-height merges may produce different binary refinements; after contracting zero-length internal edges, these refinements describe the same hierarchy of positive-length branches.
:::

::: {.proof}
Maintain the following invariant: for any two active clusters $C,D$, every original distance $d_{ij}$ with $i\in C,j\in D$ equals the stored cluster distance $d(C,D)$. This holds initially for singleton clusters.

Suppose $C,D$ are selected at minimum distance $\delta$, and let $E$ be another active cluster. Choose tips $i\in C,j\in D,k\in E$. By the invariant,
$$
d_{ij}=\delta,\qquad d_{ik}=d(C,E),\qquad d_{jk}=d(D,E).
$$
The last two distances are at least $\delta$. They must be equal: if one were larger, it would be the unique largest distance in this triple, contradicting ultrametricity. Hence $d(C,E)=d(D,E)$, and their weighted average equals that same value. All distances from $C\cup D$ to $E$ are therefore equal to the new stored distance. This proves the invariant by induction.

When two tips first enter the same cluster, they are in opposite children of a merge. Their original distance is that merge's $\delta$. The new node has height $\delta/2$; branch lengths telescope from it to each tip to give $\delta/2$. Their tree path distance is thus $\delta$, as required.

To understand ties, fix a threshold $r\geq0$. The relation $i\sim j$ when $d_{ij}\leq r$ is transitive by the ultrametric inequality, so it partitions tips into groups. These nested groups determine the hierarchy. Tied merges resolve a group at the same height in different orders, adding only zero-length internal edges; they do not change the recovered distances.
:::

This is a conditional guarantee. It does not say that noisy sequence estimates are ultrametric, or that a biological lineage evolves at a constant rate. A useful diagnostic compares the supplied distances with path distances in the output; a visibly aligned set of tips in the drawing only confirms the algorithm's imposed constraint.

### Computational cost

With $k$ active clusters, searching all pairs costs $\mathcal{O}(k^2)$ and updating distances costs $\mathcal{O}(k)$. There are $n-1$ merges, so this straightforward implementation takes
$$
\sum_{k=2}^n\mathcal{O}(k^2)=\mathcal{O}(n^3)
$$
time and $\mathcal{O}(n^2)$ space, assuming constant-time arithmetic and dictionary lookups. The code keeps old distance entries, but only $\mathcal{O}(n^2)$ pairs are ever created. Constructing the original mismatch matrix from $n$ aligned sequences of length $L$ costs $\mathcal{O}(n^2L)$; that data-processing cost is separate from tree construction.

## When evolutionary rate changes the answer

A small distance can reflect a recent ancestor, a slow rate, or both. Suppose the true tree separates the pair A,B from the pair C,D. Give A and C terminal branches of length 1, B and D terminal branches of length 8, and the edge connecting the pairs length 1, all in the same hundredths-of-a-substitution units. Summing along the paths gives

| | A | B | C | D |
|---|---:|---:|---:|---:|
| A | 0 | 9 | 3 | 10 |
| B | 9 | 0 | 10 | 17 |
| C | 3 | 10 | 0 | 9 |
| D | 10 | 17 | 9 | 0 |

A and C have the smallest distance, 3, even though they belong to different pairs in the true tree. UPGMA joins these slowly changing tips first. This changes the branching pattern, not merely its branch lengths.

```{.python .execute}
#| label: code-upgma-clock-failure
unequal_rate_distances = np.array([
    [0, 9, 3, 10],
    [9, 0, 10, 17],
    [3, 10, 0, 9],
    [10, 17, 9, 0],
])
wrong_tree, wrong_merges = upgma(labels, unequal_rate_distances)
wrong_merges
```

After merging A,C, their distances to B and D are both 9.5. The code resolves that tie by joining B next; the last merge occurs at distance 12. The inferred A-to-C path has length 3, but the inferred A-to-B path has length 9.5 instead of 9. More sequence data would not repair this example: the input distances were already exact; the failed assumption was the clock.

```{.python .execute}
#| label: fig-upgma-rate-failure
#| fig-width: 6
#| fig-height: 6
#| fig-cap: "An unequal-rate tree and UPGMA's reconstruction from its exact distances. The true tree pairs A with B and C with D; UPGMA instead joins the slowly changing A and C. Horizontal extents show branch lengths in both panels."
#| fig-alt: "Two trees with root drawn at the left. In the true tree A and C have short terminal branches and B and D long ones. UPGMA pairs A and C and places all tips at the same horizontal coordinate."
true_tree = {"name": "", "length": 0, "children": [
    {"name": "", "length": 0.5, "children": [
        {"name": "A", "length": 1}, {"name": "B", "length": 8}]},
    {"name": "", "length": 0.5, "children": [
        {"name": "C", "length": 1}, {"name": "D", "length": 8}]},
]}
fig, axes = plt.subplots(2, 1)
cf.tree(axes[0], true_tree)
cf.tree(axes[1], wrong_tree)
axes[0].set_title("True tree: root drawn at the middle of the connecting edge")
axes[1].set_title("UPGMA: closest tips are joined")
for ax in axes:
    ax.set_xlim(-0.5, 10)
plt.show()
```

The true distances still contain information about the branching pattern. For example, $d_{AB}+d_{CD}=18$, while both $d_{AC}+d_{BD}$ and $d_{AD}+d_{BC}$ equal 20: the latter paths cross the connecting edge twice. Methods for additive distances can use such relationships without imposing a clock; neighbor joining is one example and normally returns an unrooted tree [@durbin1998].

A sequence tree describes the history of the sampled region; extending its conclusions to a whole genome requires further assumptions. Distance summaries discard which columns supplied the differences. The next lesson keeps those characters and asks how many changes a candidate tree requires.

## Exercises

::: {#exr-upgma-mismatches}
Compute the mismatch matrix for the eight-column alignment near the start. Check each triangle inequality, and decide whether its distances are ultrametric.
:::

::: {#exr-upgma-path-checks}
Trace UPGMA on the four-sample matrix. List all six edge lengths and verify all six tip-to-tip distances. If the substitution rate is 0.01 substitutions per site per million years, what calendar age corresponds to the inferred root height?
:::

::: {#exr-upgma-weights}
Keep $d_{AD}=d_{BD}=16$ but change $d_{CD}$ to 20 in the four-sample example. Trace the updates, compare size-weighted averaging with an unweighted average of clusters, and identify an input distance the final tree cannot reproduce.
:::

::: {#exr-upgma-ties}
For three tips with all off-diagonal distances equal to 6, draw two possible binary outputs. Which branches have length zero? What information about branching order is missing from the matrix?
:::

::: {#exr-upgma-rate-change}
In the unequal-rate example, change both long terminal branches, to B and D, from 8 to $b\geq0$. Find the range of $b$ for which A,C is the unique first merge. Explain why knowing the true distances perfectly is not enough to justify that merge biologically.
:::

::: {#refs}
**References**
:::
