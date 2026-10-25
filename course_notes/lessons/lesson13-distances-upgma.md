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

Lessons 9 to 12 labeled regions within one genome. We now compare DNA from several organisms or viral samples to ask which are more closely related. Their common ancestors are usually unavailable to sequence. We instead compare the descendants: sequencing gives a DNA letter string from each sample, and substitutions accumulated since their ancestors can leave different letters at corresponding positions. A *phylogeny*, or evolutionary tree, describes the inferred pattern of shared ancestry. Its branches can also represent the amount of sequence change along each lineage.

During the 2003 SARS outbreak, researchers asked which animal the new coronavirus had come from. An early guess was birds, by analogy with bird flu, but the viral genome sequenced in April 2003 ruled them out. Coronavirus genomes differ in gene content and order, so instead of aligning whole genomes, researchers aligned a single gene, encoding the Spike protein the virus uses to attach to host cells, across coronaviruses from patients and from animals; viruses from palm civets were the closest relatives of the human ones [@compeau2015]. Such a comparison can only rank the viruses that have been sampled, and later sampling of bats pointed to civets as an intermediate rather than the original host. Here we develop the distance-based inference on synthetic DNA.

The sequences must be *homologous*: descended from a corresponding sequence in a common ancestor. Our input is a supplied *multiple alignment*, with one sequence per row and corresponding positions in the same columns. @fig-upgma-data shows a made-up alignment of four samples beside the kind of tree we want to infer from it.

```{.python .execute}
#| label: fig-upgma-data
#| code-fold: true
#| code-summary: "Draw the alignment and the tree it suggests"
#| fig-width: 6
#| fig-height: 2
#| fig-cap: "A made-up alignment of four sampled sequences (left), with shaded columns where the rows differ, and the rooted tree these differences suggest (right). Tips are the samples, internal nodes are unobserved ancestors, and branch extents are drawing coordinates. A and B, which differ at one column, share the most recent ancestor."
#| fig-alt: "Four rows of eight DNA letters labeled A to D, with five shaded columns. To the right, a tree with root at the left: D branches off first, then C, then A and B as sister tips."
samples = {"A": "ACGTACGT", "B": "ACGTTCGT", "C": "ATGTTCGA", "D": "ACGACCAA"}
target_tree = {"name": "", "length": 0, "children": [
    {"name": "", "length": 1, "children": [
        {"name": "", "length": 1, "children": [
            {"name": "A", "length": 1}, {"name": "B", "length": 1}]},
        {"name": "C", "length": 2}]},
    {"name": "D", "length": 3},
]}
fig, axes = plt.subplots(1, 2, width_ratios=[2.2, 1])
cf.multiple_alignment(axes[0], samples)
cf.tree(axes[1], target_tree)
axes[1].set_ylim(axes[0].get_ylim())  # put each tip level with its row
plt.show()
```

A and B differ at one of the eight columns, A and C at three, and D differs from each of the others at four. Alignment matters because comparing unrelated positions would turn errors in correspondence into apparent evolutionary changes. We use only columns containing a known DNA base in every row, so all pairwise comparisons rest on the same columns; computing the alignment itself is a separate task. We first summarize the differences in a matrix, then ask which evolutionary tree could explain that matrix.

## Observed differences and evolutionary distance

For $L$ retained columns, define the mismatch fraction
$$
f_{ij}=\frac{1}{L}\sum_{\ell=1}^L
\mathbf{1}\{x^{(i)}_\ell\ne x^{(j)}_\ell\}.
$$
The indicator contributes one at a mismatch and zero otherwise, so $f_{ij}$ is the fraction of retained columns at which samples $i$ and $j$ disagree. In @fig-upgma-data, $f_{AB}=1/8$, $f_{BC}=2/8$, $f_{AC}=3/8$, and $f_{AD}=f_{BD}=f_{CD}=4/8$. The matrix of these fractions is symmetric and nonnegative with a zero diagonal. Because every pair is compared on the same columns, it also satisfies the triangle inequality $f_{ij}\leq f_{ik}+f_{kj}$: wherever $i$ and $j$ differ, at least one of them must differ from $k$.

Mismatch fractions describe the observed letters, not the full history of substitutions. Along a lineage, A could change to G and later back to A, leaving no mismatch. A change to G followed by T leaves one mismatch after two substitutions. This loss of historical information becomes substantial for divergent sequences [@durbin1998].

An *evolutionary distance* estimates the expected number of substitutions per site along the path between two samples. Unlike a mismatch fraction, it can exceed one. The tree algorithm takes these distances as input, whether they come from a substitution model or from raw mismatch counts.

### A correction for hidden substitutions
The Jukes–Cantor model assumes equal base frequencies and equal rates for changes between each pair of different bases. Its predicted mismatch fraction after total path distance $d$ is
$$
f=\frac34\left(1-e^{-4d/3}\right).
$$
Taking this model relation as given, rearranging and taking the natural logarithm gives
$$
d=-\frac34\ln\left(1-\frac43f\right),\qquad 0\leq f<\frac34.
$$
For A and B, $f=1/8$ gives $d\approx0.137$, only slightly more than the observed fraction. For A and D, $f=1/2$ gives $d=-\tfrac34\ln\tfrac13\approx0.824$: at this level of divergence the observed mismatches miss about two-fifths of the estimated substitutions (@fig-upgma-saturation). As the lineages diverge, their bases approach independent uniform draws, so the probability of a mismatch approaches $3/4$. The observed difference then *saturates*: large additional evolutionary distances produce small changes in $f$. The correction diverges as $f$ approaches $3/4$ and is not a finite real estimate at or above that value. Near saturation, sampling error in $f$ becomes large uncertainty in $d$ [@durbin1998].

```{.python .execute}
#| label: fig-upgma-saturation
#| fig-cap: "Expected mismatch fraction under the Jukes-Cantor model (blue) and the line on which every substitution stays visible (dashed). The red arrow is the correction for samples A and D, from an observed fraction of 1/2 to about 0.82 substitutions per site. Near the 3/4 level, extra substitutions barely change the observed fraction."
#| fig-alt: "A rising concave curve of mismatch fraction against substitutions per site approaches a dotted horizontal line at 0.75; a dashed diagonal line rises above it. A red horizontal arrow at height 0.5 runs from the diagonal at 0.5 to the curve at about 0.82."
distance = np.linspace(0, 3, 301)
mismatch = 0.75 * (1 - np.exp(-4 * distance / 3))
f_ad = 0.5
d_ad = -0.75 * np.log(1 - 4 * f_ad / 3)
fig, ax = plt.subplots()
ax.plot(distance, mismatch, label="Expected mismatch fraction")
ax.plot([0, 1], [0, 1], "--", color=cf.RULE, label="Every change remains visible")
ax.axhline(0.75, color=cf.RULE, linewidth=0.8, linestyle=":")
ax.text(3, 0.77, "3/4", ha="right", va="bottom", fontsize=9, color=cf.RULE)
ax.annotate("", xy=(d_ad, f_ad), xytext=(f_ad, f_ad),
            arrowprops={"arrowstyle": "->", "color": cf.EMPHASIS, "lw": 1.5})
ax.text((f_ad + d_ad) / 2, f_ad + 0.02, "A, D", ha="center", va="bottom",
        fontsize=9, color=cf.EMPHASIS)
ax.set(xlabel="Evolutionary distance (substitutions per site)",
       ylabel="Mismatch fraction", ylim=(0, 1))
ax.legend(loc="lower right")
plt.show()
```

## What the tree must explain

We need to distinguish a branching pattern from the distances it predicts. A rooted tree has a distinguished ancestral node, the *root*. Its *leaves* or *tips* are the sampled sequences; internal nodes represent inferred common ancestors. Each edge $e$ has a nonnegative branch length $\ell_e$. For the unique path $P_T(i,j)$ between two tips in tree $T$, define
$$
d_T(i,j)=\sum_{e\in P_T(i,j)}\ell_e.
$$
These are *additive tree distances*. We want to reverse this computation [@compeau2015].

::: {#def-distance-tree-problem}
## Distance-based tree reconstruction
**Input:** an $n\times n$ distance matrix $D=(d_{ij})$ for $n$ sampled sequences.

**Output:** a tree $T$ with the samples as tips and nonnegative branch lengths that *fits* $D$: $d_T(i,j)=d_{ij}$ for every pair $i,j$.
:::

Arbitrary matrices need not admit a fitting tree; we will establish a guarantee under a stronger assumption.

Pairwise distances alone generally do not locate the root. UPGMA adds a *molecular-clock assumption*: all lineages accumulate substitutions at the same rate, and the tips are sampled at the same time. Under this model every root-to-tip path has the same length [@durbin1998]. A node's height $h_u$ is its distance down to any descendant tip. Thus tips have height zero and an edge from parent $u$ to child $v$ has length $\ell_{uv}=h_u-h_v$.

If $u$ is the most recent common ancestor of tips $i,j$, their paths to $u$ have equal length. Therefore
$$
d_{ij}=2h_u.
$$
Heights are measured in substitutions per site; converting them to calendar ages requires a calibrated substitution rate.

::: {.callout-note title="Key idea"}
Under a molecular clock, the distance between two tips is twice the height of their most recent common ancestor. The smallest distance in the matrix therefore marks the most recent split, so merging the closest pair first recovers that split.
:::

::: {#def-ultrametric}
## Ultrametric distances
A symmetric, nonnegative matrix with zero diagonal is **ultrametric** if
$$
d_{ij}\leq\max(d_{ik},d_{jk})
\quad\text{for every }i,j,k.
$$
Equivalently, among the three distances for any three tips, the largest two are equal. All three may be equal.
:::

A clock tree has this property. For three tips, either all share the same most recent common ancestor, or two meet below the ancestor they share with the third. In the latter case, the two distances to the third tip are equal and larger. This gives an observable constraint on the matrix, stronger than the ordinary triangle inequality.

## UPGMA: merge the closest clusters

UPGMA builds the tree from the bottom up as if the data came from a clock tree, which lets it return a tree even when the matrix fits none exactly [@compeau2015]. Start with one cluster per tip. After a merge, a cluster represents an inferred ancestral group. How should we measure its separation from another group? UPGMA uses the average over all original cross-cluster tip pairs:
$$
d(C,D)=\frac{1}{|C||D|}\sum_{i\in C}\sum_{j\in D}d_{ij}.
$$
Every tip pair has equal weight. The method's name is *unweighted pair group method with arithmetic mean*; equal weighting of tips requires weighting cluster updates by their sizes.

When $C$ and $D$ merge, partitioning the sum into pairs originating in $C$ and pairs originating in $D$ gives
$$
d(C\cup D,E)
=\frac{|C|d(C,E)+|D|d(D,E)}{|C|+|D|}.
$$
This update avoids revisiting every original pair.

::: {.algorithm}
**Algorithm: UPGMA.**
Initialize singleton clusters, each with size one and height zero.

1. Choose two distinct active clusters $C,D$ with smallest distance $\delta$.
2. Create their parent at height $\delta/2$. Its branches to their roots have lengths $\delta/2-h_C$ and $\delta/2-h_D$.
3. Replace $C,D$ by $C\cup D$, whose size is $|C|+|D|$. Update its distances to the other clusters by the size-weighted formula.
4. Repeat until one cluster remains; its node is the root. Break ties by a stated rule.
:::

The merge distances never decrease. If $\delta$ is the current minimum, each distance from the new cluster is a weighted average of two distances that are both at least $\delta$:
$$
d(C\cup D,E)\geq\min\{d(C,E),\,d(D,E)\}\geq\delta,
$$
and the retained distances are also at least $\delta$. Parent heights therefore never lie below child heights, so the algorithm produces nonnegative branches even when its input is not ultrametric.

## Four samples, three merges

Distances estimated from finite sequences are rarely exactly ultrametric, because sampling noise and departures from a clock perturb them. Suppose a distance estimator supplies the matrix below, in units of $0.01$ substitutions per site for convenient hand arithmetic (an entry 4 means 0.04). It is close to ultrametric but not exactly so: in the triple A, B, C the two largest distances are 7 and 9.

::: {#exm-upgma-four-samples}
## UPGMA on four samples

| | A | B | C | D |
|---|---:|---:|---:|---:|
| A | 0 | 4 | 7 | 16 |
| B | 4 | 0 | 9 | 18 |
| C | 7 | 9 | 0 | 20 |
| D | 16 | 18 | 20 | 0 |

The closest pair is A,B at distance 4. Their parent U has height 2, and both terminal branches have length 2. Averaging over tip pairs gives
$$
d(AB,C)=\frac{7+9}{2}=8,
\qquad d(AB,D)=\frac{16+18}{2}=17.
$$
The smallest remaining distance is $d(AB,C)=8$, so AB merges with C. Their parent V has height 4; the branch V→U has length $4-2=2$, while V→C has length 4. The clusters now have unequal sizes, and the weights matter:
$$
d(ABC,D)=\frac{2\,d(AB,D)+1\,d(C,D)}{3}
=\frac{2(17)+20}{3}=18.
$$
This is the average of the three original distances 16, 18, and 20 from D. An unweighted average of the two cluster distances would give $(17+20)/2=18.5$, counting the single pair C,D as heavily as the two pairs A,D and B,D together. The root R has height 9, with a branch of length 5 to V and 9 to D.
:::

Under the clock model, A and B share the most recent ancestor, C joins their lineage earlier, and D separates earlier still; the heights 0.02, 0.04, and 0.09 substitutions per site are not dates. A is not inferred to be the ancestor of B: both are observed tips descending from U.

The tree cannot reproduce every input. Its paths from A and from B to C both have length $2+2+4=8$, where the matrix has 7 and 9, and all three paths to D have length 18, where it has 16, 18, and 20. Each path distance between two merged clusters is the average of the input distances it replaces, so comparing input and path distances shows how far the data depart from a clock.

## Computing and drawing the result

The implementation stores a size, a height, and the distances for each active cluster. It returns a nested tree, with each node's name, children, and branch length from its parent, and the merge history for comparison with the hand calculation.

```{.python .execute}
#| label: code-upgma
def upgma(labels: list[str], matrix: np.ndarray) -> tuple[dict, list[tuple]]:
    """Return an average-linkage clock tree and its merge history."""
    matrix = np.asarray(matrix, dtype=float)
    n = len(labels)
    if n == 0 or len(set(labels)) != n or matrix.shape != (n, n):
        raise ValueError("Use distinct labels and a matching nonempty square matrix.")

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
    [0, 4, 7, 16],
    [4, 0, 9, 18],
    [7, 9, 0, 20],
    [16, 18, 20, 0],
])
clock_tree, merge_history = upgma(labels, distances)
assert [step[2] for step in merge_history] == [4, 8, 18]
assert [step[3] for step in merge_history] == [2, 4, 9]
merge_history
```

```{.python .execute}
#| label: fig-upgma-clock-tree
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "UPGMA tree for the matrix of @exm-upgma-four-samples. Parentheses give node heights in units of 0.01 substitutions per site; all tips have height zero, and internal branches carry their lengths. The tips line up because UPGMA imposes a clock, even though the input distances were not exactly ultrametric."
#| fig-alt: "Root at the left branches to D and a cluster containing C and the sister tips A and B; all four tips have the same horizontal coordinate."
fig, ax = plt.subplots()
# Name the internal nodes of this worked instance, leaving the algorithm general.
clock_tree["name"] = "R (9)"
clock_tree["children"][1]["name"] = "V (4)"
clock_tree["children"][1]["children"][1]["name"] = "U (2)"
cf.tree(ax, clock_tree, show_lengths=True)
plt.show()
```

## Why exact ultrametric input is recovered

Under a clock, the closest tips meet at the youngest ancestral height because their distance is twice that height. Merging them therefore recovers the most recent split. Once this group is represented by a cluster, every member has the same distance to a tip outside it: all those paths meet at the same older ancestor. Averaging preserves that distance, so the same reasoning applies at the next merge.

::: {#thm-upgma-recovery}
## Exact distance recovery
Given an ultrametric distance matrix, UPGMA constructs a rooted tree whose tip-to-tip path distances equal every input distance. Tied merges may resolve a group in different orders, which only adds internal edges of length zero.
:::

::: {.proof}
For active clusters $C,D$, maintain the invariant
$$
d_{ij}=d(C,D)\qquad\text{for every }i\in C,\ j\in D,
$$
which holds initially for singletons. Suppose $C,D$ merge at minimum distance $\delta$, and let $E$ be another active cluster. For tips $i\in C$, $j\in D$, $k\in E$, the invariant gives the triple
$$
d_{ij}=\delta,\qquad d_{ik}=d(C,E)\geq\delta,\qquad d_{jk}=d(D,E)\geq\delta.
$$
If the last two differed, the larger would be the unique largest distance in the triple, which ultrametricity forbids. Hence
$$
d(C\cup D,E)=\frac{|C|\,d(C,E)+|D|\,d(D,E)}{|C|+|D|}=d(C,E)=d(D,E),
$$
so every original distance from $C\cup D$ to $E$ equals the stored one, and the invariant survives the merge.

Two tips $i,j$ first share a cluster at a merge with distance $\delta=d_{ij}$. Its node has height $\delta/2$, and branch lengths telescope from it down to each tip, so
$$
d_T(i,j)=\frac{\delta}{2}+\frac{\delta}{2}=d_{ij}.
$$
:::

The proof also shows when the size weights matter: on ultrametric input the two distances being averaged are equal, so every weighting gives the same tree. The weights changed a distance in @exm-upgma-four-samples only because that matrix is not ultrametric.

### Computational cost

With $k$ active clusters, searching all pairs costs $\mathcal{O}(k^2)$ and updating distances costs $\mathcal{O}(k)$. Summing over the $n-1$ merges gives $\mathcal{O}(n^3)$ time and $\mathcal{O}(n^2)$ space for this straightforward implementation, assuming constant-time arithmetic and dictionary lookups.

## When evolutionary rate changes the answer

A small distance can reflect a recent ancestor, a slow rate, or both. Without a clock, the smallest entry of a matrix need not pair two tips that are neighbors in the tree that produced it [@compeau2015]. Suppose the true tree separates the pair A,B from the pair C,D. Give A and C terminal branches of length 1, B and D terminal branches of length 8, and the edge connecting the pairs length 1, all in the same hundredths-of-a-substitution units. Summing along the paths gives

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

The true distances still contain information about the branching pattern. For example, $d_{AB}+d_{CD}=18$, while both $d_{AC}+d_{BD}$ and $d_{AD}+d_{BC}$ equal 20: the latter paths cross the connecting edge twice. Methods for additive distances can use such relationships without imposing a clock; neighbor joining is one example and normally returns an unrooted tree [@durbin1998; @compeau2015].

In the SARS comparison, the smallest distances marked civet viruses as the closest relatives of the human ones; the unequal-rate example shows that this reading assumes no lineage evolved unusually slowly. Distance summaries also discard which columns supplied the differences. The next lesson keeps those characters and asks how many changes a candidate tree requires.

## Exercises

::: {#exr-upgma-mismatches}
Apply the Jukes–Cantor correction to the mismatch fractions of the alignment in @fig-upgma-data, then run UPGMA on the corrected distances. Does the result agree with the tree drawn there? Which triple of samples violates the ultrametric condition?
:::

::: {#exr-upgma-path-checks}
For the tree in @exm-upgma-four-samples, list all six edge lengths and all six tip-to-tip path distances, and compare the path distances with the input matrix. If the substitution rate is 0.01 substitutions per site per million years, what calendar age corresponds to the inferred root height?
:::

::: {#exr-upgma-weights}
In @exm-upgma-four-samples, set $d_{AC}=d_{BC}=8$ and $d_{AD}=d_{BD}=d_{CD}=18$. Trace UPGMA once with size-weighted updates and once with unweighted averages of cluster distances. Show that the traces agree, and explain using the ultrametric condition why they must.
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
