---
published: true
title: "Parsimony and Fitch's Algorithm"
subtitle: "Lesson 15 · 2027-03-17 · Week 10"
nocite: |
  @durbin1998, @compeau2015, @felsenstein1978
---

## Core question

Given aligned DNA sequences and a candidate tree, what is the fewest substitutions needed to explain the observed bases?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## From sequences to ancestral changes

Sequencing gives present-day DNA strings. To compare a gene across organisms, we align corresponding sequences: each row belongs to a sampled organism and each column represents a proposed corresponding position. A column is also called a *site*, and its base is a *character state*. The alignment is an input today; inferring the alignment and checking whether positions really correspond are separate tasks.

Lesson 14 reduced sequence differences to pairwise distances. Here we retain the bases at each site. Suppose two organisms carry A and two carry G at a corresponding position. One evolutionary explanation puts the A-bearing organisms together and needs one change between A and G. Another arrangement may need two changes. The question concerns both relationships and the changes that could have occurred along them [@compeau2015, ch. 7, pp. 35–38].

We observe the leaf sequences, not the sequences at ancestral branching points. For a candidate tree, we will assign bases to its internal nodes to minimize the number of changes across edges. This is the *fixed-tree* or *small parsimony* problem. Comparing the best scores across different trees is a second problem. Keeping them separate lets us understand what the algorithm guarantees [@durbin1998, pp. 174–176].

## The data and the objective

Use this synthetic two-site alignment throughout:

| sampled taxon | site 1 | site 2 | sequence |
|---|---|---|---|
| a | A | C | AC |
| b | A | C | AC |
| c | G | C | GC |
| d | G | T | GT |

A *taxon* here is simply one named sampled organism or group. The short rows are teaching data, not a claim about particular species. We assume every entry is an unambiguous DNA base; gaps and uncertain bases are excluded from this example.

The candidate tree pairs a with b and c with d. Write it as $ ((a,b),(c,d)) $. Let u join a and b, v join c and d, and r join u and v. This specifies the *topology*, or branching pattern. No numerical evolutionary branch lengths are supplied.

::: {#def-parsimony-score}
## Equal-cost parsimony
Let $T$ be the supplied tree and $E(T)$ its edge set. At one site, let $z_w\in\Sigma=\{A,C,G,T\}$ be the state assigned to node w, fixing each leaf to its observed base. The cost is
$$
C(z;T)=\sum_{\{v,w\}\in E(T)}\mathbf 1\{z_v\ne z_w\}.
$$
Write $x_{\ell j}$ for the observed base at leaf $\ell$ and site $j$, and let $m$ be the number of sites. The site and alignment scores are
$$
c_j(T)=\min_{z:\,z_\ell=x_{\ell j}\text{ for every leaf }\ell} C(z;T),
\qquad
c(T)=\sum_{j=1}^{m}c_j(T).
$$
:::

The indicator contributes zero for equal bases and one for different bases. Every substitution has the same cost, regardless of its type or direction. Sites can be optimized separately because their assigned states have no cross-site constraints and the objective adds their edge costs.

Fewer changes is a modeling preference, not a universal law of evolution. The minimum counts changes sufficient to explain the observed endpoints; it cannot recover repeated substitutions hidden on the same edge. An inferred ancestral base is an optimal explanation under this objective, not a sequenced ancestor or a posterior probability.

The node r gives the algorithm a direction in which to process the tree. Under this symmetric cost, inserting a root on any edge leaves the minimum score unchanged: for edge endpoints a and b, minimizing $\mathbf1\{a\ne t\}+\mathbf1\{t\ne b\}$ over an inserted state t gives $\mathbf1\{a\ne b\}$. Thus this objective alone cannot identify the biological root or date a divergence.

## Fitch's recurrence

With $N$ leaves, a rooted binary tree has $N-1$ internal nodes, so enumerating ancestral bases tries $4^{N-1}$ assignments at each site. Choosing one base too early can also lose an optimum. If a pair of leaves has A and C, either A or C at their parent costs one change; the rest of the tree may decide which choice is useful. We retain both possibilities until that context is available.

For each node v, store its subtree minimum $c_v$ and a set $S_v$ of states attaining that minimum *at the root of that subtree*. At a leaf with observed base b, set $S_v=\{b\}$ and $c_v=0$. For an internal node with children u and w, use
$$
(S_v,c_v)=
\begin{cases}
(S_u\cap S_w,\;c_u+c_w),&S_u\cap S_w\ne\varnothing,\\
(S_u\cup S_w,\;c_u+c_w+1),&S_u\cap S_w=\varnothing.
\end{cases}
$$
An intersection means both children can use the same state without adding a change. Disjoint sets require a change; their union records the parent choices that attain the new minimum. Process children before their parent, a *postorder* traversal. The final answer is $c_r$.

Why is a set enough when there are four possible bases and many ancestral assignments? What the parent needs is the cost of attaching to the solved subtree.

::: {#lem-fitch-attachment}
## Cost of attaching a subtree
Suppose a subtree rooted at v has minimum cost $c_v$ and minimizing root states $S_v$. If its parent has fixed state $t$, the minimum cost of the subtree **including its edge to that parent** is
$$
c_v+\mathbf 1\{t\notin S_v\}.
$$
:::

::: {.proof}
If $t\in S_v$, choose a minimum-cost subtree rooted at t and pay zero for the connecting edge. Otherwise, choose any root state in $S_v$ and pay one edge change, achieving $c_v+1$. A lower cost would have to be $c_v$, since costs are integers. That would require both a minimum-cost subtree and a zero-cost edge, making t a member of $S_v$, a contradiction.
:::

The connecting edge matters. For two A leaves, forcing their parent to C costs two changes inside that subtree; attaching an A-rooted subtree to a C state above it costs only one. The lemma concerns the latter quantity, which is exactly what the next recurrence step uses.

::: {#thm-fitch-correctness}
## Correctness of Fitch's algorithm
For a binary tree and the equal-cost objective, the recurrence returns the minimum number of changes at each site.
:::

::: {.proof}
Induct on subtrees. The leaf initialization fixes the observed state. Assume each child has its true minimum and its minimizing root states. By @lem-fitch-attachment, each child contributes its minimum plus zero or one according to whether the proposed parent state belongs to its set.

The two child contributions add, since their assignments can be chosen independently after fixing t. If the sets intersect, choosing t in the intersection pays no extra change; every other t pays at least one. If they are disjoint, no t can avoid both penalties. Choosing t in their union pays exactly one, while choosing outside pays two. This establishes both the minimum and its minimizing root states, completing the induction.
:::

## Two sites by hand

At site 1, a and b both have A, so $S_u=\{A\}$ and $c_u=0$. Similarly $S_v=\{G\}$ and $c_v=0$. Their sets are disjoint: $S_r=\{A,G\}$ and $c_r=1$. Assigning A to r and u, and G to v, places the one change on r–v. Either A or G at the root gives a one-change explanation.

At site 2, $S_u=\{C\}$ with cost zero. The C/T pair under v gives $S_v=\{C,T\}$ with cost one. The root intersection is $\{C\}$, so the total remains one. Assign C to all three internal nodes; only the edge to d changes. Adding sites gives alignment score two.

### Reconstructing one ancestral assignment

Choose any state in the root set. Move from parent to child: if the parent's chosen state belongs to the child set, reuse it; otherwise choose any state from that set and pay one edge change. Continue to the leaves. The proof's child contributions show that this procedure attains the score.

The sets describe minima for isolated subtrees. To see why this matters, change the leaf bases to $a=A$, $b=C$, $c=G$, $d=G$. Now $S_u=\{A,C\}$ with cost 1, $S_v=\{G\}$ with cost 0, and the full tree costs 2. Choosing $r=v=G$ and $u=A$ pays on the edges u-to-b and r-to-u. Choosing $r=u=v=G$ instead pays on u-to-a and u-to-b, still totaling 2, even though $G\notin S_u$. The extra change inside u's subtree saves a change above it. Fitch reconstruction supplies one optimum; finding all optimal ancestral assignments requires the state-conditioned costs described by weighted parsimony [@durbin1998].

## A runnable scorer and reconstruction

::: {.algorithm}
**Algorithm: Fitch scoring and one reconstruction.**

1. In postorder, initialize each leaf set to its observed base. At each internal node, use a nonempty child intersection; otherwise use the union and add one change. Sum child scores.
2. Choose a state from the root set.
3. In preorder, reuse the parent state when it belongs to the child set; otherwise choose a child-set state. Return the root score and assigned states.
:::

Represent the binary tree by a dictionary mapping each internal node to its two children; nodes absent from that dictionary are leaves. A separate dictionary supplies one observed base per leaf. The function returns a score, the Fitch sets, and one assignment. Sorting a set before choosing a state makes ties reproducible; alphabetical order has no biological meaning.

```{.python .execute}
#| label: code-fitch-site
def fitch_site(
    children: dict[str, tuple[str, str]],
    root: str,
    leaf_states: dict[str, str],
) -> tuple[int, dict[str, set[str]], dict[str, str]]:
    """Score one site on a binary tree and return one optimal assignment."""
    sets = {}

    def postorder(node: str) -> int:
        if node not in children:
            sets[node] = {leaf_states[node]}
            return 0
        left, right = children[node]
        score = postorder(left) + postorder(right)
        shared = sets[left] & sets[right]
        if shared:
            sets[node] = shared
            return score
        sets[node] = sets[left] | sets[right]
        return score + 1

    score = postorder(root)
    assigned = {root: min(sets[root])}

    def preorder(node: str) -> None:
        if node not in children:
            return
        for child in children[node]:
            if assigned[node] in sets[child]:
                assigned[child] = assigned[node]
            else:
                assigned[child] = min(sets[child])
            preorder(child)

    preorder(root)
    return score, sets, assigned


children = {"r": ("u", "v"), "u": ("a", "b"), "v": ("c", "d")}
alignment = {"a": "AC", "b": "AC", "c": "GC", "d": "GT"}
site_scores = []
for site in range(2):
    leaf_states = {}
    for taxon, sequence in alignment.items():
        leaf_states[taxon] = sequence[site]
    score, sets, assigned = fitch_site(children, "r", leaf_states)
    site_scores.append(score)
    edge_changes = 0
    for parent, daughters in children.items():
        for daughter in daughters:
            edge_changes += assigned[parent] != assigned[daughter]
    assert edge_changes == score
    assert all(assigned[leaf] == base for leaf, base in leaf_states.items())
print(site_scores, sum(site_scores))
```

The output $[1,1]$, total 2, agrees with the hand calculation. The figure draws the sets returned by the same computation.

```{.python .execute}
#| label: fig-fitch-two-sites
#| fig-cap: "Computed Fitch sets for the two alignment sites. Disjoint child sets require a change; a shared state lets the parent join the child subtrees without adding one. Branch extents are drawing coordinates, not evolutionary lengths."
#| fig-alt: "Two four-leaf trees pairing a with b and c with d. Site 1 combines A and G into root set A,G. Site 2 intersects C with C,T to give root set C."
#| fig-width: 6
#| fig-height: 5.5
def labeled_subtree(node: str, sets: dict[str, set[str]]) -> dict:
    """Build drawing labels from the computed Fitch sets."""
    label = node + ": {" + ",".join(sorted(sets[node])) + "}"
    result = {"name": label, "length": 1}
    if node in children:
        result["children"] = [labeled_subtree(child, sets) for child in children[node]]
    return result

fig, axes = plt.subplots(2, 1)
for site, ax in enumerate(axes):
    observed = {taxon: sequence[site] for taxon, sequence in alignment.items()}
    score, sets, _ = fitch_site(children, "r", observed)
    cf.tree(ax, labeled_subtree("r", sets))
    ax.set_title(f"Site {site + 1}: {score} change")
plt.show()
```

For N leaves a binary tree has $2N-1$ nodes. With the fixed four-base alphabet, set operations take bounded time, so one site costs $\mathcal O(N)$ time and space. Scoring m sites costs $\mathcal O(mN)$ time. Processing one site at a time retains $\mathcal O(N)$ working storage; storing all ancestral sequences additionally needs $\mathcal O(mN)$ space.

## Scoring trees versus searching for a tree

Four named leaves have three unrooted binary topologies, described by their splits into pairs. Score the same alignment on each:

| split | site 1 | site 2 | total |
|---|---:|---:|---:|
| ab \| cd | 1 | 1 | 2 |
| ac \| bd | 2 | 1 | 3 |
| ad \| bc | 2 | 1 | 3 |

For ac \| bd, both child pairs contain A and G at site 1. Each pair costs one; their root sets intersect, so the total is two. Site 2 costs one on every topology because T occurs only at d: putting C at every internal node and changing the edge to d always works. This site contributes changes but no preference among these splits. Site 1 drives the preference for ab \| cd.

```{.python .execute}
#| label: code-fitch-three-trees
splits = [("a", "b", "c", "d"), ("a", "c", "b", "d"), ("a", "d", "b", "c")]
for a, b, c, d in splits:
    candidate = {"r": ("u", "v"), "u": (a, b), "v": (c, d)}
    total = 0
    for site in range(2):
        observed = {}
        for taxon, sequence in alignment.items():
            observed[taxon] = sequence[site]
        score, _, _ = fitch_site(candidate, "r", observed)
        total += score
    print(f"{a}{b}|{c}{d}", total)
```

This exhaustive comparison is feasible for four taxa. For $N\ge3$ named taxa, the number of unrooted binary trees is $(2N-5)!!$, the product of odd integers from 1 through $2N-5$. Adding another leaf inserts it on any of the current tree's $2N-3$ edges, explaining the next multiplicative factor. There are $15!!=2{,}027{,}025$ trees for 10 taxa and $35!!\approx2.22\times10^{20}$ for 20 taxa. Fast scoring of one tree therefore does not make exhaustive tree search practical. Larger searches often improve a starting tree by local rearrangements, which can stop before the global minimum. Fitch's guarantee concerns the supplied tree, not such a search strategy.

## Biological limits of the answer

A score of two means two changes suffice under our alignment and equal-cost assumptions. It does not mean exactly two substitutions happened historically. Repeated change can restore an earlier base; unrelated lineages can acquire the same base independently. In particular, rapidly evolving lineages can accumulate shared-looking states through repeated or convergent change and be grouped together incorrectly, a problem called *long-branch attraction* [@felsenstein1978].

Real substitution types also have different frequencies: changes between A/G or between C/T need not behave like changes across those pairs. Weighted parsimony replaces the zero/one edge cost with an explicit cost table and keeps a minimum for each possible parent state. That changes the modeling assumptions as well as the calculation.

Alignment mistakes, missing characters, and unrepresentative sites can change the preferred topology. Our toy assumes correspondence and excludes uncertainty. A real analysis must specify how such entries are treated and examine whether the answer persists across suitable sites. Agreement with the minimum-change objective is evidence conditional on those choices, rather than proof of the true genealogy.

## Exercises

::: {#exr-fitch-hand}
Replace site 2 by A at a and c, G at b and d. Compute both site scores on all three splits above. Which splits minimize the combined score? Explain what the conflict between sites means for the answer.
:::

::: {#exr-fitch-reconstruction}
For the original site 1, start the root reconstruction with G instead of A. Give every internal state and locate the changed edge. Explain why this does not determine the direction or time of the historical substitution.
:::

::: {#exr-fitch-validation}
For a three-leaf tree, enumerate the four possible states at each of its two internal nodes. Compare that minimum with `fitch_site` for a site of your choice. Include an all-equal site as a boundary check.
:::

::: {#exr-fitch-costs}
Give A↔G changes cost 1 and all other unequal changes cost 2. Explain why the ordinary intersection/union counter is no longer sufficient, and write the state-conditioned recurrence needed for one internal node.
:::

::: {#refs}
**References**
:::
