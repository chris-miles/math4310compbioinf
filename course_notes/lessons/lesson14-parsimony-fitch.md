---
published: true
title: "Parsimony and Fitch's Algorithm"
subtitle: "Lesson 14 · 2027-03-15 · Week 10"
nocite: |
  @durbin1998, @compeau2015, @felsenstein1978
---

## Core question

Which evolutionary tree explains the observed bases with the fewest changes?

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

Sequencing gives present-day DNA strings, while the ancestral sequences that connect them are usually unknown. Suppose two sampled organisms carry A and two carry G at a corresponding position. Did one change between A and G separate two ancestral groups, or did the same change occur on different branches? The same question arises for whole traits. Compeau and Pevzner recount a 2003 study of winged and wingless stick insects: on the tree built from their DNA sequences, wings had been lost and regained several times, against the long-held view that a lost wing does not return [@compeau2015]. The answer depends on the evolutionary tree and on the states assigned to its ancestors.

To compare a gene across organisms, we use a supplied alignment of corresponding sequences: each row belongs to a sampled organism and each column represents a proposed corresponding position. A column is also called a *site*, and its base is a *character state*. Lesson 13 reduced these observations to pairwise distances, which forgets where two sequences differ and says nothing about the ancestral sequences. Here we retain the letters at each site and ask how many changes are needed to explain them on a candidate tree. This keeps the evidence about which sites support a proposed relationship and which require repeated changes, and it proposes a base for every ancestor [@compeau2015].

We observe the leaf sequences, not the sequences at ancestral branching points. For a candidate tree, we will assign bases to its internal nodes to minimize the number of changes across edges. This is the *fixed-tree* or *small parsimony* problem. Comparing the best scores across different trees is a second problem. Keeping them separate lets us understand what the algorithm guarantees [@durbin1998, pp. 174–176].

## The data and the objective

Use this synthetic two-site alignment throughout:

| sampled taxon | site 1 | site 2 | sequence |
|---|---|---|---|
| a | A | C | AC |
| b | A | C | AC |
| c | G | C | GC |
| d | G | T | GT |

A *taxon* here is simply one named sampled organism or group. We assume every entry is an unambiguous DNA base; gaps and uncertain bases are excluded from this example.

The candidate tree pairs a with b and c with d. Write it as $ ((a,b),(c,d)) $. Let u join a and b, v join c and d, and r join u and v. This specifies the *topology*, or branching pattern. No numerical evolutionary branch lengths are supplied.

@fig-fitch-competing-ancestors shows why both the tree and the ancestral assignments matter. At site 1, pairing the two A-bearing samples together and the two G-bearing samples together admits one change. Pairing each A-bearing sample with a G-bearing sample requires at least one change in each pair. We first solve the ancestral-assignment problem for a supplied tree, then use its score to compare these arrangements.

```{.python .execute}
#| label: fig-fitch-competing-ancestors
#| echo: false
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "Two candidate trees explaining site 1, where a and b carry A and c and d carry G. Red branches join different ancestral or observed bases: one change when equal bases stay together, two when they are separated across the pairs. Internal labels give one minimum-change assignment; branch extents are drawing coordinates."
#| fig-alt: "Left tree pairs a A with b A under an A ancestor and c G with d G under a G ancestor; the branch from the A root to the G ancestor is red. Right tree pairs a A with c G and b A with d G; both terminal branches to G are red."
def ancestral_tree(pairs: list[tuple[str, str]], states: list[str]) -> dict:
    """Build drawing labels for a supplied site-1 ancestral assignment."""
    observed = {"a": "A", "b": "A", "c": "G", "d": "G"}
    branches = []
    for pair, state, name in zip(pairs, states, ["u", "v"]):
        tips = []
        for taxon in pair:
            tips.append({"name": f"{taxon}: {observed[taxon]}", "length": 1,
                         "state": observed[taxon]})
        branches.append({"name": f"{name}: {state}", "length": 1,
                         "state": state, "children": tips})
    return {"name": "r: A", "length": 0, "state": "A", "children": branches}


def substitution_edges(node: dict) -> list[tuple[str, str]]:
    """Identify branches whose shown parent and child bases differ."""
    changed = []
    for child in node.get("children", []):
        if child["state"] != node["state"]:
            changed.append((node["name"], child["name"]))
        changed.extend(substitution_edges(child))
    return changed

fig, axes = plt.subplots(1, 2)
candidate_trees = [ancestral_tree([("a", "b"), ("c", "d")], ["A", "G"]),
                   ancestral_tree([("a", "c"), ("b", "d")], ["A", "A"])]
for ax, candidate in zip(axes, candidate_trees):
    cf.tree(ax, candidate, highlight_edges=substitution_edges(candidate))
axes[0].set_title("ab | cd: one change")
axes[1].set_title("ac | bd: two changes")
for ax, label in zip(axes, ["A", "B"]):
    ax.text(-0.05, 1.04, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold")
plt.show()
```

::: {#def-parsimony-score}
## Small parsimony with equal costs
**Input:** a rooted binary tree $T$ with edge set $E(T)$, and an alignment giving the observed base $x_{\ell j}$ at each leaf $\ell$ and site $j=1,\ldots,m$.

**Output:** states for the internal nodes at every site that minimize the total number of changes, and that minimum $c(T)$.

At one site, let $z_w\in\Sigma=\{A,C,G,T\}$ be the state assigned to node $w$, fixing each leaf to its observed base. The cost is
$$
C(z;T)=\sum_{\{v,w\}\in E(T)}\mathbf 1\{z_v\ne z_w\}.
$$
The site and alignment scores are
$$
c_j(T)=\min_{z:\,z_\ell=x_{\ell j}\text{ for every leaf }\ell} C(z;T),
\qquad
c(T)=\sum_{j=1}^{m}c_j(T).
$$
:::

Every substitution costs one, regardless of its type or direction. Sites can be optimized separately because their assigned states have no cross-site constraints and the objective adds their edge costs.

Fewer changes is a modeling preference, not a law of evolution, and an inferred ancestral base is an optimal explanation under this objective rather than a sequenced ancestor.

The node r gives the algorithm a direction in which to process the tree. Under this symmetric cost, inserting a root on any edge leaves the minimum score unchanged: for edge endpoint states $s$ and $s'$, minimizing $\mathbf1\{s\ne t\}+\mathbf1\{t\ne s'\}$ over an inserted state $t$ gives $\mathbf1\{s\ne s'\}$. Thus this objective alone cannot identify the biological root or date a divergence.

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

Why is a set enough when there are four possible bases and many ancestral assignments? What the parent needs is the cost of attaching to the solved subtree. Suppose a subtree rooted at v has minimum cost $c_v$ and minimizing root states $S_v$. If its parent has fixed state $t$, the minimum cost of the subtree **including its edge to that parent** is
$$
c_v+\mathbf 1\{t\notin S_v\}.
$$ {#eq-fitch-attachment}

::: {.proof}
Sketch. Two constructions give upper bounds:
$$
t\in S_v:\ \underbrace{c_v}_{\text{subtree rooted at }t}+\underbrace{0}_{\text{edge}},
\qquad
t\notin S_v:\ \underbrace{c_v}_{\text{root state in }S_v}+\underbrace{1}_{\text{edge}}.
$$
Costs are integers, so anything below $c_v+1$ equals $c_v$. That requires a minimum-cost subtree and an unchanged edge, so its root state is $t$ and $t\in S_v$.
:::

The connecting edge matters. For two A leaves, forcing their parent to C costs two changes inside that subtree; attaching an A-rooted subtree to a C state above it costs only one. @eq-fitch-attachment concerns the latter quantity, which is exactly what the next recurrence step uses.

::: {.callout-note title="Key idea"}
A parent needs only its child's pair $(S_v,c_v)$: a parent state inside $S_v$ attaches at no extra cost, and any other state costs exactly one more. One set per node therefore replaces a search over all $4^{N-1}$ ancestral assignments.
:::

::: {#thm-fitch-correctness}
## Correctness of Fitch's algorithm
For a binary tree and the equal-cost objective, the recurrence returns the minimum number of changes at each site.
:::

::: {.proof}
Sketch. Induct on subtrees; a leaf's set and zero cost are correct. Let v have children u and w with correct pairs, and fix a state $t$ at v. The two child subtrees can then be optimized independently, and @eq-fitch-attachment gives the best cost with $z_v=t$:
$$
c_v(t)=c_u+c_w+\mathbf 1\{t\notin S_u\}+\mathbf 1\{t\notin S_w\}.
$$
Minimizing over $t$,
$$
\min_t c_v(t)=\begin{cases}
c_u+c_w, & \text{attained exactly by }t\in S_u\cap S_w\text{ if it is nonempty},\\
c_u+c_w+1, & \text{attained exactly by }t\in S_u\cup S_w\text{ if }S_u\cap S_w=\varnothing,
\end{cases}
$$
since a state outside both disjoint sets pays two. This is the recurrence for both the minimum and its minimizing root states, completing the induction.
:::

## Two sites by hand

::: {#exm-fitch-hand-sites}
## Fitch sets on $((a,b),(c,d))$
At site 1, a and b both have A, so $S_u=\{A\}$ and $c_u=0$. Similarly $S_v=\{G\}$ and $c_v=0$. These sets are disjoint, so $S_r=\{A,G\}$ and $c_r=1$. Assigning A to r and u, and G to v, places the one change on r–v; G at the root works equally well.

At site 2, $S_u=\{C\}$ with cost zero. The C/T pair under v gives $S_v=\{C,T\}$ with cost one, and keeping both states now pays off. The root intersection is $\{C\}$, so the total stays at one: assign C to all three internal nodes, and only the edge to d changes. Had we fixed v to T, the best total would be two: one change on v–c and another between v and the C-bearing subtree under u. Adding the two sites gives alignment score two.
:::

### Reconstructing one ancestral assignment

Choose any state in the root set. Move from parent to child: if the parent's chosen state belongs to the child set, reuse it; otherwise choose any state from that set and pay one edge change. Continue to the leaves. By @eq-fitch-attachment, each step adds exactly the change it pays for, so this procedure attains the score.

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

The output $[1,1]$, total 2, agrees with @exm-fitch-hand-sites. @fig-fitch-two-sites draws the sets and one reconstruction returned by the same computation.

```{.python .execute}
#| label: fig-fitch-two-sites
#| fig-cap: "Computed Fitch sets for the two alignment sites, with the change placed by one reconstruction in red. Disjoint child sets require a change; a shared state lets the parent join the child subtrees without adding one. Branch extents are drawing coordinates, not evolutionary lengths."
#| fig-alt: "Two four-leaf trees pairing a with b and c with d. Site 1 combines A and G into root set A,G, with the branch from r to v in red. Site 2 intersects C with C,T to give root set C, with the branch from v to d in red."
#| fig-width: 6
#| fig-height: 5.5
def fitch_label(node: str, sets: dict[str, set[str]]) -> str:
    """Write a node name with its computed Fitch set."""
    return node + ": {" + ",".join(sorted(sets[node])) + "}"


def labeled_subtree(node: str, sets: dict[str, set[str]]) -> dict:
    """Build drawing labels from the computed Fitch sets."""
    result = {"name": fitch_label(node, sets), "length": 1}
    if node in children:
        result["children"] = [labeled_subtree(child, sets) for child in children[node]]
    return result

fig, axes = plt.subplots(2, 1)
for site, ax in enumerate(axes):
    observed = {taxon: sequence[site] for taxon, sequence in alignment.items()}
    score, sets, assigned = fitch_site(children, "r", observed)
    changed = []
    for parent, daughters in children.items():
        for daughter in daughters:
            if assigned[parent] != assigned[daughter]:
                changed.append((fitch_label(parent, sets), fitch_label(daughter, sets)))
    cf.tree(ax, labeled_subtree("r", sets), highlight_edges=changed)
    ax.set_title(f"Site {site + 1}: {score} change")
plt.show()
```

A rooted binary tree with $N$ leaves has $2N-1$ nodes. With the fixed four-base alphabet, set operations take bounded time, so one site costs $\mathcal O(N)$ time and space, and $m$ sites cost $\mathcal O(mN)$ time.

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

This exhaustive comparison is feasible for four taxa. How many trees would we have to score for $N$ named taxa? Let $U(N)$ count the unrooted binary trees. There is one tree on three leaves, and an unrooted binary tree with $N-1$ leaves has $2N-5$ edges. Attaching leaf $N$ to the middle of any one of those edges gives a different tree on $N$ leaves, and removing leaf $N$ from a tree on $N$ leaves undoes exactly one such attachment. Therefore
$$
U(N)=(2N-5)\,U(N-1),\qquad U(3)=1,
\qquad\text{so}\qquad
U(N)=1\cdot3\cdot5\cdots(2N-5)=(2N-5)!!.
$$
Each new leaf multiplies the count by a larger factor. There are $15!!=2{,}027{,}025$ trees for 10 taxa and $35!!\approx2.22\times10^{20}$ for 20 taxa, while scoring any one of them visits only $2N-1$ nodes per site (@fig-fitch-tree-count). Fast scoring of one tree therefore does not make exhaustive tree search practical. Larger searches often improve a starting tree by local rearrangements, which can stop before the global minimum; Fitch's guarantee concerns the supplied tree, not such a search.

```{.python .execute}
#| label: fig-fitch-tree-count
#| code-fold: true
#| code-summary: "Plot tree counts against the work to score one tree"
#| fig-cap: "The number of unrooted binary trees on N named taxa (red) and the number of nodes Fitch visits to score one rooted tree at one site (blue), on a logarithmic axis. Scoring a tree grows linearly in N, while the number of candidate trees grows faster than any exponential."
#| fig-alt: "Two curves against the number of taxa from 3 to 20 on a log vertical axis: a red curve rising ever more steeply to about two times ten to the twentieth, and a blue line rising slowly to 39."
taxa = list(range(3, 21))
tree_counts = []
count = 1
for n_taxa in taxa:
    if n_taxa > 3:
        count *= 2 * n_taxa - 5  # attach the new leaf to any of 2N - 5 edges
    tree_counts.append(count)
nodes_visited = [2 * n_taxa - 1 for n_taxa in taxa]
fig, ax = plt.subplots()
ax.semilogy(taxa, tree_counts, "o-", color=cf.EMPHASIS, label="unrooted trees $(2N-5)$!!")
ax.semilogy(taxa, nodes_visited, "s-", label="nodes per tree $2N-1$")
ax.set(xlabel="Number of taxa $N$", ylabel="Count", xticks=[3, 5, 10, 15, 20])
ax.legend()
plt.show()
```

## Biological limits of the answer

A score of two means that two changes suffice under our alignment and equal-cost assumptions, not that exactly two substitutions happened. The stick-insect conclusion has the same form: "regained several times" is a minimum count of changes on the inferred tree, and it is only as reliable as that tree. Repeated change can restore an earlier base, and unrelated lineages can acquire the same base independently. Rapidly evolving lineages can therefore accumulate shared-looking states and be grouped together incorrectly, a problem called *long-branch attraction* [@felsenstein1978].

Real substitution types also have different frequencies: changes between A/G or between C/T need not behave like changes across those pairs. Weighted parsimony replaces the zero/one edge cost with an explicit cost table and keeps a minimum for each possible parent state, which changes the modeling assumptions as well as the calculation. Alignment mistakes, missing characters, and unrepresentative sites can also change the preferred topology, so a real analysis states how such entries are treated and checks whether the answer persists across suitable sites.

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
