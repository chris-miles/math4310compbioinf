---
published: false
title: "Modern Frontiers in Bioinformatics"
subtitle: "Lesson 24 · 2027-04-21 · Week 15"
nocite: |
  @compeau2015, @hickey2024pangenome
---

## Core question

Which genome differences can longer DNA reads and references representing several genomes help us recover?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import math
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## Longer reads and richer references

Two earlier problems had a common cause. In assembly, reads shorter than a repeat could not say which unique sequence follows which copy. In read mapping, each read was compared with one reference genome, so a base absent from that genome counted as a mismatch even when many people carry it. The goal is to recover missing or duplicated DNA, including changes in which genes a sample carries, that short reads against one reference can overlook. Such changes can alter which proteins a cell makes.

The short reads behind most of our examples are typically under 300 bases. Long-read instruments measure single DNA molecules over far longer stretches: one platform passes a polymerase around the same circular molecule several times and combines the passes into a consensus read of more than 10,000 bases with better than 99% accuracy; nanopore instruments report reads of 10,000 to 100,000 bases, and "ultra-long" reads beyond that [@logsdon2020longread]. A long read is still an inferred base string without genomic coordinates, but it carries much more context.

Context matters because human DNA is highly repetitive. About half of the human genome derives from transposons, DNA elements that can move or copy themselves to new genome locations, and the roughly 300-base Alu element alone has about a million copies with small differences among them [@compeau2015]. About 5% of the genome lies in segmental duplications, long stretches present in two or more nearly identical copies [@logsdon2020longread]. A 150-base read lying inside an Alu copy may match many of the million copies equally well, and then it cannot say which one produced it.

Recall `ACGTCGACG` and `ACGACGTCG` from Lesson 7. They have the same 3-mer multiset, so error-free 3-base reads cannot distinguish them: in @fig-frontiers-read-layout, the same seven reads fit under both. One 4-base read changes this. `TCGA` occurs in the first candidate only, and `ACGA` in the second only.

```{.python .execute}
#| label: fig-frontiers-read-layout
#| code-fold: true
#| code-summary: "Place the same short reads under both candidates"
#| fig-cap: "The 3-base reads of each candidate at their starts, followed by one 4-base read containing a copy of CG and a base on each side. The short reads form the same collection in both panels; each 4-base read fits only its own candidate."
#| fig-alt: "Two panels, each a nine-letter candidate with seven staggered three-letter reads beneath it and a final four-letter read: TCGA under ACGTCGACG on the left, ACGA under ACGACGTCG on the right."
#| fig-width: 6
#| fig-height: 3.4
candidates = ["ACGTCGACG", "ACGACGTCG"]
bridging_reads = ["TCGA", "ACGA"]
fig, axes = plt.subplots(1, 2)
for number, (ax, candidate, bridging) in enumerate(
    zip(axes, candidates, bridging_reads), start=1
):
    short_reads = [candidate[i:i + 3] for i in range(len(candidate) - 2)]
    reads = short_reads + [bridging]
    starts = list(range(len(short_reads))) + [candidate.index(bridging)]
    labels = ["3-base"] * len(short_reads) + ["4-base"]
    cf.read_layout(ax, reads, starts=starts, reference=candidate, labels=labels)
    ax.set_title(f"Candidate {number}")
plt.show()
```

We first ask how long a read must be to settle such an ambiguity and how many reads we need. Then we turn to references that store several genomes as paths through one graph.

## Bridging a repeat

In `ACGTCGACG`, the string `CG` occurs three times. A walk through the Lesson 7 de Bruijn graph enters `CG` three times and must choose which excursion comes next. The genome fixes that choice through the letters before and after each copy: the middle copy is entered from `T` and left toward `A`. A read containing the copy and both neighbors records that pairing.

::: {#def-bridging-read}
## Bridging read
Let a string $R$ of length $\ell$ occur in a genome at positions $p,\ldots,p+\ell-1$. A read *bridges* this copy if it covers positions $p-1$ through $p+\ell$: the whole copy and at least one base on each side. A read of length $L$ can bridge a copy only if $L\geq \ell+2$.
:::

Here `TCGA` bridges the middle copy of `CG`, with $\ell=2$ and $L=4$. The next result shows that no shorter read could have helped, however many we collect.

::: {#prp-short-reads-cannot-order}
## Reads too short to bridge cannot order repeat copies
Let $R$ have length $\ell$, let $X$ and $Y$ be strings, and let $A$ and $B$ be nonempty strings. Define
$$
G_1 = X\,R\,A\,R\,B\,R\,Y, \qquad G_2 = X\,R\,B\,R\,A\,R\,Y .
$$
For every $m\leq \ell+1$, $G_1$ and $G_2$ have the same multiset of length-$m$ substrings.
:::

::: {.proof}
Sketch. Assign each length-$m$ window to the block where it starts. A window starting in $X$ is at most one letter longer than $R$, so it lies inside $XR$. Likewise a window starting in the first copy of $R$ or in $A$ lies inside $RAR$, one starting in the second copy or in $B$ lies inside $RBR$, and one starting in the last copy or in $Y$ lies inside $RY$. The windows of $G_1$ are therefore the windows of the pieces $XR$, $RAR$, $RBR$, $RY$, counted by start. $G_2$ consists of the same pieces with $RAR$ and $RBR$ exchanged, so it has the same windows.
:::

A short window containing a copy sees only one of its flanks, and the genomes differ only in how flanks pair. Any method given only such reads sees identical input from $G_1$ and $G_2$, so no better assembler can recover the order.

For the Lesson 7 strings, take $X=\texttt{A}$, $R=\texttt{CG}$, $A=\texttt{T}$, $B=\texttt{A}$, and $Y$ empty, so that $G_1=\texttt{A·CG·T·CG·A·CG}$. With $\ell=2$, substrings of lengths 1 through 3 must agree. The code prints, for each length, the substrings in excess in each candidate.

```{.python .execute}
#| label: code-frontiers-substrings
def substring_counts(sequence: str, m: int) -> Counter:
    """Multiset of length-m substrings, as a Counter."""
    return Counter(sequence[i:i + m] for i in range(len(sequence) - m + 1))


g1, g2 = "ACGTCGACG", "ACGACGTCG"
for m in range(2, 6):
    only_g1 = sorted((substring_counts(g1, m) - substring_counts(g2, m)).elements())
    only_g2 = sorted((substring_counts(g2, m) - substring_counts(g1, m)).elements())
    print(m, only_g1, only_g2)
```

The candidates separate exactly at length 4, through the bridging reads of @fig-frontiers-read-layout. A bridging read records one pairing of an entering flank with a leaving flank; it decides between candidates when that pairing occurs in only one of them, as here.

### How many reads bridge a copy

A long read can bridge a copy only if it lands in the right place. Suppose a genome of length $n$ is sequenced with $m$ reads of length $L$, giving coverage $C=mL/n$ as in Lesson 6, and assume read starts are independent and uniform, ignoring the genome's ends. A read starting at $t$ bridges the copy at positions $p,\ldots,p+\ell-1$ when $t\leq p-1$ and $t+L-1\geq p+\ell$. The good starts run from $p+\ell-L+1$ to $p-1$, which is $L-\ell-1$ positions when $L\geq\ell+2$. Each read bridges with probability $(L-\ell-1)/n$, so
$$
\Pr(\text{no read bridges the copy}) = \Bigl(1-\frac{L-\ell-1}{n}\Bigr)^{m} \leq \exp\Bigl(-\frac{m(L-\ell-1)}{n}\Bigr) = \exp\Bigl(-C\,\frac{L-\ell-1}{L}\Bigr).
$$ {#eq-frontiers-unbridged}

The inequality is $1-u\leq e^{-u}$, nearly an equality for a genome of billions of bases, and the last step substitutes $m/n=C/L$. The exponent is the expected number of bridging reads: count the $L-\ell-1$ start positions from which a read covers the repeat and both flanking bases, divide by the read length, and multiply by the coverage.

If $L\leq\ell+1$, no start bridges and the probability of a bridge is zero at every coverage. If $L$ is much longer than $\ell$, the fraction is near 1 and the chance of missing the copy is about $e^{-C}$, the same as the chance that a single base goes unsequenced.

::: {.callout-note title="Key idea"}
Coverage cannot substitute for read length. A repeat copy of length $\ell$ can be bridged only by reads of length at least $\ell+2$; once reads are that long, coverage controls how likely a bridge is.
:::

```{.python .execute}
#| label: code-frontiers-bridging
def unbridged_bound(read_length: int, repeat_length: int, coverage: float) -> float:
    """Upper bound exp(-C(L - ell - 1)/L) on missing every bridge of one copy."""
    good_starts = read_length - repeat_length - 1
    if good_starts <= 0:
        return 1.0  # no read of this length can bridge the copy
    return math.exp(-coverage * good_starts / read_length)


settings = [(150, 30), (15_000, 30), (100_000, 10)]  # (L, C)
for repeat_length in [300, 12_000, 50_000]:
    row = [unbridged_bound(L, repeat_length, C) for L, C in settings]
    print(repeat_length, [f"{value:.2g}" for value in row])
```

Each row is a repeat length, and the columns follow `settings`. For an Alu-sized copy, short reads never bridge, while 15,000-base reads expect $30\times 14{,}699/15{,}000\approx 29.4$ bridging reads and miss with probability about $2\times10^{-13}$. For a 12,000-base duplication they expect only about 6, and miss roughly one copy in 400. At $\ell=50{,}000$ only the ultra-long reads help, even at a third of the coverage.

The calculation assumes error-free reads of one length. Errors cause a separate problem: Alu copies differ at only a few positions, and a read whose error rate is comparable to that divergence cannot tell near-identical copies apart even while spanning them. Consensus over repeated passes of one molecule improves accuracy, while length places the unique flanks.

## Pangenome graphs

Long reads now give nearly complete assemblies of individual human genomes, so a reference no longer has to be one sequence. To see why that matters, take two toy genomes that differ at one site, `ACGTT` and `ACATT`, and use the first as a linear reference. With Lesson 4's fitting alignment (every read letter used, reference ends free) and match $+1$, mismatch $-1$, gap $-2$, the read `CGTT` scores 4 and the read `CATT` scores $3-1=2$. Suppose a mapper keeps only placements scoring at least 3, to reject chance matches. It keeps reads carrying `G` and discards reads carrying `A`, so even if half of a sample's copies carry `A`, the retained reads show only `G`. This systematic loss of reads that differ from the chosen reference is *reference bias*.

The Human Pangenome Reference Consortium used long reads to assemble the two inherited chromosome sets of each of 47 individuals separately. These 94 genome sequences, called haplotypes here, were combined into a graph. Compared with workflows based on the linear reference GRCh38, mapping short reads to the graph reduced small-variant discovery errors by 34% and increased the structural variants detected per haplotype by 104% [@liao2023pangenome]. @fig-frontiers-rhd-graph shows one region, where the neighboring genes *RHD* and *RHCE* have highly similar sequences. *RHD* encodes the D protein of the Rh blood group, and in people of European ancestry the D-negative blood type is mostly due to deletion of *RHD* [@gu2014rhd]. The similar pair is a repeat, so short reads from it are hard to place.

::: {#fig-frontiers-rhd-graph .course-figure}

<img src="../assets/borrowed/liao-2023-rhd-rhce-pangenome-paths.png" alt="Panel a: a curved graph with a large and a small loop labeled RHD and a long arm labeled RHCE. Panel b: four copies of the graph, each outlining one haplotype's route. Panel c: ten haplotype names with counts from 43 down to 1 and gene-arrangement diagrams." loading="lazy">

A pangenome graph of the *RHD*–*RHCE* region built from 90 haplotype assemblies: (a) gene positions on the graph; (b) four haplotypes traced as paths, including the GRCh38 reference and one that skips the *RHD* loop; (c) counts of each structural haplotype, asterisks marking newly found ones. Each genome is one route through a shared structure, and the linear reference is only one of these routes. Reproduced from @liao2023pangenome (CC BY 4.0).

:::

In panel c, 15 of the 90 haplotypes carry *RHCE* alone, and several carry hybrid genes. A linear reference containing *RHD* has no sequence for the junction left where a haplotype lacks the gene, or for a hybrid junction, so reads crossing those junctions cannot be placed correctly, although reads from unchanged sequence still map.

::: {#def-variation-graph}
## Variation graph
A *variation graph* is a directed graph whose nodes carry DNA strings. A path $v_1\to v_2\to\cdots\to v_t$ spells the concatenation of its node strings. Each genome used to build the graph is spelled by one path.
:::

Our toy genomes give nodes `AC`, `G`, `A`, `TT`, with `AC` leading to both `G` and `A` and both leading to `TT`. The two paths spell `ACGTT` and `ACATT`. Shared sequence is stored once, and the difference is a bubble, which here records known variation rather than an assembly puzzle. For alignment we split each node into one node per letter, which changes no spelled path, and write $\lambda(v)$ for the letter of node $v$. Numbering the nodes so that every edge goes from a smaller number to a larger one, a *topological order*, gives @tbl-frontiers-letter-graph.

| node | $v_1$ | $v_2$ | $v_3$ | $v_4$ | $v_5$ | $v_6$ |
|---|---|---|---|---|---|---|
| $\lambda(v)$ | A | C | G | A | T | T |
| predecessors | none | $v_1$ | $v_2$ | $v_2$ | $v_3, v_4$ | $v_5$ |

: The letter graph of `ACGTT` and `ACATT`. {#tbl-frontiers-letter-graph}

## Aligning a read to a graph

Given a read $x=x_1\cdots x_L$ and a variation graph, we want a path and a fitting alignment of $x$ to that path's string with the highest score. With two paths we can align to each: `CATT` scores 2 against `ACGTT` and 4 against `ACATT`. But $b$ two-way bubbles in a row give $2^b$ paths, so thirty variant sites give more than $10^9$ paths in a graph of about a hundred nodes. We need a table that shares work across paths.

In the tables of Lessons 3 and 4, column $j$ took its candidates from column $j-1$. Along a graph path, the letter before $\lambda(v)$ belongs to a predecessor of $v$, and there may be several.

For $0\leq i\leq L$ and each node $v$, let $Q(i,v)$ be the best score of a global alignment of $x_1\cdots x_i$ with the letters of a final segment (possibly empty) of some path ending at $v$. Leaving out the start of the path costs nothing, as in fitting alignment. Add a start column $\circ$ as the only predecessor of each node without incoming edges, and write $\operatorname{Pred}(v)$ for the predecessors of $v$. For $i\geq1$,
$$
Q(i,v)=\max\begin{cases}
\max_{u\in\operatorname{Pred}(v)} Q(i-1,u)+s\bigl(x_i,\lambda(v)\bigr), & \text{pair } x_i \text{ with } \lambda(v),\\[2pt]
\max_{u\in\operatorname{Pred}(v)} Q(i,u)+g, & \text{put } \lambda(v) \text{ opposite a gap},\\[2pt]
Q(i-1,v)+g, & \text{put } x_i \text{ opposite a gap},
\end{cases}
$$ {#eq-frontiers-graph-recurrence}
with $Q(0,v)=0$ for every node, $Q(0,\circ)=0$, $Q(i,\circ)=ig$, and answer $\max_v Q(L,v)$.

The first case takes the best alignment of the shorter read prefix ending at any predecessor of $v$, then pairs the next read letter with $v$'s letter. The zero row and the maximum over the last row let the path start and end anywhere. When every node has one predecessor, the graph is a single string and this is Lesson 4's fitting recurrence.

::: {#prp-graph-alignment}
## Alignment to an acyclic graph
If the graph has no directed cycle and each row is filled in topological order, then $\max_v Q(L,v)$ is the best fitting-alignment score of $x$ against the string of any path in the graph. Filling the $(L+1)(|V|+1)$ entries takes $\mathcal O\bigl(L(|V|+|E|)\bigr)$ time.
:::

::: {.proof}
Sketch. Remove the last column of an optimal alignment counted by $Q(i,v)$. If it pairs $x_i$ with $\lambda(v)$, the rest aligns $x_1\cdots x_{i-1}$ with a path ending at a predecessor $u$ (or $u=\circ$); if it puts $\lambda(v)$ opposite a gap, the rest aligns $x_1\cdots x_i$ with a path ending at a predecessor; if it puts $x_i$ opposite a gap, the rest aligns $x_1\cdots x_{i-1}$ with a path still ending at $v$. The rest must be optimal for its prefix and endpoint, or replacing it would improve the whole. The second case uses $Q(i,u)$ from the same row, which is ready because the topological order puts $u$ before $v$. Each node inspects its predecessors twice per row, so a row costs $\mathcal O(|V|+|E|)$.
:::

Each bubble adds a few nodes and edges, so the work grows linearly in the number of bubbles while the number of paths doubles. The saving comes from merges: later entries need only the best alignment arriving at a merge, not which branch produced it.

::: {.algorithm}
**Algorithm: fitting alignment of a read to an acyclic letter graph.**
List the nodes $v_1,\ldots,v_{|V|}$ in topological order and set row 0 to zero. For $i=1,\ldots,L$, set $Q(i,\circ)=ig$, then fill $v_1,\ldots,v_{|V|}$ in order with the maximum of the three cases in @eq-frontiers-graph-recurrence. Return $\max_v Q(L,v)$.
:::

In the code, column 0 is the start column and column $k$ is node $v_k$, so `predecessors[5] == [3, 4]`. A linear reference uses a chain of predecessors.

```{.python .execute}
#| label: code-frontiers-graph-alignment
def graph_fitting_table(read: str, letters: str, predecessors: dict[int, list[int]],
                        match: int = 1, mismatch: int = -1,
                        gap: int = -2) -> list[list[int]]:
    """Fill Q(i, v) for a read against a letter graph in topological order."""
    # row 0 is zero: any beginning of the path may be skipped at no cost
    table = [[0] * (len(letters) + 1) for _ in range(len(read) + 1)]
    for i in range(1, len(read) + 1):
        table[i][0] = i * gap
        for v in range(1, len(letters) + 1):
            pair = match if read[i - 1] == letters[v - 1] else mismatch
            best = table[i - 1][v] + gap               # read letter opposite a gap
            for u in predecessors[v]:
                best = max(best,
                           table[i - 1][u] + pair,     # pair read letter with node
                           table[i][u] + gap)          # node letter opposite a gap
            table[i][v] = best
    return table


def chain(length: int) -> dict[int, list[int]]:
    """Predecessors of a single linear reference of the given length."""
    return {v: [v - 1] for v in range(1, length + 1)}


graph_letters = "ACGATT"
graph_predecessors = {1: [0], 2: [1], 3: [2], 4: [2], 5: [3, 4], 6: [5]}
read = "CATT"
graph_table = graph_fitting_table(read, graph_letters, graph_predecessors)
linear_table = graph_fitting_table(read, "ACGTT", chain(5))
max(graph_table[-1][1:]), max(linear_table[-1][1:])
```

The graph score is 4 and the linear-reference score is 2, matching the path-by-path alignments. The informative entry is $Q(3,v_5)$, where the prefix `CAT` reaches the merge node. Its predecessors hold $Q(2,v_3)=0$, the best alignment of `CA` ending on the `G` branch, and $Q(2,v_4)=2$, which pairs `CA` with `CA` on the `A` branch. With $x_3=\texttt{T}=\lambda(v_5)$,
$$
Q(3,v_5)=\max\bigl\{\max(0,\,2)+1,\ \ \max(Q(3,v_3),\,Q(3,v_4))-2,\ \ Q(2,v_5)-2\bigr\}
=\max\{3,\ \max(-2,0)-2,\ 0-2\}=3 .
$$
Pairing wins by extending the `A` branch. In @fig-frontiers-graph-table, the two inputs to that case come from columns that are not adjacent in the table, because adjacency now comes from the graph.

```{.python .execute}
#| label: fig-frontiers-graph-table
#| code-fold: true
#| code-summary: "Draw the letter graph and its alignment table"
#| fig-cap: "Top: the letter graph of ACGTT and ACATT, with the best path for read CATT in red. Bottom: the table Q(i, v), starting with the start column; red arrows bring both predecessor values into the shaded merge cell Q(3, v5). At a merge, the table keeps the better branch and forgets which one it was."
#| fig-alt: "Upper panel: nodes A1, C2, then G3 above and A4 below, joining at T5 and T6; edges C2-A4-T5-T6 are red. Lower panel: a five-by-seven grid of integers with two red arrows from the G3 and A4 columns of row 2 into the shaded T5 cell of row 3, which holds 3."
#| fig-width: 5.4
#| fig-height: 5.6
node_names = [rf"$\mathtt{{{letter}}}_{k}$" for k, letter in enumerate(graph_letters, start=1)]
graph_edges = [(node_names[u - 1], node_names[v - 1])
               for v, us in graph_predecessors.items() for u in us if u > 0]
node_positions = dict(zip(node_names, [(0, 0), (1, 0), (2, 0.6), (2, -0.6),
                                       (3, 0), (4, 0)]))
best_path = [(node_names[1], node_names[3]), (node_names[3], node_names[4]),
             (node_names[4], node_names[5])]
fig, axes = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 1.7]})
cf.digraph(axes[0], graph_edges, pos=node_positions,
           highlight_edges=best_path, node_size=650)
cf.dp_table(axes[1], list(read), node_names, graph_table,
            highlight=[(3, 5)], arrows=[((2, 3), (3, 5)), ((2, 4), (3, 5))])
plt.show()
```

Under the error-free assumption, the exact match to `ACATT` says the read came from a copy carrying `A`. Whether the sample carries `A`, `G`, or both still depends on many reads and their base qualities, as with the error bubbles of Lesson 7; the graph keeps reads carrying `A` from being discarded before that question is asked.

Two limits matter for real graphs. The recurrence needs a topological order, but variation graphs used in practice can contain cycles or reverse-oriented nodes, for example around inversions and copy-number changes; the paths in panel b of @fig-frontiers-rhd-graph loop back through one region. And a graph allows every combination of its bubbles, including ones no input genome carried, so a high score on such a path is not evidence that the combination was ever observed (@exr-frontiers-recombinant).

## Checking a reported alignment

Graph aligners and AI assistants report results we did not compute. Suppose a tool reports that `CATT` aligns to the path $v_1\to v_2\to v_4\to v_5\to v_6$ with columns `CATT` over `CATT` and score 4. Each consecutive pair of nodes must be an edge in @tbl-frontiers-letter-graph, each alignment row without gaps must be the read or a contiguous piece of the path's string `ACATT`, and the column scores must sum to 4. Here a score of 4 already equals the most a four-letter read can earn, so no better alignment exists; in general, the filled table is what certifies that none does. Even then, an optimum is a statement about the scoring rule; a claim about the sample needs independent reads and an error model. Goh and colleagues draw the same line for AI in bioinformatics: software can speed up an analysis but cannot judge whether its result is biologically valid [@goh2026expertise].

## Exercises

::: {#exr-frontiers-bridge-count}
In `ACGTCGACG`, list every read of length 5 that bridges the middle copy of `CG`, and compare their number with $L-\ell-1$. Which of these reads occur in `ACGACGTCG`?
:::

::: {#exr-frontiers-coverage}
Using @eq-frontiers-unbridged, find the smallest coverage $C$ for which 15,000-base reads leave a given 12,000-base repeat copy unbridged with probability at most $10^{-3}$. Explain why no coverage achieves this with 12,000-base reads.
:::

::: {#exr-frontiers-recombinant}
The genomes `ACGTTCA` and `ACATTGA` differ at two sites. Draw their letter graph, count its paths, and identify the spelled strings that are neither genome. Align the read `GTTG` to the graph and explain what its best path does and does not tell us about the sample.
:::

::: {#exr-frontiers-check}
A tool reports that read `CATA` aligns to the path $v_1\to v_2\to v_4\to v_5\to v_6$ with columns `CATA` over `CATT` and score 3. Apply the checks from the last section. Which check fails, what is the correct score of this alignment, and what optimum does `graph_fitting_table` give?
:::

::: {#refs}
**References**
:::
