---
published: true
title: "Reads, k-mers, and de Bruijn Graphs"
subtitle: "Lesson 6 · 2027-02-01 · Week 4"
nocite: |
  @compeau2015, @pevzner2001
---

## Core question

How can overlapping DNA fragments reveal a genome sequence?

```{.python .execute}
#| label: setup
#| code-fold: true
from pathlib import Path
import sys
sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## Reads and coverage

Suppose we sequence DNA from a bacterium without a suitable reference genome. The instrument gives a collection of reads, each a string of base calls along a fragment. Their genomic positions are missing. Placing these reads against a reference would require the sequence we are trying to discover; instead, we must use overlaps among the reads themselves. Reconstructing sequence this way is *de novo assembly*.

For a toy example, the reads `TACG` and `ACGA` suggest the longer string `TACGA`: the final ACG of the first read matches the initial ACG of the second. @fig-reads-records-to-overlaps uses four reads to show the distinction between input records and a proposed placement. Their record names identify fragments, not genomic positions. The overlaps suggest `TACGACG` as a candidate sequence.

```{.python .execute}
#| label: fig-reads-records-to-overlaps
#| code-fold: true
#| code-summary: "Compare unordered read records with proposed overlapping placements"
#| fig-cap: "The same four read records listed independently (left) and placed beneath candidate sequence TACGACG (right). Their overlaps supply the proposed order; several reads observe the same genomic bases, so observation counts differ from genome-copy counts."
#| fig-alt: "Left panel lists r1 CGAC, r2 TACG, r3 GACG, and r4 ACGA without genomic offsets. Right panel staggers those same reads beneath TACGACG at starts two, zero, three, and one respectively."
#| fig-width: 6
#| fig-height: 3
read_records = ["CGAC", "TACG", "GACG", "ACGA"]
record_names = ["r1", "r2", "r3", "r4"]
fig, axes = plt.subplots(1, 2)
cf.read_layout(axes[0], read_records, labels=record_names)
cf.read_layout(axes[1], read_records, starts=[2, 0, 3, 1],
               reference="TACGACG", labels=record_names)
axes[0].set_title("Input: fragment records")
axes[1].set_title("Proposed overlapping placements")
plt.show()
```

A shared word supplies a possible join, but a repeated word elsewhere in the genome can supply the same apparent join. We need a representation that retains the alternatives rather than committing to the first overlap we find.

We begin with the read collection, then isolate a simpler reconstruction problem: find one string with a supplied multiset of exact short words. Those words become edges in a graph, and their overlaps determine which edges can follow one another. Comparing read-derived counts with exact genomic counts will show which assumptions make this ideal problem different from an actual assembly.

Let $G$ be a genome string. An *ideal read* is a substring of $G$ or its reverse complement. Observed reads can differ from these substrings because of sequencing errors. Sequencing data are a multiset $R=\{\!\{r_1,\ldots,r_m\}\!\}$ because identical read sequences may be observed more than once.

A multiset is an unordered collection that retains how often each item occurs. The strings record order within each read; the multiset does not record the reads' positions in the genome.

::: {.biology}
The DNA sample contains many molecules covering overlapping genomic regions; sequencing samples fragments from them, so a genomic position can be observed repeatedly. These repeated measurements help distinguish sequence from errors, but their counts also reflect uneven sampling rather than just the number of copies in the genome.
:::

If $m$ reads have length $L$ and the genome has length $n$, the nominal coverage is $C=mL/n$: the total number of measured bases divided by the number of genomic bases. Coverage is an average, not a guarantee at every position. In practice, genome length may itself need to be estimated.

For the toy examples, all reads have already been oriented to the same strand. Real reads do not arrive with that orientation, so an assembly method must handle the read and its reverse complement consistently [@compeau2015].

::: {#exm-read-coverage}
## Coverage counts observations
Let $G=\texttt{TACGACG}$ and take every length-4 substring:
$$R=\{\!\{\texttt{TACG},\texttt{ACGA},\texttt{CGAC},\texttt{GACG}\}\!\}.$$
Here $C=16/7$. In @fig-reads-records-to-overlaps, end bases occur in one read and the central base occurs in four, even though starts are perfectly regular. The word CGA at genome positions 3 through 5 lies in both r1 and r4; it occurs once in the genome but is observed twice.
:::

### Read boundaries and overlap

Breaking reads into shorter words forgets some useful order: which words came from the same read and how far apart they were. It also exposes short overlaps shared by reads with different endpoints. Those overlaps let us join evidence from many fragments into one graph. We first count the words, then turn each word occurrence into an edge whose endpoints record the overlap.

A shared k-mer is one way to find a possible overlap, but it is not required for two reads to connect in the graph we will build. For example, `TACG` and `CGAC` from @exm-read-coverage share the two-base flank `CG` and no 3-mer. At $k=3$, the last word of the first read, `ACG`, can nevertheless be followed by the first word of the second, `CGA`: their two-base suffix and prefix agree. The shared flank `CG` will be the meeting point of those two edges.

Matching words or flanks do not prove that reads came from neighboring genomic positions. Repeats can create the same matches at different locations. Small $k$ supplies less context for distinguishing those locations; large $k$ requires longer exact overlaps and can disconnect true neighbors when words are missing or contain errors. The choice of $k$ changes which sequence relationships the graph can represent [@compeau2015].

For the reads in @exm-read-coverage, k=3 retains the chain of overlaps:

| Read | 3-mers |
|---|---|
| `TACG` | `TAC`, `ACG` |
| `ACGA` | `ACG`, `CGA` |
| `CGAC` | `CGA`, `GAC` |
| `GACG` | `GAC`, `ACG` |

The table also shows why read-derived counts are not the requested genome composition. The occurrence `CGA` at one genome position is observed in two overlapping reads. Its count of two is evidence about coverage, not evidence for two genomic copies.

::: {.callout-warning title="Two multisets with different meanings"}
The multiset of genomic k-mer occurrences answers “how many copies occur in the genome?” The multiset extracted from reads answers “how many times did the sequencer observe this word?” Assembly begins with the second and tries to infer the first.
:::

## From reads to k-mers

::: {#def-kmer-composition}
## $k$-mer composition
For a length-$n$ string $x$ and integer $k\geq1$,
$$\operatorname{Comp}_k(x)=\{\!\{x_i\cdots x_{i+k-1}:1\leq i\leq n-k+1\}\!\}.$$
This is a multiset: repeated occurrences remain repeated. If $k>n$, there are no valid starting positions, so the composition is empty.
:::

::: {#exm-spectrum-multiplicity}
## Genome and read multiplicity

| 3-mer | `TAC` | `ACG` | `CGA` | `GAC` |
|---|---:|---:|---:|---:|
| copies in $G$ | 1 | 2 | 1 | 1 |
| observations in $R$ | 1 | 3 | 2 | 2 |

Here the genome supplies five 3-mer occurrences, while the four reads supply eight observations. Each length-$L$ read contributes $L-k+1$ words when $L\geq k$, and no words otherwise. This total checks whether a counting program has omitted windows or invented words across read boundaries.
:::

::: {.algorithm}
**Algorithm: count words in reads.**

1. Start an empty count dictionary.
2. Within each read, visit each starting position that leaves at least $k$ symbols.
3. Add one to the count of that length-$k$ word. Do not join separate reads.
:::

```{.python .execute}
#| label: code-reads-kmer-counts
def kmer_counts(strings: list[str], k: int) -> dict[str, int]:
    """Count every length-k substring in a collection."""
    if k < 1:
        raise ValueError("k must be positive.")
    counts = {}
    for sequence in strings:
        for start in range(len(sequence) - k + 1):
            word = sequence[start:start + k]
            counts[word] = counts.get(word, 0) + 1
    return counts

genome = "TACGACG"
reads = ["TACG", "ACGA", "CGAC", "GACG"]
kmer_counts([genome], 3), kmer_counts(reads, 3)
```

## String reconstruction

::: {#def-string-reconstruction}
## String reconstruction from composition
**Input:** an integer $k\geq2$ and a nonempty multiset $K$ of $k$-mers.

**Output:** a string $x$ with $\operatorname{Comp}_k(x)=K$, if one exists.
:::

The model assumes one linear genome, its true multiplicities, no missing $k$-mers, and no errors. Ordering
$$\texttt{TAC},\texttt{ACG},\texttt{CGA},\texttt{GAC},\texttt{ACG}$$
by $(k-1)$-symbol overlaps spells `TACGACG`. There are $5!$ orders of five labeled occurrences to try, and most fail the overlap requirement.

One graph construction puts each word occurrence at a node and draws an arrow when two words overlap in $k-1$ symbols. Reconstruction then asks for a path visiting every occurrence node once, a *Hamiltonian path*. Repeated words require separate nodes. The de Bruijn construction moves each word occurrence onto an edge instead. Its two endpoints are the prefix and suffix that neighboring words must share. A reconstruction can revisit an overlap node while still consuming each word edge once [@compeau2015].

This change matters computationally: using every edge once admits a linear-time traversal algorithm. We develop the construction here and derive that algorithm in the next lesson.

## The de Bruijn graph

::: {#def-de-bruijn-kmers}
## De Bruijn graph of a k-mer multiset
$D(K)$ has one node for every $(k-1)$-mer prefix or suffix in $K$. Each occurrence of $w\in K$ contributes an edge
$$\operatorname{prefix}(w)\longrightarrow\operatorname{suffix}(w).$$
It is a directed multigraph, so parallel edges remain separate.
:::

::: {.algorithm}
**Algorithm: construct the de Bruijn edges.**

1. For each $k$-mer occurrence, remove its last symbol to obtain the source node.
2. Remove its first symbol to obtain the target node.
3. Append that directed edge to a list, retaining repeated occurrences.
:::

```{.python .execute}
#| label: code-reads-de-bruijn
def de_bruijn_edges(kmers: list[str]) -> list[tuple[str, str]]:
    """Return one prefix-to-suffix edge per k-mer occurrence."""
    if not kmers:
        return []
    k = len(kmers[0])
    if k < 2 or any(len(word) != k for word in kmers):
        raise ValueError("k-mers need one common length of at least 2.")
    edges = []
    for word in kmers:
        source = word[:-1]
        target = word[1:]
        edges.append((source, target))
    return edges

ideal_kmers = ["TAC", "ACG", "CGA", "GAC", "ACG"]
edges = de_bruijn_edges(ideal_kmers)
edges
```

```{.python .execute}
#| label: fig-reads-de-bruijn-graph
#| fig-cap: "The de Bruijn graph for the ideal 3-mer composition of `TACGACG`. The two highlighted copies of `ACG` encode genomic multiplicity, not repeated sequencing."
#| fig-alt: "A directed graph with nodes TA, AC, CG, and GA, including two edges from AC to CG."
#| fig-width: 6
#| fig-height: 3
positions = {"TA": (0, 0), "AC": (1.6, 0), "CG": (3.2, 0), "GA": (2.4, -1.3)}
fig, ax = plt.subplots()
cf.digraph(ax, edges, pos=positions, highlight_edges=[("AC", "CG")])
plt.show()
```

### Constructing the graph by hand

Each 3-mer splits into two overlapping 2-mers:

| Edge label | Prefix node | Suffix node |
|---|---|---|
| `TAC` | `TA` | `AC` |
| `ACG` | `AC` | `CG` |
| `CGA` | `CG` | `GA` |
| `GAC` | `GA` | `AC` |
| `ACG` | `AC` | `CG` |

The two `ACG` rows are not duplicates to discard. They represent two positions in the genome and become parallel edges. If we collapsed them, the graph would contain four edges and could spell only a length-6 string, while a five-edge walk spells the required length $5+3-1=7$.

There are two useful ways to read the resulting graph. Locally, an edge says which $(k-1)$-mer can follow another after one new symbol. Globally, a walk chooses an ordering of k-mer occurrences consistent with all those local overlaps. The graph does not choose an ordering for us. A walk may use only part of the input or revisit an edge; a reconstruction must use every supplied edge occurrence exactly once. Complete walks satisfying that additional requirement represent the orderings compatible with the ideal composition.

The construction needs the supplied words, not the unknown genome. For each word occurrence in $K$, compute its prefix and suffix and add the corresponding edge. An unknown ordering of strings becomes a traversal problem in a graph built directly from the unordered multiset.

For now we supply the ideal genomic composition, so using every edge once is exactly the reconstruction requirement. The same construction on read counts gives a graph of observations; its edge counts need an additional sampling model before they can specify genomic copies.

::: {#exm-spell-de-bruijn}
## Spelling a walk
The walk
$$\texttt{TA}\to\texttt{AC}\to\texttt{CG}\to\texttt{GA}\to\texttt{AC}\to\texttt{CG}$$
starts with `TA` and appends `CGACG`, giving `TACGACG`. Its five edges use every member of the ideal multiset.
:::

::: {#prp-walk-spelling}
## Walks spell overlapping k-mers
A walk of $q$ edges in $D(K)$ spells a string of length $q+k-1$. Its $i$th $k$-mer is the $i$th edge traversed, including multiplicity.
:::

::: {.proof}
Sketch. Start with the first node, which contributes $k-1$ symbols. An edge labeled $w$ has source $w_1\cdots w_{k-1}$ and target $w_2\cdots w_k$. Appending the target's last symbol therefore makes the latest $k$ symbols exactly $w$. Each of the $q$ edges appends one symbol, giving length $k-1+q$. For example, `TA → AC` appends C to TA and records the 3-mer TAC. Neighboring node labels overlap in $k-2$ symbols, while neighboring edge words overlap in $k-1$.
:::

Construction takes $\mathcal O(|K|)$ time for fixed $k$ and stores $|K|$ edges. Reconstruction now asks for a walk using every edge occurrence once. Lesson 7 supplies that algorithm.

## What changes with k

Return to the reads of `TACGACG`, but suppose the read `ACGA` was not sampled. The remaining reads are `TACG`, `CGAC`, and `GACG`. Every genomic base still appears in at least one read. Whether all words are represented depends on $k$:

| $k$ | Distinct words observed in the remaining reads | Missing genomic words |
|---:|---|---|
| 3 | `TAC, ACG, CGA, GAC` | none |
| 4 | `TACG, CGAC, GACG` | `ACGA` |

At $k=4$, the missing edge `ACG → CGA` breaks the route through the graph. At $k=3$, the reads still supply all distinct words needed for the genomic route, although their counts remain coverage-weighted. This example distinguishes base coverage from k-mer coverage: observing every position does not guarantee an observed window of every chosen length [@compeau2015].

For a DNA model with independent, equally likely bases, a particular k-mer has probability $4^{-k}$ at any position. This rough calculation explains why longer words are less likely to match by accident. Genomes are not random, and repeats violate the calculation badly, but it gives the direction of the tradeoff. Increasing k separates more repeated contexts while requiring longer error-free coverage.

We can compare graphs at several $k$ values rather than assume one choice is always adequate. Our theorem in Lesson 7 applies to each fixed graph. Deciding which graph best represents the reads is a separate statistical question.

:::: {.callout-note collapse="true" title="Going deeper: count words in a reference genome"}
The [phiX174 reference](../data/phix174-NC_001422.1.fasta) gives us a small real sequence on which to check the same counting code. Download the repository to run this block from the lesson folder. FASTA's first line names the record; its remaining lines hold one continuous sequence.

```{.python .execute}
#| label: code-reads-reference-kmers

fasta_path = Path("../../data/phix174-NC_001422.1.fasta")
lines = fasta_path.read_text().splitlines()
reference = "".join(line for line in lines if not line.startswith(">"))
counts = kmer_counts([reference], 3)
assert len(reference) == 5386
assert sum(counts.values()) == len(reference) - 2
len(reference), sum(counts.values()), len(counts)
```

The first two numbers are sequence length and total 3-mer occurrences; the third counts distinct words. The genome is circular, but this calculation treats its FASTA representation as linear and omits the two 3-mers crossing the end/start boundary. This is the genome's composition, not the coverage-weighted counts of a sequencing experiment.
::::

A repeated $(k-1)$-mer merges distinct genomic positions into one node. That merge can leave several complete edge orders compatible with the data. The next lesson distinguishes finding an order from establishing that only one sequence is possible.

## Exercises

::: {#exr-reads-composition}
List $\operatorname{Comp}_3(\texttt{ATATG})$. Then count 3-mers in reads `ATAT` and `TATG` and explain the different multiplicities.
:::

::: {#exr-reads-build-graph}
Build $D(K)$ for $K=\{\!\{\texttt{AAG},\texttt{AGA},\texttt{GAA},\texttt{AAT}\}\!\}$. Give every in-degree and out-degree.
:::

::: {#exr-reads-parallel}
Draw the graph for `ATATAT` at $k=3$, including parallel edges. What is lost if edges are stored as a set?
:::

::: {#exr-reads-implement}
Modify `de_bruijn_edges` to return an adjacency dictionary that retains repeated destinations. Test empty input.
:::

::: {#refs}
**References**
:::
