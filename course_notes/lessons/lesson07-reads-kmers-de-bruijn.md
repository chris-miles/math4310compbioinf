---
published: true
title: "Reads, k-mers, and de Bruijn Graphs"
subtitle: "Lesson 7 · Week 4"
---

## Core question

How can a collection of short DNA fragments determine a much longer sequence?

```{.python .execute}
#| label: setup
#| code-fold: true
import sys
sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## Reads and coverage

::: {#def-read-multiset}
## Read multiset
Let $G$ be a genome string. An *ideal read* is a substring of $G$ or its reverse complement. Observed reads can differ from these substrings because of sequencing errors. Sequencing data are a multiset $R=\{\!\{r_1,\ldots,r_m\}\!\}$ because identical read sequences may be observed more than once.
:::

::: {.biology}
A short-read sequencer reports millions of fragments without their genomic positions. Sampling is uneven, and errors enter before assembly begins.
:::

If $m$ reads have length $L$ and the genome has length $n$, the nominal coverage is $C=mL/n$. Coverage is an average, not a guarantee at every position.

::: {#exm-read-coverage}
## Coverage counts observations
Let $G=\texttt{TACGACG}$ and take every length-4 substring:
$$R=\{\!\{\texttt{TACG},\texttt{ACGA},\texttt{CGAC},\texttt{GACG}\}\!\}.$$
Here $C=16/7$. End bases occur in one read and the central base occurs in four, even though starts are perfectly regular.
:::


### Read boundaries and overlap

Breaking reads into shorter words seems destructive. It forgets which k-mers came from the same read and where each k-mer sat inside that read. The benefit is that exact overlaps become easy to count and index. Instead of comparing every read against every other read, we use the k-mer itself as a key in a dictionary.

The choice is useful only when neighboring reads share k-mers. If two reads overlap in fewer than k bases, they contribute no word that certifies the overlap. If k is small, unrelated genomic positions share words by chance. If k is large, true neighbors fail to connect when coverage has a gap or an error falls inside the overlap. Thus k is part of the model rather than a harmless implementation setting.

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
For a length-$n$ string $x$,
$$\operatorname{Comp}_k(x)=\{\!\{x_i\cdots x_{i+k-1}:1\leq i\leq n-k+1\}\!\}.$$
This is a multiset: repeated occurrences remain repeated.
:::

The exact composition is an ideal object. Breaking every read into $k$-mers instead gives coverage-weighted observations. The two multiplicities differ.

::: {#exm-spectrum-multiplicity}
## Genome and read multiplicity

| 3-mer | `TAC` | `ACG` | `CGA` | `GAC` |
|---|---:|---:|---:|---:|
| copies in $G$ | 1 | 2 | 1 | 1 |
| observations in $R$ | 1 | 3 | 2 | 2 |

Coverage creates repeated observations of one genomic occurrence. An assembler must estimate genomic multiplicity from these counts.
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
by $(k-1)$-symbol overlaps spells `TACGACG`. Searching all orders hides the structure recorded by a graph.

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
    return [(word[:-1], word[1:]) for word in kmers]

ideal_kmers = ["TAC", "ACG", "CGA", "GAC", "ACG"]
edges = de_bruijn_edges(ideal_kmers)
edges
```

```{.python .execute}
#| label: fig-reads-de-bruijn-graph
#| fig-cap: "The de Bruijn graph for the ideal 3-mer composition of `TACGACG`. Two copies of `ACG` give two parallel edges."
#| fig-alt: "A directed graph with nodes TA, AC, CG, and GA, including two edges from AC to CG."
#| fig-width: 6
#| fig-height: 3
positions = {"TA": (0, 0), "AC": (1.6, 0), "CG": (3.2, 0), "GA": (2.4, -1.3)}
fig, ax = plt.subplots()
cf.digraph(ax, edges, pos=positions)
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

There are two useful ways to read the resulting graph. Locally, an edge says which $(k-1)$-mer can follow another after one new symbol. Globally, a walk chooses an ordering of k-mer occurrences consistent with all those local overlaps. The graph does not choose an ordering for us. It represents every ordering still compatible with the input.

The same construction works without knowing the genome. For each observed word, compute its prefix and suffix and add the corresponding edge. This is the key reduction: an unknown ordering of strings becomes a traversal problem in a graph built directly from the unordered multiset.

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
Sketch. Start with the first node. At each edge, the next node overlaps the current node in $k-2$ symbols, so appending its final symbol adds exactly the edge's $k$-mer.
:::

Construction takes $\mathcal O(|K|)$ time for fixed $k$ and stores $|K|$ edges. Reconstruction now asks for a walk using every edge occurrence once. Lesson 8 supplies that algorithm.


## What changes with k

Use the short genome `ACGTT`. At $k=2$, its edges are `A→C`, `C→G`, `G→T`, and `T→T`; one-letter nodes contain little context. At $k=4$, its edges are `ACG→CGT` and `CGT→GTT`; the overlap is much more specific, but there are only two chances to observe it.

For a random DNA model with equally likely bases, a particular k-mer has probability $4^{-k}$ at any position. This rough calculation explains why longer words are less likely to match by accident. Genomes are not random, and repeats violate the calculation badly, but it gives the direction of the tradeoff. Increasing k separates more repeated contexts while requiring longer error-free coverage.

We can compare graphs at several $k$ values rather than assume one choice is always adequate. Our theorem in Lesson 8 applies to each fixed graph. Deciding which graph best represents the reads is a separate statistical question.

## Limitations

1. **Exact composition assumes perfect coverage.** Missing $k$-mers remove edges; repeated sequencing adds observations without adding genomic copies.
2. **The toy model fixes one strand.** The reverse complement of `ACG` is `CGT`: reverse the order and exchange A/T and C/G. A canonical word chooses one representative of this pair, but a graph must also retain orientation to recover valid overlaps.
3. **The graph may lose order.** A repeated $(k-1)$-mer merges genomic positions and may permit several edge orders.
4. **Errors create edges.** One wrong base changes up to $k$ read-derived $k$-mers.

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

## Further reading

Compeau and Pevzner develop string reconstruction and de Bruijn graphs in Chapter 3 [@compeau2015]. Pevzner, Tang, and Waterman connect this representation to sequencing reads and errors [@pevzner2001].
