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
import numpy as np
import coursefigs as cf
cf.use_style()
```

## Reads and coverage

Lesson 1 placed each read under a reference genome and read down the columns. Suppose instead that we sequence a virus or bacterium for which no suitable reference exists. The instrument still returns reads, strings of a few hundred letters or fewer, but nothing records where each one came from and there is no reference to place them against. Reconstructing the genome from the reads alone is *de novo assembly*. A reconstructed genome lets us search for complete genes and see how they are arranged, even when individual reads cover only parts of them.

Compeau and Pevzner compare assembly to recovering a day's news from the confetti of many copies of one newspaper blown apart [@compeau2015]. Scraps are lost, so no single copy can be rebuilt, but the copies tore in different places, and a scrap from one copy can bridge a tear in another. A sequencing sample likewise holds many copies of the genome, one per cell or virus particle, each broken at different points (@fig-reads-torn-copies); the letters that overlapping reads share are the only clue to their order.

```{.python .execute}
#| label: fig-reads-torn-copies
#| code-fold: true
#| code-summary: "Draw copies of a genome broken at different places"
#| fig-cap: "Three copies of one genome, each broken at different places. The red piece from copy 2 spans a break in copy 1 (dotted line), so it overlaps the pieces on both sides of that break; an assembler uses overlaps like this to put the pieces back in order."
#| fig-alt: "Three rows labeled copy 1 to copy 3, each a black horizontal line split into three pieces at different positions. A dotted vertical line marks the first break in copy 1, and the red middle piece of copy 2 extends across it."
#| fig-width: 5
#| fig-height: 1.5
fig, ax = plt.subplots()
cf.torn_copies(ax, cuts=[(4.0, 7.5), (2.2, 5.6), (3.0, 8.4)], bridge=(1, 1))
plt.show()
```

To see what such data look like, we simulate them from a genome we already know. The bacteriophage phiX174, a virus that infects bacteria, has a circular genome of 5,386 bases; it was the first DNA genome ever sequenced, by Frederick Sanger's group in 1977 [@compeau2015], and its sequence is in the course [data folder](../data/phix174-NC_001422.1.fasta). The chunk below draws 162 synthetic reads of length 100 from it, with independent, uniformly random starts and no sequencing errors. To run it, download the repository and work from the lesson folder.

```{.python .execute}
#| label: code-reads-simulate-phix
fasta_path = Path("../../data/phix174-NC_001422.1.fasta")
lines = fasta_path.read_text().splitlines()
phix = "".join(line for line in lines if not line.startswith(">"))
n = len(phix)  # genome length: 5,386 bases
L = 100        # read length
m = 162        # number of reads
rng = np.random.default_rng(11)
# The genome is circular, so a read starting near the end wraps to the beginning.
circular = phix + phix[:L - 1]
starts = rng.integers(0, n, size=m)  # known to the simulation, hidden from the assembler
sim_reads = [circular[start:start + L] for start in starts]
depth = np.zeros(n, dtype=int)       # depth[i] = number of reads covering position i
for start in starts:
    depth[(start + np.arange(L)) % n] += 1
sim_reads[0]
```

That output is one complete record: a hundred letters and no coordinate. An assembler receives 162 such strings in no particular order. @fig-reads-phix-coverage uses the hidden starts to show what the strings leave out: where each read came from and how many reads cover each base.

```{.python .execute}
#| label: fig-reads-phix-coverage
#| code-fold: true
#| code-summary: "Draw the synthetic reads at their hidden positions"
#| fig-cap: "Synthetic, error-free reads from the phiX174 genome. Top: the 162 reads at their true positions, which the assembler never sees. Bottom: the number of reads covering each base, with the average C ≈ 3 dashed. Even under uniform sampling the depth fluctuates, and no read covers the red stretches, so no method using only these reads can recover those bases."
#| fig-alt: "Upper panel: short grey horizontal segments stacked in ten rows across positions 0 to 5,386. Lower panel: a blue step curve of read depth between 0 and 9 with a dashed horizontal line at 3. Several narrow red bands, aligned across both panels, mark positions where depth is zero."
#| fig-width: 6
#| fig-height: 3.6
fig, (top, bottom) = plt.subplots(2, 1, sharex=True,
                                  gridspec_kw={"height_ratios": [1.3, 1]})
for index, start in enumerate(np.sort(starts)):
    for left in (start, start - n):  # a wrapped read appears at both ends
        top.plot([left, left + L], [index % 10] * 2, color=cf.RULE, lw=1.5)
top.set(yticks=[], ylabel="Synthetic reads", xlim=(0, n))
top.spines["left"].set_visible(False)
bottom.plot(np.arange(n), depth, color=cf.SERIES[0], lw=1)
bottom.axhline(m * L / n, color=cf.INK, ls="--", lw=0.8)
bottom.set(xlabel="Genome position (bases)", ylabel="Depth")
for ax in (top, bottom):
    ax.fill_between(np.arange(n), 0, 1, where=depth == 0, color=cf.EMPHASIS,
                    alpha=0.3, lw=0, transform=ax.get_xaxis_transform())
plt.show()
```

Depth is uneven, with some bases read eight or nine times and a few not at all, so even a perfect method may recover the genome only in separate pieces.

To trace the reasoning by hand we need a much smaller genome. The reads `TACG` and `ACGA` suggest the longer string `TACGA`, because the final ACG of the first read matches the initial ACG of the second. @fig-reads-records-to-overlaps adds two more reads. On the left are the records as an assembler receives them, named r1 to r4 in arbitrary order; on the right is a placement in which every overlap agrees, spelling the candidate sequence `TACGACG`.

```{.python .execute}
#| label: fig-reads-records-to-overlaps
#| code-fold: true
#| code-summary: "Compare unordered read records with proposed overlapping placements"
#| fig-cap: "The same four read records listed independently (left) and placed beneath candidate sequence TACGACG (right). Their overlaps supply the proposed order, and several reads observe the same genomic bases."
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

A shared word supplies a possible join, but a repeated word elsewhere in the genome can supply the same apparent join. Joining reads greedily, by the first overlap found, can then take a wrong branch that shows up only much later, as reads left over that fit nowhere [@compeau2015]. We need a representation that keeps the alternatives. We first isolate a simpler problem: find one string whose short words, counted with multiplicity, match a supplied collection.

Let $G$ be the genome string. An *ideal read* is a substring of $G$ or of its reverse complement; observed reads can differ from these substrings because of sequencing errors. The data form a multiset $R=\{\!\{r_1,\ldots,r_m\}\!\}$, an unordered collection that retains how often each item occurs, because the same read sequence may be observed more than once.

Because fragments are drawn from many copies of the genome, a genomic position can be read many times or not at all. If $m$ reads have length $L$ and the genome has length $n$, the *coverage*
$$
C=\frac{mL}{n}
$$
is the total number of sequenced bases divided by the number of genomic bases, which is the average number of reads covering a position. In @fig-reads-phix-coverage, $C=162\cdot100/5386\approx3.0$, the dashed line.

An average of three does not mean every base is read three times. How much of the genome should we expect to miss? Model each read start as uniform on the circular genome and independent of the other reads, as in the simulation. A fixed position $i$ is covered by a read exactly when the read starts at one of the $L$ positions $i-L+1,\ldots,i$, so a single read misses $i$ with probability $1-L/n$. Independence lets us multiply this probability over the $m$ reads:
$$
\Pr(\text{no read covers position } i)=\Bigl(1-\frac{L}{n}\Bigr)^{m}\approx e^{-mL/n}=e^{-C}.
$$
The approximation uses $\ln(1-L/n)\approx-L/n$, which is accurate when reads are much shorter than the genome. Every position has the same probability of being missed, so $e^{-C}$ is also the expected fraction of the genome that no read covers. Counting the uncovered positions in the simulation checks the formula:

```{.python .execute}
#| label: code-reads-coverage-check
C = m * L / n
fraction_uncovered = np.mean(depth == 0)
round(C, 2), round(float(fraction_uncovered), 3), round(float(np.exp(-C)), 3)
```

The simulated reads leave 4.3% of phiX174 unread, close to the predicted 4.9%. At $C=3$, roughly one base in twenty lies in a gap, and each gap splits the genome into pieces that no overlap can join. The expected number of unread bases is $ne^{-C}$, which falls below one only when $C>\ln n$: about 8.6 for phiX174 and about 15 for a five-million-base bacterial genome. Real libraries sample some regions less often than this uniform model assumes, which raises the coverage needed further.

For the toy examples, all reads have already been oriented to the same strand; real reads are not, so an assembler must handle each read and its reverse complement consistently [@compeau2015].

::: {#exm-read-coverage}
## Coverage counts observations
Let $G=\texttt{TACGACG}$ and take every length-4 substring:
$$R=\{\!\{\texttt{TACG},\texttt{ACGA},\texttt{CGAC},\texttt{GACG}\}\!\}.$$
Four reads of length 4 on a genome of length 7 give $C=16/7\approx2.3$. Reading down the columns of @fig-reads-records-to-overlaps, the depth along the genome is $1,2,3,4,3,2,1$: even with perfectly regular starts, the end bases are observed once and the central base four times.
:::

### Read boundaries and overlap

Breaking reads into shorter words forgets some useful order: which words came from the same read and how far apart they were. In exchange it exposes short overlaps shared by reads with different endpoints, and those overlaps let us join evidence from many fragments into one graph.

Two reads need not share a k-mer to connect in the graph we will build. `TACG` and `CGAC` from @exm-read-coverage share the two-base flank `CG` and no 3-mer, yet at $k=3$ the last word of the first read, `ACG`, can be followed by the first word of the second, `CGA`, because their two-base suffix and prefix agree. The flank `CG` will be the node where those two edges meet.

For the reads in @exm-read-coverage, $k=3$ retains the chain of overlaps:

| Read | 3-mers |
|---|---|
| `TACG` | `TAC`, `ACG` |
| `ACGA` | `ACG`, `CGA` |
| `CGAC` | `CGA`, `GAC` |
| `GACG` | `GAC`, `ACG` |

::: {.callout-warning title="Two multisets with different meanings"}
The genome contains `CGA` once, at positions 3 through 5, but reads `ACGA` and `CGAC` both contain it, so counting words in the reads gives two. The multiset of genomic k-mer occurrences answers “how many copies occur in the genome?” The multiset extracted from reads answers “how many times did the sequencer observe this word?” Assembly begins with the second and tries to infer the first.
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

The genome supplies $7-3+1=5$ occurrences of 3-mers, while the reads supply
$$
\sum_{j=1}^{m}\max\bigl(0,\;|r_j|-k+1\bigr)=4\cdot(4-3+1)=8
$$
observations, since a read shorter than $k$ contributes no words.
:::

The total of 8 is a useful check on a counting program: a different total means it has omitted windows or invented words across read boundaries.

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

The most direct graph, the *overlap graph* on the left of @fig-reads-de-bruijn-graph, puts each word occurrence at a node, with an arrow from $u$ to $v$ when the last $k-1$ symbols of $u$ equal the first $k-1$ symbols of $v$. The two occurrences of `ACG` need separate nodes because a reconstruction must use both. A reconstruction is then a path visiting every node exactly once, a *Hamiltonian path*, shown in red. No local rule says which node comes next: from `TAC` both ACG nodes are possible, and only a complete path shows which choices work. No efficient algorithm is known for finding Hamiltonian paths in general graphs [@compeau2015].

The de Bruijn construction on the right moves each word onto an edge, running from the word's prefix to its suffix, the $(k-1)$-mers that neighboring words must share. To picture it, write the genome as a path of five edges, $\texttt{TA}\to\texttt{AC}\to\texttt{CG}\to\texttt{GA}\to\texttt{AC}\to\texttt{CG}$, each labeled by the 3-mer it spells, and glue together nodes with the same label, so the two AC nodes become one and the two CG nodes become one [@compeau2015]. Repeated overlaps collapse while each word occurrence keeps its own edge, and the genome's path survives the gluing as a walk that may revisit nodes but uses every edge exactly once. Visiting every node once has become using every edge once. Whether an edge-covering walk exists can be decided by counting edges at each node, and one can then be built in linear time, as the next lesson shows; no comparable method is known for node-covering paths.

```{.python .execute}
#| label: fig-reads-de-bruijn-graph
#| code-fold: true
#| code-summary: "Draw the overlap graph and the de Bruijn graph for the same words"
#| fig-cap: "The five 3-mer occurrences of `TACGACG` drawn two ways. Left: the overlap graph, one node per occurrence, with the Hamiltonian path that spells `TACGACG` in red. Right: the de Bruijn graph, one edge per occurrence from its 2-mer prefix to its 2-mer suffix, with two parallel edges for ACG. Moving words from nodes to edges turns “visit every word once” into “use every edge once.”"
#| fig-alt: "Left: five nodes TAC, ACG1, ACG2, CGA, GAC with seven arrows; a red path runs TAC, ACG1, CGA, GAC, ACG2. Right: four nodes TA, AC, CG, GA with edges labeled TAC, ACG times 2 (two parallel arcs from AC to CG), CGA, and GAC."
#| fig-width: 6
#| fig-height: 3
# an arrow wherever the last two letters of one word are the first two of another
overlap_edges = [("TAC", "ACG₁"), ("TAC", "ACG₂"), ("ACG₁", "CGA"), ("ACG₂", "CGA"),
                 ("CGA", "GAC"), ("GAC", "ACG₁"), ("GAC", "ACG₂")]
hamiltonian = [("TAC", "ACG₁"), ("ACG₁", "CGA"), ("CGA", "GAC"), ("GAC", "ACG₂")]
overlap_positions = {"TAC": (0, 0), "ACG₁": (1.3, 0.9), "ACG₂": (1.3, -0.9),
                     "CGA": (2.6, 0.9), "GAC": (2.6, -0.9)}
word_edges = [("TA", "AC"), ("AC", "CG"), ("CG", "GA"), ("GA", "AC"), ("AC", "CG")]
word_positions = {"TA": (0, 0), "AC": (1.6, 0), "CG": (3.2, 0), "GA": (2.4, -1.2)}
word_labels = {("TA", "AC", 0): "TAC", ("AC", "CG", 0): "ACG ×2",
               ("CG", "GA", 0): "CGA", ("GA", "AC", 0): "GAC"}
fig, axes = plt.subplots(1, 2)
cf.digraph(axes[0], overlap_edges, pos=overlap_positions,
           highlight_edges=hamiltonian)
cf.digraph(axes[1], word_edges, pos=word_positions, edge_labels=word_labels)
axes[0].set_title("Overlap graph: words on nodes")
axes[1].set_title("de Bruijn graph: words on edges")
axes[1].set_ylim(-1.75, 0.55)  # leave room for the lowest node
plt.show()
```

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

### Constructing the graph by hand

Each 3-mer splits into two overlapping 2-mers, giving the edges in the right panel of @fig-reads-de-bruijn-graph:

| Edge label | Prefix node | Suffix node |
|---|---|---|
| `TAC` | `TA` | `AC` |
| `ACG` | `AC` | `CG` |
| `CGA` | `CG` | `GA` |
| `GAC` | `GA` | `AC` |
| `ACG` | `AC` | `CG` |

The two `ACG` rows represent two positions in the genome, so they become parallel edges rather than duplicates to discard. If we collapsed them, the graph would contain four edges and could spell only a length-6 string, while a five-edge walk spells the required length $5+3-1=7$.

Locally, an edge says which $(k-1)$-mer can follow another after one new symbol. Globally, a walk that uses every edge occurrence exactly once is an ordering of the words compatible with the composition. The gluing picture started from the genome, but the construction does not need it: gluing identical labels in five separate edges, one per word, gives the same graph, and that is what the algorithm above does.

::: {.callout-note title="Key idea"}
Placing each $k$-mer occurrence on an edge from its prefix to its suffix turns “in what order do the words occur?” into “which walk uses every edge exactly once?” The graph is built from the unordered words alone.
:::

::: {#exm-spell-de-bruijn}
## Spelling a walk
The walk
$$\texttt{TA}\to\texttt{AC}\to\texttt{CG}\to\texttt{GA}\to\texttt{AC}\to\texttt{CG}$$
starts with `TA` and appends `CGACG`, giving `TACGACG`. Its five edges use every member of the ideal multiset.
:::

The example generalizes: a walk of $q$ edges in $D(K)$ spells a string of length $q+k-1$, and the $i$th $k$-mer of that string is the label of the $i$th edge traversed, repeats included.

::: {.proof}
Sketch. An edge labeled $w=w_1\cdots w_k$ joins the two $(k-1)$-mers inside $w$:
$$
\underbrace{w_1\cdots w_{k-1}}_{\text{source}}\;\longrightarrow\;\underbrace{w_2\cdots w_k}_{\text{target}}.
$$
Spelling starts with the $k-1$ letters of the first node, and each edge appends the last letter of its target. Just after edge $w$, the newest $k$ letters are its source followed by $w_k$, which is $w$ itself; for example, `TA → AC` appends C to TA and records TAC. Counting letters,
$$
\underbrace{(k-1)}_{\text{first node}}+\underbrace{1+\cdots+1}_{q\text{ edges}}=q+k-1.
$$
:::

Construction takes $\mathcal O(|K|)$ time for fixed $k$ and stores $|K|$ edges. What remains is to find a walk that uses every edge occurrence once, which is the subject of the next lesson.

## What changes with k

Return to the reads of `TACGACG`, but suppose the read `ACGA` was not sampled. The remaining reads are `TACG`, `CGAC`, and `GACG`. Every genomic base still appears in at least one read. Whether all words are represented depends on $k$:

| $k$ | Distinct words observed in the remaining reads | Missing genomic words |
|---:|---|---|
| 3 | `TAC, ACG, CGA, GAC` | none |
| 4 | `TACG, CGAC, GACG` | `ACGA` |

At $k=4$, the missing edge `ACG → CGA` breaks the route through the graph. At $k=3$, the reads still supply every distinct word on the genomic route. Observing every base does not guarantee observing every window of a chosen length, which is why assemblers break reads into words shorter than the reads: shorter words are more likely to be fully represented in the reads [@compeau2015].

Why not simply take $k$ large? Longer words are less likely to repeat by accident. In a model DNA sequence with independent, equally likely bases, a particular $k$-mer appears at a given position with probability $4^{-k}$. Adding this probability over the $n-k+1$ starting positions gives its expected number of occurrences in a random genome of length $n$:
$$
(n-k+1)\,4^{-k}.
$$
For $n=5386$, the length of phiX174, this is about 1.3 at $k=6$ and about 0.02 at $k=9$. Short words recur throughout a genome by chance, so their $(k-1)$-mer nodes merge unrelated positions, while words longer than about $\log_4 n$ rarely recur unless the genome contains a real repeat. Real genomes are far from random: about half of the human genome consists of repeats, including roughly a million copies of the 300-base Alu element, similar enough that an assembler cannot tell many of them apart [@compeau2015]. The calculation still gives the direction of the tradeoff: increasing $k$ separates more repeated contexts while requiring longer stretches of error-free coverage, as the missing `ACGA` showed.

:::: {.callout-note collapse="true" title="Going deeper: repeated words in phiX174"}
To compare the random-genome calculation with a real sequence, count the distinct $k$-mers that occur more than once in the phiX174 sequence `phix` loaded at the start of the lesson.

```{.python .execute}
#| label: code-reads-reference-kmers
repeated = {}
for k in [3, 6, 9, 12, 13]:
    k_counts = kmer_counts([phix], k)
    repeated[k] = sum(1 for count in k_counts.values() if count > 1)
repeated
```

All 64 possible 3-mers occur more than once, but the number of repeated words falls quickly once $k$ passes $\log_4 5386\approx6.2$, and no 13-mer occurs twice. The count treats the circular genome's FASTA sequence as linear, so it omits words that cross the end/start boundary. Larger genomes also contain biological repeats far longer than any practical $k$ [@compeau2015].
::::

A repeated $(k-1)$-mer merges distinct genomic positions into one node. That merge can leave several complete edge orders compatible with the data. The next lesson distinguishes finding an order from establishing that only one sequence is possible.

## Exercises

::: {#exr-reads-composition}
List $\operatorname{Comp}_3(\texttt{ATATG})$. Then count 3-mers in reads `ATAT` and `TATG` and explain the different multiplicities.
:::

::: {#exr-reads-coverage}
A bacterial genome has about 5,000,000 bases and reads have length 150. Under the uniform sampling model, how many reads are needed before the expected number of bases covered by no read falls below one? What coverage $C$ does that correspond to?
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
