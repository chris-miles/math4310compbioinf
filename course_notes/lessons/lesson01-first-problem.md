---
published: true
title: "Sequences, Distances, and Alignments"
subtitle: "Lesson 1 · 2027-01-11 · Week 1"
nocite: |
  @compeau2015, @durbin1998
---

## Core question

How can we compare a DNA fragment with a known sequence when their letters differ?

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

## Why compare sequences

Suppose we want to find the DNA differences between a bacterial sample and a previously studied strain of the same species. *DNA sequencing* estimates the order of bases in DNA, but not in one piece. The sample's DNA is broken into many fragments, an instrument measures a signal along each fragment (light as a copied strand grows, or electrical current as DNA passes through a pore), and *basecalling* software converts that signal into a string of A, C, G, and T called a *read* [@illumina-sbs; @nanopore-technology]. A short-read instrument reports millions of reads, each a few hundred letters or fewer, and no read says where in the genome it came from. What we do have is a *reference sequence*: the previously assembled genome of the studied strain.

The plan is to place each read at the stretch of reference it resembles most and then read down each column. @fig-first-problem-pileup shows the result for a few made-up reads. Where every read covering a position disagrees with the reference in the same way, the sample probably carries a real difference. A disagreement that appears in one read and is contradicted by its neighbors is more likely a measurement error. Both conclusions rest on a step the figure takes for granted: deciding where each read belongs and which of its letters correspond to which reference letters.

```{.python .execute}
#| label: fig-first-problem-pileup
#| code-fold: true
#| code-summary: "Draw the made-up reads under the reference"
#| fig-cap: "Six made-up reads placed under a reference segment. Grey letters agree with the reference and bold letters do not. All four reads covering the red column report T where the reference has C, so the sample probably differs there; the lone bold G in the last read is contradicted by the read above it and is more likely a sequencing error."
#| fig-alt: "A reference string of 22 bold letters above six staggered rows of read letters. Most read letters are grey. Four reads have a red-filled T in the same column; one read has a single bold G elsewhere."
#| fig-width: 6
#| fig-height: 2.6
reference = "TTGACGTTAGCCATGCAAGTCA"
reads = ["TTGACGTTA", "ACGTTAGTC", "GTTAGTCAT", "TAGTCATGC", "GTCATGCAA", "TGCGAGTCA"]
starts = [0, 3, 5, 7, 9, 13]
fig, ax = plt.subplots()
cf.read_pileup(ax, reference, reads, starts, highlight_column=10)
plt.show()
```

Exact matching is too strict for placing reads, because the reads we most want to place are the ones that differ from the reference. Suppose a read is `ACGTT` and one candidate reference segment is `ACTT`. No position-by-position comparison makes them agree, yet pairing G with a gap lets the final TT bases line up (@fig-first-problem-read-reference). We need a rule for deciding which letters correspond and how much their differences count. This lesson compares rules that allow different kinds of changes.

```{.python .execute}
#| label: fig-first-problem-read-reference
#| code-fold: true
#| code-summary: "Draw the read/reference correspondence"
#| fig-cap: "A proposed correspondence between the toy read ACGTT and reference segment ACTT. The gap preserves the pairing of the final TT bases and identifies an extra G in the read relative to this segment."
#| fig-alt: "ACGTT in the upper row and AC-TT in the lower row, with matching A, C, T, and T columns and a gap opposite G."
#| fig-width: 4
#| fig-height: 1.6
fig, ax = plt.subplots()
cf.alignment(ax, "ACGTT", "AC-TT", labels=("read", "ref"))
plt.show()
```

The gap proposes an extra G in this read relative to the candidate segment. Whether the sample really carries that G is the column-by-column question of @fig-first-problem-pileup, answered with many reads. First we need to make the comparison of one read with one segment precise.

## Sequences as strings

::: {#def-alphabet}
## Alphabet and sequence
An *alphabet* $\Sigma$ is a finite set of symbols. A *sequence* over $\Sigma$ is a finite string $x = x_1 x_2 \ldots x_n$ with each $x_i \in \Sigma$; its length is $|x| = n$.
:::

For DNA, $\Sigma = \{A, C, G, T\}$; for proteins, the twenty amino-acid letters. The string keeps the order of the bases and discards everything else: the second strand, the three-dimensional structure, which stretches are genes.

A read can come from either strand, so read mapping also tries its *reverse complement*: reverse the string and exchange A with T and C with G, so GAG becomes CTC. Our examples assume the orientation has already been chosen.

## Hamming distance

To compare one read with one candidate segment, start with the simplest rule. Let $x=\texttt{ACGTT}$ be a read, and consider two candidate reference segments in the same orientation. The first, $y=\texttt{ACGCT}$, differs from $x$ at one base, as a true variant or a basecalling error would make it. The second, $z=\texttt{CGTTA}$, shares the four letters CGTT with $x$ but has them one position earlier, as happens when the read carries an extra base at its start relative to this segment. Which segment is closer to the read? The answer depends on which differences we allow. If the bases are already known to correspond position by position, counting disagreeing positions gives a first comparison.

```{.python .execute}
#| label: code-first-problem-strings
x = "ACGTT"
y = "ACGCT"
z = "CGTTA"
```

::: {#def-hamming-distance}
## Hamming distance
For sequences $x$ and $y$ with $|x| = |y| = n$,
$$
d_H(x, y) = \bigl|\{\, i \in \{1, \ldots, n\} : x_i \ne y_i \,\}\bigr| .
$$
:::

$d_H(x, y)$ is the number of substitutions that turn $x$ into $y$ when substitution is the only allowed change.

::: {.algorithm}
**Algorithm: count fixed-position mismatches.**
Require equal-length strings $x$ and $y$. Start a counter at zero. For each position $i$, add one if $x_i\ne y_i$. Return the counter.
:::

```{.python .execute}
#| label: code-first-problem-hamming
def hamming_distance(x: str, y: str) -> int:
    """Count positions at which two equal-length strings differ."""
    if len(x) != len(y):
        raise ValueError("Hamming distance needs equal-length strings.")
    distance = 0
    for a, b in zip(x, y):
        if a != b:
            distance += 1
    return distance


hamming_distance(x, y), hamming_distance(x, z)
```

::: {#exm-hamming}
## Hamming distance on $x$, $y$, $z$
$$
\begin{array}{ll}
x = \texttt{A C G T T} \qquad & x = \texttt{A C G T T} \\
y = \texttt{A C G C T} \qquad & z = \texttt{C G T T A} \\
\phantom{y = }\texttt{. . . * .} \qquad & \phantom{z = }\texttt{* * * . *}
\end{array}
$$
$d_H(x, y) = 1$ and $d_H(x, z) = 4$.
:::

Hamming distance already gives a way to place the reads of @fig-first-problem-pileup. Slide a read along the reference, compute $d_H$ between the read and the reference letters beneath it at each start position, and choose the start with the smallest value.

```{.python .execute}
#| label: code-first-problem-hamming-scan
def hamming_scan(read: str, reference: str) -> list[int]:
    """Hamming distance from the read to the reference window at each start."""
    scores = []
    for start in range(len(reference) - len(read) + 1):
        window = reference[start:start + len(read)]
        scores.append(hamming_distance(read, window))
    return scores


substitution_read = "ACGTTAGTC"   # second read of the pileup
deletion_read = "ACGTAGCCA"       # reference from position 4 with one T deleted
hamming_scan(substitution_read, reference), hamming_scan(deletion_read, reference)
```

For the second read, `ACGTTAGTC`, the scan has a sharp minimum: one mismatch at position 4, where @fig-first-problem-pileup places it, and at least three everywhere else. Now delete one T from the reference segment that begins at position 4. The result, `ACGTAGCCA`, is the read a sample with a one-base deletion would produce. Its best Hamming score is 3, at position 5, one place to the right of where the read belongs; the correct start scores 4 (@fig-first-problem-offset-scan). One deletion explains this read completely, but fixed-position comparison cannot see it.

```{.python .execute}
#| label: fig-first-problem-offset-scan
#| code-fold: true
#| code-summary: "Plot the two scans"
#| fig-cap: "Hamming distance between a read and the reference window at each start position. The read with one substitution (blue) has a sharp minimum at its true start, position 4 (dotted line). The read with one deleted base (orange) has no good placement: its best score, 3, falls at position 5, and the true start scores 4."
#| fig-alt: "Two jagged curves of Hamming distance against start positions 1 to 14. The blue curve drops to 1 at position 4 and stays at 3 or more elsewhere. The orange curve has its lowest value, 3, at position 5 and a value of 4 at position 4. A dotted vertical line marks position 4."
positions_sub = range(1, len(reference) - len(substitution_read) + 2)
positions_del = range(1, len(reference) - len(deletion_read) + 2)
fig, ax = plt.subplots()
ax.plot(positions_sub, hamming_scan(substitution_read, reference), "o-",
        label="one substitution: " + substitution_read)
ax.plot(positions_del, hamming_scan(deletion_read, reference), "s-",
        label="one deletion: " + deletion_read)
ax.axvline(4, color=cf.RULE, linestyle=":", linewidth=1)
ax.set(xlabel="Start position in the reference", ylabel="Hamming distance $d_H$",
       xticks=[1, 4, 7, 10, 13], ylim=(0, 10))
ax.legend(frameon=False, loc="lower right")
plt.show()
```

The same failure appears on the small example. $z$ can be obtained from $x$ by deleting the first A and appending A at the end. The surviving CGTT stays in order, but moves one position to the left. Fixed-position comparison now pairs C with A, G with C, and T with G. One deletion can therefore create many mismatched positions without changing those surviving bases. In @fig-first-problem-shift, gaps restore their correspondence.

```{.python .execute}
#| label: fig-first-problem-shift
#| fig-cap: "The read x and candidate z paired at fixed positions (left) and with a shifted correspondence (right). Four fixed-position mismatches become four matches and two gap columns."
#| fig-alt: "Two pairs of sequence rows: ACGTT over CGTTA has four mismatch dots; ACGTT- over -CGTTA pairs the middle CGTT letters."
#| fig-width: 6
#| fig-height: 2
fig, axes = plt.subplots(1, 2)
cf.alignment(axes[0], x, z, labels=("x", "z"))
cf.alignment(axes[1], x + "-", "-" + z, labels=("x", "z"))
plt.show()
```

## Edit distance

To compare strings whose corresponding bases may have shifted, we allow insertions and deletions as well as substitutions. We then minimize over possible transformations, because the observed strings do not tell us which changes produced their differences.

::: {#def-edit-distance}
## Edit distance
An *edit* is a substitution of one symbol for another, an insertion of a symbol at any position, or a deletion of a symbol. The *edit distance* $d_E(x, y)$ is the minimum number of edits that transforms $x$ into $y$.
:::

::: {#exm-edit-distance}
## Edit distance on $x$, $y$, $z$

| Pair | $d_H$ | $d_E$ | Edits achieving $d_E$ |
|---|---:|---:|---|
| $x, y$ | 1 | 1 | substitute $T \to C$ at position 4 |
| $x, z$ | 4 | 2 | delete $x_1$; append $A$ |
| $y, z$ | 5 | 3 | delete $y_1$; substitute $C \to T$ at position 3; append $A$ |

:::

::: {#prp-edit-vs-hamming}
## Edit distance is at most Hamming distance
For sequences of equal length, $d_E(x, y) \le d_H(x, y)$, and the inequality can be strict.
:::

::: {.proof}
Sketch. Every substitution is an edit, so the $d_H(x, y)$ substitutions form a valid edit sequence and the minimum is no larger. For strictness take $x = ACGTT$, $z = CGTTA$: $d_H = 4$ and two edits suffice. Two are necessary: zero edits would mean $x = z$; a single substitution changes one position, but $d_H = 4$; a single insertion or deletion changes the length.
:::

The proof uses a two-part pattern for establishing a minimum: exhibit a solution (upper bound), then rule out anything smaller (lower bound). One lower bound comes free. Every insertion or deletion changes the length by one and a substitution does not, so
$$
d_E(x,y)\geq\bigl||x|-|y|\bigr|.
$$
Changing $\texttt{ACGTT}$ to $\texttt{ACT}$ therefore needs at least two edits, and deleting G and one T achieves two, so the distance is exactly two. For equal-length strings this bound is zero, which is why the proof above needed its own argument to exclude a single edit.

Unit costs treat every change alike: a likely substitution and an unlikely one, a sequencing error and a true variant. Lesson 2 asks how much a particular pair of letters should count.

## Counting alignments

Because $d_H$ fixes the correspondence in advance, it needs only $n$ comparisons. Edit distance instead ranges over every correspondence, and an alignment is the record of one of them, with its cost visible as mismatched letter pairs and gap columns.

::: {#def-alignment}
## Alignment
An *alignment* of $x$ and $y$ writes the two strings one above the other with gap symbols $-$ inserted so that every column is a letter over a letter, a letter over a gap, or a gap over a letter, and reading each row without its gaps gives back the original string.
:::

Write $\operatorname{cost}(\mathcal A)$ for the number of mismatched letter pairs plus the number of gap columns in an alignment $\mathcal A$.

::: {#prp-edit-alignment}
## Edit distance is a minimum over alignments
For any sequences $x$ and $y$,
$$
d_E(x,y)=\min_{\mathcal A}\ \operatorname{cost}(\mathcal A),
$$
where the minimum runs over every alignment $\mathcal A$ of $x$ and $y$.
:::

::: {.proof}
Sketch. Reading an alignment from left to right gives an edit sequence of the same cost: keep a match, substitute a mismatch, delete a top-row letter opposite a gap, and insert a bottom-row letter opposite a gap. Hence $d_E(x,y)$ is at most the minimum cost. Conversely, take a shortest edit sequence. In it, inserted letters are chosen in their final form, deleted letters are never substituted first, and each surviving letter changes at most once; otherwise an operation could be removed. Surviving letters keep their order, so pairing each with its final letter and placing deleted and inserted letters opposite gaps gives an alignment of the same cost.
:::

On our example, the edit sequence "delete $x_1$; append A" becomes the left alignment below: the deleted A sits opposite a gap, the surviving CGTT pairs with itself, and the appended A sits opposite a gap in the top row. The single substitution turning $x$ into $y$ becomes one mismatched column on the right.
$$
\begin{array}{c}
\texttt{A C G T T -} \\
\texttt{- C G T T A}
\end{array}
\qquad\qquad
\begin{array}{c}
\texttt{A C G T T} \\
\texttt{A C G C T}
\end{array}
$$

Edit distance is therefore a search over alignments, and the size of that search decides whether we can carry it out directly. How many alignments are there? If the number were small, we could simply try them all.

Let $A(n, m)$ count alignments of strings of lengths $n$ and $m$, with $n,m\geq0$. This counts gap placements, regardless of whether the letters match. For $n,m\geq1$, partition the alignments by their last column. A letter pair consumes one letter from each string; a letter over a gap consumes only a letter from $x$; a gap over a letter consumes only a letter from $y$. Removing the last column leaves an alignment of the corresponding shorter prefixes. Conversely, appending that column extends every such shorter alignment in exactly one way. The three cases are disjoint and exhaustive, giving
$$
A(n, m) = A(n-1, m-1) + A(n-1, m) + A(n, m-1), \qquad A(n, 0) = A(0, m) = 1 .
$$

If one string is empty, there is exactly one alignment: every letter of the other string is paired with a gap. In particular, $A(0,0)=1$ counts the empty alignment. For two one-letter strings, the recurrence gives $A(1,1)=1+1+1=3$: pair the letters, put the $x$ letter before the $y$ letter in separate gap columns, or reverse that order. At the shaded cell in @fig-first-problem-count-table, $A(2,2)=3+5+5=13$.

::: {.callout-note title="Key idea"}
Every alignment ends in one of three kinds of column, and removing that column leaves an alignment of shorter prefixes. A table indexed by prefix lengths can therefore summarize all alignments without listing any of them.
:::

::: {.algorithm}
**Algorithm: count alignments by prefix lengths.**
Create a table indexed by $i=0,\ldots,n$ and $j=0,\ldots,m$. Set its boundary row and column to one. Fill the remaining entries row by row, adding the entries above, to the left, and diagonally above-left. Return $A(n,m)$.
:::

The three required entries have already been filled in this order. Storing their counts avoids reconstructing all the alignments they represent.

```{.python .execute}
#| label: code-first-problem-count-alignments
def count_alignments(n: int, m: int) -> int:
    """Number of alignments of a length-n string with a length-m string."""
    # table[i][j] holds A(i, j); the boundary row and column are all ones
    table = [[1] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            table[i][j] = table[i - 1][j - 1] + table[i - 1][j] + table[i][j - 1]
    return table[n][m]


[count_alignments(k, k) for k in range(1, 9)]
```

```{.python .execute}
#| label: fig-first-problem-count-table
#| fig-cap: "The table `count_alignments` fills for two three-letter strings. Each entry counts the alignments of the prefixes ending at that row and column; the arrows show the three entries the shaded cell is computed from. Lesson 3 fills a table of this shape with a different rule."
#| fig-alt: "A four-by-four grid of integers with a boundary row and column of ones; three red arrows point into the shaded middle cell from its upper, left, and upper-left neighbors."
#| fig-width: 4
#| fig-height: 4
table = [[count_alignments(i, j) for j in range(4)] for i in range(4)]

fig, ax = plt.subplots()
cf.dp_table(ax, ["$x_1$", "$x_2$", "$x_3$"], ["$y_1$", "$y_2$", "$y_3$"], table,
            highlight=[(2, 2)],
            arrows=[((1, 1), (2, 2)), ((1, 2), (2, 2)), ((2, 1), (2, 2))])
plt.show()
```

The first output is 3, agreeing with our one-letter enumeration. Two five-letter strings already have $1{,}683$ alignments; two of length 20 have more than $10^{14}$. Computing a count does not enumerate those objects: the table has only $(n+1)(m+1)$ entries and uses three previous counts per interior entry. @fig-first-problem-count-growth puts the two numbers side by side.

```{.python .execute}
#| label: fig-first-problem-count-growth
#| code-fold: true
#| code-summary: "Plot alignment counts against table size"
#| fig-cap: "The number of alignments of two length-n strings (red) and the number of table entries used to count them (blue), on a logarithmic axis. The count grows exponentially while the table grows only quadratically, so any method that examines alignments one at a time is hopeless for reads of realistic length."
#| fig-alt: "Two curves against string length from 1 to 20 on a log vertical axis: a red curve rising steeply to above ten to the fourteenth, and a blue curve rising slowly to about four hundred."
lengths = list(range(1, 21))
alignment_counts = [count_alignments(n, n) for n in lengths]
table_entries = [(n + 1) ** 2 for n in lengths]
fig, ax = plt.subplots()
ax.semilogy(lengths, alignment_counts, "o-", color=cf.EMPHASIS, label="alignments $A(n,n)$")
ax.semilogy(lengths, table_entries, "s-", label="table entries $(n+1)^2$")
ax.set(xlabel="String length $n$", ylabel="Count", xticks=[1, 5, 10, 15, 20])
ax.legend(frameon=False)
plt.show()
```

Lesson 3 uses the same prefix decomposition to maximize an alignment score, adding a score for each final column instead of counting each extension. It computes that optimum in $\mathcal{O}(nm)$ time under a unit-cost arithmetic assumption; an optional extension applies the table to edit distance.

## k-mer composition

Suppose a long reference contains many possible places for a read. Before computing an alignment at each place, we can ask which segments contain short words also present in the read. A preserved word can move to a new position after an insertion or deletion, so comparing words offers a different notion of similarity. This comparison is useful for selecting candidates; it does not yet identify which individual bases correspond.

::: {#def-kmer}
## $k$-mer and $k$-mer spectrum
For an integer $1\leq k\leq n$, a *$k$-mer* of $x$ is a substring $x_i x_{i+1} \ldots x_{i+k-1}$ of length $k$. The *$k$-mer spectrum* of $x$ is the multiset of its $n - k + 1$ $k$-mers, usually stored as a table of counts.
:::

For $k=3$, the windows of $x=\texttt{ACGTT}$ start at positions 1, 2, and 3: ACG, CGT, and GTT. They overlap because advancing a window by one preserves its last two letters. The code stores multiplicities, then compares only the distinct words present in both strings.

::: {.algorithm}
**Algorithm: count and compare short words.**
For a chosen positive $k$, slide a length-$k$ window along each string and increment the count of each word encountered. Return the number of distinct words with a positive count in both strings. A string shorter than $k$ contributes no words.
:::

```{.python .execute}
#| label: code-first-problem-kmers
def count_kmers(sequence: str, k: int) -> dict[str, int]:
    """Count the length-k substrings of a sequence."""
    if k < 1:
        raise ValueError("k must be positive.")
    counts = {}
    for start in range(len(sequence) - k + 1):
        kmer = sequence[start:start + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def shared_kmers(x: str, y: str, k: int) -> int:
    """Number of distinct k-mers that occur in both strings."""
    return len(count_kmers(x, k).keys() & count_kmers(y, k).keys())


shared_kmers(x, y, 3), shared_kmers(x, z, 3)
```

::: {#exm-shared-kmers}
## Shared 3-mers of $x$, $y$, $z$

| String | 3-mers | Shared with $x$ |
|---|---|---:|
| $x = ACGTT$ | ACG, CGT, GTT | |
| $y = ACGCT$ | ACG, CGC, GCT | 1 |
| $z = CGTTA$ | CGT, GTT, TTA | 2 |

By shared 3-mers, $z$ is closer to $x$ than $y$ is. Hamming distance and edit distance said the opposite.
:::

The three definitions disagree because they measure different things: $d_H$ counts mismatched positions, $d_E$ minimizes the number of allowed edits, and shared $k$-mers count preserved words wherever they sit. In $y$, changing position 4 destroys CGT and GTT, leaving only ACG shared with $x$. In $z$, the shift preserves CGT and GTT, although ACG is lost and TTA is gained at the boundaries. Thus a shift can preserve most words without preserving their positions. A substitution affects up to $k$ windows, depending on its position; choosing $k$ changes how much local agreement a shared word demands.

For fixed $k$, with constant-time dictionary operations, the comparison takes $\mathcal{O}(n+m)$ time: one pass over each string. That is cheap enough to screen a long reference for candidate places before computing any alignment.

The spectrum discards the order of its words, so distinct strings can have identical spectra (@exr-first-problem-kmers). Counting only distinct shared words also discards multiplicities. These losses matter when reconstructing a sequence from its pieces: the available words need not specify one unique reconstruction.

## Exercises

::: {#exr-first-problem-distances}
Let $u = GATTACA$ and $v = ATTACAG$. Compute $d_H(u, v)$ and $d_E(u, v)$ by hand, and prove your value of $d_E$ is minimal using the upper-bound, lower-bound pattern from @prp-edit-vs-hamming.
:::

::: {#exr-first-problem-construct}
Construct two strings of length 6 with $d_H = 6$ and $d_E = 2$. State in one sentence which mutation events your construction imitates.
:::

::: {#exr-first-problem-count}
Show from the recurrence that $A(n, 1) = 2n + 1$, and describe the $2n + 1$ alignments.
:::

::: {#exr-first-problem-kmers}
Find two different strings of length 5 with the same 2-mer spectrum. Then decide whether two different strings of length 5 can have the same 4-mer spectrum, and explain.
:::

::: {#exr-first-problem-object}
Name one kind of biological data that is not naturally a single string. Say what mathematical object you would use for it and one thing that object throws away.
:::


::: {#refs}
**References**
:::
