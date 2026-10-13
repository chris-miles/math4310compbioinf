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

Suppose we want to find DNA differences between a bacterial sample and a previously studied bacterium. *DNA sequencing* estimates the order of bases in a DNA molecule. An instrument measures signals along DNA fragments, and software converts those signals into ordered A, C, G, T letters, called *reads*. Each read records part of a DNA molecule without its position in the genome. A *reference sequence* is a previously recorded genome that supplies candidate locations for comparison.

Suppose a read is `ACGTT` and one candidate reference segment is `ACTT`. Pairing G with a gap lets the final TT bases line up (@fig-first-problem-read-reference). We need a rule for deciding which letters correspond and how much their differences count. This lesson compares rules that allow different kinds of changes.

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

The gap proposes an extra G in this read relative to the candidate segment. Other reads and their measurement quality help determine whether that difference is real. First we need to make the comparison itself precise.

## Sequences as strings

::: {#def-alphabet}
## Alphabet and sequence
An *alphabet* $\Sigma$ is a finite set of symbols. A *sequence* over $\Sigma$ is a finite string $x = x_1 x_2 \ldots x_n$ with each $x_i \in \Sigma$; its length is $|x| = n$.
:::

For DNA, $\Sigma = \{A, C, G, T\}$; for proteins, the twenty amino-acid letters. The string keeps the order of the bases and discards everything else: the second strand, the three-dimensional structure, which stretches are genes.

The conversion from instrument signals to letters is called *basecalling*. Some instruments measure light as a copied strand grows, and others measure electrical current as DNA passes through a pore [@illumina-sbs; @nanopore-technology]. Both give an inferred base string. Exact matching to a reference can fail because of a true sequence difference or a basecalling error. An *alignment* proposes a correspondence while allowing mismatches and gaps; we define it precisely after comparing two ways to count changes.

The read/reference comparison in @fig-first-problem-read-reference is one candidate pairing. Before choosing among candidates, we need a precise rule for counting their differences.

DNA reads can come from either strand. Read mapping therefore also considers the *reverse complement*, formed by reversing the string and exchanging A with T and C with G; for example, GAG becomes CTC. Our examples assume the orientation has already been chosen.

## Hamming distance

```{.python .execute}
#| label: code-first-problem-strings
x = "ACGTT"
y = "ACGCT"
z = "CGTTA"
```

Suppose $x$ is a read and $y$ and $z$ are candidate reference segments in the same orientation. Which segment is closer to the read? The answer depends on the allowed differences. If the bases are already known to correspond position by position, counting disagreeing positions gives a useful first comparison.

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

$z$ can be obtained from $x$ by deleting the first A and appending A at the end. The surviving CGTT stays in order, but moves one position to the left. Fixed-position comparison now pairs C with A, G with C, and T with G. One deletion can therefore create many mismatched positions without changing those surviving bases. In @fig-first-problem-shift, gaps restore their correspondence.

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

Every allowed edit costs one in this definition. The minimum does not distinguish a likely substitution from an unlikely one, or a sequencing error from a true variant. Lesson 2 derives letter-pair scores by comparing related-sequence and background models; Lesson 5 treats gap runs separately.

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

The proof uses a useful two-part pattern for establishing a minimum: exhibit a solution (upper bound), then rule out anything smaller (lower bound).

### Bounds before computation

Every insertion or deletion changes the length by one, while a substitution leaves it unchanged. Therefore $d_E(x,y)\geq\bigl||x|-|y|\bigr|$. A proposed edit distance smaller than the length difference is impossible, regardless of the letters. This bound can certify an answer when a construction attains it.

For example, changing $\texttt{ACGTT}$ to $\texttt{ACT}$ requires at least two edits because the lengths differ by two. Deleting G and one T achieves two, so the distance is exactly two. For equal-length strings the length bound is zero and may tell us very little; the earlier example needed a separate argument to exclude one substitution.

## Counting alignments

$d_H$ uses $n$ comparisons because it fixes the correspondence in advance. Edit distance allows many correspondences. An alignment records one of them and makes its cost visible: count the mismatched letter pairs and the gap columns. We first count the possible alignments to see why searching each separately is impractical.

::: {#def-alignment}
## Alignment
An *alignment* of $x$ and $y$ writes the two strings one above the other with gap symbols $-$ inserted so that every column is a letter over a letter, a letter over a gap, or a gap over a letter, and reading each row without its gaps gives back the original string.
:::

Reading an alignment from left to right gives a transformation: preserve a match, substitute a mismatch, delete a top-row letter opposite a gap, or insert a bottom-row letter opposite a gap. Thus every alignment supplies an edit sequence with the same cost.

The converse needs a little care. In a shortest edit sequence, inserted letters are chosen in their final form, deleted originals need no substitutions, and surviving originals change at most once. Otherwise an operation could be removed or several substitutions replaced by one. Surviving original letters also retain their order. Pair each with its final letter, put deleted originals opposite gaps, and put inserted letters opposite gaps in the other row. This constructs an alignment with the same edit cost. Minimizing over alignments therefore gives exactly $d_E$.

Here are alignments realizing $d_E(x, z) = 2$ and $d_E(x, y) = 1$:
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

Let $A(n, m)$ count alignments of strings of lengths $n$ and $m$, with $n,m\geq0$. This counts gap placements, regardless of whether the letters match. For $n,m\geq1$, partition the alignments by their last column. A letter pair consumes one letter from each string; a letter over a gap consumes only a letter from $x$; a gap over a letter consumes only a letter from $y$. Removing the last column leaves an alignment of the corresponding shorter prefixes. Conversely, appending that column extends every such shorter alignment in exactly one way. The three cases are disjoint and exhaustive, giving
$$
A(n, m) = A(n-1, m-1) + A(n-1, m) + A(n, m-1), \qquad A(n, 0) = A(0, m) = 1 .
$$

If one string is empty, there is exactly one alignment: every letter of the other string is paired with a gap. In particular, $A(0,0)=1$ counts the empty alignment. For two one-letter strings, the recurrence gives $A(1,1)=1+1+1=3$: pair the letters, put the $x$ letter before the $y$ letter in separate gap columns, or reverse that order. At the shaded cell in @fig-first-problem-count-table, $A(2,2)=3+5+5=13$.

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

The first output is 3, agreeing with our one-letter enumeration. Two five-letter strings already have $1{,}683$ alignments; two of length 20 have more than $10^{14}$. Computing a count does not enumerate those objects: the table has only $(n+1)(m+1)$ entries and uses three previous counts per interior entry. Lesson 3 uses the same prefix decomposition to maximize an alignment score, adding a score for each final column instead of counting each extension. It computes that optimum in $\mathcal{O}(nm)$ time under a unit-cost arithmetic assumption; an optional extension applies the table to edit distance.

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

For fixed $k$ and expected constant-time dictionary operations, this comparison costs $\mathcal{O}(n+m)$. If $k$ varies, copying and hashing each Python substring costs $\mathcal{O}(k)$, giving $\mathcal{O}(k(n+m))$. Short-word matches can narrow a search before alignment. A shared word alone does not establish common ancestry or a unique genomic placement.

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
