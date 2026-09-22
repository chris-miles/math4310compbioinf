---
published: true
title: "Sequences, Distances, and Alignments"
subtitle: "Lesson 1 · 2027-01-11 · Week 1"
---

## Core question

When are two DNA sequences similar, and what does the answer depend on?

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

In 1983 Russell Doolittle compared the protein sequence of *v-sis*, the cancer-causing gene of a monkey virus, against every protein sequence then known. It matched platelet-derived growth factor, a normal human protein that signals cells to divide [@doolittle1983; @waterfield1983]. The viral gene caused tumors because it was a growth signal that never switched off. No experiment had connected the two; the connection came from sequence similarity alone.

This is the daily routine of the field. A newly sequenced gene of unknown function is compared against databases of known sequences, and a strong match is taken as evidence of a shared ancestor and, usually, a shared function. The same comparison finds the mutation that separates a patient from a reference, locates a gene in a genome, and, in the phylogenetics lessons, reconstructs which species are related. All of it depends on two things: a precise definition of "similar," and a way to compute it fast enough to search a database. Today is about the first; Lessons 3 to 6 are about the second.

## Sequences as strings

::: {#def-alphabet}
## Alphabet and sequence
An *alphabet* $\Sigma$ is a finite set of symbols. A *sequence* over $\Sigma$ is a finite string $x = x_1 x_2 \ldots x_n$ with each $x_i \in \Sigma$; its length is $|x| = n$.
:::

For DNA, $\Sigma = \{A, C, G, T\}$; for proteins, the twenty amino-acid letters. The string keeps the order of the bases and discards everything else: the second strand, the three-dimensional structure, which stretches are genes.

::: {.biology}
A sequencer reads the order of bases along a DNA fragment and writes it as text, one letter per base. One run produces millions of fragments (reads), a few hundred letters each on the common instruments and tens of thousands on long-read machines, with errors, and with no read spanning a chromosome.
:::

The two strands of DNA pair $A$ with $T$ and $C$ with $G$ and run in opposite directions. The *reverse complement* of a string is the other strand's reading: the reverse complement of $GAG$ is $CTC$. A sequencer does not report which strand it read, so a string and its reverse complement are the same molecule.

## Hamming distance

```{.python .execute}
#| label: code-first-problem-strings
x = "ACGTT"
y = "ACGCT"
z = "CGTTA"
```

Is $y$ or $z$ closer to $x$? The question has no answer until "closer" is defined. The simplest definition counts disagreeing positions.

::: {#def-hamming-distance}
## Hamming distance
For sequences $x$ and $y$ with $|x| = |y| = n$,
$$
d_H(x, y) = \bigl|\{\, i \in \{1, \ldots, n\} : x_i \ne y_i \,\}\bigr| .
$$
:::

$d_H(x, y)$ is the number of substitutions that turn $x$ into $y$ when substitution is the only allowed change.

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

$z$ is $x$ with the first letter deleted and one letter appended: two mutation events. Hamming distance reports four because it assumes position $i$ of one string corresponds to position $i$ of the other; one deletion shifts every later position and breaks the assumption. The 4 is a fact about the assumption, not about the DNA.

## Edit distance

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

The proof has the two-part shape used for every optimization claim in the course: exhibit a solution (upper bound), then rule out anything smaller (lower bound).

## Counting alignments

$d_H$ costs one pass, $n$ comparisons. $d_E$ by the definition means a search over edit sequences, and the size of that search is the number of alignments.

::: {#def-alignment}
## Alignment
An *alignment* of $x$ and $y$ writes the two strings one above the other with gap symbols $-$ inserted so that every column is a letter over a letter, a letter over a gap, or a gap over a letter, and reading each row without its gaps gives back the original string.
:::

Each edit sequence is an alignment: a match or substitution is a letter over a letter, a deletion is a letter over a gap, an insertion is a gap over a letter. The alignments realizing $d_E(x, z) = 2$ and $d_E(x, y) = 1$:
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

Let $A(n, m)$ be the number of alignments of a length-$n$ string with a length-$m$ string. The last column is one of the three kinds, so
$$
A(n, m) = A(n-1, m-1) + A(n-1, m) + A(n, m-1), \qquad A(n, 0) = A(0, m) = 1 .
$$

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

Two five-letter strings have $1{,}683$ alignments; two of length 20 have more than $10^{14}$; two genes of ordinary length have more than there are atoms in the observable universe. Enumeration is not an option. Lesson 3 computes $d_E$ in $nm$ steps with the same recurrence, a minimum in place of the sum.

## k-mer composition

A third definition of similarity ignores positions entirely and compares the short words two strings contain.

::: {#def-kmer}
## $k$-mer and $k$-mer spectrum
A *$k$-mer* of $x$ is a substring $x_i x_{i+1} \ldots x_{i+k-1}$ of length $k$. The *$k$-mer spectrum* of $x$ is the multiset of its $n - k + 1$ $k$-mers, usually stored as a table of counts.
:::

```{.python .execute}
#| label: code-first-problem-kmers
def count_kmers(sequence: str, k: int) -> dict[str, int]:
    """Count the length-k substrings of a sequence."""
    counts: dict[str, int] = {}
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

The three definitions disagree because they measure different things: $d_H$ counts substitutions at fixed positions, $d_E$ counts events of any kind, and shared $k$-mers count preserved words wherever they sit. A single substitution destroys up to $k$ $k$-mers, so $k$-mer similarity is sensitive to substitutions and blind to shifts, the reverse of Hamming distance.

The cost is one pass over each string with a dictionary, $\mathcal{O}(n + m)$, against $\mathcal{O}(nm)$ for edit distance in Lesson 3. That gap is why $k$-mers are the first step of database search (Lesson 6) and the whole basis of genome assembly (Lesson 7).

## Limitations

1. **Hamming distance assumes no shifts.** One insertion or deletion makes $d_H$ large between strings that differ by a single event. It is the right distance only when positions are known to correspond, as for reads already placed at the same reference position.
2. **Edit distance assumes every edit costs one.** In real genomes an $A \leftrightarrow G$ substitution is more common than $A \leftrightarrow T$, and insertions and deletions are rarer than substitutions. Lesson 2 replaces the count with scores derived from how often each event occurs.
3. **$k$-mer composition ignores order.** Two different strings can have identical $k$-mer spectra (@exr-first-problem-kmers), and comparing only distinct shared words discards their multiplicities. Both facts are exactly what makes assembly hard (Lesson 8).
4. **All three ignore the other strand.** $GAG$ and $CTC$ are the same molecule. Any comparison of DNA strings has to decide whether to compare against the reverse complement as well (Lesson 7).

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
