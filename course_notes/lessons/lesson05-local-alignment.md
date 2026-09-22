---
published: true
title: "Local Alignment"
subtitle: "Lesson 5 · 2027-01-27 · Week 3"
nocite: |
  @smith1981, @durbin1998, @compeau2015
---

## Core question

How do we find a strong shared region without paying for unrelated flanks?

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

## Substrings instead of whole strings

A protein can share a functional domain with another protein while having different sequence on either side. Comparing both proteins end to end can bury the shared region under penalties for unrelated flanks. We want the strongest matching pair of substrings, with their locations in the original records.

Our input remains two strings and column scores. The change is which letters the output must use: we may now leave prefixes and suffixes outside the alignment. The Smith-Waterman algorithm makes that change precise with one extra recurrence candidate and new traceback rules.

Global alignment uses every letter. Local alignment can isolate a shared protein domain. Mapping an entire read inside a chromosome often instead uses a free-end alignment: it consumes the read while leaving reference flanks uncharged. Which letters must appear in the answer is part of the problem specification.

::: {#def-local-alignment}
## Local alignment problem
Given strings $x$ and $y$ and additive column scores, find substrings of $x$ and $y$ whose alignment has maximum score. The empty alignment, with score zero, is allowed.
:::

Allowing the empty alignment means the answer is never negative. More importantly, it lets a poor prefix be discarded before a promising region begins.

Let $H(i,j)$ be the best score of an alignment of suffixes ending at positions $i$ and $j$, with zero available for an empty alignment. The last column can be a letter pair or a gap column.

## The Smith-Waterman recurrence

::: {#def-local-recurrence}
## Local alignment recurrence
With substitution score $s$ and linear gap score $g<0$,
$$
H(i,j)=\max\begin{cases}
0,\\
H(i-1,j-1)+s(x_i,y_j),\\
H(i-1,j)+g,\\
H(i,j-1)+g,
\end{cases}
$$
with $H(i,0)=H(0,j)=0$.
:::

The global and local algorithms differ in three places:

| Rule | Global | Local |
|---|---|---|
| Border | cumulative gap scores | zero |
| Answer | bottom-right cell | largest cell anywhere |
| Traceback stop | $(0,0)$ | first zero |

: Three changes from global to local alignment. {#tbl-global-local}

The zero is a fourth case: discard the current prefix and start a new alignment later.

::: {.algorithm}
**Algorithm: find one optimal local alignment.**

1. Set the first row and column of $H$ to zero.
2. Fill each cell with @def-local-recurrence and remember a maximizing predecessor.
3. Start traceback at a largest cell in the table.
4. Follow pointers until reaching a cell with value zero.
:::

## A shared interior region

Take $x=\texttt{TACG}$ and $y=\texttt{ACGT}$, with match $+1$, mismatch $-1$, and gap $-2$.

::: {#exm-local-reset}
## The effect of zero
At $H(1,1)$, T meets A:
$$
H(1,1)=\max\{0,-1,-2,-2\}=0.
$$
The mismatch is discarded. The diagonal run A/A, C/C, G/G then reaches $H(4,3)=3$.
:::

The best local alignment is
$$
\begin{array}{c}\texttt{ACG}\\\texttt{ACG}\end{array}
$$
with score 3. A global alignment must also account for the leading T in $x$ and trailing T in $y$; its best score is $-1$.

```{.python .execute}
#| label: fig-local-traceback
#| fig-cap: "Local-alignment table for TACG and ACGT, with the optimal traceback in red. The path begins at an interior maximum and stops at zero, leaving the terminal T symbols unaligned."
#| fig-alt: "A five-by-five table with zero borders and a red diagonal path from score three at row G, column G back to a zero at row T, empty column."
#| fig-width: 4.5
#| fig-height: 4.5
local_scores = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1],
    [0, 1, 0, 0, 0],
    [0, 0, 2, 0, 0],
    [0, 0, 0, 3, 1],
]
fig, ax = plt.subplots()
cf.dp_table(ax, list("TACG"), list("ACGT"), local_scores,
            path=[(1, 0), (2, 1), (3, 2), (4, 3)])
plt.show()
```

```{.python .execute}
#| label: code-local-alignment
def local_alignment(x: str, y: str, match: int = 1,
                    mismatch: int = -1, gap: int = -2) -> tuple[int, str, str]:
    """Return a local-alignment score and one optimal alignment."""
    table = [[0] * (len(y) + 1) for _ in range(len(x) + 1)]
    pointer = [["stop"] * (len(y) + 1) for _ in range(len(x) + 1)]
    best_score, best_cell = 0, (0, 0)

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            pair = match if x[i - 1] == y[j - 1] else mismatch
            candidates = [
                (0, "stop"),
                (table[i - 1][j - 1] + pair, "diag"),
                (table[i - 1][j] + gap, "up"),
                (table[i][j - 1] + gap, "left"),
            ]
            table[i][j], pointer[i][j] = max(candidates, key=lambda item: item[0])
            if table[i][j] > best_score:
                best_score, best_cell = table[i][j], (i, j)

    top, bottom = [], []
    i, j = best_cell
    while table[i][j] > 0:
        move = pointer[i][j]
        if move == "diag":
            top.append(x[i - 1])
            bottom.append(y[j - 1])
            i -= 1
            j -= 1
        elif move == "up":
            top.append(x[i - 1])
            bottom.append("-")
            i -= 1
        else:
            top.append("-")
            bottom.append(y[j - 1])
            j -= 1
    return best_score, "".join(reversed(top)), "".join(reversed(bottom))

local_alignment("TACG", "ACGT")
```

::: {#exm-local-no-match}
## No positive region
For $x=\texttt{AAA}$ and $y=\texttt{CCC}$ with these scores, every nonzero candidate is negative. The whole table is zero, so the reported score is 0 and the optimal local alignment is empty.
:::

### Local-alignment coordinates

Traceback identifies more than the aligned letters. If it starts at a maximum cell $(i_1,j_1)$ and stops at a zero cell $(i_0,j_0)$, then the selected substrings begin just after the stop coordinates and end at the start coordinates, respectively. For the TACG/ACGT example, traceback starts at $(4,3)$ and stops at $(1,0)$. It therefore selects $x_2\cdots x_4=\texttt{ACG}$ and $y_1\cdots y_3=\texttt{ACG}$.

This coordinate information matters in applications. A local protein alignment reports where a shared domain lies in each protein. A mapped read reports where its alignment begins and ends on the reference. Returning only the aligned text discards those locations.

::: {#exm-local-two-islands}
## Two separated matches
Compare $x=\texttt{AAACCCAAA}$ with $y=\texttt{AAAGGGAAA}$ using match $+1$, mismatch $-1$, and gap $-2$. Each AAA block scores 3. Joining the blocks requires crossing three mismatches, so the combined ungapped score is also 3. The table therefore has multiple optima: either block alone and the full nine-column alignment all tie. Making the mismatch score $-2$ prevents joining the blocks in an optimum. Each of the two AAA blocks in $x$ can match either AAA block in $y$, giving four tied pairs of locations.
:::

This example shows that "the best local region" need not have a unique location or length. A reporting program needs a tie rule, just as global traceback does.

## Correctness and cost

::: {#thm-local-correctness}
## Correctness of local alignment
$H(i,j)$ is the maximum score among the empty alignment and all alignments of a suffix of $x_1\cdots x_i$ with a suffix of $y_1\cdots y_j$ that end at position $(i,j)$. Therefore the largest table entry is the best score over every pair of substrings.
:::

::: {.proof}
Use induction on $i+j$. A nonempty alignment ending at $(i,j)$ has one of the same three final columns as a global alignment. Removing it leaves an alignment ending at the matching predecessor, so the three recurrence candidates exhaust the nonempty cases. The zero candidate represents choosing an empty alignment after the current prefixes. Every pair of substrings has some endpoint $(i,j)$, so maximizing over all cells considers them all.
:::

The algorithm still fills $(n+1)(m+1)$ cells and traces at most $n+m$ moves. Its time and memory costs are $\mathcal O(nm)$.

### Free-end alignment variants

The recurrence, borders, answer cell, and stopping rule act together. Start from the three-way global recurrence, without a zero candidate in interior cells. Zero borders and a required bottom-right answer allow free leading ends but charge trailing ends. Global borders and an answer taken from the last row or last column instead charge leading ends and allow free trailing ends. These are free-end variants; neither permits the arbitrary internal restart used by local alignment.

A program that reports coordinates should check them as well as the score. The teaching function returns only the score and aligned strings; it could also return the traceback's start and stop cells. On TACG/ACGT, an answer of 3 must select ACG/ACG. A program that reports the same score while including a terminal mismatch has inconsistent traceback.

## Scoring and region length

Local alignment is meaningful when unrelated sequence tends to decrease the score. If random columns have positive average score, extending through unrelated material is often rewarded and the selected region expands. A zero mismatch score does not by itself force the whole sequence into the answer, but it can make neutral columns free and produce long, poorly localized optima.

::: {.callout-warning title="Similarity is not evidence by itself"}
A high local score can arise in a low-complexity repeat such as ACACAC. Database tools mask or downweight such regions because the scoring model alone does not know that many similar-looking hits occur by chance.
:::

Smith-Waterman finds the best region under its score. It does not supply a probability that the sequences are homologous. Significance also depends on sequence lengths, composition, and how many database records were searched.

## Limitations

1. **Only the best region is returned.** Two disjoint domains require repeated searches or methods designed for multiple hits.
2. **The score controls length.** Weak mismatch or gap penalties allow unrelated columns to join nearby good regions.
3. **Low complexity creates strong accidental matches.** Repetitive composition violates the simple background model.
4. **Quadratic time remains.** Exact local alignment against every database record is too costly at database scale. Lesson 6 introduces seed-and-extend search.

## Exercises

::: {#exr-local-table}
Fill the local table for $x=\texttt{GAC}$ and $y=\texttt{TGA}$ using match $+2$, mismatch $-1$, and gap $-2$. Give the start and end coordinates of one best region.
:::

::: {#exr-local-global}
Construct two strings whose best local score is positive while their best global score is negative under match $+1$, mismatch $-1$, gap $-2$. Explain the flanks' effect.
:::

::: {#exr-local-zero}
Remove the zero candidate from @def-local-recurrence while leaving the zero borders. State which part of @thm-local-correctness fails.
:::

::: {#exr-local-tie}
Find an example with two disjoint local alignments tied for best score. What additional output should a program provide if both matter biologically?
:::

::: {#refs}
**References**
:::
