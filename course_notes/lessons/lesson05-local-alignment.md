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

After translating a predicted coding sequence, we may compare the resulting amino-acid string with a characterized protein. A *domain* is a region of a protein that often folds as a unit and carries a particular activity, such as binding another molecule [@ebi-domains]. Two proteins may share a domain while having different regions around it. Requiring an end-to-end alignment can bury the domain match under penalties for those flanks; local alignment lets the shared region determine its own endpoints.

Our input remains two strings and column scores. The change is which letters the output must use: we may now leave prefixes and suffixes outside the alignment. The Smith-Waterman algorithm makes that change precise with one extra recurrence candidate and new traceback rules.

Which letters must appear in the answer depends on the biological task:

| Input and intended comparison | Letters required | Alignment choice |
|---|---|---|
| Two complete versions of a corresponding gene | both full records | global |
| Two proteins sharing a domain | the shared regions | local |
| A read placed within a reference chromosome | the whole read, a reference interval | fitting (free reference ends) |

: Biological tasks determine alignment boundaries. {#tbl-alignment-biological-tasks}

The fitting model assumes the whole read belongs in one reference interval. A read with an adapter (a laboratory-added sequence used during sequencing) or an unreliable end may need trimming or an alignment that leaves those bases unaligned. We will keep the simpler whole-read requirement when deriving the fitting recurrence.

::: {#def-local-alignment}
## Local alignment problem
Given strings $x$ and $y$ and additive column scores, find substrings of $x$ and $y$ whose alignment has maximum score. The empty alignment, with score zero, is allowed.
:::

Allowing the empty alignment means the answer is never negative. More importantly, it lets a poor prefix be discarded before a promising region begins.

At grid coordinate $(i,j)$, we have consumed the prefixes $x_1\cdots x_i$ and $y_1\cdots y_j$. Let $H(i,j)$ be the best alignment score among suffixes of those prefixes, allowing a suffix to be empty. A selected region may therefore begin after any earlier coordinate, but must end at $(i,j)$. Zero is always available by choosing both suffixes empty. As before, a nonempty alignment ends in a letter pair or a gap column.

## The Smith-Waterman recurrence

::: {#def-local-recurrence}
## Local alignment recurrence
For $i,j\geq1$, with substitution score $s$ and linear gap score $g<0$,
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

The global and local algorithms differ in four connected rules:

| Rule | Global | Local |
|---|---|---|
| Interior recurrence | three predecessor extensions | those extensions plus zero |
| Border | cumulative gap scores | zero |
| Answer | bottom-right cell | largest cell anywhere |
| Traceback stop | $(0,0)$ | first zero |

: Four changes from global to local alignment. {#tbl-global-local}

The zero candidate places an empty alignment at the current coordinate. If every extension loses score, starting there costs less than retaining the poor region. Later cells can extend this restart by pairing letters, without paying for the discarded prefixes. The largest cell selects the other endpoint.

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
with score 3. This is an independent check on the expected output: ACG supplies three matches, and there cannot be four matches in order because the terminal T symbols lie on opposite sides of ACG. Any mismatch or gap decreases the score, so 3 is an upper bound as well as an attained score. A global alignment must also account for the leading T in $x$ and trailing T in $y$; its best score is $-1$.

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

The function returns the score, the two gapped strings, and the start and end grid coordinates of the selected region. If the latter are $(i_0,j_0)$ and $(i_1,j_1)$, the ungapped strings must equal Python slices `x[i_0:i_1]` and `y[j_0:j_1]`. These are zero-based half-open intervals: the start is included and the end excluded. In our one-based mathematical notation, the first selected letters are $x_{i_0+1}$ and $y_{j_0+1}$. An empty result uses start and end $(0,0)$.

```{.python .execute}
#| label: code-local-alignment
def local_alignment(x: str, y: str, match: int = 1,
                    mismatch: int = -1, gap: int = -2
                    ) -> tuple[int, str, str, tuple[int, int], tuple[int, int]]:
    """Return an optimal local score, alignment, and region boundary coordinates."""
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
    aligned_x = "".join(reversed(top))
    aligned_y = "".join(reversed(bottom))
    return best_score, aligned_x, aligned_y, (i, j), best_cell

x, y = "TACG", "ACGT"
result = local_alignment(x, y)
score, aligned_x, aligned_y, start, end = result
assert aligned_x.replace("-", "") == x[start[0]:end[0]]
assert aligned_y.replace("-", "") == y[start[1]:end[1]]
result
```

The output score 3 agrees with the hand bound. Python's `max` keeps the first candidate in a tie, so placing `stop` first discards a zero-score prefix. Updating `best_cell` only for a strictly larger score keeps the first maximizing endpoint in row-major order. These conventions choose one optimum; they do not enumerate every tied region.

::: {#exm-local-no-match}
## No positive region
For $x=\texttt{AAA}$ and $y=\texttt{CCC}$ with these scores, every nonzero candidate is negative. The whole table is zero, so the reported score is 0 and the optimal local alignment is empty.
:::

### Local-alignment coordinates

Traceback identifies more than the aligned letters. If it starts at a maximum cell $(i_1,j_1)$ and stops at a zero cell $(i_0,j_0)$, then the selected substrings begin just after the stop coordinates and end at the start coordinates, respectively. For the TACG/ACGT example, traceback starts at $(4,3)$ and stops at $(1,0)$. It therefore selects $x_2\cdots x_4=\texttt{ACG}$ and $y_1\cdots y_3=\texttt{ACG}$.

This coordinate information matters in applications. A local protein alignment proposes where a shared domain lies in each protein, and a read alignment identifies a reference interval. The returned boundaries let us verify the selected substrings as well as the score.

::: {#exm-local-two-islands}
## Two separated matches
Compare $x=\texttt{AAACCCAAA}$ with $y=\texttt{AAAGGGAAA}$ using match $+1$, mismatch $-1$, and gap $-2$. Each AAA block scores 3. Joining the corresponding blocks crosses three mismatches, so the full ungapped score is $3-3+3=3$ as well. All four pairs of AAA locations and the full nine-column alignment attain the same optimum.

Along the diagonal, scores rise to $H(3,3)=3$, fall to $H(6,6)=0$, and rise again to $H(9,9)=3$. Tracing from $(9,9)$ and stopping at its first zero returns the final block, removing the neutral prefix AAACCC/AAAGGG. The full nine-column alignment is optimal under the problem's definition, but this stop convention does not report it. Our function instead chooses the first maximum, $(3,3)$, and returns the first block.

Making the mismatch score $-2$ prevents the full alignment from tying. Now the diagonal score reaches zero before the mismatches finish, and the final block begins from a fresh restart. Each AAA block can still match either block in the other string, giving four tied location pairs.
:::

This example shows that "the best local region" need not have a unique location or length. A reporting program needs a tie rule, just as global traceback does.

## Correctness and cost

::: {#thm-local-correctness}
## Correctness of local alignment
$H(i,j)$ is the maximum score among all alignments of suffixes of $x_1\cdots x_i$ and $y_1\cdots y_j$, including empty suffixes. Therefore the largest table entry is the best score over every pair of substrings.
:::

::: {.proof}
Use induction on $i+j$. On a boundary, any nonempty alignment consists only of negative-score gap columns, so the empty alignment attains the maximum zero.

For an interior cell, a nonempty optimal alignment has one of the three final columns. Removing it leaves suffixes at the corresponding predecessor. Their alignment must be optimal there: replacing it with a better one and appending the same final column would improve the original. Thus its score is bounded by one predecessor candidate. The empty case has score zero. Conversely, each candidate is attained by extending an optimal predecessor alignment, or by choosing the empty alignment. Their maximum is therefore exactly $H(i,j)$.

Every pair of substrings has an endpoint $(i,j)$, so the largest cell considers them all. A negative-score initial stretch can be removed to increase an alignment's score; a zero-score initial stretch can be removed without changing it. Consequently some optimum starts at a zero cell, and stopping traceback at the first zero recovers one optimum even when it omits a tied neutral prefix.
:::

The algorithm fills $(n+1)(m+1)$ cells with constant work per cell and traces at most $n+m$ moves. For positive $n,m$ and constant-time arithmetic, its time and memory costs are $\mathcal O(nm)$.

### Free-end alignment variants

The recurrence, borders, answer cell, and stopping rule act together. Start from the three-way global recurrence, without a zero candidate in interior cells. Zero borders and a required bottom-right answer allow free leading ends but charge trailing ends. Global borders and an answer taken from the last row or last column instead charge leading ends and allow free trailing ends. These are free-end variants; neither permits the arbitrary internal restart used by local alignment.

### Fitting an entire read inside a reference

Suppose $x$ is a read of length $n$ and $y$ is a longer reference of length $m$. We want every read base represented, while reference bases before and after the placement are free. Let $F(i,j)$ use the global three-way recurrence, with rows for read prefixes and columns for reference prefixes. Choose
$$
F(0,j)=0,\qquad F(i,0)=ig,\qquad
\text{answer}=\max_{0\leq j\leq m}F(n,j).
$$
The zero top row lets us skip a reference prefix. The penalized first column prevents skipping a read prefix for free. Taking the answer anywhere in the last row leaves the unused reference suffix uncharged but requires consuming the whole read. Traceback stops on the top row, where $i=0$; it does not stop at an arbitrary interior zero.

For read $\texttt{AC}$ and reference $\texttt{TTACGG}$, using match $+1$, mismatch $-1$, and gap $-2$, the fitting alignment scores 2 and places the read at reference positions 3 and 4. Two is also an upper bound: each of the two read letters can contribute at most one match reward, and gaps have negative scores. Global alignment instead must account for four additional reference letters; local alignment is allowed to omit part of the read.

These boundary choices express a biological requirement in the table. Before modifying an alignment routine, state which input letters the output must consume, then choose its initialization and endpoint to enforce that requirement.

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
