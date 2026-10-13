---
published: true
title: "Local Alignment"
subtitle: "Lesson 4 · 2027-01-25 · Week 3"
nocite: |
  @smith1981, @durbin1998, @compeau2015
---

## Core question

Which regions of two proteins have similar sequences?

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

Suppose a DNA assembly contains a predicted gene whose protein function is unknown. Reading its coding sequence in three-base groups gives the predicted protein's amino-acid string, which we can compare with a protein whose function has been studied. A *domain* is a protein region that often folds as a unit and carries a particular activity, such as binding another molecule [@ebi-domains]. Finding a shared domain can suggest an activity to test, even when the proteins have different regions around it.

Global alignment requires every letter in both records. Its penalties for those surrounding regions can outweigh the evidence in the shared domain. For a short calculation, use DNA strings $x=\texttt{TTACGCC}$ and $y=\texttt{GGACGAA}$ with match $+1$, mismatch $-1$, and gap $-2$. The whole-string alignment scores $-1$, while the interior ACG/ACG scores 3. The same issue arises with amino-acid strings: the useful correspondence may cover only part of each record.

```{.python .execute}
#| label: fig-local-whole-versus-region
#| code-fold: true
#| code-summary: "Compare the whole strings with their shared region"
#| fig-cap: "The toy DNA strings paired end to end (left) and restricted to their shared ACG region (right). Four flank mismatches cancel the full alignment's three match rewards; selecting the interior retains the shared-region score of 3."
#| fig-alt: "TTACGCC over GGACGAA has mismatches at the first two and last two positions. A second pair of rows contains only ACG over ACG, with three match marks."
#| fig-width: 6
#| fig-height: 1.8
fig, axes = plt.subplots(1, 2)
cf.alignment(axes[0], "TTACGCC", "GGACGAA")
cf.alignment(axes[1], "ACG", "ACG")
plt.show()
```

Our input remains two strings and column scores. We now seek both the best correspondence and the start and end of its region in each string, allowing unused prefixes and suffixes. This is *local alignment*. The Smith-Waterman algorithm lets a promising region start after a poor prefix and end before poor trailing columns; its table finds those endpoints along with the score.

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

The start of the shared region is unknown. At a fixed endpoint $(i,j)$, its part in $x$ could begin at any position up to $i$, and likewise in $y$. Those choices are precisely the suffixes of $x_1\cdots x_i$ and $y_1\cdots y_j$. Let $H(i,j)$ store the best alignment score over all such suffix pairs, including the empty alignment. One entry therefore summarizes many possible starts. Taking the largest entry will also choose the endpoint.

Recall the global traceback: choose a maximizing predecessor, emit its column, and reverse the collected columns. Here the correspondence must also tell us where a shared region starts and ends. The restart candidate, best endpoint, and stopping rule change together.

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

The zero candidate means "discard the alignment built so far and start after this coordinate." Suppose a prefix has score $-2$ and a later run would add 3. Keeping the prefix gives 1; restarting gives 3. Because a local alignment is allowed to omit that prefix, keeping its loss can never help. A restart replaces all negative alternatives by zero, while a positive extension keeps accumulating evidence.

::: {.algorithm}
**Algorithm: find one optimal local alignment.**

1. Set the first row and column of $H$ to zero.
2. Fill each cell with @def-local-recurrence and remember a maximizing predecessor.
3. Start traceback at a largest cell in the table.
4. Follow pointers until reaching a cell with value zero.
:::

## A shared interior region

Use a toy DNA comparison with a shared core and differing flanks: $x=\texttt{TTACGCC}$ and $y=\texttt{GGACGAA}$, with match $+1$, mismatch $-1$, and gap $-2$. The letters make the calculation short; the same boundary choice lets a protein-domain match omit the surrounding protein regions.

::: {#exm-local-reset}
## Discarding a poor prefix
At $H(1,1)$, T meets G:
$$
H(1,1)=\max\{0,-1,-2,-2\}=0.
$$
The next diagonal T/G mismatch is also discarded. Pairing the following A/A, C/C, and G/G gives $H(3,3)=1$, $H(4,4)=2$, and $H(5,5)=3$. Continuing through the two C/A mismatches lowers the diagonal scores to 2 and then 1. The largest score, rather than the final cell, selects the end of the core.
:::

The local alignment is ACG/ACG with score 3. The global alignment TTACGCC/GGACGAA scores $-2+3-2=-1$: its four flank mismatches cancel more than the core earns. Local alignment omits the flanks and identifies a region worth investigating. Its coordinates will tell us where that region lies in both input records.

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

x, y = "TTACGCC", "GGACGAA"
result = local_alignment(x, y)
score, aligned_x, aligned_y, start, end = result
assert aligned_x.replace("-", "") == x[start[0]:end[0]]
assert aligned_y.replace("-", "") == y[start[1]:end[1]]
result
```

The output is `(3, "ACG", "ACG", (2, 2), (5, 5))`, agreeing with the hand calculation. Placing `stop` first in a score tie discards a neutral prefix; updating `best_cell` only when the score increases chooses the first maximizing endpoint.


```{.python .execute}
#| label: fig-local-traceback
#| fig-cap: "Computed local-alignment table for TTACGCC and GGACGAA, with the ACG traceback in red. The zero restart drops the leading mismatches, and the interior maximum excludes the trailing mismatches."
#| fig-alt: "An eight-by-eight table with zero borders. Red diagonal cells run from a zero at (2,2) through scores one, two, and three at (5,5)."
#| fig-width: 6
#| fig-height: 5.5
local_scores = [[0] * (len(y) + 1) for _ in range(len(x) + 1)]
for i in range(1, len(x) + 1):
    for j in range(1, len(y) + 1):
        pair = 1 if x[i - 1] == y[j - 1] else -1
        local_scores[i][j] = max(0, local_scores[i - 1][j - 1] + pair,
                                 local_scores[i - 1][j] - 2,
                                 local_scores[i][j - 1] - 2)
fig, ax = plt.subplots()
cf.dp_table(ax, list(x), list(y), local_scores,
            path=[(2, 2), (3, 3), (4, 4), (5, 5)])
plt.show()
```

::: {#exm-local-no-match}
## No positive region
For $x=\texttt{AAA}$ and $y=\texttt{CCC}$ with these scores, every nonzero candidate is negative. The whole table is zero, so the reported score is 0 and the optimal local alignment is empty.
:::

### Local-alignment coordinates

Traceback identifies more than the aligned letters. If it starts at a maximum cell $(i_1,j_1)$ and stops at a zero cell $(i_0,j_0)$, then the selected substrings begin just after the stop coordinates and end at the start coordinates, respectively. For TTACGCC/GGACGAA, traceback starts at $(5,5)$ and stops at $(2,2)$. It selects $x_3\cdots x_5=y_3\cdots y_5=\texttt{ACG}$.

This coordinate information matters in applications. A local protein alignment proposes where a shared domain lies in each protein, and a read alignment identifies a reference interval. The returned boundaries let us verify the selected substrings as well as the score.

::: {#exm-local-two-islands}
## Two separated matches
For AAACCCAAA and AAAGGGAAA, each AAA block scores 3. Joining the corresponding blocks crosses three $-1$ mismatches, so the full ungapped alignment also scores $3-3+3=3$. The diagonal table reaches 3, drops to zero, and reaches 3 again. Stopping at that zero returns the final block; retaining the neutral prefix would give a longer tied optimum. Our code selects the first maximum and returns the first block. A best score therefore need not specify a unique region or length.
:::

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

A length-$n$ string has $n(n+1)/2$ nonempty substrings, so separately aligning every substring pair would consider on the order of $n^2m^2$ pairs. Restarting within a shared table avoids those separate searches. The algorithm fills $(n+1)(m+1)$ cells with constant work per cell and traces at most $n+m$ moves. For positive $n,m$ and constant-time arithmetic, its time and memory costs are $\mathcal O(nm)$.

:::: {.callout-note collapse="true" title="Going deeper: free ends and fitting an entire read"}

### Free-end alignment variants

The recurrence, borders, answer cell, and stopping rule act together. Start from the three-way global recurrence, without a zero candidate in interior cells. Zero borders and a required bottom-right answer allow free leading ends but charge trailing ends. Global borders and an answer taken from the last row or last column instead charge leading ends and allow free trailing ends. These are free-end variants; neither permits the arbitrary internal restart used by local alignment.

### Fitting an entire read inside a reference

Suppose $x$ is a read of length $n$ and $y$ is a longer reference of length $m$. We want every read base represented, while reference bases before and after the placement are free. Let $Q(i,j)$ use the global three-way recurrence, with rows for read prefixes and columns for reference prefixes. Choose
$$
Q(0,j)=0,\qquad Q(i,0)=ig,\qquad
\text{answer}=\max_{0\leq j\leq m}Q(n,j).
$$
The zero top row lets us skip a reference prefix. The penalized first column prevents skipping a read prefix for free. Taking the answer anywhere in the last row leaves the unused reference suffix uncharged but requires consuming the whole read. Traceback stops on the top row, where $i=0$; it does not stop at an arbitrary interior zero.

For read $\texttt{AC}$ and reference $\texttt{TTACGG}$, using match $+1$, mismatch $-1$, and gap $-2$, the fitting alignment scores 2 and places the read at reference positions 3 and 4. Two is also an upper bound: each of the two read letters can contribute at most one match reward, and gaps have negative scores. Global alignment instead must account for four additional reference letters; local alignment is allowed to omit part of the read.

These boundary choices express a biological requirement in the table. Before modifying an alignment routine, state which input letters the output must consume, then choose its initialization and endpoint to enforce that requirement.

::::

## Scoring and region length

The restart rule relies on the scoring choice from Lesson 2. A calibrated log-odds table has negative expected score under its background model, so extending an unselected alignment through unrelated sequence tends to decrease its score. This makes the edges of a shared region detectable. If background columns instead have positive average score, the algorithm is rewarded for extending into unrelated material and the selected region expands. A zero mismatch score does not by itself force the whole sequence into the answer, but it can make neutral columns free and produce long, poorly localized optima.

::: {.callout-warning title="Similarity is not evidence by itself"}
A high local score can arise in a low-complexity repeat such as ACACAC. Database tools mask or downweight such regions because the scoring model alone does not know that many similar-looking hits occur by chance.
:::

Smith-Waterman finds the best region under its score. It does not supply a probability that the sequences are homologous. Significance also depends on sequence lengths, composition, and how many database records were searched.

## Limitations

1. **Only the best region is returned.** Two disjoint domains require repeated searches or methods designed for multiple hits.
2. **The score controls length.** Weak mismatch or gap penalties allow unrelated columns to join nearby good regions.
3. **Low complexity creates strong accidental matches.** Repetitive composition violates the simple background model.
4. **Quadratic time remains.** Exact local alignment against every database record is too costly at database scale. Lesson 5 introduces seed-and-extend search.

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
