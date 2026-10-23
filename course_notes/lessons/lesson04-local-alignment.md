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
import numpy as np
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

The input is still two strings and column scores, but the output must now say where the matching region starts and ends in each string, leaving the prefixes and suffixes outside it unused. This problem is *local alignment*; the Smith-Waterman algorithm solves it with one table that finds the endpoints along with the score.

Which letters must appear in the answer depends on the biological task:

| Input and intended comparison | Letters required | Alignment choice |
|---|---|---|
| Two complete versions of a corresponding gene | both full records | global |
| Two proteins sharing a domain | the shared regions | local |
| A read placed within a reference chromosome | the whole read, a reference interval | fitting (free reference ends) |

: Biological tasks determine alignment boundaries. {#tbl-alignment-biological-tasks}

::: {#def-local-alignment}
## Local alignment problem
Given strings $x$ and $y$ and additive column scores, find substrings of $x$ and $y$ whose alignment has maximum score. The empty alignment, with score zero, is allowed.
:::

Neither end of the shared region is known in advance. Fix a candidate endpoint $(i,j)$ for the moment, so that the region ends at $x_i$ and $y_j$. Its start could be any position up to $i$ in $x$, and likewise in $y$; those choices are the suffixes of $x_1\cdots x_i$ and $y_1\cdots y_j$. Write $\operatorname{opt}(u,v)$ for the best global-alignment score of two strings $u$ and $v$, with $\operatorname{opt}=0$ when both are empty. The table entry we want is
$$
H(i,j)=\max_{0\le i_0\le i,\ 0\le j_0\le j}\ \operatorname{opt}\bigl(x_{i_0+1}\cdots x_i,\ y_{j_0+1}\cdots y_j\bigr).
$$
Read it from the inside out: choose how many letters $i_0$ and $j_0$ to skip at the start of each prefix, align what remains end to end, and keep the largest score over all $(i+1)(j+1)$ choices of starts. Choosing $i_0=i$ and $j_0=j$ selects the empty alignment, so $H(i,j)\geq0$. Maximizing over the endpoint as well gives the best local score, $\max_{i,j}H(i,j)$. Computing each $\operatorname{opt}$ separately would repeat the global algorithm for every start; the recurrence below folds the choice of start into the table itself.

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

::: {.callout-note title="Key idea"}
The zero candidate starts a new alignment wherever every way of continuing the old one has a negative score. One table therefore maximizes over all starting points at once, and its largest entry chooses where the shared region ends.
:::

::: {.algorithm}
**Algorithm: find one optimal local alignment.**

1. Set the first row and column of $H$ to zero.
2. Fill each cell with @def-local-recurrence and remember a maximizing predecessor.
3. Start traceback at a largest cell in the table.
4. Follow pointers until reaching a cell with value zero.
:::

## A shared interior region

Return to $x=\texttt{TTACGCC}$ and $y=\texttt{GGACGAA}$ with match $+1$, mismatch $-1$, and gap $-2$, a shared core between differing flanks.

::: {#exm-local-reset}
## Discarding a poor prefix
At $H(1,1)$, T meets G:
$$
H(1,1)=\max\{0,-1,-2,-2\}=0.
$$
The next diagonal T/G mismatch is also discarded. Pairing the following A/A, C/C, and G/G gives $H(3,3)=1$, $H(4,4)=2$, and $H(5,5)=3$. Continuing through the two C/A mismatches lowers the diagonal scores to 2 and then 1. The largest score, rather than the final cell, selects the end of the core.
:::

The local alignment is ACG/ACG with score 3. The global alignment TTACGCC/GGACGAA scores $-2+3-2=-1$, because its four flank mismatches cancel more than the core earns. Here the core occupies positions 3 to 5 in both strings, so global alignment still pairs the right letters; it fails only by burying a perfect three-letter agreement inside a negative total. When the shared region sits at different positions in the two records, global alignment can also choose the wrong correspondence.

::: {#exm-local-global-shift}
## A shared core that global alignment misplaces
Suppose the domain lies near the start of one record and near the end of the other. Removing the first two letters of $x$ and the last two of $y$ gives $x=\texttt{ACGCC}$ and $y=\texttt{GGACG}$, so the core is now $x_1x_2x_3$ and $y_3y_4y_5$. A global alignment must use every letter, and pairing the two copies of ACG then needs two gap columns on each side. Compare it with the ungapped alignment:
$$
\begin{array}{ccc}
\text{alignment} & \text{columns} & \text{score}\\[2pt]
\begin{array}{c}\texttt{--ACGCC}\\\texttt{GGACG--}\end{array} & \text{3 matches, 4 gaps} & 3-8=-5\\[10pt]
\begin{array}{c}\texttt{ACGCC}\\\texttt{GGACG}\end{array} & \text{1 match, 4 mismatches} & 1-4=-3
\end{array}
$$
The ungapped alignment is the unique global optimum. It pairs A with G and G with A, and its one match is a coincidental C/C in column 4. Local alignment pays for none of those gaps: $H$ climbs along the diagonal from $(0,2)$ to $H(3,5)=3$, and traceback selects $x_1x_2x_3=y_3y_4y_5=\texttt{ACG}$.
:::

In the protein comparison, global alignment would report a negative score and line up the domain of one protein with unrelated residues of the other. Local alignment reports the shared region and its coordinates in both records, which is what we need to choose a region whose activity to test.

The function returns the score, the two gapped strings, and the grid coordinates $(i_0,j_0)$ where traceback stops and $(i_1,j_1)$ where it starts. The selected letters are $x_{i_0+1}\cdots x_{i_1}$ and $y_{j_0+1}\cdots y_{j_1}$, which are exactly the Python slices `x[i_0:i_1]` and `y[j_0:j_1]` because Python slices include the start index and exclude the end.

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
result, local_alignment("ACGCC", "GGACG")
```

The first output is `(3, "ACG", "ACG", (2, 2), (5, 5))`, agreeing with the hand calculation. For the shifted pair of @exm-local-global-shift, the coordinates `(0, 2)` and `(3, 5)` locate the same core at the start of one string and the end of the other. Placing `stop` first in a score tie discards a neutral prefix; updating `best_cell` only when the score increases chooses the first maximizing endpoint.


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

For TTACGCC/GGACGAA, traceback starts at $(5,5)$ and stops at $(2,2)$, selecting $x_3x_4x_5=y_3y_4y_5=\texttt{ACG}$. These coordinates are the biological output: they propose where a shared domain lies in each protein. They need not be unique, even when the score is.

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

The recurrence, borders, answer cell, and stopping rule act together. Keeping the three-way global recurrence but changing only the borders or the answer cell gives free-end variants: zero borders make leading end gaps free, and taking the answer anywhere in the last row or column makes trailing end gaps free. Neither permits the arbitrary internal restart used by local alignment.

### Fitting an entire read inside a reference

Suppose $x$ is a read of length $n$ and $y$ is a longer reference of length $m$. We want every read base represented, while reference bases before and after the placement are free. Let $Q(i,j)$ use the global three-way recurrence, with rows for read prefixes and columns for reference prefixes. Choose
$$
Q(0,j)=0,\qquad Q(i,0)=ig,\qquad
\text{answer}=\max_{0\leq j\leq m}Q(n,j).
$$
The zero top row lets us skip a reference prefix. The penalized first column prevents skipping a read prefix for free. Taking the answer anywhere in the last row leaves the unused reference suffix uncharged but requires consuming the whole read. Traceback stops on the top row, where $i=0$; it does not stop at an arbitrary interior zero.

For read $\texttt{AC}$ and reference $\texttt{TTACGG}$, using match $+1$, mismatch $-1$, and gap $-2$, the fitting alignment scores 2 and places the read at reference positions 3 and 4. Global alignment instead must account for four additional reference letters; local alignment is allowed to omit part of the read.

These boundary choices express a biological requirement in the table. Before modifying an alignment routine, state which input letters the output must consume, then choose its initialization and endpoint to enforce that requirement.

::::

## Scoring and region length

The restart rule finds the edges of a shared region only if unrelated sequence tends to lower the score. Lesson 2 measured this tendency with the expected score of one background column. If the two letters of an unrelated column are drawn independently with base probabilities $q_a$ and $q_b$, then
$$
\mathbb E_R[s]=\sum_{a,b}q_aq_b\,s(a,b).
$$
The sum weights each letter pair's score by how often unrelated sequence produces that pair, then adds. For uniform DNA with match $+1$ and mismatch $-1$, a column matches with probability $1/4$, so $\mathbb E_R[s]=\tfrac14(1)+\tfrac34(-1)=-\tfrac12$, and $\ell$ unrelated columns cost about $\ell/2$ points on average. A calibrated log-odds table always has negative background expectation.

To see what this does to the recurrence, follow one diagonal of a longer simulated pair and ignore gaps. With $s_k$ the score of column $k$, the running total and the restarted score are
$$
C_k=C_{k-1}+s_k,\qquad H_k=\max\{0,\ H_{k-1}+s_k\},\qquad C_0=H_0=0.
$$
Both add the new column's score; $H_k$ then replaces a negative result by zero, as the local recurrence does along diagonal moves. In @fig-local-score-drift, columns 81 to 110 come from a shared region with 90% identity and the flanks are unrelated. The running total drifts down through the flanks, while the restarted score waits near zero, climbs through the shared region, and peaks where it ends. Chance runs of matches still raise small bumps in the flanks, so a small positive local score is weak evidence.

```{.python .execute}
#| label: fig-local-score-drift
#| code-fold: true
#| code-summary: "Simulate one diagonal through unrelated and shared columns"
#| fig-cap: "Running total C_k (blue) and restarted score H_k (red) along one simulated ungapped diagonal, using match +1 and mismatch -1. Columns 81 to 110 (shaded) match with probability 0.9 and the rest with probability 1/4. Negative background expectation pulls the running total down in both flanks, while the restart holds the local score near zero until the shared region, so its maximum marks the region's end."
#| fig-alt: "Two score curves against column number from 1 to 190. A blue line falls steadily to about minus 50, rises across a shaded band from 81 to 110, and falls again to about minus 75. A red line stays near zero with small bumps, rises across the shaded band to a peak of 28 at column 110, and declines back to zero."
rng = np.random.default_rng(4)
flank_length, core_length = 80, 30
# Unrelated columns match with probability 1/4; shared-region columns with 0.9.
match_probability = ([0.25] * flank_length + [0.9] * core_length
                     + [0.25] * flank_length)
column_scores = np.where(rng.random(len(match_probability)) < match_probability, 1, -1)

running_total, restarted = [], []
c_k, h_k = 0, 0
for s_k in column_scores:
    c_k = c_k + s_k
    h_k = max(0, h_k + s_k)
    running_total.append(c_k)
    restarted.append(h_k)

columns = range(1, len(column_scores) + 1)
fig, ax = plt.subplots()
ax.axvspan(flank_length + 0.5, flank_length + core_length + 0.5,
           color=cf.LIGHT, label="shared region")
ax.axhline(0, color=cf.RULE, linewidth=0.8)
ax.plot(columns, running_total, label="running total $C_k$")
ax.plot(columns, restarted, color=cf.EMPHASIS, label="restarted score $H_k$")
ax.set(xlabel="Column $k$ along the diagonal", ylabel="Score")
ax.legend(loc="lower left")
plt.show()
```

If background columns instead had positive average score, both curves would climb through the flanks, and the maximum would move toward the end of the sequence; the selected region would expand into unrelated material. A zero mismatch score does not by itself force the whole sequence into the answer, but it can make neutral columns free and produce long, poorly localized optima.

::: {.callout-warning title="Similarity is not evidence by itself"}
A high local score can arise in a low-complexity repeat such as ACACAC. Database tools mask or downweight such regions because the scoring model alone does not know that many similar-looking hits occur by chance.
:::

Smith-Waterman finds the best region under its score. It does not supply a probability that the sequences are homologous. Significance also depends on sequence lengths, composition, and how many database records were searched.

## Limitations

1. **Only the best region is returned.** Two disjoint domains require repeated searches or methods designed for multiple hits.
2. **The score controls length.** Weak mismatch or gap penalties allow unrelated columns to join nearby good regions.
3. **Quadratic time remains.** Exact local alignment against every database record is too costly at database scale. Lesson 5 introduces seed-and-extend search.

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
