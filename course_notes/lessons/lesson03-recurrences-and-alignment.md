---
published: true
title: "Global Alignment: Dynamic Programming"
subtitle: "Lesson 3 · 2027-01-20 · Week 2"
nocite: |
  @compeau2015, @durbin1998, @needleman1970
---

## Core question

How can we find the best end-to-end comparison of two gene sequences?

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

## Alignment as a grid path

Suppose we sequence two bacterial isolates, assemble their reads into longer DNA strings, and extract the same gene from each assembly. We want to identify which bases correspond and where the gene sequences differ. Because both records cover the complete gene, we expect a comparison from beginning to end. An insertion or deletion can shift the internal positions, so equal-numbered bases need not be the corresponding ones.

For a toy version, take $x=\texttt{ACGT}$ and $y=\texttt{AGT}$. Pairing C with a gap allows the later G and T bases to match. Longer strings admit many competing gap placements. Given scores for paired letters and gaps, we seek the highest-scoring correspondence that uses every letter of both records, called a *global alignment*. Trying each correspondence separately would repeat much of the same work. A dynamic program stores the best score for each pair of prefixes, then recovers an alignment attaining the full score.

Let $x=x_1\cdots x_n$ and $y=y_1\cdots y_m$. After reading some alignment columns, record $(i,j)$, the number of letters consumed from $x$ and $y$. An alignment consumes letters in order. A letter over a letter consumes one symbol from each string; a gap column consumes a symbol from one string. These counts turn a correspondence between letters into a path on a grid.

::: {#def-alignment-path}
## Alignment path
An *alignment path* runs from $(0,0)$ to $(n,m)$ using three moves: diagonal $(i-1,j-1)\to(i,j)$ for $x_i/y_j$, vertical $(i-1,j)\to(i,j)$ for $x_i/-$, and horizontal $(i,j-1)\to(i,j)$ for $-/y_j$.
:::

Start at $(0,0)$ and read the columns to obtain the path. Conversely, read the moves and emit the indicated letters or gap to recover the alignment. Thus the correspondence is one to one. In our tables, $i$ increases down the rows and $j$ across the columns; the empty prefixes supply row zero and column zero.

::: {#exm-alignment-path}
## A path for ACGT and AGT
The alignment
$$
\begin{array}{c}\texttt{ACGT}\\\texttt{A-GT}\end{array}
$$
uses a diagonal move for A/A, a vertical move for C/-, then diagonals for G/G and T/T. Its coordinates are
$$
(0,0)\to(1,1)\to(2,1)\to(3,2)\to(4,3).
$$
At $(2,1)$ we have consumed AC from $x$ and A from $y$. The path ends at $(4,3)$ because every letter has been used.
:::

```{.python .execute}
#| label: fig-global-prefix-path
#| code-fold: true
#| code-summary: "Draw the alignment and its grid path"
#| fig-cap: "The proposed alignment ACGT/A-GT (left) and its path through consumed-prefix coordinates (right). The vertical step records the C opposite a gap, allowing the remaining G and T bases to correspond."
#| fig-alt: "Two sequence rows ACGT and A-GT beside a five-by-four grid. Red grid cells follow (0,0), (1,1), (2,1), (3,2), and (4,3)."
#| fig-width: 6
#| fig-height: 3.5
path_cells = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3)]
empty_grid = [[None] * 4 for _ in range(5)]
fig, axes = plt.subplots(1, 2)
cf.alignment(axes[0], "ACGT", "A-GT")
cf.dp_table(axes[1], list("ACGT"), list("AGT"), empty_grid, path=path_cells)
plt.show()
```

Different paths can reach the same prefix pair. They then leave the same unused suffixes, so a table can retain one summary per pair. A *recurrence* computes that summary from smaller prefixes; row-major order makes the row above and the cell to the left available before each new cell. We will retain the best score, rather than enumerate every path.

## The global alignment objective

The grid represents every allowed correspondence, but the biological comparison calls for a preferred one. Let $s(a,b)$ score a letter pair, and let $g<0$ score a letter paired with a gap. For an alignment $A$ with rows $u,v$ of common length $L$, define the column score $c(a,b)=s(a,b)$ when both symbols are letters and $c(a,b)=g$ when exactly one is a gap. Thus
$$
\operatorname{score}(A)=\sum_{\ell=1}^{L}c(u_\ell,v_\ell),\qquad
A^*\in\operatorname*{arg\,max}_{A\in\mathcal A(x,y)}\operatorname{score}(A),
$$
where $\mathcal A(x,y)$ is the set of global alignments of $x$ and $y$. We seek an alignment attaining the largest total, using the same three moves with these scores attached.

::: {#def-global-alignment}
## Global alignment problem
Given strings $x,y$, substitution score $s$, and linear gap score $g$, find a maximum-score alignment that uses every letter of both strings.
:::

::: {#def-global-subproblem}
## Global alignment subproblem
For $0\le i\le n$ and $0\le j\le m$, let $F(i,j)$ be the maximum score of a global alignment of prefixes $x_1\cdots x_i$ and $y_1\cdots y_j$.
:::

$F(i,j)$ summarizes every possible alignment of those prefixes by its best score. In particular, $F(n,m)$ answers the full comparison. The optional path-counting extension applies a different summary to the same grid.

### The three-way recurrence

Suppose two paths reach the same predecessor cell, with scores 5 and 3. Appending the same final column adds the same amount to both, so the second path cannot overtake the first. The same is true for any common continuation because scores add by columns. We may therefore discard every lower-scoring path to that cell while retaining its best score.

An optimal alignment ending at $(i,j)$ must use one of the three final columns from @def-alignment-path. Remove that column. If the remaining prefix were worse than another path to the same predecessor, replacing it would improve the whole alignment. Conversely, appending the indicated column to an optimal predecessor gives a valid candidate. For $i,j\geq1$, compare all three:
$$
F(i,j)=\max\begin{cases}
F(i-1,j-1)+s(x_i,y_j), & \text{pair the final letters},\\
F(i-1,j)+g, & \text{pair }x_i\text{ with a gap},\\
F(i,j-1)+g, & \text{pair a gap with }y_j.
\end{cases}
$$ {#eq-global-recurrence}

The boundaries are
$$
F(0,0)=0,\qquad F(i,0)=ig,\qquad F(0,j)=jg.
$$
A nonempty prefix aligned to an empty string needs one gap per letter. A zero border makes end letters free and answers a different problem.

These scores can be computed in row-major order. Every predecessor has either a smaller row index or the same row index and a smaller column index. No complete alignment is enumerated; the maximum compares three best prefix scores.

::: {.algorithm}
**Algorithm: fill the global-alignment table.**

1. Create an $(n+1)\times(m+1)$ table.
2. Set the boundary values.
3. Fill interior cells in row or column order using @eq-global-recurrence.
4. Return $F(n,m)$.
:::

::: {#thm-global-correctness}
## Correctness of the global-alignment recurrence
For every $0\le i\le n$ and $0\le j\le m$, $F(i,j)$ equals the maximum score among all global alignments of $x_1\cdots x_i$ and $y_1\cdots y_j$.
:::

::: {.proof}
We use induction on $i+j$. The claim holds on the boundaries because there is only one alignment of a nonempty prefix with an empty string: every letter is paired with a gap.

For $i,j>0$, take an optimal alignment $A$. Its final column is $x_i/y_j$, $x_i/-$, or $-/y_j$. Removing that column leaves a valid alignment for the corresponding predecessor cell. It must be optimal for those prefixes. Otherwise, replacing it with a better prefix alignment would improve $A$, a contradiction. Thus $A$ has the score of one candidate in the recurrence, so the optimum is no larger than the table maximum.

Conversely, append the appropriate final column to an optimal alignment from each predecessor. Each candidate is attained by a valid alignment. Their maximum is therefore attainable and cannot be larger than the optimum. Together with the preceding inequality, this proves equality.
:::

The proof uses two facts: the final-column cases are exhaustive, and an optimal solution cannot contain a suboptimal prefix. The second fact is called *optimal substructure*.

## A complete table

Use $x=\texttt{ACGT}$, $y=\texttt{AGT}$, match $+1$, mismatch $-1$, and gap $-2$. These are illustrative scores, rather than the calibrated log-odds table from Lesson 2. The algorithm applies to either choice of additive scores.

The earlier path placed C opposite a gap. The table must decide whether that choice beats pairing C with G and repairing the correspondence later. Each prefix gets its own optimum; a choice that wins for AC/AG may lose in the complete comparison.

::: {#exm-global-cell}
## Computing one cell
At $F(2,2)$, C meets G. The candidates are
$$
F(1,1)-1=0,\qquad F(1,2)-2=-3,\qquad F(2,1)-2=-3.
$$
Here $F(1,1)=1$ comes from A/A. The predecessors $F(1,2)=F(2,1)=-1$ each pair A/A and add one gap. Thus $F(2,2)=0$, attained by AC/AG with one match and one mismatch. This is the best alignment for these prefixes; it need not be the prefix used by the best complete alignment.
:::

In Python, `table[i][j]` stores $F(i,j)$, so its indices are prefix lengths. Their final letters are `x[i - 1]` and `y[j - 1]` because Python strings start at index zero. The function returns the whole table for inspection; its bottom-right entry is the optimum score.

```{.python .execute}
#| label: code-global-score
def global_score(x: str, y: str, match: int = 1,
                 mismatch: int = -1, gap: int = -2) -> list[list[int]]:
    """Return the global-alignment score table."""
    table = [[0] * (len(y) + 1) for _ in range(len(x) + 1)]
    for i in range(1, len(x) + 1):
        table[i][0] = i * gap
    for j in range(1, len(y) + 1):
        table[0][j] = j * gap
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            pair = match if x[i - 1] == y[j - 1] else mismatch
            table[i][j] = max(table[i - 1][j - 1] + pair,
                              table[i - 1][j] + gap,
                              table[i][j - 1] + gap)
    return table

x, y = "ACGT", "AGT"
scores = global_score(x, y)
scores
```

```{.python .execute}
#| label: fig-global-table
#| fig-cap: "Global-alignment table for ACGT and AGT with match +1, mismatch -1, and gap -2. Three arrows show the candidates for the shaded AC/AG cell. The red optimal full path bypasses that cell, pairing C with a gap."
#| fig-alt: "A five-by-four score table with three arrows entering cell (2,2), shaded gray, and red path cells (0,0), (1,1), (2,1), (3,2), (4,3)."
#| fig-width: 4.5
#| fig-height: 4
fig, ax = plt.subplots()
cf.dp_table(ax, list(x), list(y), scores,
            path=path_cells, highlight=[(2, 2)],
            arrows=[((1, 1), (2, 2)), ((1, 2), (2, 2)), ((2, 1), (2, 2))])
plt.show()
```

In @fig-global-table, the borders decrease by two each step and the A/A cell is 1. No interior cell can exceed the largest predecessor by more than the match reward, 1.

The final entry is $F(4,3)=1$. At $(3,2)$, pairing G/G extends the prefix AC/A, whose best alignment puts C opposite a gap. The winning full path therefore bypasses the shaded cell $(2,2)$: its mismatch was optimal for AC/AG, but those are not the prefixes used by the full optimum.

We can now check the final answer independently. No alignment has more than three letter pairs, each worth at most $+1$, and the length difference forces at least one $-2$ gap. The score is at most $3-2=1$, attained by ACGT/A-GT. For the gene comparison, this optimum proposes a correspondence under our chosen scores. To identify the paired bases, we must recover the path that attains this score.

## Traceback

The best score does not yet identify which bases correspond. Record the maximizing predecessor for each cell so that its final column can be recovered.

::: {#def-backpointer}
## Backpointer
A *backpointer* from $(i,j)$ records a predecessor whose candidate attains the maximum defining $F(i,j)$. A cell can have more than one backpointer when candidates tie.
:::


::: {.algorithm}
**Algorithm: recover one optimal global alignment.**

1. Start at $(n,m)$.
2. Follow one backpointer to a predecessor.
3. Write the alignment column represented by that move.
4. Stop at $(0,0)$ and reverse the collected columns.
:::

For $x=\texttt{ACGT}$ and $y=\texttt{AGT}$, with match $+1$, mismatch $-1$, and gap $-2$, the unique path yields
$$
\begin{array}{c}
\texttt{ACGT}\\
\texttt{A-GT}
\end{array}
$$
with score $1-2+1+1=1=F(4,3)$.

```{.python .execute}
#| label: fig-global-traceback
#| fig-cap: "Traceback through the ACGT/AGT table. The red cells form the unique optimal path from the bottom-right entry to the origin."
#| fig-alt: "A five-by-four dynamic-programming table with a red path passing diagonally except for one vertical step."
#| fig-width: 4.5
#| fig-height: 4
x, y = "ACGT", "AGT"
scores = [[0, -2, -4, -6], [-2, 1, -1, -3], [-4, -1, 0, -2],
          [-6, -3, 0, -1], [-8, -5, -2, 1]]
fig, ax = plt.subplots()
cf.dp_table(ax, list(x), list(y), scores,
            path=[(0, 0), (1, 1), (2, 1), (3, 2), (4, 3)])
plt.show()
```

::: {#prp-traceback-score}
## Traceback preserves the table score
Any path that follows maximizing backpointers from $(n,m)$ to $(0,0)$ produces an alignment with score $F(n,m)$.
:::

::: {.proof}
Each pointer satisfies $F(i,j)=F(i',j')+c$, where $(i',j')$ is its predecessor and $c$ is the emitted column score. Summing along the path cancels the intermediate table values, leaving $F(n,m)=F(0,0)+\sum c=\sum c$.
:::

### Traceback by hand

At $(4,3)$ in the ACGT/AGT table, T meets T. The diagonal candidate is $F(3,2)+1=1$; the vertical and horizontal candidates are $-3$ and $-4$. We move diagonally and emit T/T. The same calculation at $(3,2)$ emits G/G.

At $(2,1)$, C meets A. The candidates are $F(1,0)-1=-3$, $F(1,1)-2=-1$, and $F(2,0)-2=-6$. The vertical pointer wins, so we emit C/-. The remaining diagonal emits A/A. Traceback generates columns from right to left, so an implementation must reverse both collected rows.

Every move decreases $i$, $j$, or both. The process reaches $(0,0)$ after at most $n+m$ moves. A pointer that leaves both indices unchanged creates an infinite loop; one that increases an index was stored in the wrong direction.

::: {.callout-tip title="Check the recovered alignment"}
Sum its printed column scores and compare with $F(n,m)$. Removing gaps must recover both input strings; matching the score alone cannot detect a wrong letter.
:::

For AA against A, either AA/A- or AA/-A scores $1-2=-1$. Both placements are optimal: the scores do not identify which A corresponds to the shorter record. A *tie policy* chooses a maximizing predecessor when candidates agree, for example preferring diagonal, then vertical, then horizontal moves. It prints one optimal alignment without establishing uniqueness.

::: {#prp-global-fill-cost}
## Cost of filling the table
Here $\mathcal O(nm)$ means that work or storage grows at most proportionally to the product of the sequence lengths, apart from a fixed constant. For $n,m\geq1$ and constant-time score lookup and arithmetic, the score table takes $\mathcal O(nm)$ time and $\mathcal O(nm)$ memory.
:::

::: {.proof}
There are $(n+1)(m+1)$ cells. Each interior cell evaluates three candidates, and each boundary cell is set once.
:::

Each entry reads only the previous row and the current row. If only the score is needed, keeping these two rows reduces storage to $\mathcal O(m)$; reconstructing an alignment also requires retaining or recomputing the choices.


The table and traceback give an exact optimum for the supplied additive scores. For our gene comparison, the gap proposes a difference between the records; without an ancestor, it cannot establish which lineage gained or lost C. Lesson 4 changes the end-to-end requirement to find a shared region. Lesson 5 changes the per-character gap model.

:::: {.callout-note collapse="true" title="Going deeper: other summaries of the alignment grid"}

The same table shape can minimize a cost. If a match costs 0 and a substitution, insertion, or deletion costs 1, define $E(i,j)$ as the minimum edit cost between the two prefixes:
$$
E(i,j)=\min\begin{cases}
E(i-1,j-1)+[x_i\ne y_j],\\
E(i-1,j)+1,\\
E(i,j-1)+1.
\end{cases}
$$
Here $[x_i\ne y_j]$ is 1 when the letters differ and 0 when they match. The boundaries are $E(i,0)=i$ and $E(0,j)=j$.

For ACGT and AGT, the minimum edit cost is 1, realized by deleting C. The maximum-score model above also chooses that alignment, but the numerical answers, 1 and 1, happen to agree only for this instance. One is a similarity score to maximize; the other is a cost to minimize.

::: {#exm-score-cost}
## Two conventions, one path
For ACGT/AGT, the path with A/A, C/-, G/G, T/T has global score $1-2+1+1=1$ and edit cost $0+1+0+0=1$. For the exact match ACG/ACG, the best global score is 3 and the edit distance is 0. The two numerical scales have different meanings even when the chosen alignment agrees.
:::

Negating every column score converts maximization to an equivalent minimization problem. Changing the relative costs of substitutions and gaps changes the model. In either case, the proof considers the same exhaustive final-column cases.

### Counting paths

Different paths can consume the same two prefixes. Once they reach $(i,j)$, they have the same unused suffixes, even if their earlier gap placements differ. This shared continuation is what makes a table useful. Lesson 1 counted all ways of reaching each prefix pair using $A(i,j)$; we now interpret those counts as numbers of grid paths.

For $i,j\geq1$, the counts satisfy
$$
A(i,j)=A(i-1,j-1)+A(i-1,j)+A(i,j-1),
$$
with $A(i,0)=A(0,j)=1$. Removing the final move leaves a path to one of the three predecessors. Every predecessor path extends in exactly one way for that move, so the three counts add without omission or double counting. On a boundary there is only one path, consisting entirely of horizontal or vertical moves.

::: {#exm-delannoy-count}
## Counting before optimizing
The recurrence gives $A(2,2)=13$, $A(3,3)=63$, and $A(10,10)=8{,}097{,}453$. Two length-10 strings have more than eight million alignments, but their table has only $11\times11=121$ cells.
:::

The counting table also has $(n+1)(m+1)$ cells, though its integers grow large. A direct recursive program repeatedly requests the same prefix count; storing each count once avoids that repetition.

::: {#exr-global-count}
Fill $A(i,j)$ through $A(3,3)$. Group the 13 paths to $(2,2)$ by final move.
:::

::::

:::: {.callout-note collapse="true" title="Going deeper: implementation and alignment ambiguity"}

### Reconstruction in Python

```{.python .execute}
#| label: code-global-traceback
def global_alignment(x: str, y: str, match: int = 1,
                     mismatch: int = -1, gap: int = -2) -> tuple[int, str, str]:
    """Return a global-alignment score and one optimal alignment."""
    n, m = len(x), len(y)
    table = [[0] * (m + 1) for _ in range(n + 1)]
    pointer = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        table[i][0], pointer[i][0] = i * gap, "up"
    for j in range(1, m + 1):
        table[0][j], pointer[0][j] = j * gap, "left"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = match if x[i - 1] == y[j - 1] else mismatch
            candidates = {
                "diag": table[i - 1][j - 1] + pair,
                "up": table[i - 1][j] + gap,
                "left": table[i][j - 1] + gap,
            }
            pointer[i][j] = max(candidates, key=candidates.get)
            table[i][j] = candidates[pointer[i][j]]
    top, bottom = [], []
    i, j = n, m
    while i > 0 or j > 0:
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
    return table[n][m], "".join(reversed(top)), "".join(reversed(bottom))

global_alignment("ACGT", "AGT")
```

The dictionary insertion order supplies the tie policy: diagonal, then up, then left. That choice affects which optimum is printed, not its score.

### Checking the reported score

A traceback implementation can return the correct table value beside the wrong alignment. We therefore score the printed columns independently.

::: {.algorithm}
**Algorithm: check an alignment score.**

Scan the displayed columns, reject a gap-gap column, and add the specified score for each letter pair or gap. Compare the total with the reported optimum and check that removing gaps recovers both inputs.
:::

```{.python .execute}
#| label: code-score-printed-alignment
def score_alignment(top: str, bottom: str, match: int = 1,
                    mismatch: int = -1, gap: int = -2) -> int:
    """Score an already constructed alignment."""
    if len(top) != len(bottom):
        raise ValueError("Aligned rows must have equal length.")
    total = 0
    for a, b in zip(top, bottom):
        if a == "-" and b == "-":
            raise ValueError("A column cannot contain two gaps.")
        if a == "-" or b == "-":
            total += gap
        elif a == b:
            total += match
        else:
            total += mismatch
    return total

reported, top, bottom = global_alignment("ACGT", "AGT")
assert score_alignment(top, bottom) == reported
assert top.replace("-", "") == "ACGT"
assert bottom.replace("-", "") == "AGT"
```

Rescoring is independent of the table recurrence and can detect an omitted gap or a misplaced column. Agreement alone is insufficient: `ATGT/A-GT` also scores 1, but it does not align the input `ACGT`. The two gap-removal checks detect that wrong letter even though its score agrees. Score agreement and input reconstruction check different requirements.

::: {.callout-tip title="Check: unconsumed prefixes"}
Traceback must continue while either index is positive. Stopping as soon as one index reaches zero drops the other string's remaining prefix. AA against A exposes this bug; an exact equal-length match may not.
:::


### Ties and optimal alignments

Take $x=\texttt{ACGT}$ and $y=\texttt{AGCT}$ with match $+1$, mismatch $-1$, and gap $-1$. The optimum score is 1, with two optimal alignments:
$$
\begin{array}{cc}
\texttt{ACG-T} & \texttt{A-CGT}\\
\texttt{A-GCT} & \texttt{AGC-T}
\end{array}
$$

::: {#def-tie-policy}
## Tie policy
A *tie policy* is a deterministic rule for choosing among maximizing predecessors. It selects an optimum but does not make that optimum unique.
:::

Different programs can print different optimal alignments while agreeing on the score. Enumerating all optimal alignments can take exponential time because there can be exponentially many tied paths.

### Sensitivity to score parameters

An optimality proof applies to the specified scores. We can inspect how a conclusion depends on them by scoring a few competing alignments symbolically.

For ACGT and AGCT, keep the match reward at 1, write the mismatch penalty as $p>0$, and the per-character gap penalty as $d>0$. The ungapped alignment has two matches and two mismatches, giving $2-2p$. Either displayed gapped alignment has three matches and two gap characters, giving $3-2d$. The gapped candidates beat the ungapped candidate exactly when
$$
3-2d>2-2p,\qquad\text{or equivalently}\qquad p>d-\tfrac12.
$$
With $p=d=1$, they win by one point. Keeping $p=1$ but raising $d$ to 2 reverses this comparison: the ungapped candidate scores 0 and the gapped candidates score $-1$. At $p=d-1/2$ they tie.

This calculation compares these candidates; a complete claim of optimality still requires checking all admissible alignments, which the dynamic program does. It gives a useful prediction before rerunning the code. If a small parameter change moves the chosen gap, that movement can reflect the scoring assumptions even when both implementations are correct. Counting ties at one parameter setting describes exact ambiguity; varying the scores asks whether the answer remains stable under nearby modeling choices.

### Counting optimal alignments

A tie can leave the score settled while the correspondence remains ambiguous. To count alternatives without printing them, let $C(i,j)$ count optimal paths to $(i,j)$ and let $\operatorname{Pred}^*(i,j)$ be the predecessors whose candidate scores attain $F(i,j)$. Sum their counts:
$$
C(i,j)=\sum_{(i',j')\in\operatorname{Pred}^*(i,j)}C(i',j'),\qquad
C(i,0)=C(0,j)=1.
$$
The boundary rule includes $C(0,0)=1$, counting the empty alignment. Each optimal path has one final predecessor, so the sum counts tied paths without double counting.

For AA against A with match $+1$ and gap $-2$, the score is $-1$ and the count is 2: AA/A- and AA/-A. The count tells us how many optimal paths exist, while the pointer graph shows where they diverge.

At $(2,1)$ in this example, the diagonal candidate is $F(1,0)+1=-1$ and the vertical candidate is $F(1,1)-2=-1$. The horizontal candidate is $F(2,0)-2=-6$ and does not contribute to the optimal count. Both maximizing predecessors have count 1, so
$$
C(2,1)=C(1,0)+C(1,1)=1+1=2.
$$

```{.python .execute}
#| label: fig-traceback-optimal-counts
#| fig-cap: "Optimal-path counts for AA against A. The two arrowed predecessors attain the final score, so their counts add instead of selecting one pointer."
#| fig-alt: "A three-by-two table of path counts. The final cell contains 2, with arrows entering from the cell above and the cell diagonally above-left; all other counts are 1."
#| fig-width: 3.5
#| fig-height: 3
optimal_counts = [[1, 1], [1, 1], [1, 2]]
fig, ax = plt.subplots()
cf.dp_table(
    ax, list("AA"), list("A"), optimal_counts,
    highlight=[(2, 1)],
    arrows=[((1, 0), (2, 1)), ((1, 1), (2, 1))],
)
plt.show()
```

These two alignments disagree about which A in the longer string corresponds to the A in the shorter one. Their tie does not establish equally probable evolutionary histories; it says the supplied scores do not distinguish those correspondences.

::: {#prp-count-optima}
## Counting optimal paths
If $C(i,j)$ sums the counts of exactly those predecessors that attain $F(i,j)$, then $C(i,j)$ equals the number of optimal alignment paths to $(i,j)$.
:::

::: {.proof}
Induct on $i+j$. Every optimal path reaches $(i,j)$ through one maximizing predecessor, and appending that final move to every optimal predecessor path produces every optimal path to $(i,j)$. Distinct final moves make the groups disjoint.
:::

Reporting a single traceback can hide this ambiguity. A large count means many correspondences attain the same optimum; additional biological evidence or a different scoring model would be needed to distinguish them.

::: {#exr-traceback-bug}
Change the traceback condition from "or" to "and". Predict the output for AA against A, then explain why its printed alignment cannot have the reported score.
:::

::: {#exr-correctness-template}
Suppose diagonal moves are forbidden. Write the recurrence and adapt the proof of @thm-global-correctness. Which sentence changes?
:::

::::

## Exercises

::: {#exr-global-hand-table}
Fill the table for $x=\texttt{CAT}$ and $y=\texttt{CT}$ with match $+2$, mismatch $-1$, gap $-2$. Show the candidates at $F(2,2)$.
:::

::: {#exr-global-boundary}
Initialize the borders to zero and score AC against C. Explain which letter is made free.
:::

::: {#exr-global-order}
Give two valid cell-filling orders. State the dependency every valid order respects.
:::

::: {#exr-traceback-hand}
Using the table for ACGT and AGT, mark every maximizing predecessor and perform the traceback. Verify the score by summing column scores.
:::

::: {#exr-traceback-ties}
Build the table for AA and A with match $+1$, mismatch $-1$, and gap $-2$. Identify both optimal alignments.
:::

::: {#refs}
**References**
:::
