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

Suppose we sequence two bacterial isolates, assemble their reads into longer DNA strings, and extract the same gene from each assembly. We want to identify which bases correspond and where the gene sequences differ. Locating those differences gives us specific DNA changes to test for effects on the gene's activity. Because both records cover the complete gene, we expect a comparison from beginning to end. An insertion or deletion can shift the internal positions, so equal-numbered bases need not be the corresponding ones.

For a toy version, take $x=\texttt{ACGT}$ and $y=\texttt{AGT}$: pairing C with a gap lets the later G and T bases match. Longer strings admit many competing gap placements, and their number grows exponentially with length. Given scores for paired letters and gaps, we seek the highest-scoring correspondence that uses every letter of both records, called a *global alignment*. Because many correspondences share their beginnings, trying each one separately would repeat much of the same work; a dynamic program instead stores the best score for each pair of prefixes and then recovers an alignment attaining the full score.

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

Different paths can reach the same prefix pair, and from there they face the same unused suffixes. A table indexed by prefix pairs can therefore keep one summary per cell, computed from smaller prefixes by a *recurrence*. Here the summary will be the best score.

## The global alignment objective

The grid represents every allowed correspondence, but the biological comparison calls for a preferred one. Let $s(a,b)$ score a letter pair, and let $g<0$ score a letter paired with a gap. For an alignment $\mathcal A$ with rows $u,v$ of common length $L$, define the column score $c(a,b)=s(a,b)$ when both symbols are letters and $c(a,b)=g$ when exactly one is a gap. Thus
$$
\operatorname{score}(\mathcal A)=\sum_{\ell=1}^{L}c(u_\ell,v_\ell).
$$
We seek an alignment attaining $\max_{\mathcal A}\operatorname{score}(\mathcal A)$, where the maximum runs over every alignment $\mathcal A$ of $x$ and $y$.

::: {#def-global-alignment}
## Global alignment problem
**Input:** strings $x$ and $y$, a substitution score $s$, and a linear gap score $g<0$.

**Output:** an alignment of $x$ and $y$, using every letter of both strings, with maximum score.
:::

::: {#def-global-subproblem}
## Global alignment subproblem
For $0\le i\le n$ and $0\le j\le m$, let $F(i,j)$ be the maximum score of a global alignment of prefixes $x_1\cdots x_i$ and $y_1\cdots y_j$.
:::

$F(i,j)$ summarizes every possible alignment of those prefixes by its best score. In particular, $F(n,m)$ answers the full comparison.

### The three-way recurrence

Many paths reach each cell, yet the table stores a single number there. Why is the best score enough? Take the five paths into $(2,1)$, which align AC with A, and score them with match $+1$, mismatch $-1$, and gap $-2$, the scores we use for the full table below. They score $-1$, $-3$, $-6$, $-6$, and $-6$ (@fig-global-path-bundle). Any continuation $C$ from $(2,1)$ to $(4,3)$ consumes the same remaining letters, and because scores add by columns it adds the same amount to every path $P$ it extends:
$$
\operatorname{score}(P\,C)=\operatorname{score}(P)+\operatorname{score}(C).
$$
Adding the same number to every candidate cannot change their order. Along the dashed G/G, T/T continuation, the five totals become $1,-1,-4,-4,-4$, and the path that was best at $(2,1)$ is still best. We may therefore discard every lower-scoring path into a cell and keep only its best score.

```{.python .execute}
#| label: fig-global-path-bundle
#| code-fold: true
#| code-summary: "Draw the five paths into one cell"
#| fig-cap: "The five alignment paths of AC with A, all ending at cell (2,1), with their scores under match +1, mismatch −1, gap −2. Every path continues along the same dashed G/G, T/T columns, which add +2 to each, so the red path that is best at (2,1) stays best and the grey paths can be discarded."
#| fig-alt: "A five-by-four grid of dots. Grey and red line segments run from the top-left dot to the dot in row C, column A; a dashed black diagonal continues to the bottom-right dot. A legend lists AC/A- with score −1, AC/-A with −3, and three other paths with −6."
#| fig-width: 5
#| fig-height: 3.4
paths_into_cell = [
    [(0, 0), (1, 1), (2, 1)],
    [(0, 0), (1, 0), (2, 1)],
    [(0, 0), (1, 0), (2, 0), (2, 1)],
    [(0, 0), (1, 0), (1, 1), (2, 1)],
    [(0, 0), (0, 1), (1, 1), (2, 1)],
]
path_labels = ["AC/A-  −1", "AC/-A  −3", "three other paths  −6", None, None]
fig, ax = plt.subplots()
cf.path_bundle(ax, list("ACGT"), list("AGT"), paths_into_cell, keep=0,
               continuation=[(2, 1), (3, 2), (4, 3)], labels=path_labels)
plt.show()
```

An optimal alignment ending at $(i,j)$ must use one of the three final columns from @def-alignment-path. Remove that column. If the remaining prefix were worse than another path to the same predecessor, replacing it would improve the whole alignment. Conversely, appending the indicated column to an optimal predecessor gives a valid candidate. For $i,j\geq1$, compare all three:
$$
F(i,j)=\max\begin{cases}
F(i-1,j-1)+s(x_i,y_j), & \text{pair the final letters},\\
F(i-1,j)+g, & \text{pair }x_i\text{ with a gap},\\
F(i,j-1)+g, & \text{pair a gap with }y_j.
\end{cases}
$$ {#eq-global-recurrence}

Read the recurrence as three lookups and a comparison: take the best scores already stored for the three neighboring prefix pairs, add the score of the column that connects each one to $(i,j)$, and keep the largest. The boundaries are
$$
F(0,0)=0,\qquad F(i,0)=ig,\qquad F(0,j)=jg.
$$
A nonempty prefix aligned to an empty string needs one gap per letter. A zero border makes end letters free and answers a different problem.

Filling the table row by row works because every predecessor lies in an earlier row or earlier in the same row.

::: {.callout-note title="Key idea"}
Every continuation from a cell adds the same score to all paths that reach it, so a best path into $(i,j)$ extends a best path into one of its three predecessors. One number per cell therefore summarizes exponentially many alignments.
:::

::: {.algorithm}
**Algorithm: fill the global-alignment table.**

1. Create an $(n+1)\times(m+1)$ table.
2. Set the boundary values.
3. Fill interior cells in row or column order using @eq-global-recurrence.
4. Return $F(n,m)$.
:::

::: {#thm-global-correctness}
## Correctness of the global-alignment recurrence
The optimal prefix scores $F(i,j)$ of @def-global-subproblem take the boundary values above and, for $i,j\geq1$, satisfy @eq-global-recurrence. Filling the table in row order with these rules therefore computes every $F(i,j)$, including the answer $F(n,m)$.
:::

::: {.proof}
On a boundary, a nonempty prefix has only one alignment with an empty string, every letter opposite a gap, so $F(i,0)=ig$ and $F(0,j)=jg$.

For $i,j\geq1$, let $R(i,j)$ be the maximum on the right side of @eq-global-recurrence. Take an optimal alignment $\mathcal A$ of the two prefixes. Its final column is one of the three kinds in @def-alignment-path, with score $c$, and removing it leaves an alignment $\mathcal A'$ of the prefixes at the matching predecessor $(i',j')$. That shorter alignment must be optimal, or replacing it with a better one would improve $\mathcal A$. Hence
$$
F(i,j)=\operatorname{score}(\mathcal A')+c=F(i',j')+c\leq R(i,j).
$$
Conversely, each candidate in $R(i,j)$ is the score of an actual alignment of the two prefixes, an optimal predecessor alignment with the indicated column appended, so
$$
R(i,j)\leq F(i,j).
$$
Thus $F(i,j)=R(i,j)$. Row order computes every predecessor before the cell that reads it, so by induction each table entry equals $F(i,j)$.
:::

The proof uses two facts: the final-column cases are exhaustive, and an optimal solution cannot contain a suboptimal prefix. The second fact is called *optimal substructure*.

## A complete table

Use $x=\texttt{ACGT}$, $y=\texttt{AGT}$, match $+1$, mismatch $-1$, and gap $-2$. On a uniform background, match $+1$ and mismatch $-1$ are log-odds scores for related sequences about 75% identical; the gap score $-2$ is a separate modeling choice. The recurrence works for any additive scores.

The earlier path placed C opposite a gap, which costs more at that step than pairing C with G. A greedy rule that always appends the best-scoring next column takes the cheaper step and pays later: it pairs C/G and then G/T for $-1$ each and must end with T opposite a gap, so ACGT/AGT- scores $1-1-1-2=-3$. A tourist on a street grid who always picks the more interesting next block can miss the best route in the same way [@compeau2015]. The table instead keeps an optimum for each prefix pair; a choice that wins for AC/AG may lose in the complete comparison.

::: {#exm-global-cell}
## Computing one cell
At $F(2,2)$, C meets G. The candidates are
$$
F(1,1)-1=0,\qquad F(1,2)-2=-3,\qquad F(2,1)-2=-3.
$$
Here $F(1,1)=1$ comes from A/A. The predecessors $F(1,2)=F(2,1)=-1$ each pair A/A and add one gap. Thus $F(2,2)=0$, attained by AC/AG with one match and one mismatch.
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

In @fig-global-table, the borders decrease by two each step and the final entry is $F(4,3)=1$. At $(3,2)$, pairing G/G extends the prefix AC/A, whose best alignment puts C opposite a gap. The winning full path therefore bypasses the shaded cell $(2,2)$: its mismatch was optimal for AC/AG, but those are not the prefixes used by the full optimum.

A short bound confirms this value: no alignment has more than three letter pairs, each worth at most $+1$, and the length difference forces at least one $-2$ gap. The score is at most $3-2=1$, attained by ACGT/A-GT. For the gene comparison, this optimum proposes a correspondence under our chosen scores.

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

For $x=\texttt{ACGT}$ and $y=\texttt{AGT}$, with match $+1$, mismatch $-1$, and gap $-2$, the unique path (the red cells in @fig-global-table) yields
$$
\begin{array}{c}
\texttt{ACGT}\\
\texttt{A-GT}
\end{array}
$$
with score $1-2+1+1=1=F(4,3)$. In general, any path that follows maximizing backpointers from $(n,m)$ to $(0,0)$ produces an alignment with score $F(n,m)$.

::: {.proof}
Sketch. Write $v_L=(n,m),v_{L-1},\ldots,v_0=(0,0)$ for the cells on the path and $c_\ell$ for the score of the column emitted on the move from $v_\ell$ to $v_{\ell-1}$. Each backpointer satisfies $F(v_\ell)=F(v_{\ell-1})+c_\ell$, so the intermediate table values cancel:
$$
F(n,m)=c_L+F(v_{L-1})=c_L+c_{L-1}+F(v_{L-2})=\cdots=\sum_{\ell=1}^{L}c_\ell+F(0,0)=\sum_{\ell=1}^{L}c_\ell.
$$
:::

### Traceback by hand

At $(4,3)$ in the ACGT/AGT table, T meets T. The diagonal candidate is $F(3,2)+1=1$; the vertical and horizontal candidates are $-3$ and $-4$. We move diagonally and emit T/T. The same calculation at $(3,2)$ emits G/G.

At $(2,1)$, C meets A. The candidates are $F(1,0)-1=-3$, $F(1,1)-2=-1$, and $F(2,0)-2=-6$. The vertical pointer wins, so we emit C/-. The remaining diagonal emits A/A. Traceback generates columns from right to left, so an implementation must reverse both collected rows.

Every move decreases $i$, $j$, or both, so the process reaches $(0,0)$ after at most $n+m$ moves.

For AA against A, either AA/A- or AA/-A scores $1-2=-1$. Both placements are optimal: the scores do not identify which A corresponds to the shorter record. A *tie policy* chooses a maximizing predecessor when candidates agree, for example preferring diagonal, then vertical, then horizontal moves. It prints one optimal alignment without establishing uniqueness.

### Cost of the table

Here $\mathcal O(nm)$ means that work or storage grows at most proportionally to the product of the sequence lengths, apart from a fixed constant. For $n,m\geq1$ and constant-time score lookup and arithmetic, filling the score table takes $\mathcal O(nm)$ time and $\mathcal O(nm)$ memory.

::: {.proof}
Sketch. There are $(n+1)(m+1)$ cells; each interior cell compares three candidates and each boundary cell is set once. Since $n+1\leq2n$ and $m+1\leq2m$,
$$
\text{work}\leq3\,(n+1)(m+1)\leq12\,nm .
$$
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
Here $[x_i\ne y_j]$ is 1 when the letters differ and 0 when they match. The boundaries are $E(i,0)=i$ and $E(0,j)=j$. Because edit distance is the minimum cost over alignments, $E(n,m)=d_E(x,y)$ by @thm-edit-alignment.

::: {#exm-score-cost}
## Two conventions, one path
For ACGT/AGT, the path with A/A, C/-, G/G, T/T has global score $1-2+1+1=1$ and edit cost $0+1+0+0=1$. For the exact match ACG/ACG, the best global score is 3 and the edit distance is 0. The two numerical scales have different meanings even when the chosen alignment agrees.
:::

Negating every column score converts maximization to an equivalent minimization problem. Changing the relative costs of substitutions and gaps changes the model. In either case, the proof considers the same exhaustive final-column cases.

### Counting paths

Lesson 1 counted all ways of reaching each prefix pair using $A(i,j)$; we now interpret those counts as numbers of grid paths.

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

Rescoring is independent of the table recurrence and can detect an omitted gap or a misplaced column. Agreement alone is insufficient: `ATGT/A-GT` also scores 1, but it does not align the input `ACGT`. The two gap-removal checks detect that wrong letter even though its score agrees.

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

A tie policy, a deterministic rule for choosing among maximizing predecessors, prints one of these optima but does not make it unique. Enumerating all optimal alignments can take exponential time because there can be exponentially many tied paths.

### Sensitivity to score parameters

An optimality proof applies to the specified scores. We can inspect how a conclusion depends on them by scoring a few competing alignments symbolically.

For ACGT and AGCT, keep the match reward at 1, write the mismatch penalty as $p>0$, and the per-character gap penalty as $d>0$. The ungapped alignment has two matches and two mismatches, giving $2-2p$. Either displayed gapped alignment has three matches and two gap characters, giving $3-2d$. The gapped candidates beat the ungapped candidate exactly when
$$
3-2d>2-2p,\qquad\text{or equivalently}\qquad p>d-\tfrac12.
$$
With $p=d=1$, they win by one point. Keeping $p=1$ but raising $d$ to 2 reverses this comparison: the ungapped candidate scores 0 and the gapped candidates score $-1$. At $p=d-1/2$ they tie.

This calculation compares only these candidates; the dynamic program is what certifies an optimum over all alignments. Still, it shows before any code is rerun that a small change in $p$ or $d$ can move the chosen gap, and when that happens the change reflects the scoring assumptions rather than a bug.

### Counting optimal alignments

A tie can leave the score settled while the correspondence remains ambiguous. To count alternatives without printing them, let $C(i,j)$ count optimal paths to $(i,j)$ and let $\operatorname{Pred}^*(i,j)$ be the predecessors whose candidate scores attain $F(i,j)$. Sum their counts:
$$
C(i,j)=\sum_{(i',j')\in\operatorname{Pred}^*(i,j)}C(i',j'),\qquad
C(i,0)=C(0,j)=1.
$$
The boundary rule includes $C(0,0)=1$, counting the empty alignment. By induction on $i+j$, $C(i,j)$ is exactly the number of optimal paths to $(i,j)$.

::: {.proof}
Sketch. An optimal path to $(i,j)$ arrives from one maximizing predecessor, and its part up to that predecessor is optimal there. Conversely, extending any optimal path to a maximizing predecessor by the final move gives an optimal path to $(i,j)$. Grouping paths by their final predecessor gives disjoint groups, so
$$
\#\{\text{optimal paths to }(i,j)\}=\sum_{(i',j')\in\operatorname{Pred}^*(i,j)}\#\{\text{optimal paths to }(i',j')\}=\sum_{(i',j')\in\operatorname{Pred}^*(i,j)}C(i',j').
$$
:::

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

These two alignments disagree about which A in the longer string corresponds to the A in the shorter one, and the supplied scores cannot distinguish them. Reporting a single traceback can hide this ambiguity. A large count means many correspondences attain the same optimum; additional biological evidence or a different scoring model would be needed to distinguish them.

::: {#exr-traceback-bug}
Change the traceback condition from "or" to "and". Predict the output for AA against A, then explain why its printed alignment cannot have the reported score.
:::

::: {#exr-correctness-template}
Suppose diagonal moves are forbidden. Write the recurrence and adapt the proof of @thm-global-correctness. Which step changes?
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
