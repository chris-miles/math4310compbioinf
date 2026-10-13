---
published: true
title: "Global Alignment: Dynamic Programming"
subtitle: "Lesson 3 · 2027-01-20 · Week 2"
nocite: |
  @compeau2015, @durbin1998, @needleman1970
---

## Core question

How can a table represent every alignment without listing them?

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

Suppose we have assembled DNA from two bacterial isolates and extracted the same gene from each assembly. The inputs are now two gene sequences, rather than two arbitrary raw reads. We expect the records to correspond from beginning to end, but an insertion or deletion can shift their internal coordinates. Global alignment makes that end-to-end expectation part of the problem: every letter in both records must appear in the answer.

We have scores for matches, mismatches, and gaps, but trying every gap placement is already impractical for short strings. A dynamic program will summarize the competing alignments in a table and compute the best score.

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
#| fig-cap: "The alignment ACGT/A-GT as a path through consumed-prefix coordinates. The vertical step consumes C from x without consuming a letter from y."
#| fig-alt: "A five-by-four grid labeled A C G T on rows and A G T on columns, with red cells at (0,0), (1,1), (2,1), (3,2), and (4,3)."
#| fig-width: 4.5
#| fig-height: 3.5
path_cells = [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3)]
empty_grid = [[None] * 4 for _ in range(5)]
fig, ax = plt.subplots()
cf.dp_table(ax, list("ACGT"), list("AGT"), empty_grid, path=path_cells)
plt.show()
```

## Counting paths

Different paths can consume the same two prefixes. Once they reach $(i,j)$, they have the same unused suffixes, even if their earlier gap placements differ. This shared continuation is what makes a table useful. Lesson 1 counted all ways of reaching each prefix pair using $A(i,j)$; we now interpret those counts as numbers of grid paths.

A *recurrence* defines a quantity from values on smaller inputs. Boundary conditions give the smallest values directly.

For $i,j\geq1$, the counts satisfy
$$
A(i,j)=A(i-1,j-1)+A(i-1,j)+A(i,j-1),
$$
with $A(i,0)=A(0,j)=1$. Removing the final move leaves a path to one of the three predecessors. Every predecessor path extends in exactly one way for that move, so the three counts add without omission or double counting. On a boundary there is only one path, consisting entirely of horizontal or vertical moves.

::: {#exm-delannoy-count}
## Counting before optimizing
The recurrence gives $A(2,2)=13$, $A(3,3)=63$, and $A(10,10)=8{,}097{,}453$. Two length-10 strings have more than eight million alignments, but their table has only $11\times11=121$ cells.
:::

We write $f(n)\in\mathcal O(g(n))$ if constants $C>0$ and $n_0$ exist such that $f(n)\le Cg(n)$ for all $n\ge n_0$.

A table for lengths $n,m$ has $(n+1)(m+1)$ cells. For $n,m\geq1$, constant work per cell costs $\mathcal O(nm)$ time. Here we count integer addition as one operation; the counts themselves grow large. Big-O compares growth rates; it does not say fixed costs are irrelevant. See the [Big-O reference](../appendices/computing-math-reference.html#big-o-notation) for an example and graph.

### Repeated recursion and tabulation

A direct recursive program mirrors the counting formula, but it recomputes the same quantities. To find $A(10,10)$, both $A(10,9)$ and $A(9,10)$ request $A(9,9)$, and the duplication continues below it. Memoization stores a value the first time it is requested. Bottom-up dynamic programming chooses an order in which every dependency is already available.

For this grid, row-major order works because each cell depends only on the row above and the cell to its left. Column-major order works too. An order that attempts $(i,j)$ before $(i-1,j)$ does not.

The table stores one count per prefix pair instead of recomputing it in many recursion branches. To optimize an alignment, we retain the same index pairs but store a different summary: the best score reaching each pair.

## The global alignment objective

Counting tells us how many alignments exist, but the biological comparison calls for a preferred correspondence. Let $s(a,b)$ score a letter pair, and let $g<0$ score a letter paired with a gap. For an alignment $A$ with rows $u,v$ of common length $L$, define the column score $c(a,b)=s(a,b)$ when both symbols are letters and $c(a,b)=g$ when exactly one is a gap. Thus
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

$F(i,j)$ summarizes every possible alignment of those prefixes by its best score. In particular, $F(n,m)$ answers the full comparison. The count $A(i,j)$ and optimum $F(i,j)$ describe the same set of paths but retain different information about it.

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

These scores can be computed in the same row-major order as the counts. Every predecessor has either a smaller row index or the same row index and a smaller column index. No complete alignment is enumerated; the maximum compares three best prefix scores.

::: {.algorithm}
**Algorithm: fill the global-alignment table.**

1. Create an $(n+1)\times(m+1)$ table.
2. Set the boundary values.
3. Fill interior cells in row or column order using @eq-global-recurrence.
4. Return $F(n,m)$.
:::

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

We can now check the final answer independently. No alignment has more than three letter pairs, each worth at most $+1$, and the length difference forces at least one $-2$ gap. The score is at most $3-2=1$, attained by ACGT/A-GT. For the gene comparison, this optimum proposes a correspondence under our chosen scores. Lesson 4 recovers that correspondence from the table and proves the recurrence correct.

::: {#prp-global-fill-cost}
## Cost of filling the table
For $n,m\geq1$ and constant-time score lookup and arithmetic, the score table takes $\mathcal O(nm)$ time and $\mathcal O(nm)$ memory.
:::

::: {.proof}
There are $(n+1)(m+1)$ cells. Each interior cell evaluates three candidates, and each boundary cell is set once.
:::

Each entry reads only the previous row and the current row. If only the score is needed, keeping these two rows reduces storage to $\mathcal O(m)$; reconstructing an alignment also requires retaining or recomputing the choices.

:::: {.callout-note collapse="true" title="Going deeper: the same table for edit distance"}

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
::::

The recurrence finds an exact optimum for its supplied scores. The output here is a score; recovering the corresponding letters requires traceback. The end-to-end requirement and per-character gap rule are modeling choices that later lessons will change.

## Exercises

::: {#exr-global-count}
Fill $A(i,j)$ through $A(3,3)$. Group the 13 paths to $(2,2)$ by final move.
:::

::: {#exr-global-hand-table}
Fill the table for $x=\texttt{CAT}$ and $y=\texttt{CT}$ with match $+2$, mismatch $-1$, gap $-2$. Show the candidates at $F(2,2)$.
:::

::: {#exr-global-boundary}
Initialize the borders to zero and score AC against C. Explain which letter is made free.
:::

::: {#exr-global-order}
Give two valid cell-filling orders. State the dependency every valid order respects.
:::

::: {#refs}
**References**
:::
