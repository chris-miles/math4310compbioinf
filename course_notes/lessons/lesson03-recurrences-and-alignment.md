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

We have scores for matches, mismatches, and gaps, but trying every gap placement is already impractical for short strings. Today's dynamic program computes the best score; the next lesson recovers the aligned strings. If the input were instead a short read and an entire chromosome, requiring both full strings would be inappropriate. Lesson 5 changes the boundary rules to express that different task.

Let $x=x_1\cdots x_n$ and $y=y_1\cdots y_m$. An alignment consumes letters in order. A letter over a letter consumes one symbol from each string; a gap column consumes a symbol from one string.

::: {#def-alignment-path}
## Alignment path
An *alignment path* runs from $(0,0)$ to $(n,m)$ using three moves: diagonal $(i-1,j-1)\to(i,j)$ for $x_i/y_j$, vertical $(i-1,j)\to(i,j)$ for $x_i/-$, and horizontal $(i,j-1)\to(i,j)$ for $-/y_j$.
:::

Every alignment gives one path, and every path gives one alignment.

::: {#exm-alignment-path}
## A path for ACGT and AGT
The alignment
$$
\begin{array}{c}\texttt{ACGT}\\\texttt{A-GT}\end{array}
$$
uses a diagonal move for A/A, a vertical move for C/-, then diagonals for G/G and T/T.
:::

## Counting paths

::: {#def-recurrence}
## Recurrence
A *recurrence* defines a quantity from values on smaller inputs. Boundary conditions give the smallest values directly.
:::

If $D(i,j)$ counts alignment paths to $(i,j)$, then
$$
D(i,j)=D(i-1,j-1)+D(i-1,j)+D(i,j-1),
$$
with $D(i,0)=D(0,j)=1$. The terms classify paths by their final move, so no path is omitted or counted twice.

::: {#exm-delannoy-count}
## Counting before optimizing
The recurrence gives $D(2,2)=13$, $D(3,3)=63$, and $D(10,10)=8{,}097{,}453$. Two length-10 strings have more than eight million alignments, but their table has only $11\times11=121$ cells.
:::

::: {#def-big-o}
## Asymptotic upper bound
We write $f(n)\in\mathcal O(g(n))$ if constants $C>0$ and $n_0$ exist such that $f(n)\le Cg(n)$ for all $n\ge n_0$.
:::

A table for lengths $n,m$ has $(n+1)(m+1)$ cells. Constant work per cell costs $\mathcal O(nm)$ time. Big-O compares growth rates; it does not say fixed costs are irrelevant. See the [Big-O reference](../appendices/computing-math-reference.html#big-o-notation) for an example and graph.

### Repeated recursion and tabulation

A direct recursive program mirrors the counting formula, but it recomputes the same quantities. To find $D(10,10)$, both $D(10,9)$ and $D(9,10)$ request $D(9,9)$, and the duplication continues below it. Memoization stores a value the first time it is requested. Bottom-up dynamic programming goes one step further: it chooses an order in which every dependency is already available.

For this grid, row-major order works because each cell depends only on the row above and the cell to its left. Column-major order works too. An order that attempts $(i,j)$ before $(i-1,j)$ does not.

::: {#prp-dp-reuse}
## One value per subproblem
If each pair $(i,j)$ is computed once and each computation uses constant time after its predecessors are available, the entire recurrence table takes $\mathcal O(nm)$ time.
:::

::: {.proof}
There are $(n+1)(m+1)$ distinct index pairs. Computing each once replaces the repeated branches of the recursion tree with one table entry.
:::

Dynamic programming replaces repeated calculations with stored subproblem values. We still reason about all complete paths, but paths that share a prefix also share the stored answer for that prefix.

## The global alignment objective

Let $s(a,b)$ score a letter pair, and let $g<0$ score a letter paired with a gap. Column scores add.

::: {#def-global-alignment}
## Global alignment problem
Given strings $x,y$, substitution score $s$, and linear gap score $g$, find a maximum-score alignment that uses every letter of both strings.
:::

Global alignment fits two records believed to correspond end to end.

::: {#def-global-subproblem}
## Global alignment subproblem
For $0\le i\le n$ and $0\le j\le m$, let $F(i,j)$ be the maximum score of a global alignment of prefixes $x_1\cdots x_i$ and $y_1\cdots y_j$.
:::

Defining a cell before writing its recurrence prevents many indexing errors.

### The three-way recurrence

The final column has exactly one form from @def-alignment-path. Removing it leaves a smaller prefix alignment:
$$
F(i,j)=\max\begin{cases}
F(i-1,j-1)+s(x_i,y_j),\\
F(i-1,j)+g,\\
F(i,j-1)+g.
\end{cases}
$$ {#eq-global-recurrence}

The boundaries are
$$
F(0,0)=0,\qquad F(i,0)=ig,\qquad F(0,j)=jg.
$$
A nonempty prefix aligned to an empty string needs one gap per letter. A zero border makes end letters free and answers a different problem.

::: {.algorithm}
**Algorithm: fill the global-alignment table.**

1. Create an $(n+1)\times(m+1)$ table.
2. Set the boundary values.
3. Fill interior cells in row or column order using @eq-global-recurrence.
4. Return $F(n,m)$.
:::

## A complete table

Use $x=\texttt{ACGT}$, $y=\texttt{AGT}$, match $+1$, mismatch $-1$, and gap $-2$.

::: {#exm-global-cell}
## Computing one cell
At $F(2,2)$, C meets G. The candidates are
$$
F(1,1)-1=0,\qquad F(1,2)-2=-3,\qquad F(2,1)-2=-3.
$$
Thus $F(2,2)=0$.
:::

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
#| fig-cap: "Global-alignment table for ACGT and AGT. Three arrows show the candidates for the shaded cell."
#| fig-alt: "A five-by-four score table with three red arrows entering one shaded cell."
#| fig-width: 4.5
#| fig-height: 4
fig, ax = plt.subplots()
cf.dp_table(ax, list(x), list(y), scores,
            highlight=[(2, 2)],
            arrows=[((1, 1), (2, 2)), ((1, 2), (2, 2)), ((2, 1), (2, 2))])
plt.show()
```

The completed table is
$$
\begin{array}{c|rrrr}
 & - & A & G & T\\ \hline
- & 0 & -2 & -4 & -6\\
A & -2 & 1 & -1 & -3\\
C & -4 & -1 & 0 & -2\\
G & -6 & -3 & 0 & -1\\
T & -8 & -5 & -2 & 1
\end{array}
$$

Reading a table should include quick checks. The borders decrease by two each step. The A/A cell is 1. No interior cell can exceed the largest predecessor by more than the match reward, 1. These checks do not prove correctness, but each catches a common arithmetic or initialization error.

The final entry is $F(4,3)=1$. Lesson 4 recovers the alignment and proves the recurrence cannot miss a better one.

::: {#prp-global-fill-cost}
## Cost of filling the table
The score table takes $\mathcal O(nm)$ time and $\mathcal O(nm)$ memory.
:::

::: {.proof}
There are $(n+1)(m+1)$ cells. Each interior cell evaluates three candidates, and each boundary cell is set once.
:::

### Computing with two rows

To compute row $i$, the recurrence reads row $i-1$ and the entries already computed in row $i$. Older rows have no further role in the score calculation. We can therefore retain a previous row and a current row, replacing the previous row after each pass.

In the ACGT/AGT example, the row for prefix AC is $(-4,-1,0,-2)$. To compute the row for ACG, initialize its boundary to $-6$. The first interior cell compares $-4-1$, $-1-2$, and $-6-2$, giving $-3$. The next cell compares $-1+1$, $0-2$, and $-3-2$, giving 0. The final cell is $-1$. We have recovered $(-6,-3,0,-1)$ without reading the rows for the empty prefix or A.

The storage is now $2(m+1)$ numbers, or $\mathcal O(m)$, while the number of cell updates stays $nm$. This distinction matters: a memory improvement need not be a time improvement. For two strings of length 10,000, a full table has roughly $10^8$ cells; two rows have roughly $2\times10^4$ entries. Python containers have additional overhead, so cell counts are more portable than a byte estimate.

The discarded rows matter if we later want the alignment itself. A rolling score calculation supplies the optimum value but loses the intermediate choices needed for ordinary traceback. Lesson 4 keeps those choices explicitly.

### Maximization and edit distance

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

## Limitations

1. **Scores remain assumptions.** The recurrence optimizes supplied scores; it does not validate them.
2. **Both ends are charged.** Unrelated flanks can overwhelm a shared region. Lesson 5 changes this rule.
3. **Linear gaps forget history.** A length-$k$ gap contributes $kg$. Lesson 6 adds state for gap opening.
4. **The result is only a score.** Traceback in Lesson 4 recovers an alignment.

## Exercises

::: {#exr-global-count}
Fill $D(i,j)$ through $D(3,3)$. Group the 13 paths to $(2,2)$ by final move.
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
