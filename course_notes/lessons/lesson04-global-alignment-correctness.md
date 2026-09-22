---
published: true
title: "Global Alignment: Traceback and Correctness"
subtitle: "Lesson 4 · 2027-01-25 · Week 3"
nocite: |
  @durbin1998, @compeau2015, @needleman1970
---

## Core question

Why does the recurrence find the best score, and how do we recover an alignment that achieves it?

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

## The table as a directed graph

A score of 1 does not tell a biologist where a deletion occurred. To inspect a predicted sequence change, we need the aligned letters and gaps that produced the score. We also need to know whether several equally good explanations place that gap differently.

The data and scoring model are unchanged from Lesson 3. Today's output is one optimal alignment together with checks that its printed columns recover the input strings and reproduce the reported score. The proof explains why the table finds an optimum; the checks help us catch code that fails to implement that proof.

Keep the model from Lesson 3. For strings $x=x_1\cdots x_n$ and $y=y_1\cdots y_m$, $F(i,j)$ is the best score of a global alignment of their length-$i$ and length-$j$ prefixes. Each cell is a vertex. Its three incoming edges represent the three possible final columns.

::: {#def-backpointer}
## Backpointer
A *backpointer* from $(i,j)$ records a predecessor whose candidate attains the maximum defining $F(i,j)$. A cell can have more than one backpointer when candidates tie.
:::

The score table gives an optimum value. A path of backpointers gives an alignment that attains it.

## Correctness of the recurrence

::: {#thm-global-correctness}
## Correctness of the global-alignment recurrence
For every $0\le i\le n$ and $0\le j\le m$, $F(i,j)$ equals the maximum score among all global alignments of $x_1\cdots x_i$ and $y_1\cdots y_j$.
:::

::: {.proof}
We use induction on $i+j$. The claim holds on the boundaries because there is only one alignment of a nonempty prefix with an empty string: every letter is paired with a gap.

For $i,j>0$, take an optimal alignment $A$. Its final column is $x_i/y_j$, $x_i/-$, or $-/y_j$. Removing that column leaves a valid alignment for the corresponding predecessor cell. It must be optimal for those prefixes. Otherwise, replacing it with a better prefix alignment would improve $A$, a contradiction. Thus $A$ has the score of one candidate in the recurrence.

Conversely, append the appropriate final column to an optimal alignment from each predecessor. Each candidate is attained by a valid alignment. Their maximum is therefore attainable and cannot be smaller than the optimum.
:::

The proof uses two facts: the final-column cases are exhaustive, and an optimal solution cannot contain a suboptimal prefix. The second fact is called *optimal substructure*.

## Traceback

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

At $(2,1)$, C meets A. The candidates are $F(1,0)-1=-3$, $F(1,1)-2=-1$, and $F(2,0)-2=-6$. The vertical pointer wins, so we emit C/-. The remaining diagonal emits A/A. Traceback generates columns from right to left, which is why the implementation reverses both lists.

Every move decreases $i$, $j$, or both. The process reaches $(0,0)$ after at most $n+m$ moves. A pointer that leaves both indices unchanged creates an infinite loop; one that increases an index was stored in the wrong direction.

## Ties and optimal alignments

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

## Implementation

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
```

This check is independent of the table recurrence. It catches a traceback that emits the wrong letter, forgets a leading gap, or reverses only one row. We also remove gaps from each row and confirm that the original strings are recovered.

### Boundary and traceback tests

::: {.callout-tip title="Check: boundaries, indices, and traceback"}
Test "A" against the empty string to check boundary penalties. Test "A" against "A" to catch an off-by-one error between table index $i$ and string index $i-1$. Test "AA" against "A" to make sure traceback continues until both indices reach zero.
:::

```{.python .execute}
#| label: code-global-tests
assert global_alignment("A", "")[0] == -2
assert global_alignment("A", "A")[0] == 1
score, top, bottom = global_alignment("AA", "A")
assert score == -1
assert top.replace("-", "") == "AA"
assert bottom.replace("-", "") == "A"

score, top, bottom = global_alignment("ACGT", "AGT")
assert score_alignment(top, bottom) == score
assert top.replace("-", "") == "ACGT"
assert bottom.replace("-", "") == "AGT"
```

Stopping when $i=0$ *or* $j=0$ drops the unconsumed prefix. Initializing the borders to zero silently changes global alignment into a free-end variant. Both bugs can pass equal-length examples, so the tests include empty and unequal-length strings.

### Counting optimal alignments

If biological ambiguity matters, storing one pointer is insufficient. A second table can count optimal paths. Set $C(0,0)=1$. At each cell, sum $C$ over every predecessor whose candidate equals the optimal score. Boundary cells have count 1.

For AA against A with match $+1$ and gap $-2$, the score is $-1$ and the count is 2: AA/A- and AA/-A. The count tells us how many optimal paths exist, while the pointer graph shows where they diverge.

::: {#prp-count-optima}
## Counting optimal paths
If $C(i,j)$ sums the counts of exactly those predecessors that attain $F(i,j)$, then $C(i,j)$ equals the number of optimal alignment paths to $(i,j)$.
:::

::: {.proof}
Induct on $i+j$. Every optimal path reaches $(i,j)$ through one maximizing predecessor, and appending that final move to every optimal predecessor path produces every optimal path to $(i,j)$. Distinct final moves make the groups disjoint.
:::

A large count does not mean the score is uncertain. It means the scoring model cannot distinguish many explanations. Reporting one alignment without the count can hide that ambiguity.

### Cost and memory

The fill performs constant work in $(n+1)(m+1)$ cells, and traceback takes at most $n+m$ moves. The total is $\mathcal O(nm)$ time and $\mathcal O(nm)$ memory. If only the score is needed, two rows suffice, reducing memory to $\mathcal O(m)$. Hirschberg's algorithm recovers an alignment in linear memory by recomputing selected table slices.

## Limitations

1. **One traceback does not establish uniqueness.** A tie policy hides alternative optima unless the program stores or counts all maximizing predecessors.
2. **Correctness is conditional on additivity.** The proof removes one final column whose score does not depend on earlier columns. Affine gaps require extra state.
3. **Quadratic work remains expensive.** One table is routine; one table for each record in a large database is not. Lesson 6 trades guaranteed optimality for speed.

## Exercises

::: {#exr-traceback-hand}
Using the Lesson 3 table for ACGT and AGT, mark every maximizing predecessor and perform the traceback. Verify the score by summing column scores.
:::

::: {#exr-traceback-ties}
Build the table for AA and A with match $+1$, mismatch $-1$, and gap $-2$. Identify both optimal alignments.
:::

::: {#exr-traceback-bug}
Change the traceback condition from "or" to "and". Predict the output for AA against A, then explain why its printed alignment cannot have the reported score.
:::

::: {#exr-correctness-template}
Suppose diagonal moves are forbidden. Write the recurrence and adapt the proof of @thm-global-correctness. Which sentence changes?
:::

::: {#refs}
**References**
:::
