---
published: false
title: "Problem Set 2: Global Alignment"
subtitle: "Lesson 3 · About 80 minutes"
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Any language is fine, and you may adapt the Lesson 3 code and use AI to help write or debug it. Disclose that help briefly and say how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies. Unless a problem says otherwise, score a match $+1$, a mismatch $-1$, and each gap column $g=-2$, as in Lesson 3. Make a first attempt at the hand calculations without tools; that is exam practice.

## 1. A table with a tie (25 points, about 20 minutes)

Let $x=\texttt{GAAT}$ and $y=\texttt{GAT}$.

a. Fill the global-alignment table $F(i,j)$ by hand, including the boundary row and column, and report $F(4,3)$.

b. Trace back from $(4,3)$ to $(0,0)$. At every cell where more than one predecessor attains the maximum, record all of them. Write out every optimal alignment you find.

c. The optimal alignments disagree about which A in $x$ has no partner in $y$. Explain in two or three sentences why the scoring rule cannot decide between them, and name one kind of additional information that could.

**Evidence to submit:** the filled table with the tied predecessors marked, every optimal alignment, and your explanation.

## 2. Why a table, not a greedy walk (20 points, about 15 minutes)

A classmate proposes a faster method: start at $(0,0)$ and repeatedly take whichever of the three moves (pair the next letters, gap in $y$, gap in $x$) has the best *immediate* column score, preferring a pair when scores tie, until both strings are used up.

a. Run the greedy walk on $x=\texttt{ACGT}$, $y=\texttt{AGT}$ and report its alignment and score. Compare with the optimal score from Lesson 3.

b. The classmate's AI assistant replies that greedy is optimal whenever $|x|=|y|$, because then no gaps are needed. Refute this with two equal-length strings of length at most 4 on which greedy is not optimal. Give both scores and an optimal alignment as your certificate.

c. In two or three sentences, explain which step of the Lesson 3 correctness argument the greedy walk skips. Refer to what a table entry $F(i,j)$ stores.

**Evidence to submit:** both greedy runs written as alignments with column scores, the counterexample with its certificate, and the explanation.

## 3. A substitution-score function and an independent check (35 points, about 30 minutes)

Transitions (A↔G, C↔T) are more common than transversions in many genomes. Use the scores
$$
s(a,b)=\begin{cases}+2 & a=b,\\ -1 & \text{transition},\\ -2 & \text{transversion},\end{cases}
\qquad g=-3 .
$$

a. Modify `global_score` from Lesson 3 so that it takes a score function `s(a, b)` and a gap score instead of fixed match and mismatch values. The traceback below is supplied and works with your table; you need not write your own.

```python
# --- supplied for Problem Set 2, Problem 3 ---
def traceback(table, x, y, s, gap):
    """Recover one optimal global alignment from a filled score table."""
    i, j, top, bottom = len(x), len(y), [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + s(x[i - 1], y[j - 1]):
            top.append(x[i - 1]); bottom.append(y[j - 1]); i -= 1; j -= 1
        elif i > 0 and table[i][j] == table[i - 1][j] + gap:
            top.append(x[i - 1]); bottom.append("-"); i -= 1
        else:
            top.append("-"); bottom.append(y[j - 1]); j -= 1
    return "".join(reversed(top)), "".join(reversed(bottom))
```

b. Write a function `check_alignment(top, bottom, x, y, s, gap)` that does not use the table. It should confirm that removing gaps from `top` and `bottom` gives back `x` and `y`, that no column is a gap over a gap, and that the alignment's column scores sum to a reported score. State what each of the three checks would catch.

c. Run the table, the traceback, and your check on these three cases and report each score and alignment:
- $x=\texttt{GAAT}$, $y=\texttt{GAT}$;
- $x$ empty, $y=\texttt{ACG}$;
- $x=\texttt{ACGTTGCA}$, $y=\texttt{ATGTCA}$.

d. Before running the third case, predict whether its optimal alignment will pair any transversion. Then compare with the output. If your traceback returns a different alignment from a classmate's with the same score, explain why both can be correct.

**Evidence to submit:** your modified function, the three outputs with your check's verdicts, the prediction made before running, and the comparison.

## 4. When the gap score decides the story (20 points, about 15 minutes)

Let $x=\texttt{CTAGGT}$ and $y=\texttt{CTGGAT}$.

a. With $g=-1$ and with $g=-2$, compute the optimal score and one optimal alignment (by hand or with your code).

b. Write the score of each of your two alignments as a function of $g$, and find the value of $g$ at which they tie. Explain why the optimum switches there.

c. In biological terms, describe the history each alignment proposes, in counts of substitutions and insertions or deletions. A colleague says the $g=-1$ alignment "shows the two sequences differ by one deletion and one insertion." State what the alignment does support, and what it cannot establish on its own.

**Evidence to submit:** both alignments with their scores, the two linear functions of $g$ and the crossing value, and the interpretation.

## Submission check

Include enough output to reproduce each conclusion without opening your code. End with two or three sentences naming any code, tools, AI assistance, or people you used, and one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: correct tables and tracebacks, certificates that settle a claim, checks that would catch a specific mistake, and conclusions limited to what the scores support.
