---
published: true
title: "Problem Set 2: Global Alignment"
subtitle: "Lesson 3 · About 80 minutes"
nocite: |
  @compeau2015
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Any language is fine, and you may adapt the Lesson 3 code and use AI to help write or debug it. Disclose that help briefly and say how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies. Unless a problem says otherwise, score a match $+1$, a mismatch $-1$, and each gap column $g=-2$, as in Lesson 3. Make a first attempt at Problems 1, 2, and 4 without tools; that is exam practice.

## 1. A table with a tie (25 points, about 20 minutes)

Let $x=\texttt{GAAT}$ and $y=\texttt{GAT}$.

a. Fill the global-alignment table $F(i,j)$ by hand, including the boundary row and column, and report $F(4,3)$.

b. Trace back from $(4,3)$ to $(0,0)$. At every cell where more than one predecessor attains the maximum, record all of them. Write out every optimal alignment you find.

c. The optimal alignments disagree about which A in $x$ has no partner in $y$. Explain in two or three sentences why the scoring rule cannot decide between them, and name one kind of additional information that could.

**Evidence to submit:** the filled table with the tied predecessors marked, every optimal alignment, and your explanation.

## 2. Why a table, not a greedy walk (15 points, about 12 minutes)

A classmate proposes a faster method: start at $(0,0)$ and repeatedly take whichever of the three moves (pair the next letters, gap in $y$, gap in $x$) has the best *immediate* column score, preferring a pair when scores tie, until both strings are used up.

a. Run the greedy walk on $x=\texttt{ACGT}$, $y=\texttt{AGT}$ and report its alignment and score. Compare with the optimal score from Lesson 3.

b. The classmate's AI assistant replies that greedy is optimal whenever $|x|=|y|$, because then no gaps are needed. Refute this with two equal-length strings of length at most 4 on which greedy is not optimal. Give both scores and an optimal alignment as your certificate.

c. In two or three sentences, explain which step of the Lesson 3 correctness argument the greedy walk skips. Refer to what a table entry $F(i,j)$ stores.

Adapted from @compeau2015, Ch. 5.

**Evidence to submit:** both greedy runs written as alignments with column scores, the counterexample with its certificate, and the explanation.

## 3. A substitution-score function and an independent check (25 points, about 22 minutes)

Transitions (A↔G, C↔T) are more common than transversions in many genomes. Use the scores
$$
s(a,b)=\begin{cases}+2 & a=b,\\ -1 & \text{transition},\\ -2 & \text{transversion},\end{cases}
\qquad g=-3 .
$$

a. Modify `global_score` from Lesson 3 so that it takes a score function `s(a, b)` and a gap score instead of fixed match and mismatch values. The traceback below is supplied and works with your table; you need not write your own.

```python
# --- supplied for Problem Set 2, Problems 3 and 5 ---
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

b. Write a function `check_alignment(top, bottom, x, y, s, gap, reported)` that does not use the table. It should confirm that removing gaps from `top` and `bottom` gives back `x` and `y`, that no column is a gap over a gap, and that the alignment's column scores sum to `reported`, the score your table gives. State what each of the three checks would catch.

c. Run the table, the traceback, and your check on these three cases and report each score and alignment:
- $x=\texttt{GAAT}$, $y=\texttt{GAT}$;
- $x$ empty, $y=\texttt{ACG}$;
- $x=\texttt{ACGTTGCA}$, $y=\texttt{ATGTCA}$.

d. Before running the third case, predict whether its optimal alignment will pair any transversion. Then compare with the output. If your traceback returns a different alignment from a classmate's with the same score, explain why both can be correct.

**Evidence to submit:** your modified function, the three outputs with your check's verdicts, the prediction made before running, and the comparison.

## 4. When the gap score decides the story (15 points, about 10 minutes)

Let $x=\texttt{CTAGGT}$ and $y=\texttt{CTGGAT}$.

a. With $g=-1$ and with $g=-2$, compute the optimal score and one optimal alignment.

b. Write the score of each of your two alignments as a function of $g$, and find the value of $g$ at which they tie. Explain why the optimum switches there.

c. In counts of substitutions and of insertions or deletions, describe the history each alignment proposes. The scoring rule picks one of these histories; name one kind of evidence beyond the two strings that could support one history over the other.

**Evidence to submit:** both alignments with their scores, the two linear functions of $g$ and the crossing value, and the interpretation.

## 5. Two phage genes for the same protein (20 points, about 16 minutes)

phiX174 and G4 are related phages of *E. coli* with small single-stranded DNA genomes. Each has a gene D that encodes a scaffolding protein used to assemble the capsid. NCBI annotates the gene D coding sequence at positions 390–848 of the [phiX174 genome](../data/phix174-NC_001422.1.fasta) and at positions 1976–2434 of the [G4 genome](../data/g4-NC_001420.2.fasta). Positions start at 1 and include both ends, so in Python the genes are `phix[389:848]` and `g4[1975:2434]`. Both coding sequences have 459 bases: 152 codons followed by a stop codon. Score matches $+1$ and mismatches $-1$ throughout. Read the files with `read_fasta` from Problem Set 1. The helper below summarizes two aligned rows; mismatches are its transitions plus transversions.

```python
# --- supplied for Problem Set 2, Problem 5 ---
PURINES = set("AG")

def summarize_alignment(top: str, bottom: str) -> dict:
    """Count matches, mismatch types, and gaps in two aligned rows."""
    summary = {"matches": 0, "transitions": 0, "transversions": 0, "gaps": []}
    for column, (a, b) in enumerate(zip(top, bottom), start=1):
        if a == "-" or b == "-":
            summary["gaps"].append((column, "top" if a == "-" else "bottom"))
        elif a == b:
            summary["matches"] += 1
        elif (a in PURINES) == (b in PURINES):
            summary["transitions"] += 1
        else:
            summary["transversions"] += 1
    return summary
```

a. Align the two genes with your Problem 3 table and the supplied traceback at $g=-2$, and run `check_alignment` on the result. Report the score, the numbers of matches, mismatches, and gap columns, and for each gap its row and its 1-based column in the alignment (columns counted including gaps). Also compute the score of the ungapped alignment, which pairs position $i$ of one gene with position $i$ of the other.

b. Show that any alignment of two strings of equal length has as many gaps in the top row as in the bottom row. Your gaps should come in pairs of single-base gaps, one in each row a few columns apart. In a coding gene, what happens to the codons that lie between the two gaps of such a pair?

c. Write the score of your $g=-2$ alignment as a linear function of $g$, and use it to predict the gap score below which the ungapped alignment should become optimal. Then run $g=-3$ and $g=-4$, and report the optimal score and the number of gap columns for each. Compare these outputs with your prediction and explain any disagreement.

d. In the ungapped alignment, count the mismatches that are transitions and the mismatches that are transversions. If each substitution were equally likely to produce any of the other three bases, what fraction would be transitions? Compare with your count, and give one reason the number of mismatched columns can undercount the substitutions that occurred since the two genes shared an ancestor.

**Evidence to submit:** the alignment summary with your check's verdict, the ungapped score, the counting argument, the score line with your prediction, the outputs at $g=-3$ and $g=-4$, the transition and transversion counts, and the interpretation.

## Submission check

Include enough output to reproduce each conclusion without opening your code. End with two or three sentences naming any code, tools, AI assistance, or people you used, and one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: correct tables and tracebacks, certificates that settle a claim, checks that would catch a specific mistake, predictions compared honestly with outputs, and conclusions limited to what the scores support.

::: {#refs}
**References**
:::
