---
published: false
title: "Problem Set 6: Hidden Regions, Decoding, and Uncertainty"
subtitle: "Lessons 10–12 · About 80 minutes"
nocite: |
  @durbin1998, @mit6047
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Python, R, Julia, MATLAB, or another language is fine. Working Python for Viterbi decoding and forward–backward appears in Lessons 11 and 12; you may copy, translate, or adapt it. You are not asked to write Viterbi, forward, and backward from scratch. Each coding task asks for one small change or one check, and the credit goes to evidence that does not depend on the code being tested. You may use AI to help write or debug code. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies.

All problems use the B/I hidden-region model from Lessons 10 to 12 unless a part says otherwise. Observations are $S$ for C or G and $W$ for A or T; hidden states are background $B$ and island-like $I$:
$$
p_B=p_I=0.5,
\qquad
A=\begin{array}{c|cc}
 & B&I\\\hline
B&0.8&0.2\\
I&0.4&0.6
\end{array},
\qquad
E=\begin{array}{c|cc}
 & S&W\\\hline
B&0.3&0.7\\
I&0.8&0.2
\end{array}.
$$
Rows of $A$ are departure states: $a_{\ell k}$ is the probability of moving from $\ell$ to $k$.

## 1. Two S symbols or three (25 points, about 20 minutes)

The DNA fragment `ACGT` becomes $x=\texttt{WSSW}$, and `ACGCT` becomes $\texttt{WSSSW}$. Both contain a short GC-rich run between AT bases. Does the model call that run an island?

a. Before filling any table, predict the Viterbi path for each observation. Lesson 11 compared labeling a run of $r$ consecutive $S$ symbols as $I$ against labeling it $B$, with both flanking positions fixed to $B$; the ratio was $2^{r-1}/3$. Use that comparison and give one sentence of reasoning for each prediction.

b. Fill the Viterbi table for $\texttt{WSSW}$ by hand, in ordinary probabilities, as in @exm-viterbi-table. Record a backpointer for every entry at positions 2 to 4, then trace back the path. Verify its joint probability by writing its initial, transition, and emission factors in sequence order and multiplying them. At one position the larger table entry is in row $I$, yet the decoded path is in $B$ there. Identify that position and explain how the final $W$ changes the comparison.

c. The first three columns of your table also apply to $\texttt{WSSSW}$. Add its last two columns, by hand or with Lesson 11's `viterbi` function, and report the decoded path. Compare both decoded paths with your predictions in part (a).

d. For $\texttt{WSSW}$, compute $\gamma_I(2)=P(\pi_2=I\mid x)$ and $\gamma_I(3)$ by enumerating all 16 paths: add the joint probabilities of the paths with $I$ at that position and divide by the total. Confirm both values with Lesson 12's `forward_backward` code. What labels does posterior decoding give at positions 2 and 3? In two or three sentences, explain how the most probable path can use $B$ at positions 2 and 3 while $I$ is the more probable label at each of those positions. Name the two paths that contribute most to $\gamma_I(2)$.

**Evidence to submit:** both predictions with reasons, the hand table with backpointers, the factor product, the enumeration sums behind each posterior, and the explanation.

## 2. Diagnose a forward–backward translation (25 points, about 20 minutes)

A classmate rewrote Lesson 12's forward–backward code and reports that it passes the lesson's checks. Here is the function:

```python
def forward_backward(observations, states, initial, transition, emission):
    """Return forward and backward tables and P(x); this version has a defect."""
    forward = [{k: initial[k] * emission[k][observations[0]] for k in states}]
    for symbol in observations[1:]:
        column = {}
        for current in states:
            total = 0.0
            for previous in states:
                total += forward[-1][previous] * transition[current][previous]
            column[current] = emission[current][symbol] * total
        forward.append(column)

    backward = [{k: 1.0 for k in states} for _ in observations]
    for i in range(len(observations) - 2, -1, -1):
        for current in states:
            total = 0.0
            for following in states:
                total += (transition[following][current]
                          * emission[following][observations[i + 1]]
                          * backward[i + 1][following])
            backward[i][current] = total
    probability = sum(forward[-1].values())
    if probability == 0:
        raise ValueError("Impossible observations or numerical underflow.")
    return forward, backward, probability
```

a. Run the function on $\texttt{SSWS}$ and apply Lesson 12's two checks: every forward–backward cut gives the same $P(x)$, and every posterior column sums to one. Report $P(x)$ and $\gamma_I(i)$ at each position, and compare them with @exm-forward-ssws and @exm-posterior-ssws. Both checks pass even though the numbers differ from the lesson's. What do these two checks actually test?

b. Locate the defect. Quote the lines involved, state which convention for $a_{\ell k}$ they violate, and give corrected lines.

c. Write an independent check that uses only the joint-probability formula from Lesson 10: enumerate every path and add. Run the classmate's function and your check on all observations of length 1, 2, and 3. List the observations for which the returned $P(x)$ is correct and those for which the returned posteriors are correct. Explain why no observation of length 1 can expose this defect.

d. The observations with a correct $P(x)$ share a pattern. Compare the returned $P(\texttt{SW})$ with the true $P(\texttt{SW})$ and the true $P(\texttt{WS})$. State what the defective code computes when $p_B=p_I$, and prove it by reindexing the sum over paths. Would checking that the returned $P(x)$ sums to one over all $2^n$ observations of length $n$ catch this defect in the B/I model? Name one change to the start probabilities that would let that check catch it. For which transition matrices does the defective code return correct answers on every input?

**Evidence to submit:** the outputs from part (a), the corrected lines, the table of length-1 to length-3 comparisons, and the reindexing argument.

## 3. Log-space forward on a phage genome (30 points, about 25 minutes)

Lesson 12 warned that the probability-space code is only for short sequences. Here you make the small change that removes that limit and apply it to a real genome: the 5,386-base phiX174 reference in the course [data folder](https://github.com/chris-miles/math4310compbioinf/tree/main/data). The B/I parameters are illustrative and were not fitted to phiX174, and the phage is not a CpG-island benchmark, so the outputs below test the model and the code rather than annotate the phage.

This code reads the FASTA file and recodes it. Adjust the path to wherever you saved the file.

```python
from pathlib import Path

lines = Path("data/phix174-NC_001422.1.fasta").read_text().splitlines()
genome = "".join(line.strip() for line in lines if not line.startswith(">"))
observations = "".join("S" if base in "CG" else "W" for base in genome)
print(len(observations), observations.count("S"))
```

It should print `5386 2411`.

a. Start from Lesson 11's log-space `viterbi` function. Replace the maximum over previous states with log-sum-exp from Lesson 12, and drop the backpointers, so the function returns $\ln P(x)$. Submit your log-sum-exp helper and the lines you changed. In one sentence, explain why this replacement turns the Viterbi recurrence into the forward recurrence.

b. Check the change without relying on it.
   (i) On $\texttt{SSWS}$, the result should equal $\ln 0.060576$.
   (ii) Choose an observation of length 8 and compare your result with the logarithm of the sum of all 256 enumerated path probabilities.
   (iii) Explain why $\ln P(x)$ can never be smaller than the Viterbi log score, and confirm the inequality on both examples.

c. Run your function and Lesson 11's `viterbi` on the recoded phiX174 sequence. Report $\ln P(x)$, the average $\ln P(x)/n$, and the Viterbi log joint probability. Before running any probability-space code, use the average and the smallest positive double, about $4.9\times10^{-324}$, to predict roughly the position at which a probability-space forward column first becomes all zeros. Then find that position by running the probability-space forward recurrence column by column, and reconcile it with your prediction.

d. Assess each claim in two or three sentences.
   (i) "The Viterbi path holds only a negligible fraction of the conditional probability, so the Viterbi annotation of phiX174 is meaningless." Use your output to state the fraction.
   (ii) "The B/I model gives phiX174 a higher likelihood than a fair coin that emits $S$ and $W$ with probability $1/2$ each, so the phage contains island-like regions." Support your assessment with one more calculation: $\ln P(x)$ under a one-state model that emits $S$ with probability equal to phiX174's overall $S$ fraction.

**Evidence to submit:** the changed lines, the three checks with numbers, the phiX174 outputs, the underflow prediction and actual position, and the two assessments.

## 4. Labels that no path can produce (20 points, about 15 minutes)

Adapted from @durbin1998 and @mit6047.

Region and gene models often forbid some transitions. For example, a model might require every island to begin in a separate boundary state. Posterior decoding (@def-posterior-decoding) chooses each label separately, so nothing forces the reported labels to form a path the model allows.

a. Suppose a two-state HMM has $a_{k\ell}=0$ for some states $k$ and $\ell$ (possibly $k=\ell$). Prove that $\gamma_k(i)>1/2$ and $\gamma_\ell(i+1)>1/2$ cannot both hold. Conclude that when no posterior column contains a tie, two-state posterior decoding never places a forbidden transition between adjacent positions. Hint: what is $P(\pi_i=k,\pi_{i+1}=\ell\mid x)$?

b. Construct a three-state HMM over $\{S,W\}$ with at least one zero transition probability, together with an observation of length at most 4, for which $P(x)>0$ but posterior decoding returns a labeling with joint probability zero. One possible design adds an island-entry state that must be followed by the island state. Certify your instance: give all parameters, the posterior columns, the decoded labels, and an enumeration of all $3^n$ paths showing that the decoded labeling has probability zero. Also report the Viterbi path.

c. In two or three sentences, explain why the argument in part (a) fails for three states, and why your result is not a programming error. If an annotation must be a valid path, what could it report instead?

**Evidence to submit:** the proof, your model parameters, the posterior table, the enumeration, and the explanation.

## Optional enrichment (not graded)

The Going deeper note in Lesson 12 gives the posterior transition probability $\xi_{k\ell}(i)$. For $\texttt{WSSW}$, compute the expected number of uses of each transition, $\sum_i\xi_{k\ell}(i)$, check that the four expected counts sum to $n-1=3$, and normalize each row to obtain one Baum–Welch update of $A$.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: predictions made before computing, hand calculations, checks whose answers do not come from the code being tested, and conclusions limited to what the model supports.

::: {#refs}
**References**
:::
