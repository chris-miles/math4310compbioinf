---
published: true
title: "CpG Islands and Likelihood Ratios"
subtitle: "Lesson 10 · 2027-02-17 · Week 6"
nocite: |
  @durbin1998
---

## Core question

How can two sequence models classify a genomic window?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports"
import math
```

## CpG islands

Sequencing measurements are converted to A, C, G, T base calls and assembled into genomic sequence. Our input is a contiguous interval of that DNA string, called a *window*, with its position in the genome. We want to identify windows whose sequence resembles a CpG island, a region retaining many adjacent CG pairs. Many such regions lie near *promoters*, DNA regions where transcription can begin [@durbin1998]. Candidate locations can guide further study of gene regulation.

CpG denotes C followed by G on the same strand; it is a *dinucleotide*, or pair of neighboring bases. The p denotes their phosphate connection, distinguishing this pair from C paired with G across the two strands. In vertebrate DNA, cytosine in a CpG is often *methylated*, meaning that a methyl group is chemically attached to it. Methylated cytosine can undergo deamination, a chemical change that converts it to thymine, contributing to CpG depletion over evolutionary time. Regions around many promoters are less depleted and often remain unmethylated [@durbin1998].

An ordinary DNA letter string does not record attached methyl groups or gene activity. We predict a sequence feature, rather than measure methylation or establish that a promoter is active. Today's output is a window score and a threshold decision. It does not locate an island's exact boundaries.

A CpG island is not simply a GC-rich string. The order matters: `GCTA` and `CGTA` have the same base composition, but only the second contains `CG`. A first-order chain is the smallest model that can assign different probabilities to these two strings.

## Two models and one window

For training, suppose a separate annotation lists genomic intervals labeled island or background. Extracting their DNA strings gives two collections of labeled sequence. We fit a first-order Markov chain to each collection, then compare the island model $M_+$ and background model $M_-$ on a new window $x=x_1\ldots x_n$. The annotation defines the target we learn; it is not a methylation measurement supplied by the DNA letters.

As in Lesson 9, each chain has an initial distribution $q^\pm$ and a transition matrix $A^\pm$. Initial probabilities sum to one, and each transition row sums to one. We assume every parameter is strictly positive so both window likelihoods and their logarithms are defined.

::: {#def-cpg-classification-problem}
## CpG window classification
Given $x$, two initial distributions $q^+,q^-$, two transition matrices $A^+,A^-$, and a threshold $c$, classify the window as island when
$$
\log_2\frac{P(x\mid M_+)}{P(x\mid M_-)}>c.
$$
:::

The ratio asks which model makes this exact window less surprising. At threshold zero, the classifier chooses the model with larger likelihood.

::: {#def-log-likelihood-ratio}
## Log-likelihood ratio
The *log-likelihood ratio* of $M_+$ against $M_-$ is
$$
S(x)=\log_2\frac{P(x\mid M_+)}{P(x\mid M_-)}.
$$
Positive scores favor $M_+$; negative scores favor $M_-$; a score of $r$ bits means a likelihood ratio of $2^r$.
:::

::: {#prp-markov-ratio-adds}
## Markov log ratios are additive
For two first-order Markov chains with positive parameters, the score is the initial-base log ratio plus one fixed score for each adjacent pair:
$$
S(x)=\log_2\frac{q^+_{x_1}}{q^-_{x_1}}
+\sum_{i=2}^{n}\log_2
\frac{a^+_{x_{i-1}x_i}}{a^-_{x_{i-1}x_i}}.
$$
:::

::: {.proof}
Each likelihood factors as $P(x\mid M_\pm)=q^\pm_{x_1}\prod_{i=2}^n a^\pm_{x_{i-1}x_i}$. Dividing gives one initial ratio and a product of transition ratios. Taking the logarithm turns that product into the displayed sum.
:::

This uses the same product-to-sum rule as Lesson 2, but adjacent pairs now overlap. We do not assume the pairs themselves are independent; the Markov assumption supplies the factorization. A CG pair contributes evidence according to how often G follows C in the two fitted classes.

## A synthetic CpG classifier

Use the illustrative models below. The island model has the same row for every current base. The background model makes $G$ rare after $C$.

$$
A^+=
\begin{array}{c|rrrr}
 & A&C&G&T\\\hline
A&.2&.3&.3&.2\\
C&.2&.3&.3&.2\\
G&.2&.3&.3&.2\\
T&.2&.3&.3&.2
\end{array}
\qquad
A^-=
\begin{array}{c|rrrr}
 & A&C&G&T\\\hline
A&.3&.2&.2&.3\\
C&.3&.3&.1&.3\\
G&.3&.2&.2&.3\\
T&.3&.2&.2&.3
\end{array}.
$$

We take $q_a^+=q_a^-=1/4$ for every base, so the initial term cancels. Write $\beta(a,b)=\log_2(a^+_{ab}/a^-_{ab})$ for a pair contribution. The resulting table is

| pair begins with | $A$ | $C$ | $G$ | $T$ |
|---|---:|---:|---:|---:|
| $A$ | -0.585 | 0.585 | 0.585 | -0.585 |
| $C$ | -0.585 | 0.000 | 1.585 | -0.585 |
| $G$ | -0.585 | 0.585 | 0.585 | -0.585 |
| $T$ | -0.585 | 0.585 | 0.585 | -0.585 |

Rows denote the first base of a pair and columns the second. For example,
$$
\beta(C,G)=\log_2(0.3/0.1)=\log_2 3\approx1.585.
$$

::: {#exm-cpg-two-windows}
## Equal composition, opposite classifications
For `GCTA`, the adjacent pairs are `GC`, `CT`, and `TA`:
$$
S(\texttt{GCTA})=0.585-0.585-0.585=-0.585.
$$
For `CGTA`, the pairs are `CG`, `GT`, and `TA`:
$$
S(\texttt{CGTA})=1.585-0.585-0.585=0.415.
$$
At threshold zero, the first window is background-like and the second island-like under these models. Base composition alone cannot make this distinction. Neither four-base decision identifies a real CpG island or reports whether these cytosines are methylated.
:::

::: {#exm-cpg-threshold}
## The threshold changes the decision
The likelihood ratio for `CGTA` is $2^{0.415}\approx1.33$. Threshold zero calls it an island. Threshold one bit requires the island model to be more than twice as likely, so the same window is classified as background.
:::

The parameters and four-base windows are teaching examples. Real estimates come from much longer labeled island and flanking sequences.

## Scoring windows

::: {.algorithm}
**Algorithm: score a window by log odds.**

1. Set the score to the initial-base log ratio.
2. For every adjacent pair $x_{i-1}x_i$, add its transition log ratio.
3. Compare the total with threshold $c$.
:::

```{.python .execute}
#| label: code-cpg-likelihood-score
def markov_log_ratio(
    sequence: str,
    initial_plus: dict[str, float],
    initial_minus: dict[str, float],
    transition_plus: dict[str, dict[str, float]],
    transition_minus: dict[str, dict[str, float]],
) -> float:
    """Return the base-2 log likelihood ratio of two chains."""
    if not sequence:
        return 0.0
    first = sequence[0]
    score = math.log2(initial_plus[first] / initial_minus[first])
    for previous, current in zip(sequence, sequence[1:]):
        numerator = transition_plus[previous][current]
        denominator = transition_minus[previous][current]
        score += math.log2(numerator / denominator)
    return score


plus_row = {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2}
minus_row = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
transition_plus = {base: plus_row.copy() for base in "ACGT"}
transition_minus = {base: minus_row.copy() for base in "ACGT"}
transition_minus["C"] = {"A": 0.3, "C": 0.3, "G": 0.1, "T": 0.3}
initial_plus = {base: 0.25 for base in "ACGT"}
initial_minus = initial_plus.copy()

[
    markov_log_ratio(
        sequence,
        initial_plus,
        initial_minus,
        transition_plus,
        transition_minus,
    )
    for sequence in ["GCTA", "CGTA"]
]
```

One length-$w$ window costs $\mathcal{O}(w)$ time. For a genome of length $N$ and $1\leq w\leq N$, separately scoring its $N-w+1$ windows costs $\mathcal{O}(Nw)$. A rolling calculation can reuse the contributions from overlapping pairs.

On the string `GCTAA`, shift the four-base window from `GCTA` to `CTAA`. Its pairs change from GC, CT, TA to CT, TA, AA. With equal initial distributions, only GC leaves and AA enters:
$$
S(\texttt{CTAA})=S(\texttt{GCTA})-\beta(G,C)+\beta(A,A)
\approx-0.585-0.585-0.585=-1.755.
$$
The shift removes a pair favoring the island model and adds one favoring background, lowering the evidence for an island-like region.

::: {.algorithm}
**Algorithm: update an overlapping window score.**
For $w\geq2$, subtract the outgoing pair's score and add the incoming pair's score. Retain the shared pair contributions. If the initial distributions differ, also replace the old initial-base log ratio by the new one.
:::

```{.python .execute}
#| label: code-cpg-rolling-check
old_score = markov_log_ratio(
    "GCTA", initial_plus, initial_minus, transition_plus, transition_minus
)
outgoing = math.log2(transition_plus["G"]["C"] / transition_minus["G"]["C"])
incoming = math.log2(transition_plus["A"]["A"] / transition_minus["A"]["A"])
rolling_score = old_score - outgoing + incoming
direct_score = markov_log_ratio(
    "CTAA", initial_plus, initial_minus, transition_plus, transition_minus
)
assert math.isclose(rolling_score, direct_score)
round(rolling_score, 6), round(direct_score, 6)
```

After scoring the first window, each shift requires constant work, giving $\mathcal{O}(N)$ total time under constant-time arithmetic and lookups. Width-one windows have no pairs and use only the initial-base score.

## Likelihood ratios and prior odds

A prior probability describes the class frequency before examining the window. A posterior probability describes it after using the observed sequence. Bayes' rule gives
$$
\frac{P(M_+\mid x)}{P(M_-\mid x)}
=\frac{P(x\mid M_+)}{P(x\mid M_-)}
\frac{P(M_+)}{P(M_-)}.
$$
The common factor $P(x)$ cancels. The prior therefore matters even when we can compute both likelihoods exactly.

A threshold also allows unequal class prevalence. Under equal costs for the two kinds of classification error, choosing the more probable class by Bayes' rule favors the island class when
$$
S(x)>\log_2\frac{P(M_-)}{P(M_+)}.
$$
If islands are rare, the right side is positive. A window needs stronger sequence evidence to overcome the prior odds. In practice, the threshold can instead be selected on validation data to meet a desired false-positive rate.

For example, if the prior island probability is 0.1, the prior odds are $1/9$. The likelihood ratio $4/3$ for `CGTA` changes the odds to $4/27$, giving posterior probability $4/31\approx0.129$. Evidence in favor of the island model has raised its probability, but background remains the more probable class.

## Estimation and evaluation

Training counts must be kept separate by class. With pseudocount $\alpha>0$,
$$
\widehat a^{+}_{ab}
=\frac{n^{+}_{ab}+\alpha}
{\sum_c n^{+}_{ac}+4\alpha},
$$
and similarly for the background. Here $n^+_{ab}$ counts occurrences of adjacent pair $ab$ in island-labeled training fragments; the denominator counts all departures from $a$. Do not count pairs across separate fragments, where adjacency is unknown. Pseudocounts give unseen transitions positive probability; initial distributions likewise need positive estimates to keep the full score finite.

Evaluate the full procedure on held-out chromosomes or long held-out genomic regions. Randomly splitting overlapping windows leaks nearly identical sequence into training and test sets and makes accuracy look better than it is. A confusion matrix counts true positives, false positives, true negatives, and false negatives. Sensitivity is the fraction of labeled islands detected, $TP/(TP+FN)$; specificity is the fraction of labeled background windows rejected, $TN/(TN+FP)$. Report both at the chosen threshold. If the threshold was selected on the test set, the test set is no longer an independent evaluation.

::: {#exm-cpg-confusion}
## Evaluating a chosen threshold
Suppose a held-out set contains 20 labeled island windows and 80 background windows. At one threshold, the classifier detects 16 islands, misses 4, and incorrectly calls 8 background windows islands. Then
$$
\text{sensitivity}=\frac{16}{20}=0.80,
\qquad
\text{specificity}=\frac{72}{80}=0.90.
$$
Raising the threshold can remove false positives, but it can also turn true positives into false negatives. There is no threshold-free statement that the classifier has 90% accuracy in every application: performance also changes with window length, genome, and the frequency of islands in the evaluated set.
:::

A useful comparison shuffles each window many times while preserving its base counts. This changes dinucleotide order while keeping GC content fixed. Comparing original scores with the shuffled score distribution measures how much the ordering contributes beyond composition. Similar classification accuracy alone is inconclusive: a threshold can conceal substantial changes in individual scores.

::: {.callout-note collapse="true" title="Going deeper: precision and class prevalence"}

Sensitivity and specificity condition on the true class. A user examining a list of predicted islands asks a different question: what fraction of that list is correct? This fraction is called *precision*, $TP/(TP+FP)$ when at least one positive prediction is made.

In @exm-cpg-confusion, precision is $16/(16+8)=2/3$, even though specificity is 0.90. To see the role of prevalence, suppose sensitivity remains 0.80 and specificity remains 0.90 in a population of 1,000 windows containing only 10 islands. The expected counts are 8 true positives and 99 false positives. Their ratio gives a precision of about $8/107=0.075$ in the corresponding population calculation.

More generally, let $p$ denote island prevalence, $s$ sensitivity, and $t$ specificity. A randomly selected window is a correctly predicted island with probability $sp$, and a falsely predicted island with probability $(1-t)(1-p)$. Therefore
$$
P(\text{island}\mid\text{positive prediction})
=\frac{sp}{sp+(1-t)(1-p)}.
$$
This calculation assumes the class-conditional performance transfers to the new population. It explains why evaluating on a balanced sample can misrepresent the usefulness of a genome-wide candidate list. To interpret that list, we need the class frequency as well as performance within each class. Raising the threshold may improve precision, but its new sensitivity and specificity must be measured rather than assumed.
:::

## Limitations

1. **Fixed windows blur boundaries.** The classifier labels a whole window and cannot say where an island begins. Lesson 11 makes the region label a state that changes along the sequence.
2. **Window length changes the evidence.** A score of 0.1 bit per pair becomes about 20 bits over 200 bases but is weak evidence over ten bases. Raw scores from different lengths are not directly comparable.
3. **Training labels define the target.** A chain trained on one annotation rule or species may not transfer to another genome with different composition.
4. **Related windows leak information.** Overlapping windows and homologous regions must stay in the same data split. Otherwise held-out performance measures memorized neighbors.
5. **A first-order chain sees only dinucleotides.** It cannot represent longer motifs, distance to a promoter, or experimental methylation measurements.

## Exercises

::: {#exr-cpg-table}
Recompute the four entries in row $C$ of the log-odds table. Explain why $CC$ contributes zero bits.
:::

::: {#exr-cpg-score}
Compute $S(\texttt{AAAA})$ and classify it at thresholds $-2$, $0$, and $1$ bit.
:::

::: {#exr-cpg-prior}
Suppose $P(M_+)=0.1$ and $P(M_-)=0.9$. Derive the Bayes threshold in bits and classify `CGTA`.
:::

::: {#exr-cpg-leakage}
A 500-base region is converted to all overlapping 100-base windows and the windows are randomly split. Explain why this is data leakage and propose a valid split.
:::

::: {#refs}
**References**
:::
