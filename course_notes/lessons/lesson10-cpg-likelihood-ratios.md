---
published: true
title: "CpG Islands as Likelihood Ratios"
subtitle: "Lesson 10 · Week 6"
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

::: {.biology}
In vertebrate genomes, cytosine in a `CG` dinucleotide is often methylated. Methylated cytosine can deaminate to thymine, so mutation depletes `CG` over evolutionary time. Regions near many promoters often remain unmethylated and retain more `CG` dinucleotides. These regions are called CpG islands; "p" denotes the phosphate joining adjacent bases on one strand [@durbin1998].
:::

A CpG island is not simply a GC-rich string. The order matters: `GCTA` and `CGTA` have the same base composition, but only the second contains `CG`. A first-order chain is the smallest model that can assign different probabilities to these two strings.

## Two models and one window

The data are a DNA window $x=x_1\ldots x_n$. We compare an island model $M_+$ with a background model $M_-$. Each is a first-order Markov chain trained from labeled sequence.

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

For a Markov chain, the score separates into local terms:
$$
S(x)=\log_2\frac{q^+_{x_1}}{q^-_{x_1}}
+\sum_{i=2}^{n}\log_2
\frac{a^+_{x_{i-1}x_i}}{a^-_{x_{i-1}x_i}}.
$$

This is the same reason alignment scores add in Lesson 2. Each adjacent pair contributes one log-odds score.

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

We take equal initial distributions so the initial term cancels. The resulting pair-score table is

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
At threshold zero, the first window is background and the second is island. Base composition alone cannot make this distinction.
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

One length-$w$ window costs $\mathcal{O}(w)$ time. Scoring every overlapping window separately costs $\mathcal{O}(Nw)$ on a genome of length $N$. Because the score is a sum of pair contributions, when the two models have the same initial distribution, a rolling implementation can subtract the pair leaving the window and add the pair entering it, reducing the scan to $\mathcal{O}(N)$. With unequal initial distributions, it must also update the initial-base term.

## Likelihood ratios and prior odds

::: {#prp-markov-ratio-adds}
## Markov log ratios are additive
For two first-order Markov chains with positive parameters, the log-likelihood ratio is the initial-base log ratio plus one fixed score for each adjacent pair.
:::

::: {.proof}
Divide the two Markov factorizations. Products in the numerator and denominator pair term by term. Taking a logarithm changes the resulting product of ratios into a sum.
:::

A threshold also allows unequal class prevalence. Under equal costs for the two kinds of classification error, choosing the more probable class by Bayes' rule favors the island class when
$$
S(x)>\log_2\frac{P(M_-)}{P(M_+)}.
$$
If islands are rare, the right side is positive. A window needs stronger sequence evidence to overcome the prior odds. In practice, the threshold can instead be selected on validation data to meet a desired false-positive rate.

## Estimation and evaluation

Training counts must be kept separate by class. With pseudocount $\alpha>0$,
$$
\widehat a^{+}_{ab}
=\frac{n^{+}_{ab}+\alpha}
{\sum_c n^{+}_{ac}+4\alpha},
$$
and similarly for the background. Pseudocounts keep every score finite.

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

## Further reading

Durbin et al. use two Markov chains to motivate CpG island classification and hidden states in Section 3.1 [@durbin1998].
