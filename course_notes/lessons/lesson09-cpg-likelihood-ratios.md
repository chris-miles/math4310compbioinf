---
published: true
title: "CpG Islands and Likelihood Ratios"
subtitle: "Lesson 9 · 2027-02-10 · Week 5"
nocite: |
  @durbin1998
---

## Core question

Which stretches of DNA retain unusually many neighboring CG pairs?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import math
import random
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## CpG islands

Many genes have nearby *promoters*, DNA regions where the cell's machinery can begin copying DNA into RNA. Some promoters lie in *CpG islands*, stretches of DNA that retain many adjacent CG pairs compared with the surrounding genome [@durbin1998]. If we have the genome sequence but no map of these regions, where should we look?

Sequencing measurements are converted to A, C, G, T base calls and assembled into genomic sequence. Our input is a contiguous interval of that DNA string, called a *window*, with its position in the genome. A record might contain `CGTA` and the coordinates of its four bases. We want to score how closely such a window resembles an island rather than surrounding sequence, then mark candidates for further study of gene regulation.

CpG means C followed by G on the same strand, a neighboring pair rather than a base pair across the two DNA strands. Why is this pair unusual? In vertebrate DNA, C at a CpG often carries a chemical modification called *methylation*. The modified C is more prone to conversion to T, depleting CG pairs over evolutionary time. Many promoter regions escape this depletion [@durbin1998].

An ordinary DNA letter string does not record attached methyl groups or gene activity. We predict a sequence feature, rather than measure methylation or establish that a promoter is active. Today's output is a window score and a threshold decision. It does not locate an island's exact boundaries.

A CpG island is not simply a GC-rich string. The order matters: `GCTA` and `CGTA` have the same base composition, but only the second contains `CG`. An independent-base model gives these strings the same probability. A first-order Markov chain can distinguish them by using the previous base when assigning the next-base probability. We fit one chain to island sequence and another to background sequence, then compare their explanations of the observed window.

## Two models and one window

For training, suppose a separate annotation lists genomic intervals labeled island or background. Extracting their DNA strings gives two collections of labeled sequence. We fit a first-order Markov chain to each collection, then compare the island model $M_+$ and background model $M_-$ on a new window $x=x_1\ldots x_n$. The annotation defines the target we learn; it is not a methylation measurement supplied by the DNA letters.

As in Lesson 8, each chain has an initial distribution $q^\pm$ and a transition matrix $A^\pm$. Initial probabilities sum to one, and each transition row sums to one. We assume every parameter is strictly positive so both window likelihoods and their logarithms are defined.

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

Use the illustrative models below. Both models have neighboring-base dependence: the island model favors G especially after C, while the background model suppresses that transition. They also differ in their overall preferences for C and G. The numbers are chosen for hand calculation, not fitted to a genome.

$$
A^+=
\begin{array}{c|rrrr}
 & A&C&G&T\\\hline
A&.2&.3&.3&.2\\
C&.1&.3&.5&.1\\
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
| $C$ | -1.585 | 0.000 | 2.322 | -1.585 |
| $G$ | -0.585 | 0.585 | 0.585 | -0.585 |
| $T$ | -0.585 | 0.585 | 0.585 | -0.585 |

Rows denote the first base of a pair and columns the second. A pair receives a positive contribution when it is more likely in the island model and a negative contribution when it is more likely in background. For example,
$$
\beta(C,G)=\log_2(0.5/0.1)=\log_2 5\approx2.322.
$$

The remaining C-row entries explain how that preference changes the whole model. CC contributes zero because both models give it probability 0.3. CA and CT each contribute $\log_2(0.1/0.3)\approx-1.585$ bits: these exits from C favor background. Outside the C row, the $\pm0.585$ scores reflect the models' differing base preferences. The score combines composition and order, with an additional effect specifically after C.

::: {#exm-cpg-two-windows}
## Equal composition, opposite classifications
For `GCTA`, the adjacent pairs are `GC`, `CT`, and `TA`:
$$
S(\texttt{GCTA})=0.585-1.585-0.585=-1.585.
$$
For `CGTA`, the pairs are `CG`, `GT`, and `TA`:
$$
S(\texttt{CGTA})=2.322-0.585-0.585=1.152.
$$
At threshold zero, the first window is background-like and the second island-like under these models. Their 2.737-bit score difference comes entirely from order because their base counts agree. The classifications describe which synthetic model better fits each window.
:::

@fig-cpg-order-evidence separates the contributions of the three neighboring pairs. Both windows contain one of each base, but putting G after C supplies a large positive contribution in the second window. The model therefore distinguishes a neighboring-base pattern that a base-count comparison misses.

```{.python .execute}
#| label: fig-cpg-order-evidence
#| echo: false
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "Pair contributions for GCTA and CGTA under the synthetic island and background chains. The red CG bar contributes log2(5) bits; changing order reverses the total score even though the two windows have identical base counts."
#| fig-alt: "Two bar plots with a shared bit scale. GCTA has contributions 0.585, minus1.585, minus0.585 for GC, CT, TA; CGTA has 2.322, minus0.585, minus0.585 for CG, GT, TA."
pair_ratios = {"GC": 1.5, "CT": 1 / 3, "TA": 2 / 3,
               "CG": 5, "GT": 2 / 3}
fig, axes = plt.subplots(1, 2, sharey=True)
for ax, window in zip(axes, ["GCTA", "CGTA"]):
    pairs = [window[i:i + 2] for i in range(len(window) - 1)]
    contributions = [math.log2(pair_ratios[pair]) for pair in pairs]
    colors = [cf.EMPHASIS if pair == "CG" else "C0" for pair in pairs]
    ax.bar(pairs, contributions, color=colors)
    ax.axhline(0, color=cf.RULE, linewidth=0.8)
    ax.set(xlabel="Neighboring pair", ylim=(-2, 2.7))
    ax.set_title(f"{window}: total {sum(contributions):.3f} bits")
axes[0].set_ylabel("Island-versus-background evidence (bits)")
for ax, label in zip(axes, ["A", "B"]):
    ax.text(-0.12, 1.04, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold")
plt.show()
```

::: {#exm-cpg-threshold}
## The threshold changes the decision
The likelihood ratio for `CGTA` is $2^{1.152}\approx2.22$. Threshold zero calls it an island. Threshold two bits requires the island model to be more than four times as likely, so the same window is classified as background.
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
transition_plus["C"] = {"A": 0.1, "C": 0.3, "G": 0.5, "T": 0.1}
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

One length-$w$ window costs $\mathcal{O}(w)$ time: it has one initial-base contribution and $w-1$ pair contributions. Every pair can strengthen or weaken the accumulated evidence.

## Likelihood ratios and prior odds

Recall the prior-odds calculation from Lesson 2. Here the prior $P(M_+)$ is the fraction of windows expected to be island-like before inspecting their letters, for the window-sampling scheme we plan to use. The posterior $P(M_+\mid x)$ is the class probability after inspecting this window. Bayes' rule gives
$$
\frac{P(M_+\mid x)}{P(M_-\mid x)}
=\frac{P(x\mid M_+)}{P(x\mid M_-)}
\frac{P(M_+)}{P(M_-)}.
$$
Both posterior probabilities divide by $P(x)$, the overall probability of observing the window after allowing either class. Taking their ratio cancels that common denominator. The sequence likelihood ratio therefore multiplies the starting class odds.

A threshold also allows unequal class prevalence. Under equal costs for the two kinds of classification error, choosing the more probable class by Bayes' rule favors the island class when
$$
S(x)>\log_2\frac{P(M_-)}{P(M_+)}.
$$
If islands are rare, the right side is positive. A window needs stronger sequence evidence to overcome the prior odds. In practice, the threshold can instead be selected on validation data to meet a desired false-positive rate.

For example, if the prior island probability is 0.1, the prior odds are $1/9$. The exact likelihood ratio $20/9$ for `CGTA` changes the odds to $20/81$, giving posterior probability $20/101\approx0.198$. The letters nearly double the island probability, but background remains the more probable class.

## Evidence accumulates along the sequence

A four-base window provides only three pair contributions. Longer windows have more opportunities to distinguish the chains. If a window has $w-1$ pairs with average contribution $\bar\beta$, then, with our equal initial distributions,
$$
S(x)=(w-1)\bar\beta.
$$
For example, an average of 0.1 bit per pair yields 0.9 bit in ten bases but 19.9 bits in 200 bases. Those scores mean likelihood ratios of about 1.87 and nearly one million. They describe two whole-window hypotheses, so the gain in evidence assumes the whole window belongs to the same regime.

We can see the accumulation by simulating a sequence with background, then an island, then background again. At each new position, use the row of the current region's chain. Both chains are the synthetic models already scored by hand.

::: {.algorithm}
**Algorithm: simulate a changing region and accumulate evidence.**
Draw the first base from the shared initial distribution. For each later position, select the island or background transition matrix according to its known simulation label and draw from the row of the previous base. Append the base and add its pair log ratio to the cumulative score.
:::

```{.python .execute}
#| label: fig-cpg-cumulative-evidence
#| code-fold: true
#| fig-cap: "Cumulative island-versus-background log ratio for a simulated sequence; the shaded interval uses the island chain. Upward and downward trends show where successive pairs favor each model."
#| fig-alt: "A cumulative score curve falls through background, rises within the shaded island interval from bases 201 to 400, and falls again in background."
rng = random.Random(0)
alphabet = list(initial_plus)
sequence = rng.choices(alphabet, weights=list(initial_plus.values()), k=1)
cumulative_score = [0.0]
for position in range(1, 800):
    region_matrix = (transition_plus if 200 <= position < 400
                     else transition_minus)
    previous = sequence[-1]
    row = region_matrix[previous]
    current = rng.choices(alphabet, weights=[row[base] for base in alphabet], k=1)[0]
    sequence.append(current)
    contribution = math.log2(transition_plus[previous][current]
                              / transition_minus[previous][current])
    cumulative_score.append(cumulative_score[-1] + contribution)

fig, ax = plt.subplots()
ax.axvspan(200.5, 400.5, color=cf.LIGHT, label="Simulated island")
ax.plot(range(1, 801), cumulative_score)
ax.set_xlabel("Base position")
ax.set_ylabel("Cumulative log ratio (bits)")
ax.legend()
plt.show()
```

The curve need not rise at every island base: CA and CT still favor background there. Its trend reflects the balance of all pair contributions. The height at position $i$ scores the whole prefix $x_1\ldots x_i$ under two homogeneous chains; the change in slope shows where local evidence changes. A positive prefix score is therefore not a label for its final base.

Window width creates a related compromise. Short windows respond quickly to a change of region but contain few pairs, so individual contributions can reverse a call. Long windows accumulate stronger evidence within a region but mix the two regimes near a boundary. For instance, a 100-base window straddling a boundary could contain 50 background bases and 50 island bases; either whole-window label assigns one class to bases from the other. Allowing the region label to change inside the sequence is the problem taken up by hidden Markov models.

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
Raising the threshold can remove false positives, but it can also turn true positives into false negatives. The 0.90 specificity describes background windows at this threshold. It does not describe all windows: 88 of the 100 windows were classified correctly, giving accuracy 0.88.
:::

The labels used for fitting define the target, and the model sees only adjacent pairs. Longer motifs, promoter context, and experimental methylation measurements would require additional information.

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
