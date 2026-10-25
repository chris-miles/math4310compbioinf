---
published: true
title: "CpG Islands and Likelihood Ratios"
subtitle: "Lesson 9 · 2027-02-10 · Week 5"
nocite: |
  @durbin1998, @compeau2015
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

Many genes have nearby *promoters*, DNA regions where the cell's machinery can begin copying DNA into RNA. Some promoters lie in *CpG islands*, stretches of DNA that retain many adjacent CG pairs compared with the surrounding genome [@durbin1998]. CpG means C followed by G on the same strand, a neighboring pair rather than a base pair across the two DNA strands. In vertebrate DNA, the C of a CpG often carries a chemical modification called *methylation*, and a methylated C is prone to mutate to T. Over evolutionary time this depletes CG pairs from most of the genome, in human DNA to about a quarter of the frequency that base composition predicts, while many promoter regions escape the depletion [@compeau2015]. Islands are typically a few hundred to a few thousand bases long [@durbin1998]. Because many islands sit at gene starts, they serve as landmarks for where genes may lie in a mammalian genome that has no annotation yet [@compeau2015].

@fig-cpg-island-track shows what an island looks like along a stretch of DNA. The 1,000 bases are simulated from the two models used later in this lesson, with one 300-base island placed around the start of a gene, and each tick marks a CG pair. Outside the shaded island, CG pairs are scattered tens of bases apart; inside it, they occur every few bases, and a count of CG pairs in a sliding 100-base window rises and falls with the island. The printed excerpts show the same contrast in the letters themselves. If we have the genome sequence but no map of these regions, where should we look?

```{.python .execute}
#| label: fig-cpg-island-track
#| code-fold: true
#| code-summary: "Simulate a stretch of DNA with one island and mark its CG pairs"
#| fig-cap: "A simulated 1,000-base stretch with an island at positions 301 to 600 (shaded). Top: each tick is one CG pair, and the arrow marks an illustrative gene start. Bottom: number of CG pairs in a 100-base window centered at each position. CG pairs cluster inside the island and are sparse outside it."
#| fig-alt: "Two stacked panels sharing a position axis from 1 to 1000. The top panel has sparse vertical ticks outside a shaded band from 301 to 600 and dense ticks inside it, with a right-pointing arrow labeled gene start inside the band. The bottom panel shows a curve near 1 or 2 outside the band rising to about 20 inside it."
#| fig-width: 6
#| fig-height: 3.4
# The two synthetic chains introduced in "A synthetic CpG classifier" below.
plus_row = {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2}
minus_row = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
transition_plus = {base: plus_row.copy() for base in "ACGT"}
transition_plus["C"] = {"A": 0.1, "C": 0.3, "G": 0.5, "T": 0.1}
transition_minus = {base: minus_row.copy() for base in "ACGT"}
transition_minus["C"] = {"A": 0.3, "C": 0.3, "G": 0.1, "T": 0.3}

# Draw each base from the previous base's row of the island chain inside
# positions 301 to 600 and of the background chain elsewhere.
length, island_start, island_end = 1000, 301, 600
rng = random.Random(1)
track = rng.choices("ACGT", k=1)
for position in range(2, length + 1):
    if island_start <= position <= island_end:
        row = transition_plus[track[-1]]
    else:
        row = transition_minus[track[-1]]
    track += rng.choices("ACGT", weights=[row[base] for base in "ACGT"], k=1)
track = "".join(track)

cg_starts = [i + 1 for i in range(length - 1) if track[i:i + 2] == "CG"]
centers = list(range(51, length - 49, 5))
window_counts = [track[c - 51:c + 49].count("CG") for c in centers]

fig, axes = plt.subplots(2, 1, sharex=True, height_ratios=[1, 2])
for ax in axes:
    ax.axvspan(island_start - 0.5, island_end + 0.5, color=cf.LIGHT, zorder=0)
axes[0].vlines(cg_starts, 0, 1, color=cf.INK, linewidth=0.8)
axes[0].annotate("", xy=(560, 1.3), xytext=(470, 1.3),
                 arrowprops=dict(arrowstyle="-|>", color=cf.INK, lw=1))
axes[0].text(465, 1.3, "gene start", ha="right", va="center", fontsize=9)
axes[0].set(yticks=[], ylim=(0, 1.6))
axes[0].spines["left"].set_visible(False)
axes[1].plot(centers, window_counts)
axes[1].set(xlabel="Base position", ylabel="CG pairs per 100 bases")
plt.show()
print("background, positions 101-140:", track[100:140])
print("island,     positions 401-440:", track[400:440])
```

Our input is a contiguous interval of assembled genomic DNA, called a *window*, together with its position in the genome. We want to score how closely a window resembles an island rather than surrounding sequence, then mark candidates for further study of gene regulation. The letters record neither methyl groups nor gene activity, so the score predicts a sequence feature; it does not measure methylation or show that a promoter is active.

A CpG island is not simply a GC-rich string. The order matters: `GCTA` and `CGTA` have the same base composition, but only the second contains `CG`, so an independent-base model gives the two strings the same probability. A first-order Markov chain distinguishes them because each base's probability depends on the base before it. We fit one chain to island sequence and another to background sequence, then ask which chain explains the observed window better.

## Two models and one window

For training, suppose a separate annotation lists genomic intervals labeled island or background. We fit a first-order Markov chain to the DNA strings of each class, giving an island model $M_+$ and a background model $M_-$, and compare them on a new window $x=x_1\ldots x_n$.

As in Lesson 8, each chain has an initial distribution $q^\pm$ and a transition matrix $A^\pm$. Initial probabilities sum to one, and each transition row sums to one. We assume every parameter is strictly positive so both window likelihoods and their logarithms are defined.

::: {#def-log-likelihood-ratio}
## Log-likelihood ratio
The *log-likelihood ratio* of $M_+$ against $M_-$ for a window $x$ is
$$
S(x)=\log_2\frac{P(x\mid M_+)}{P(x\mid M_-)}.
$$
:::

Positive scores favor $M_+$ and negative scores favor $M_-$; a score of $r$ bits means a likelihood ratio of $2^r$. The ratio asks which model makes this exact window less surprising.

::: {#def-cpg-classification-problem}
## CpG window classification
**Input:** a window $x=x_1\ldots x_n$; island and background chains with initial distributions $q^+,q^-$ and transition matrices $A^+,A^-$; and a threshold $c$.

**Output:** the label island if $S(x)>c$, and background otherwise.
:::

At threshold zero, the classifier chooses the model with larger likelihood. For two first-order chains with positive parameters, the score is the initial-base log ratio plus one fixed score for each adjacent pair:
$$
S(x)=\log_2\frac{q^+_{x_1}}{q^-_{x_1}}
+\sum_{i=2}^{n}\log_2
\frac{a^+_{x_{i-1}x_i}}{a^-_{x_{i-1}x_i}}.
$$ {#eq-markov-ratio-adds}

::: {.proof}
Sketch. Each likelihood is a first-base probability times one transition probability per adjacent pair, so their ratio splits factor by factor:
$$
\frac{P(x\mid M_+)}{P(x\mid M_-)}
=\frac{q^+_{x_1}\prod_{i=2}^n a^+_{x_{i-1}x_i}}{q^-_{x_1}\prod_{i=2}^n a^-_{x_{i-1}x_i}}
=\frac{q^+_{x_1}}{q^-_{x_1}}\prod_{i=2}^{n}\frac{a^+_{x_{i-1}x_i}}{a^-_{x_{i-1}x_i}}.
$$
Taking $\log_2$ turns the product into the sum.
:::

This uses the same product-to-sum rule as Lesson 2, but adjacent pairs now overlap. We do not assume the pairs themselves are independent; the Markov assumption supplies the factorization. A CG pair contributes evidence according to how often G follows C in the two fitted classes.

## A synthetic CpG classifier

Use the illustrative models below. Both models have neighboring-base dependence: the island model favors G especially after C, while the background model suppresses that transition. The background chain is the CG-depleted chain simulated in Lesson 8. The two models also differ in their overall preferences for C and G. The numbers are chosen for hand calculation, not fitted to a genome. For comparison, chains fitted to about 60,000 bases of labeled human sequence put G after C with probability 0.27 in islands and 0.08 elsewhere, about 1.8 bits of evidence per CG pair [@durbin1998, pp. 51–52].

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
At threshold zero, the first window is background-like and the second island-like under these models. Their 2.737-bit score difference comes entirely from order because their base counts agree.
:::

@fig-cpg-order-evidence separates the contributions of the three neighboring pairs. Both windows contain one of each base, but putting G after C supplies a large positive contribution in the second window.

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

One length-$w$ window costs $\mathcal{O}(w)$ time: it has one initial-base contribution and $w-1$ pair contributions.

## Likelihood ratios and prior odds

Recall the prior-odds calculation from Lesson 2. Here the prior $P(M_+)$ is the fraction of windows expected to be island-like before inspecting their letters, for the window-sampling scheme we plan to use. The posterior $P(M_+\mid x)$ is the class probability after inspecting this window. Bayes' rule gives
$$
\frac{P(M_+\mid x)}{P(M_-\mid x)}
=\frac{P(x\mid M_+)}{P(x\mid M_-)}
\frac{P(M_+)}{P(M_-)}.
$$
Both posterior probabilities divide by $P(x)$, the overall probability of the window under either class, so that factor cancels and the sequence likelihood ratio multiplies the starting class odds.

A threshold also allows unequal class prevalence. Under equal costs for the two kinds of classification error, choosing the more probable class by Bayes' rule favors the island class when
$$
S(x)>\log_2\frac{P(M_-)}{P(M_+)}.
$$
If islands are rare, the right side is positive, and a window needs stronger sequence evidence to overcome the prior odds. In practice, the threshold can instead be selected on validation data to meet a desired false-positive rate.

For example, if the prior island probability is 0.1, the prior odds are $1/9$. The exact likelihood ratio $20/9$ for `CGTA` changes the odds to $20/81$, giving posterior probability $20/101\approx0.198$. The letters nearly double the island probability, but background remains the more probable class.

## Evidence accumulates along the sequence

A four-base window provides only three pair contributions. To see how evidence builds up over a longer stretch, score the simulated sequence of @fig-cpg-island-track one pair at a time. Write $S_i=S(x_1\ldots x_i)$ for the score of the prefix that ends at position $i$. Because our initial distributions are equal, @eq-markov-ratio-adds gives
$$
S_1=0,
\qquad
S_i=S_{i-1}+\beta(x_{i-1},x_i)\quad(i\geq2).
$$
Each step adds the table entry for the newest adjacent pair, so plotting $S_i$ against $i$ draws a walk whose steps are the pair contributions.

::: {.algorithm}
**Algorithm: accumulate evidence along a sequence.**
Set $S_1$ to the initial-base log ratio. For $i=2,\ldots,n$, add $\beta(x_{i-1},x_i)$ to $S_{i-1}$ to obtain $S_i$. Return $S_1,\ldots,S_n$.
:::

```{.python .execute}
#| label: fig-cpg-cumulative-evidence
#| code-fold: true
#| code-summary: "Accumulate pair contributions along the simulated stretch"
#| fig-cap: "Prefix score S_i along the simulated stretch of the opening figure; the shaded interval is the island. The walk falls through background, climbs inside the island, and falls again, so its slope shows which chain explains the local pairs better."
#| fig-alt: "A jagged line against base position from 1 to 1000. It falls through the first 300 positions, rises inside a shaded band from 301 to 600, then falls again to position 1000."
prefix_scores = [math.log2(initial_plus[track[0]] / initial_minus[track[0]])]
for previous, current in zip(track, track[1:]):
    contribution = math.log2(transition_plus[previous][current]
                             / transition_minus[previous][current])
    prefix_scores.append(prefix_scores[-1] + contribution)

fig, ax = plt.subplots()
ax.axvspan(island_start - 0.5, island_end + 0.5, color=cf.LIGHT, label="Island")
ax.plot(range(1, length + 1), prefix_scores)
ax.axhline(0, color=cf.RULE, linewidth=0.8)
ax.set(xlabel="Base position i", ylabel="Prefix score $S_i$ (bits)")
ax.legend()
plt.show()
for label, first, last in [("background", 1, island_start - 1),
                           ("island", island_start, island_end),
                           ("background", island_end, length)]:
    slope = (prefix_scores[last - 1] - prefix_scores[first - 1]) / (last - first)
    print(f"{label} positions {first}-{last}: {slope:+.3f} bits per base")
```

Why does the walk fall in background and climb in the island? Suppose the previous base is $a$ and the next base is drawn from row $a$ of the background chain. The expected contribution of the new pair is
$$
\sum_{b}a^-_{ab}\,\beta(a,b)
=\sum_{b}a^-_{ab}\log_2\frac{a^+_{ab}}{a^-_{ab}}
\leq 0 .
$$
This is the inequality of @thm-background-log-odds applied to two rows instead of two pair tables: a log ratio averaged under the distribution in its denominator cannot be positive, and it equals zero only when the two rows agree. Exchanging the roles of the chains shows that the average is nonnegative when the next base comes from the island row. For our tables, the background averages are $-0.117$ bit after A, G, or T and $-0.719$ bit after C; the island averages are $+0.117$ and $+0.844$. Most of the evidence therefore arrives right after a C, where the two chains disagree most about what comes next. Weighting these averages by how often each base occurs over a long stretch gives slopes of about $-0.25$ bit per base in background and $+0.34$ in the island, close to the printed slopes of the simulated walk.

::: {.callout-note title="Key idea"}
The log-likelihood ratio of two Markov chains is a walk with one step per adjacent pair. It drifts upward where the island chain explains the pairs better and downward where the background chain does, and a window's score is the walk's rise across that window.
:::

The height $S_i$ scores the whole prefix $x_1\ldots x_i$, which may span several regions, so a positive height is not a label for base $i$; location is carried by the slope. The walk also need not rise at every island base, because CA and CT still favor background there. Only the trend reflects the balance of evidence.

### Window length

A window score is a difference of two heights of the walk. The pairs inside a window $x_i\ldots x_j$ are exactly the steps taken between positions $i$ and $j$, so, with our equal initial distributions,
$$
S(x_i\ldots x_j)=S_j-S_i .
$$
A window is called island when the walk climbs more than $c$ bits across it. A window of $w$ bases inside one region contains $w-1$ steps, so its expected score is about $w-1$ times that region's slope: roughly $+0.34(w-1)$ bits in an island and $-0.25(w-1)$ in background. The expected separation therefore grows in proportion to $w$. Individual steps, however, range from $-1.585$ to $+2.322$ bits, so a short window can land on the wrong side of the threshold by chance.

@fig-cpg-window-calls applies the threshold $c=0$ to every window of the simulated stretch, drawing each call at its window's center. With $w=20$, calls flicker: 47 background windows, each containing at least one CG, are called island, and 18 island windows are called background. With $w=100$, the chance calls disappear and only 16 windows are miscalled, all centered within ten bases of a boundary. These errors come from mixing: every window centered within 50 bases of a boundary contains bases from both regions, so where the called interval begins and ends depends on how the evidence in a mixed window happens to balance.

```{.python .execute}
#| label: fig-cpg-window-calls
#| code-fold: true
#| code-summary: "Call every window of two lengths from differences of the prefix scores"
#| fig-cap: "Island calls at threshold 0 for every window of 20 and of 100 bases along the simulated stretch, each drawn at its window's center, beneath the true island. Short windows produce scattered false calls; long windows remove them but cannot place the island boundaries exactly."
#| fig-alt: "Three horizontal strips over positions 1 to 1000. The top strip, labeled true island, is filled from 301 to 600. The middle strip, labeled 20-base windows, is filled mostly across the island with a few gaps and has several thin slivers outside it. The bottom strip, labeled 100-base windows, is one filled interval close to 301 to 600."
#| fig-width: 6
#| fig-height: 2.2
def window_calls(width: int, threshold: float) -> tuple[list[float], list[bool]]:
    """Return the center and island call of every window of the given width."""
    centers, calls = [], []
    for i in range(1, length - width + 2):
        j = i + width - 1
        window_score = prefix_scores[j - 1] - prefix_scores[i - 1]
        centers.append((i + j) / 2)
        calls.append(window_score > threshold)
    return centers, calls

fig, ax = plt.subplots()
positions = list(range(1, length + 1))
truth = [island_start <= i <= island_end for i in positions]
ax.fill_between(positions, -0.35, 0.35, where=truth, color=cf.INK, step="mid")
labels = ["true island"]
for row, width in enumerate([20, 100], start=1):
    centers, calls = window_calls(width, threshold=0)
    ax.fill_between(centers, -row - 0.35, -row + 0.35, where=calls,
                    color="C0", step="mid")
    labels.append(f"{width}-base windows")
    false_calls, missed = 0, 0
    for center, is_island in zip(centers, calls):
        inside = island_start <= center <= island_end
        false_calls += is_island and not inside
        missed += inside and not is_island
    print(f"w = {width}: {false_calls} background windows called island, "
          f"{missed} island windows called background")
ax.set_yticks([0, -1, -2], labels)
ax.set(xlabel="Base position", xlim=(1, length))
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
plt.show()
```

Window length is therefore a compromise between noise and resolution. A window must be long enough that the expected separation outgrows the pair-to-pair variation; for these synthetic chains, a few dozen bases suffice. It must also be shorter than the islands we hope to find, and near a boundary it mixes the two regimes: a 100-base window straddling a boundary may contain 50 background and 50 island bases, and either whole-window label assigns one class to bases from the other. Real islands range from a few hundred to a few thousand bases, so no single width serves every island. Overlapping windows can also give one base both labels [@compeau2015]. Assigning a label to each base, one that may change inside the sequence, is the problem taken up by hidden Markov models.

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
Of the 100 windows, 88 were classified correctly, giving accuracy 0.88; the specificity of 0.90 describes only the background windows.
:::

Raising the threshold can remove false positives, but it can also turn true positives into false negatives.

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

::: {#exr-cpg-window-boundary}
A 100-base window straddles an island boundary, with a fraction $f$ of its pairs inside the island. Using the long-run slopes $+0.34$ and $-0.25$ bit per base, find the value of $f$ at which the expected window score is zero. Where, relative to the true boundary, will the centers of the first island calls tend to fall?
:::

::: {#exr-cpg-leakage}
A 500-base region is converted to all overlapping 100-base windows and the windows are randomly split. Explain why this is data leakage and propose a valid split.
:::

::: {#refs}
**References**
:::
