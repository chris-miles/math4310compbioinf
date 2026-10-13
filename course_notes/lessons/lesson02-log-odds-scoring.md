---
published: true
title: "Log-Odds Scoring"
subtitle: "Lesson 2 · 2027-01-13 · Week 1"
nocite: |
  @durbin1998
---

## Core question

How much evidence does a match provide when matching bases can occur by chance?

```{.python .execute}
#| label: setup
#| code-fold: true
import math
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import numpy as np
import coursefigs as cf

cf.use_style()
```

## Matches can occur by chance

Suppose we have DNA sequences from two organisms and want to know whether their agreement supports descent from a common ancestral sequence. Sequencing and assembly supply the base strings. For now, take two equal-length segments whose positions have already been paired, with no gaps. A match count summarizes that pairing, but it does not say how much agreement unrelated sequences would produce.

For example, suppose ten paired positions contain eight A/A matches and two mismatches. If both sequences come from A-rich DNA, many A/A pairs can occur by chance. @fig-log-odds-chance-pair shows this for one column: each axis represents an independent base draw, and the red area contains pairs where both draws give A. With A frequency $1/4$ in each draw, that area is $(1/4)^2=1/16$; with frequency $0.8$, it is $0.8^2=0.64$. The same matching letters therefore need different scores in those two settings.

```{.python .execute}
#| label: fig-log-odds-chance-pair
#| code-fold: true
#| code-summary: "Show how independent base draws create chance matches"
#| fig-cap: "Each square represents all outcomes of two independent base draws; red marks A/A outcomes. Raising the A frequency from one quarter to 0.8 enlarges the chance-match area from one sixteenth to 0.64, so an A/A match carries less evidence against an A-rich background."
#| fig-alt: "Two unit squares divided into A and other outcomes along each axis. The red lower-left A/A area occupies one sixteenth of the left square and 64 percent of the right square."
#| fig-width: 6
#| fig-height: 3
fig, axes = plt.subplots(1, 2)
cf.independent_pair_area(axes[0], 0.25)
cf.independent_pair_area(axes[1], 0.8)
axes[0].set_title("A frequency 1/4: chance A/A = 1/16")
axes[1].set_title("A frequency 0.8: chance A/A = 0.64")
plt.show()
```

We will assign probabilities to each possible paired column under two explanations: shared ancestry and independent background sequence. Comparing those probabilities gives a score that measures which explanation better predicts the observed pairs. The input is a fixed alignment and supplied model probabilities; the output is its evidence score. Searching for the best gap placement is a separate computation in Lesson 3.

::: {#def-aligned-pair-data}
## Aligned pair data
For strings $x=x_1\ldots x_n$ and $y=y_1\ldots y_n$, the data are the ordered pairs $(x_i,y_i)$ for $i=1,\ldots,n$. A pair-count table $N_{ab}$ records how often pair $(a,b)$ occurs.
:::

For $\texttt{ACGT}$ aligned to $\texttt{AGGT}$, the pairs are $AA,CG,GG,TT$. The count table discards their order. That will be sufficient under the independence assumption below.

## Two probability models

To measure evidence, we need to predict the observed letter pairs under two possible origins: related sequences whose corresponding positions descend from an ancestor, and unrelated sequences drawn from a background population. Each model assigns a probability to every possible pair. We compare those probabilities for the pairs actually observed.

::: {#def-pair-models}
## Background and related-pair models
The background model $R$ draws the two letters independently with base probabilities $q_a>0$, where $\sum_aq_a=1$, so $r_{ab}=q_aq_b$. The related-pair model $M$ assigns each ordered pair a probability $p_{ab}>0$, with $\sum_{a,b}p_{ab}=1$.

Both models assume that distinct alignment columns are independent and use the same pair distribution at every column.
:::

A probability specifies a sampling frequency. If we repeatedly generate aligned pairs under $M$, about a fraction $p_{ab}$ are the pair $(a,b)$. It describes both letters together, so the entire table sums to one. A probability conditional on the first letter being A would instead describe frequencies only among pairs whose first letter is A. Those conditional probabilities would sum to one across that row.

One way to estimate $p_{ab}$ is to count aligned pairs in representative sequences already known to be related and normalize the counts. Unobserved pairs need a small positive allowance to satisfy our strictly positive model. Background letter frequencies provide $q_a$. In this lesson, we treat the probabilities as calibrated separately and supplied as inputs, rather than estimating them from the candidate pair we want to judge. The choice of related sequences matters: close relatives generally preserve more identities than distant relatives [@durbin1998]. With the probabilities fixed, we can examine the scoring rule without first solving an estimation problem.

A read-mapping model would need to account for sequencing errors as well as sample-to-reference variation; the same mathematical score need not use the same probabilities in both applications.

Consider a synthetic DNA model with $q_a=1/4$ for every base. In the related model, each of the four matching pairs has probability $1/8$, and each of the twelve mismatching pairs has probability $1/24$.

| Pair type | Number of ordered pairs | $p_{ab}$ per pair | $r_{ab}$ per pair |
|---|---:|---:|---:|
| Match | 4 | $1/8$ | $1/16$ |
| Mismatch | 12 | $1/24$ | $1/16$ |

Both tables are normalized. Under $M$, the total probability of a match is $4/8=1/2$. Under $R$, it is $4/16=1/4$. These values are chosen for arithmetic, not fitted to genomic data.

For example, AA has probability $1/8$ under $M$ and $1/16$ under $R$: we expect twice as many AA pairs when sampling from $M$. The ratio 2 expresses that relative evidence. Equal probability differences can represent different relative evidence: increasing from 0.01 to 0.02 doubles the probability, whereas increasing from 0.50 to 0.51 gives a ratio of only 1.02. A CG pair has ratio $2/3$, favoring $R$. Independent columns contribute multiplying ratios; logarithms turn the product into a sum that we can compute one column at a time.

::: {#def-pair-log-odds}
## Log-odds score
The score of aligned pair $(a,b)$, in bits, is
$$
s(a,b)=\log_2\frac{p_{ab}}{q_aq_b}.
$$
The total score of the fixed alignment is $S(x,y)=\sum_{i=1}^n s(x_i,y_i)$.
:::

A positive entry means that this pair is more likely under $M$ than under $R$. It need not be a match: a sufficiently common substitution in related proteins can also receive a positive score.

## Adding evidence across columns

::: {#prp-pair-additivity}
## Additive scores are log likelihood ratios
Under @def-pair-models,
$$
S(x,y)=\log_2\frac{P(x,y\mid M)}{P(x,y\mid R)}.
$$
:::

::: {.proof}
Independence across columns gives $P(x,y\mid M)=\prod_{i=1}^n p_{x_i y_i}$ and $P(x,y\mid R)=\prod_{i=1}^n q_{x_i}q_{y_i}$. Therefore
$$
\log_2\frac{P(x,y\mid M)}{P(x,y\mid R)}
=\log_2\prod_{i=1}^n\frac{p_{x_i y_i}}{q_{x_i}q_{y_i}}
=\sum_{i=1}^n\log_2\frac{p_{x_i y_i}}{q_{x_i}q_{y_i}}
=S(x,y).
$$
The second equality uses the logarithm-of-a-product rule. Each observed pair contributes exactly its column score.
:::

A score of 3 bits means the observed alignment is eight times as likely under $M$ as under $R$. It does not give a probability that $M$ is true. For illustration, suppose $M$ and $R$ are the only possibilities. Before seeing these letters, we assign $M$ probability $1/100$ and $R$ probability $99/100$. These *prior probabilities* express how often we expect each origin among candidate pairs before using their alignment evidence. Bayes' rule combines those starting odds with the likelihood ratio:
$$
\frac{P(M\mid x,y)}{P(R\mid x,y)}
=\frac{P(x,y\mid M)P(M)}{P(x,y\mid R)P(R)}
=8\times\frac{1}{99}=\frac{8}{99}.
$$
The common normalizing probability of the data cancels in this ratio. The posterior probability of $M$ is then $8/(8+99)\approx0.075$. Positive evidence has increased its probability, but a rare starting hypothesis can remain unlikely.

::: {#exm-pair-log-odds}
## Three matches and one mismatch
For the synthetic table,
$$
s_{\mathrm{match}}=\log_2 2=1,\qquad
s_{\mathrm{mismatch}}=\log_2(2/3)\approx-0.585.
$$
Thus $\texttt{ACGT}/\texttt{AGGT}$ scores $3-\!0.585\approx2.415$ bits. Directly,
$$
\frac{(1/8)^3(1/24)}{(1/16)^4}=\frac{16}{3}\approx5.333=2^{S}.
$$
:::

The product calculation and score calculation are independent ways to check the same result. A score of $+1$ for a match and $-1$ for a mismatch is a convenient teaching convention, but it is not the exact log-odds table derived from these probabilities.

Because the model treats columns identically and independently, the count table from @def-aligned-pair-data also suffices: $S=\sum_{a,b}N_{ab}s(a,b)$. Rearranging the columns leaves this score unchanged. A model that uses neighboring bases would need their order as well.

::: {.algorithm}
**Algorithm: score a fixed ungapped alignment.**

1. Check that the strings have equal length.
2. For each aligned pair, divide its related-model probability by its background probability.
3. Add the base-2 logarithms of these ratios.
:::

```{.python .execute}
#| label: code-log-odds-fixed-alignment
def alignment_log_odds(
    x: str, y: str, pair: dict[tuple[str, str], float],
    background: dict[str, float],
) -> float:
    """Score an ungapped alignment using positive model probabilities."""
    if len(x) != len(y):
        raise ValueError("An ungapped alignment needs equal-length strings.")
    score = 0.0
    for a, b in zip(x, y):
        score += math.log2(pair[a, b] / (background[a] * background[b]))
    return score


background = {base: 0.25 for base in "ACGT"}
pair = {}
for a in "ACGT":
    for b in "ACGT":
        pair[a, b] = 1 / 8 if a == b else 1 / 24

score = alignment_log_odds("ACGT", "AGGT", pair, background)
assert math.isclose(sum(pair.values()), 1.0)
assert math.isclose(2 ** score, 16 / 3)
score, 2 ** score
```

There are $n$ pair lookups and additions. With a fixed DNA alphabet, scoring takes $\mathcal O(n)$ time and constant extra space. It evaluates one supplied alignment; searching all possible alignments is a different task.

## Expected scores and chance matches

Our four-column example favors $M$, but a chance background alignment could produce the same letters. To understand what scores are typical before selecting any candidate, we calculate the average score under $R$. A sensible background-calibrated score should not reward unrelated sequence merely for being long.

::: {#def-background-expected-score}
## Expected background score
For pair probabilities $r_{ab}$ and scores $s(a,b)$, the expected score of one background column is
$$
\mathbb E_R[s]=\sum_{a,b}r_{ab}s(a,b).
$$
:::

Expectation is a probability-weighted average. It predicts the average over repeated draws, not the score of every individual draw. For our table,
$$
\mathbb E_R[s]=\tfrac14(1)+\tfrac34\log_2(2/3)\approx-0.189.
$$
One hundred unselected background columns therefore have expected total score about $-18.9$ bits. Averaging the same scores under the uniform related model gives $\mathbb E_M[s]=\tfrac12(1)+\tfrac12\log_2(2/3)\approx0.208$: the two origins predict opposite average directions. Individual chance alignments can still score positively, especially when we search many candidates and retain the best one.

::: {#exm-score-background-change}
## Changing the background
Use the simpler score $+1$ for a match and $-1$ for a mismatch. Uniform independent DNA has match probability $1/4$ and expected score $-1/2$.

Now let $q_A=0.8$ and $q_C=q_G=q_T=1/15$. Two independent letters match with probability
$$
\sum_a q_a^2=0.8^2+3(1/15)^2=\frac{49}{75}.
$$
The same score has expectation $49/75-26/75=23/75>0$. Unrelated A-rich strings accumulate positive scores on average.
:::

We can repair the score by building a related model for this same A-rich background. Suppose $M$ generates half its pairs by two independent draws from $q$, and half by drawing one letter from $q$ and copying it. Then
$$
p_{ab}=\tfrac12q_aq_b+\tfrac12q_a[a=b],
$$
where $[a=b]$ equals 1 for a match and 0 otherwise. This model preserves the background frequencies in each sequence but favors agreement between them. Its log-odds scores are
$$
s(A,A)=\log_2(9/8)\approx0.170,\qquad
s(C,C)=s(G,G)=s(T,T)=3,\qquad
s(a,b)=-1\quad(a\ne b).
$$
An A/A match now supplies little evidence because it is already common by chance. A match of a rare base supplies much more. The opening pair with eight A/A matches and two mismatches scores $8\log_2(9/8)-2\approx-0.641$ bits, favoring the background model even though the simple $+1/-1$ rule awards 6 points.

Averaging this recalibrated table under $R$ gives
$$
\mathbb E_R[s]=\frac{16}{25}\log_2(9/8)+\frac{1}{75}(3)-\frac{26}{75}
\approx-0.198.
$$
Unselected A-rich background columns now decrease the score on average. The sign is a consequence of using two normalized probability models, as the next proposition shows.

::: {#prp-background-log-odds}
## Background expectation of a log-odds score
For strictly positive normalized pair distributions $p$ and $r$,
$$
\sum_{a,b}r_{ab}\log_2(p_{ab}/r_{ab})\leq0,
$$
with equality only when $p=r$.
:::

::: {.proof}
Sketch. For $u>0$, the derivative of $u-1-\ln u$ is $1-1/u$, negative below 1 and positive above 1. Its minimum is therefore zero at $u=1$, proving $\ln u\leq u-1$. Apply this inequality to $u=p_{ab}/r_{ab}$, multiply by $r_{ab}$, and sum. The upper bound is $\sum p_{ab}-\sum r_{ab}=0$. Dividing by $\ln2$ converts to bits. Equality requires $p_{ab}/r_{ab}=1$ for every pair, since every weight $r_{ab}$ is positive.
:::

Thus a properly calibrated log-odds table has negative background expectation unless its two models coincide. This is an average statement for unselected columns. It does not prevent the best candidate in a large search from receiving a positive score.

The expectation describes repeated sampling. In @fig-log-odds-composition-walk, the same simulated pairs receive two scores, so the change in behavior comes entirely from calibration.

::: {.algorithm}
**Algorithm: compare score accumulation under the two models.**
Draw independent ordered pairs from each supplied model. Score each stream first by match $+1$/mismatch $-1$, then by its calibrated log-odds table. Plot the running sums against the number of observed pairs.
:::

```{.python .execute}
#| label: fig-log-odds-composition-walk
#| code-fold: true
#| code-summary: "Simulate and compare accumulated scores"
#| fig-cap: "Simulated pair streams under the A-rich background R and related model M, scored in two ways. Ordinary match rewards rise even for chance pairs; calibrated log-odds scores accumulate evidence in opposite directions."
#| fig-alt: "Two line plots against number of paired columns. Both cumulative point scores rise on the left; the background log-odds falls and the related-model log-odds rises on the right."
#| fig-width: 6
#| fig-height: 3.5
np.random.seed(2)
rich_background = {"A": 0.8, "C": 1 / 15, "G": 1 / 15, "T": 1 / 15}
words, random_probabilities, related_probabilities = [], [], []
for a in "ACGT":
    for b in "ACGT":
        random_probability = rich_background[a] * rich_background[b]
        copy_probability = rich_background[a] if a == b else 0
        words.append((a, b))
        random_probabilities.append(random_probability)
        related_probabilities.append(0.5 * random_probability + 0.5 * copy_probability)

fig, axes = plt.subplots(1, 2)
for label, probabilities, line_style in [
    ("Background R", random_probabilities, "-"),
    ("Related M", related_probabilities, "--"),
]:
    draws = np.random.choice(len(words), size=300, p=probabilities)
    raw_scores, calibrated_scores = [], []
    for word_index in draws:
        a, b = words[word_index]
        raw_scores.append(1 if a == b else -1)
        calibrated_scores.append(math.log2(
            related_probabilities[word_index] / random_probabilities[word_index]))
    columns = np.arange(1, len(draws) + 1)
    axes[0].plot(columns, np.cumsum(raw_scores), line_style, label=label)
    axes[1].plot(columns, np.cumsum(calibrated_scores), line_style, label=label)
for ax in axes:
    ax.set_xlabel("Paired columns")
    ax.axhline(0, color=cf.RULE, linewidth=0.8)
axes[0].set_ylabel("Cumulative points")
axes[1].set_ylabel("Cumulative log-odds (bits)")
axes[1].legend()
plt.show()
```

The derivation scores a supplied ungapped alignment. Adding gap scores creates an alignment objective, but a normalized model for gapped strings needs further assumptions. Searching for the highest-scoring alignment also changes the distribution of the reported score: it selects the best among many candidates.

## Exercises

::: {#exr-log-odds-expected-score}
Assume each ordered DNA pair is equally likely under the background. Compute the expected score for match $+1$ and mismatch $-1$.
:::

::: {#exr-log-odds-hand}
Using the uniform-background pair model, score $\texttt{ACGA}/\texttt{ATGT}$ in bits and check the likelihood ratio by multiplication.
:::

::: {#exr-log-odds-score-assumption}
Give one biological reason a protein scoring matrix might assign different scores to two kinds of mismatch.
:::

::: {#exr-log-odds-match-fraction}
For the uniform-background pair model, write the score of $n$ paired columns with $m$ matches. Find the smallest match fraction $m/n$ that gives positive evidence, and explain why the same fraction supplies more evidence as $n$ grows.
:::

::: {#exr-log-odds-units}
Prove that multiplying all column scores by a positive constant preserves every maximizing alignment. To see why adding a constant need not, score A/A and A-/-A using match $+1$ and gap $-2$, then add 6 to every column score, including gaps.
:::

::: {#refs}
**References**
:::
