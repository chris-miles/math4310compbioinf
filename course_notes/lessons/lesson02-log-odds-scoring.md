---
published: true
title: "Log-Odds Scoring"
subtitle: "Lesson 2 · 2027-01-13 · Week 1"
nocite: |
  @durbin1998
---

## Core question

How can an alignment score measure evidence for relatedness rather than just count matches?

```{.python .execute}
#| label: setup
#| code-fold: true
import math
```

## Matches and background composition

Two DNA fragments agree at eight of ten positions. That sounds persuasive until we learn that both fragments consist mostly of A. Even unrelated A-rich fragments agree often. To judge a match, we need a comparison: how often would these aligned pairs occur in unrelated sequence?

The data today are two equal-length strings with their positions already paired. We want a numerical score for that fixed, ungapped alignment. Lesson 1 supplied several ways to count differences; here we ask which differences should carry the most evidence. Finding the best placement of gaps comes in Lesson 3.

::: {.biology}
A protein sequence can be obtained by translating an identified DNA coding sequence into amino-acid letters; comparing it with known proteins can suggest a function to investigate. Amino-acid replacements occur at different frequencies in related proteins, motivating a score table that distinguishes kinds of mismatch [@ebi-pairwise].
:::

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

For this biological question, $p_{ab}$ describes pairs in corresponding positions of related sequences. These are joint probabilities, not probabilities of the second letter conditional on the first: the entire table sums to one. Later, a Markov transition table will instead have rows that sum to one.

One way to estimate $p_{ab}$ is to count aligned pairs in representative sequences already known to be related and normalize the counts. Unobserved pairs need a small positive allowance to satisfy our strictly positive model. Background letter frequencies provide $q_a$. In this lesson, we treat the probabilities as calibrated separately and supplied as inputs, rather than estimating them from the candidate pair we want to judge. The choice of related sequences matters: close relatives generally preserve more identities than distant relatives [@durbin1998]. With the probabilities fixed, we can examine the scoring rule without first solving an estimation problem.

A read-mapping model would need to account for sequencing errors as well as sample-to-reference variation; the same mathematical score need not use the same probabilities in both applications.

Consider a synthetic DNA model with $q_a=1/4$ for every base. In the related model, each of the four matching pairs has probability $1/8$, and each of the twelve mismatching pairs has probability $1/24$.

| Pair type | Number of ordered pairs | $p_{ab}$ per pair | $r_{ab}$ per pair |
|---|---:|---:|---:|
| Match | 4 | $1/8$ | $1/16$ |
| Mismatch | 12 | $1/24$ | $1/16$ |

Both tables are normalized. Under $M$, the total probability of a match is $4/8=1/2$. Under $R$, it is $4/16=1/4$. These values are chosen for arithmetic, not fitted to genomic data.

For example, the pair AA has probability $1/8$ under $M$ and $1/16$ under $R$. Its ratio is 2, favoring $M$ by a factor of two. A CG pair instead has ratio $2/3$, favoring $R$. Multiplying ratios accumulates evidence across independent columns. Taking logarithms turns this multiplication into addition, which lets us score one column at a time.

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

A score of 3 bits means the observed alignment is eight times as likely under $M$ as under $R$. It does not give a probability that $M$ is true. For illustration, suppose $M$ and $R$ are the only possibilities, with prior probabilities $1/100$ and $99/100$. Bayes' rule gives
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
One hundred unselected background columns therefore have expected total score about $-18.9$ bits. Some chance alignments still score positively, especially when we search many candidates and retain the best one.

::: {#exm-score-background-change}
## Changing the background
Use the simpler score $+1$ for a match and $-1$ for a mismatch. Uniform independent DNA has match probability $1/4$ and expected score $-1/2$.

Now let $q_A=0.8$ and $q_C=q_G=q_T=1/15$. Two independent letters match with probability
$$
\sum_a q_a^2=0.8^2+3(1/15)^2=\frac{49}{75}.
$$
The same score has expectation $49/75-26/75=23/75>0$. Unrelated A-rich strings accumulate positive scores on average.
:::

This example changes the background while keeping the score fixed. A genuine log-odds table must be recalculated when its background changes.

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

::: {.callout-note collapse="true" title="Going deeper: match fraction and accumulated evidence"}

For the synthetic pair model, let $m$ be the number of matches in $n$ fixed columns. The remaining $n-m$ columns are mismatches, so
$$
S=m+(n-m)\log_2(2/3).
$$
Writing $f=m/n$, a positive score requires
$$
f>\frac{-\log_2(2/3)}{1-\log_2(2/3)}\approx0.369.
$$
Thus a majority of matches is unnecessary for positive evidence under these particular models. The threshold lies between the background match probability $1/4$ and the related-model match probability $1/2$. Changing either model changes this threshold.

For ten fixed columns with four matches, the score is about $0.490$ bits and the likelihood ratio is only $1.405$. For one hundred columns with forty matches, the score is about $4.902$ bits and the ratio is about $29.9$. The same match fraction supplies different amounts of evidence because the second calculation contains ten times as many observations. This interpretation uses the independent-column assumption; copying the same ten observed columns ten times does not create new independent evidence.

We can also compute the expected score under the related model:
$$
\mathbb E_M[s]=\tfrac12(1)+\tfrac12\log_2(2/3)\approx0.208.
$$
Together with the negative background expectation, this explains the direction of the score: under repeated independent sampling, related-model columns accumulate positive evidence on average and background columns accumulate negative evidence. Individual alignments can still favor the wrong model.
:::

## Score scales and gap costs

The total score for a gapped alignment adds substitution scores and gap scores. A linear gap model contributes $g<0$ for each gap character. We have derived the substitution terms; this chosen gap penalty is an additional assumption about insertions and deletions. Lesson 6 will distinguish opening a gap from extending it.

Multiplying every substitution and gap score by the same positive constant preserves the ordering of all alignments. Changing units from bits to natural logarithms does exactly this. Multiplying only the substitution scores changes their weight relative to gaps and can change the chosen alignment.

::: {#exm-score-shift}
## Adding a constant can change the optimum
Align $\texttt{A}$ with $\texttt{A}$, using match $+1$ and gap $-2$. The single-column alignment scores 1; the two-column alignment $\texttt{A-}/\texttt{-A}$ scores $-4$.

Add 6 to every column score, including gaps. The first alignment now scores 7 and the second scores 8. A per-column shift rewards extra columns, whose number can differ between alignments.
:::

Gap-gap columns are never allowed. The two gap columns in this example are valid under Lesson 1's alignment definition, though some later scoring conventions forbid an immediate switch between opposite gap types.

The likelihood-ratio derivation applies to a supplied ungapped alignment under our two column models. Adding chosen gap costs gives an optimization objective, but does not by itself supply a normalized probability model for gapped strings. Lesson 3 finds the alignment maximizing that objective. Interpreting a score after such a search still requires accounting for how the alignment was selected.

## Exercises

::: {#exr-lesson03-expected-score}
Assume each ordered DNA pair is equally likely under the background. Compute the expected score for match $+1$ and mismatch $-1$.
:::

::: {#exr-log-odds-hand}
Using the synthetic pair model, score $\texttt{ACGA}/\texttt{ATGT}$ in bits and check the likelihood ratio by multiplication.
:::

::: {#exr-lesson03-score-assumption}
Give one biological reason a protein scoring matrix might assign different scores to two kinds of mismatch.
:::

::: {#exr-log-odds-units}
Prove that multiplying all column scores by a positive constant preserves every maximizing alignment. Explain why adding a constant need not.
:::

::: {#refs}
**References**
:::
