---
published: true
title: "Probability and Markov Chains for DNA"
subtitle: "Lesson 9 · 2027-02-10 · Week 5"
nocite: |
  @durbin1998
---

## Core question

What is the probability of a sequence under a stated model?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports"
import math
import random
```

## Sequences as random variables

A stretch of DNA contains several CG pairs. Is that unusual, or just what we should expect from a region rich in C and G? Sequence alone does not define surprise: we need a background model that says which strings are common and which are rare.

The input today is one DNA string, together with probabilities estimated from representative sequence. We will compute its probability under two models, first treating bases independently and then allowing the previous base to affect the next. Simulation gives another way to examine what each model assumes. These calculations prepare a classifier for island-like regions in Lesson 10.

A genome file contains one observed string. A probability model describes a collection of strings that could have been observed and assigns a probability to each one. The model lets us ask whether an observed pattern is ordinary or surprising under specific assumptions.

::: {#def-random-sequence}
## Random sequence and realization
A *random sequence* $X=X_1X_2\ldots X_n$ is a sequence of random variables taking values in an alphabet $\Sigma$. A *realization* $x=x_1x_2\ldots x_n$ is one possible observed value of $X$.
:::

For DNA, $\Sigma=\{A,C,G,T\}$. Uppercase $X$ denotes the random object and lowercase $x$ the string in the file.

::: {#def-probability-model}
## Probability model
A probability model on length-$n$ strings assigns each $x\in\Sigma^n$ a number $P(X=x)$ such that
$$
P(X=x)\geq 0
\qquad\text{and}\qquad
\sum_{x\in\Sigma^n}P(X=x)=1.
$$
:::

We abbreviate $P(X=x)$ as $P(x)$. The probability refers to a model, not uncertainty inside the letters already read by the sequencer.

::: {#def-conditional-probability}
## Conditional probability
For events $A$ and $B$ with $P(B)>0$, the *conditional probability* of $A$ given $B$ is
$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)}.
$$
Equivalently, $P(A\cap B)=P(A\mid B)P(B)$. Repeated application gives the chain rule
$$
P(x_1,\ldots,x_n)=P(x_1)\prod_{i=2}^n P(x_i\mid x_1,\ldots,x_{i-1}).
$$
:::

Conditional probability restricts attention to outcomes consistent with what follows the conditioning bar. Sequence models become tractable by simplifying the conditional probabilities in the chain rule.

## Independent and identically distributed bases

The simplest sequence model uses one probability $q_b$ for each base $b$ and makes the positions independent.

::: {#def-iid-sequence-model}
## iid sequence model
In an *independent and identically distributed* (iid) model, every position has the same base probabilities $q_b=P(X_i=b)$ and the positions are independent. Therefore
$$
P(x)=\prod_{i=1}^{n}q_{x_i}.
$$
:::

The word independent gives the product. The words identically distributed say that the same four $q_b$ values apply at every position.

::: {#exm-iid-sequence}
## Likelihood under an iid model
Consider the synthetic model $q_A=q_T=0.3$ and $q_C=q_G=0.2$. For $x=\texttt{ACGCG}$,
$$
P(x)=0.3(0.2)(0.2)(0.2)(0.2)=0.00048.
$$
Changing the order to `AGCCG` leaves the probability unchanged because the model sees only base counts.
:::

::: {#def-likelihood}
## Likelihood
For fixed observed data $x$, the *likelihood* of model parameters $\theta$ is $L(\theta;x)=P_\theta(x)$.
:::

Probability and likelihood use the same number with different quantities held fixed. A probability model fixes $\theta$ and compares possible strings. A likelihood calculation fixes $x$ and compares parameter values or models. Likelihoods need not sum to one across models.

## First-order Markov chains

The iid model cannot distinguish strings with the same base counts. A first-order chain remembers the current base when generating the next one.

::: {#def-first-order-markov-chain}
## First-order Markov chain
A *time-homogeneous first-order Markov chain* uses the same transition probabilities at every position and satisfies
$$
P(X_i=b\mid X_1=x_1,\ldots,X_{i-1}=a)
=P(X_i=b\mid X_{i-1}=a)=a_{ab}.
$$
The transition matrix $A=(a_{ab})$ has nonnegative entries and each row sums to one. With initial probabilities $q_b=P(X_1=b)$,
$$
P(x)=q_{x_1}\prod_{i=2}^{n}a_{x_{i-1}x_i}.
$$
:::

"First-order" means one position of memory. It does not assert strong dependence. Equal rows recover the iid model when the initial distribution also equals that common row.

::: {#exm-markov-sequence}
## One sequence under two models
Use the synthetic transition matrix below. Rows give the current base and columns the next base.

| current | $A$ | $C$ | $G$ | $T$ |
|---|---:|---:|---:|---:|
| $A$ | 0.3 | 0.2 | 0.3 | 0.2 |
| $C$ | 0.2 | 0.3 | 0.4 | 0.1 |
| $G$ | 0.2 | 0.3 | 0.3 | 0.2 |
| $T$ | 0.2 | 0.2 | 0.3 | 0.3 |

For $x=\texttt{ACGCG}$,
$$
P_{\mathrm{MC}}(x)=0.3(0.2)(0.4)(0.3)(0.4)=0.00288.
$$
The chain assigns six times the iid probability because both CG transitions and the GC transition are favored. These values illustrate the calculation; they are not estimates from a genome.
:::

## Log-likelihoods

Sequence probabilities become too small for ordinary floating-point arithmetic. Even $0.25^{1000}\approx 10^{-602}$ is below the smallest positive number represented by a standard Python `float`.

::: {#def-log-likelihood}
## Log-likelihood
The *log-likelihood* is $\ell(\theta;x)=\log P_\theta(x)$. For a Markov chain,
$$
\log P(x)=\log q_{x_1}+\sum_{i=2}^{n}\log a_{x_{i-1}x_i}.
$$
:::

The logarithm converts products into sums and preserves order. Base 2 reports evidence in bits; the natural logarithm is common in software. For the example, $\log_2(0.00288/0.00048)=\log_2 6\approx2.585$ bits.

::: {.algorithm}
**Algorithm: score a sequence under a Markov chain.**

1. Start with $\log q_{x_1}$.
2. For each adjacent pair $x_{i-1}x_i$, add $\log a_{x_{i-1}x_i}$.
3. Return the sum.
:::

```{.python .execute}
#| label: code-probability-markov-score
def markov_log_likelihood(
    sequence: str,
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
) -> float:
    """Return the natural-log likelihood under a Markov chain."""
    if not sequence:
        return 0.0
    if initial[sequence[0]] == 0:
        return -math.inf
    score = math.log(initial[sequence[0]])
    for previous, current in zip(sequence, sequence[1:]):
        probability = transition[previous][current]
        if probability == 0:
            return -math.inf
        score += math.log(probability)
    return score


initial = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
transition = {
    "A": {"A": 0.3, "C": 0.2, "G": 0.3, "T": 0.2},
    "C": {"A": 0.2, "C": 0.3, "G": 0.4, "T": 0.1},
    "G": {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2},
    "T": {"A": 0.2, "C": 0.2, "G": 0.3, "T": 0.3},
}

markov_log_likelihood("ACGCG", initial, transition)
```

## Estimating a chain

A training sequence is a sequence used to choose the parameters. A held-out sequence is set aside until those choices are finished, so it can test prediction on data the model did not fit. Our first estimate is simply the fraction of departures from each base that go to each possible next base.

For a training sequence from the population being modeled, count adjacent pairs. Let $n_{ab}$ be the number of times base $a$ is followed by $b$.

::: {#prp-markov-mle}
## Maximum-likelihood transition estimates
For every row with at least one observed transition, the maximum-likelihood estimate is
$$
\widehat a_{ab}=\frac{n_{ab}}{\sum_{c\in\Sigma}n_{ac}}.
$$
:::

::: {.proof}
Sketch. Consider one transition row, and let $r_b=n_{ab}/N$ be its observed proportions, where $N=\sum_b n_{ab}>0$. For a proposed row $u$, moving probability to an unobserved category cannot improve the fit. On observed categories, the change in log-likelihood relative to $r$ is $N\sum_b r_b\ln(u_b/r_b)$. The inequality $\ln t\le t-1$ from Lesson 2 bounds this above by $N(\sum_b u_b-1)=0$. A zero proposed probability for an observed category gives likelihood zero. Thus the observed proportions maximize the row likelihood. A row with no observations is undetermined.
:::

For `ACGCGT`, the transitions are `AC`, `CG`, `GC`, `CG`, `GT`. Both transitions out of $C$ go to $G$, so the raw estimate has $\widehat a_{CG}=1$ and zeros elsewhere in that row. A new `CA` would have probability zero. Adding one pseudocount to every cell gives $\widehat a_{CG}=(2+1)/(2+4)=0.5$ and reserves probability for unobserved transitions.

## Simulation

::: {.algorithm}
**Algorithm: simulate a Markov chain.**

1. Draw $X_1$ from the initial distribution $q$.
2. For $i=2,\ldots,n$, draw $X_i$ from row $X_{i-1}$ of $A$.
3. Return $X_1\ldots X_n$.
:::

```{.python .execute}
#| label: code-probability-markov-simulate
def simulate_markov_chain(
    n: int,
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    rng: random.Random,
) -> str:
    """Simulate a length-n sequence from a first-order chain."""
    if n < 0:
        raise ValueError("Length must be nonnegative.")
    if n == 0:
        return ""
    alphabet = list(initial)
    sequence = rng.choices(alphabet, weights=[initial[base] for base in alphabet], k=1)
    while len(sequence) < n:
        row = transition[sequence[-1]]
        sequence.extend(rng.choices(alphabet, weights=[row[base] for base in alphabet], k=1))
    return "".join(sequence)


rng = random.Random(0)
simulate_markov_chain(30, initial, transition, rng)
```

Scoring and simulation take $\mathcal{O}(n)$ time. A DNA chain stores four initial probabilities and 16 transitions, reduced by one normalization constraint for each distribution.

### Training fit and prediction

A likelihood comparison is meaningful only when both models score the same observed object. A longer sequence usually has a smaller raw probability than a shorter one, even when both fit equally well per base. Compare two models on the same held-out strings. When summarizing predictive performance across strings of different lengths, an average log-likelihood per base can be useful, but it is not a posterior model probability. Modeling variable sequence length requires an additional length distribution.

The number of fitted parameters also matters. A first-order DNA chain fits four conditional distributions, whereas an iid model fits one base distribution. On the training sequence, the more flexible model cannot fit worse: it includes the iid model as the special case with identical rows. That does not prove it predicts new sequence better. Fit parameters on one region and compare log-likelihoods on a separate region. A chain that learns a short training string perfectly can assign zero or poor probability to new transitions.

::: {#exm-markov-held-out}
## Training fit and prediction
The raw estimate from $\texttt{ACGCGT}$ gives $\widehat a_{CG}=1$. It fits both observed $CG$ transitions perfectly. It also assigns probability zero to a held-out $CA$ transition. With one pseudocount per cell, $\widehat a_{CG}=0.5$ and $\widehat a_{CA}=1/6$. The smoothed model gives up some training likelihood to make a finite prediction for an event absent from six training bases.
:::

## Limitations

1. **One base of memory misses longer structure.** Codon periodicity depends on position modulo three, and regulatory elements can join bases separated by many positions.
2. **The same chain applies everywhere.** One matrix cannot represent alternating regions such as CpG islands and background. Lessons 10 and 11 introduce region-specific models.
3. **Zero counts create impossible events.** Pseudocounts encode how much probability to reserve for transitions absent from finite training data.
4. **Likelihood does not include class prevalence.** A model can fit a string well and still describe a rare class. Prior probabilities enter posterior comparisons through Bayes' rule.

## Exercises

::: {#exr-markov-hand-score}
Compute the iid and Markov probabilities of `CGCG` using @exm-markov-sequence. Give the likelihood ratio and its base-2 logarithm.
:::

::: {#exr-markov-estimation}
Estimate a transition matrix from `ACACGTC` using one pseudocount per cell. Show the four entries in the $C$ row and verify that they sum to one.
:::

::: {#exr-iid-as-markov}
Prove that an iid model is a Markov chain whose rows are identical. Then show that a Markov chain with identical rows gives the iid sequence probability when its initial distribution equals that common row.
:::

::: {#exr-markov-validation}
State two tests for `markov_log_likelihood`: one using a length-one sequence and one using a chain with identical rows. Give the expected result of each.
:::

::: {#refs}
**References**
:::
