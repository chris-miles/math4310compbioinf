---
published: true
title: "Probability and Markov Chains for DNA"
subtitle: "Lesson 8 · 2027-02-08 · Week 5"
nocite: |
  @durbin1998
---

## Core question

When do neighboring DNA bases matter for a sequence model?

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

## Sequences as random variables

Sequencing and assembly give DNA strings whose adjacent letters represent neighboring bases along one strand. We can use that order to look for regions with distinctive sequence patterns. In many vertebrate genomic regions, C followed by G occurs less often than base composition alone would suggest. Chemical modification of C and its subsequent conversion to T contribute to that depletion; some regions near gene starts retain more CG neighbors [@durbin1998, pp. 48–49].

Does a region contain many CG neighbors simply because it contains many C and G bases, or does their order supply additional evidence? Consider the toy strings `ACGCG` and `AGCCG`. Each contains one A, two C bases, and two G bases, but the first has two CG pairs and the second has one. A model that uses only individual base frequencies assigns them the same probability. A model that uses the preceding base can distinguish them.

```{.python .execute}
#| label: fig-probability-same-counts-neighbors
#| code-fold: true
#| code-summary: "Mark neighboring CG pairs"
#| fig-cap: "Two strings have identical base counts but different numbers of adjacent CG pairs, marked by underlines. Base composition alone loses the order that a neighboring-base model needs."
#| fig-alt: "ACGCG has underlines below positions two to three and four to five; AGCCG has one underline below positions four to five."
#| fig-width: 6
#| fig-height: 1.8
fig, axes = plt.subplots(1, 2)
cf.sequence_windows(axes[0], "ACGCG", word="CG", label="x")
cf.sequence_windows(axes[1], "AGCCG", word="CG", label="y")
plt.show()
```

Our input is contiguous assembled DNA. Joining separate reads would invent neighbors that the experiment did not observe. We want to assign a probability to an observed string under explicit assumptions, then ask which model reproduces its patterns. Independent-base models supply a starting point; first-order Markov chains retain dependence on the preceding base. Simulation will help us see the difference before the next lesson uses two fitted chains to classify genomic regions.

A genome file contains one observed string. A probability model describes a collection of strings that could have been observed and assigns a probability to each one. The model lets us ask whether an observed pattern is ordinary or surprising under specific assumptions.

::: {#def-random-sequence}
## Random sequence and realization
A *random sequence* $X=X_1X_2\ldots X_n$ records the letters produced by a random model. Each $X_i$ is a random variable: its value is a letter in $\Sigma$, determined when we draw the sequence. A *realization* $x=x_1x_2\ldots x_n$ is one resulting string.
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
An *event* is a condition on the random outcome, such as "the second base is G." For events $A$ and $B$ with $P(B)>0$, the *conditional probability* of $A$ given $B$ is
$$
P(A\mid B)=\frac{P(A\cap B)}{P(B)}.
$$
Equivalently, $P(A\cap B)=P(A\mid B)P(B)$. Repeated application gives the chain rule
$$
P(x_1,\ldots,x_n)=P(x_1)\prod_{i=2}^n P(x_i\mid x_1,\ldots,x_{i-1}).
$$
:::

For example, $P(X_2=G\mid X_1=C)$ is the fraction of model-generated strings beginning with C whose second base is G. The denominator restricts the comparison to strings satisfying the condition.

The chain rule describes a way to assemble the probability of one whole string: first ask for its first letter, then its second letter given the first, then its third given the first two, and continue. This identity holds before making any independence or memory assumption. The difficulty is that the full past can have $4^{i-1}$ possible values at position $i$. A tractable model must specify which parts of that past matter.

## Independent and identically distributed bases

The simplest sequence model uses one probability $q_b$ for each base $b$ and makes the positions independent.

::: {#def-iid-sequence-model}
## iid sequence model
In an *independent and identically distributed* (iid) model, every position has the same base probabilities $q_b=P(X_i=b)$ and the positions are independent. Therefore
$$
P(x)=\prod_{i=1}^{n}q_{x_i}.
$$
:::

Independence means that learning any earlier letters with positive probability does not change the probabilities at the next position: $P(X_i=b\mid X_1=x_1,\ldots,X_{i-1}=x_{i-1})=q_b$. Substituting these unchanged probabilities into the chain rule gives the product. "Identically distributed" says that the same four $q_b$ values apply at every position.

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

The iid model cannot distinguish strings with the same base counts. In particular, it cannot represent a preference for G specifically after C. A first-order chain remembers the current base when generating the next one. This represents neighboring-base dependence without claiming that a genome is literally generated one base at a time.

::: {#def-first-order-markov-chain}
## First-order Markov chain
A *time-homogeneous first-order Markov chain* uses the same transition probabilities at every position and satisfies, for conditioning events of positive probability,
$$
P(X_i=b\mid X_1=x_1,\ldots,X_{i-1}=a)
=P(X_i=b\mid X_{i-1}=a)=a_{ab}.
$$
The transition matrix $A=(a_{ab})$ has nonnegative entries and each row sums to one. With initial probabilities $q_b=P(X_1=b)$,
$$
P(x)=q_{x_1}\prod_{i=2}^{n}a_{x_{i-1}x_i}.
$$
:::

The simplification is in the conditioning: all prefixes ending in C use the same C row, regardless of the letters before C. We replace a probability table for every possible history with four rows. Equal rows recover iid when the initial distribution also equals that common row.

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

:::: {.callout-note collapse="true" title="Going deeper: why normalized counts maximize likelihood"}
::: {.proof}
Sketch. Consider one transition row, with $N=\sum_b n_{ab}>0$ and observed proportions $r_b=n_{ab}/N$. Let $S=\{b:n_{ab}>0\}$ be the bases actually observed after $a$. For a proposed probability row $u$ positive on $S$, the change in log-likelihood relative to $r$ is
$$
N\sum_{b\in S}r_b\ln(u_b/r_b)
\le N\left(\sum_{b\in S}u_b-1\right)\le 0.
$$
The first inequality uses $\ln t\le t-1$ from Lesson 2; the second uses that $u$ has total probability one, possibly with some mass outside $S$. A zero proposed probability for an observed category gives likelihood zero. Thus the observed proportions maximize the row likelihood. A row with no observations is undetermined.
:::
::::

For `ACGCGT`, the transitions are `AC`, `CG`, `GC`, `CG`, `GT`. Both transitions out of $C$ go to $G$, so the raw estimate has $\widehat a_{CG}=1$ and zeros elsewhere in that row. A new `CA` would have probability zero. Adding one pseudocount to every cell gives $\widehat a_{CG}=(2+1)/(2+4)=0.5$ and reserves probability for unobserved transitions.

### Transition counts and sequence boundaries

For a fixed initial distribution, collect equal transition factors in the sequence likelihood:
$$
P(x)=q_{x_1}\prod_{a,b\in\Sigma}a_{ab}^{n_{ab}},\qquad
\log P(x)=\log q_{x_1}+\sum_{a,b}n_{ab}\log a_{ab}.
$$
Here the count form uses positive parameters. With zeros, an unobserved transition contributes no term, while an observed impossible transition makes the log-likelihood $-\infty$. Two strings with the same first base and transition counts receive the same likelihood. A first-order chain retains more order information than an iid model, but still discards distinctions beyond these counts.

When fitting separate sequence records, keep their boundaries: records `ACG` and `TAC` supply AC twice, CG once, and TA once. Concatenating them invents a GT pair between unrelated endpoints. Each record also needs its own initial-base factor. For nonempty records of lengths $L_1,\ldots,L_r$, the transition total is $\sum_j(L_j-1)$.

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
simulated = simulate_markov_chain(10000, initial, transition, rng)
simulated[:30]
```

Scoring and simulation take $\mathcal{O}(n)$ time because each new base requires one transition lookup or draw.

### Does a simulation reproduce neighboring-base structure?

A GC-rich sequence can contain many CG pairs even if its positions are independent. To isolate order, let $f(C)$ and $f(G)$ be observed base fractions and let $f(CG)$ be the fraction of adjacent pairs equal to CG. Under iid, a pair is CG with probability $q_Cq_G$. Thus the ratio
$$
R(CG)=\frac{f(CG)}{f(C)f(G)}
$$
compares observed CG frequency with the frequency predicted from composition alone. A ratio near one is consistent with independence for this pair; above one indicates enrichment and below one depletion. The estimate is undefined if either base is absent.

Our teaching chain favors C-to-G transitions rather than suppressing them. We compare its simulated ratio with that of an iid sequence having the same estimated base composition.

```{.python .execute}
#| label: fig-probability-neighbor-ratio
#| fig-cap: "CG frequency relative to the frequency predicted from base composition in each simulation. The chain preserves an excess of CG neighbors that an iid model with matching composition misses."
#| fig-alt: "Two bars labeled iid and chain show the observed-to-expected CG frequency, with a horizontal reference at one."

base_fractions = {}
for base in initial:
    base_fractions[base] = simulated.count(base) / len(simulated)
iid_rows = {base: base_fractions.copy() for base in initial}
iid_sequence = simulate_markov_chain(10000, base_fractions, iid_rows, rng)
ratios = []
for sequence in [iid_sequence, simulated]:
    cg_count = sum(a == "C" and b == "G"
                   for a, b in zip(sequence, sequence[1:]))
    observed = cg_count / (len(sequence) - 1)
    expected = (sequence.count("C") / len(sequence)
                * sequence.count("G") / len(sequence))
    ratios.append(observed / expected)

fig, ax = plt.subplots()
ax.bar(["iid", "chain"], ratios)
ax.axhline(1, color=cf.RULE, linestyle=":", linewidth=1)
ax.set_ylabel("Observed / expected CG frequency")
plt.show()
print("CG ratios:", [round(ratio, 3) for ratio in ratios])
```

The simulated chain supplies both composition and neighbor preferences. A depleted-CG background model would need a smaller C-to-G probability; the same diagnostic would then ask whether CG occurs less often than its composition predicts.

### Training fit and prediction

Compare both models on the same held-out sequence. Raw probability decreases with sequence length even under a good model, so log-likelihood per base is useful when summarizing across lengths.

A chain has four transition rows to estimate, while iid has one base distribution. Unrestricted maximum likelihood cannot give the chain a worse training fit because identical rows recover iid. Separate smoothing rules need not preserve that guarantee. Better training fit can still mean worse prediction: with little data, extra rows can learn chance patterns instead of reproducible neighboring-base structure. The held-out calculation makes that distinction visible.

::: {#exm-markov-held-out}
## Fit both models, then predict
Use the synthetic training string `ACGCGT` and set aside `AGCCGA` before fitting. Add one pseudocount to each base count for the iid model and to each transition count for the chain. The base counts are $(1,2,2,1)$ in A, C, G, T order, giving $\widehat q=(0.2,0.3,0.3,0.2)$. The smoothed transition rows are:

| current | $A$ | $C$ | $G$ | $T$ |
|---|---:|---:|---:|---:|
| $A$ | $1/5$ | $2/5$ | $1/5$ | $1/5$ |
| $C$ | $1/6$ | $1/6$ | $3/6$ | $1/6$ |
| $G$ | $1/6$ | $2/6$ | $1/6$ | $2/6$ |
| $T$ | $1/4$ | $1/4$ | $1/4$ | $1/4$ |

No departure from T was observed, so its row is determined entirely by the pseudocounts. Use $\widehat q$ as the initial distribution for both models to isolate the effect of transition probabilities. For example, the held-out chain probability is
$$
P_{\mathrm{MC}}(\texttt{AGCCGA})
=\frac15\frac15\frac26\frac16\frac36\frac16
=\frac1{5400}.
$$
An iid model can use the same scoring function by repeating $\widehat q$ in every transition row.

```{.python .execute}
#| label: code-markov-held-out
fitted_initial = {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2}
fitted_transition = {
    "A": {"A": 1/5, "C": 2/5, "G": 1/5, "T": 1/5},
    "C": {"A": 1/6, "C": 1/6, "G": 3/6, "T": 1/6},
    "G": {"A": 1/6, "C": 2/6, "G": 1/6, "T": 2/6},
    "T": {"A": 1/4, "C": 1/4, "G": 1/4, "T": 1/4},
}
iid_transition = {}
for base in fitted_initial:
    iid_transition[base] = fitted_initial.copy()

for label, sequence in [("training", "ACGCGT"), ("held-out", "AGCCGA")]:
    iid_score = markov_log_likelihood(sequence, fitted_initial, iid_transition)
    chain_score = markov_log_likelihood(sequence, fitted_initial, fitted_transition)
    print(label, round(iid_score / len(sequence), 3),
          round(chain_score / len(sequence), 3))
```

The columns report natural-log likelihood per base for iid and chain models:

| sequence | iid | chain |
|---|---:|---:|
| training `ACGCGT` | −1.339 | −1.018 |
| held-out `AGCCGA` | −1.339 | −1.432 |

Larger means more probability assigned to the observed string. The chain favors the CG-rich training pattern but predicts this held-out string less well than iid. The extra transition rows learned the training pattern but gave too little probability to the held-out CC pair. This is the cost of estimating more probabilities from only five transitions. With genomic data, we need representative training regions and separate evaluation regions to learn whether neighboring-base structure improves prediction.
:::

The first-order assumption retains transition counts but loses longer patterns, such as three-base coding units or a motif with separated important bases. The homogeneous assumption also uses one matrix everywhere. CpG islands and surrounding background require different neighboring-base models, which motivates the next lesson.

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
