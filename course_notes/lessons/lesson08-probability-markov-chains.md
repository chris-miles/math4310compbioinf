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

Lessons 6 and 7 recovered a genome's sequence from its reads. The next question is what that sequence contains, and the first step is to say which patterns of letters chance alone would produce. Sequencing and assembly turn a genome into a long string of letters. Here are the first 60 bases of the bacteriophage phiX174 genome as stored in the [course FASTA file](../data/phix174-NC_001422.1.fasta):

```text
GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAA
```

Adjacent letters are neighboring bases along one strand, so the file records the order of the bases as well as their counts, and that order can matter biologically. In many vertebrate genomic regions, C followed by G occurs less often than base composition alone would suggest. Chemical modification of C and its subsequent conversion to T contribute to that depletion, while some GC-rich regions, called *CpG islands*, retain more CG neighbors. CpG means C followed by G on one strand. These regions often lie near gene starts, so locating them can help identify DNA involved in controlling gene activity [@durbin1998, pp. 48–49].

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

Our input is contiguous assembled DNA, because joining separate reads would invent neighbors that the experiment never observed. We will compare two kinds of model for such a string, independent bases and a first-order Markov chain, and simulate from both to see which can produce a CG deficit.

A genome file contains one observed string, whereas a probability model describes a whole collection of strings that could have been observed and assigns each a probability. That lets us ask whether an observed pattern is ordinary or surprising under specific assumptions. We write a *random sequence* as $X=X_1X_2\ldots X_n$: each $X_i$ is a *random variable* whose value is a letter of $\Sigma=\{A,C,G,T\}$, determined when the model draws a sequence. Lowercase $x=x_1x_2\ldots x_n$ denotes one resulting string, such as the string in the file.

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

Conditional probabilities have a counting interpretation. Generate a very long sequence from a model that uses the same rule at every position and can move from any letter to any other, and look only at the positions whose preceding base is C. Among those positions, the fraction holding G approaches $P(X_i=G\mid X_{i-1}=C)$:
$$
P(X_i=G\mid X_{i-1}=C)\approx
\frac{\#\{i:\ x_{i-1}=C \text{ and } x_i=G\}}{\#\{i:\ x_{i-1}=C\}}.
$$
The numerator counts CG pairs and the denominator counts all pairs that begin with C. Dividing by the denominator restricts attention to positions satisfying the condition, just as dividing by $P(B)$ does in the definition.

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

Independence means that learning the earlier letters does not change the probabilities at the next position:
$$
P(X_i=b\mid X_1=x_1,\ldots,X_{i-1}=x_{i-1})=q_b
$$
for every history with positive probability. In counting terms, the fraction of G among positions that follow a C equals the fraction of G among all positions. Substituting these unchanged probabilities into the chain rule gives the product, and "identically distributed" says that the same four $q_b$ values apply at every position. For two neighbors, independence predicts
$$
P(X_{i-1}=a,\ X_i=b)=q_aq_b,
$$
a prediction about pair frequencies that we can test against sequence later in this lesson.

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

The iid model cannot distinguish strings with the same base counts. In particular, it cannot represent a preference for G specifically after C. A first-order chain remembers the current base when generating the next one.

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

The simplification is in the conditioning. For a three-letter string, the chain rule and the Markov assumption give
$$
P(\texttt{ACG})
=P(X_1=A)\,P(X_2=C\mid X_1=A)\,P(X_3=G\mid X_1=A,X_2=C)
=q_A\,a_{AC}\,a_{CG}.
$$
The last factor forgets the A because only the current base C matters. Every prefix ending in C uses the same C row, so four rows replace a probability table for every possible history. Equal rows recover iid when the initial distribution also equals that common row.

::: {.callout-note title="Key idea"}
The chain rule writes the probability of any string as a product of next-base probabilities. A first-order Markov chain keeps only the previous base in each condition, so a $4\times4$ transition table and the first-base probabilities determine the probability of every DNA string.
:::

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
The chain assigns six times the iid probability because both CG transitions and the GC transition are favored.
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

A *training sequence* is used to choose the parameters; a *held-out sequence* is set aside until those choices are finished, so it can test prediction on data the model did not fit. For a training sequence from the population being modeled, let $n_{ab}$ be the number of times base $a$ is followed by $b$. The natural estimate of $a_{ab}$ is the counting ratio from the start of the lesson: the fraction of departures from $a$ that go to $b$.

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

For a fixed initial distribution and positive transition probabilities, collect equal transition factors in the sequence likelihood:
$$
P(x)=q_{x_1}\prod_{a,b\in\Sigma}a_{ab}^{n_{ab}},\qquad
\log P(x)=\log q_{x_1}+\sum_{a,b}n_{ab}\log a_{ab}.
$$
Two strings with the same first base and transition counts receive the same likelihood. A first-order chain retains more order information than an iid model, but still discards distinctions beyond these counts.

When fitting separate sequence records, keep their boundaries: records `ACG` and `TAC` supply AC twice, CG once, and TA once. Concatenating them invents a GT pair between unrelated endpoints. Each record also needs its own initial-base factor.

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

Scoring and simulation take $\mathcal{O}(n)$ time because each new base requires one transition lookup or draw.

### Does a simulation reproduce neighboring-base structure?

A GC-rich sequence can contain many CG pairs even if its positions are independent. To separate order from composition, let $f(a)$ be the fraction of positions holding base $a$ and $f(ab)$ the fraction of adjacent pairs equal to $ab$. Independence predicts that a pair is $ab$ with probability $q_aq_b$, so the ratio
$$
R(ab)=\frac{f(ab)}{f(a)f(b)}
$$
compares the observed frequency of the pair with the frequency predicted from composition alone. A ratio near one is consistent with independence for this pair; above one indicates enrichment and below one depletion.

To see the ratio detect depletion, start from the iid model of @exm-iid-sequence, with $q_A=q_T=0.3$ and $q_C=q_G=0.2$, and change only what happens after C. The new C row is $(0.3,0.3,0.1,0.3)$ in A, C, G, T order: G after C becomes half as likely, and the removed probability goes to C. Every other row still equals $q$. We simulate 10,000 bases from each model and compute all sixteen ratios.

```{.python .execute}
#| label: fig-probability-neighbor-ratio
#| fig-cap: "Observed-to-expected ratio R(ab) for all sixteen adjacent pairs in 10,000 simulated bases from the iid model and from the chain with a depleted C row. Only the chain dips at CG (red); dividing by composition removes the effect of how common C and G are, so the dip measures order."
#| fig-alt: "Grouped bars for sixteen pairs AA through TT with a dotted reference line at one. Blue iid bars all lie near one. Orange chain bars lie near one except a red CG bar near 0.5 and a CC bar near 1.4."
#| fig-width: 6
#| fig-height: 3
def pair_ratios(sequence: str) -> dict[str, float]:
    """Return f(ab) / (f(a) f(b)) for every ordered pair of DNA bases."""
    n = len(sequence)
    pair_counts = {}
    for previous, current in zip(sequence, sequence[1:]):
        pair = previous + current
        pair_counts[pair] = pair_counts.get(pair, 0) + 1
    ratios = {}
    for a in "ACGT":
        for b in "ACGT":
            observed = pair_counts.get(a + b, 0) / (n - 1)
            expected = (sequence.count(a) / n) * (sequence.count(b) / n)
            ratios[a + b] = observed / expected
    return ratios

iid_row = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
iid_rows = {base: iid_row.copy() for base in "ACGT"}
depleted_rows = {base: iid_row.copy() for base in "ACGT"}
depleted_rows["C"] = {"A": 0.3, "C": 0.3, "G": 0.1, "T": 0.3}

rng = random.Random(1)
iid_ratios = pair_ratios(simulate_markov_chain(10000, iid_row, iid_rows, rng))
chain_ratios = pair_ratios(simulate_markov_chain(10000, iid_row, depleted_rows, rng))

pairs = list(iid_ratios)
chain_colors = [cf.EMPHASIS if pair == "CG" else "C1" for pair in pairs]
fig, ax = plt.subplots()
ax.bar([i - 0.2 for i in range(16)], iid_ratios.values(), 0.4, label="iid model")
ax.bar([i + 0.2 for i in range(16)], chain_ratios.values(), 0.4,
       color=chain_colors, label="depleted C row")
ax.axhline(1, color=cf.RULE, linestyle=":", linewidth=1)
ax.set_xticks(range(16), pairs, fontfamily=cf.MONO)
ax.set(xlabel="Adjacent pair ab", ylabel="Observed / expected, R(ab)")
ax.legend()
plt.show()
print("R(CG): iid", round(iid_ratios["CG"], 2), "chain", round(chain_ratios["CG"], 2))
```

In the iid sequence every ratio stays near one. In the chain, CG drops to about one half, the only pair far below one, and CC rises because it received the probability taken from G. The other ratios move by roughly ten percent, because the change makes G slightly rarer and C slightly more common overall. Since $R$ divides out composition, a GC-rich region with independent positions would still give $R(CG)$ near one; a dip means that G follows C less often than the frequencies of C and G predict. This CG-depleted chain will reappear as the background model when we classify genomic regions in the next lesson.

### Training fit and prediction

Compare both models on the same held-out sequence. Raw probability decreases with sequence length even under a good model, so log-likelihood per base is useful when summarizing across lengths.

A chain has four transition rows to estimate, while iid has one base distribution. Unrestricted maximum likelihood cannot give the chain a worse training fit, because choosing identical rows recovers iid. Better training fit can still mean worse prediction: with little data, the extra rows can learn chance patterns instead of reproducible neighboring-base structure, and only a held-out sequence reveals the difference.

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

Larger means more probability assigned to the observed string. The chain fits the CG-rich training string better but predicts the held-out string less well than iid, because its C row, learned from only two departures from C, gives the held-out CC pair too little probability. That is the cost of estimating sixteen probabilities from five transitions. With genomic data, we need representative training regions and separate evaluation regions to learn whether neighboring-base structure improves prediction.
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

::: {#exr-markov-phix-ratios}
Read the phiX174 genome from [`data/phix174-NC_001422.1.fasta`](../data/phix174-NC_001422.1.fasta) and compute $R(ab)$ for all sixteen pairs. Is CG depleted in this bacteriophage? Which pairs depart most from one, and what do those departures say about an iid model for this genome?
:::

::: {#refs}
**References**
:::
