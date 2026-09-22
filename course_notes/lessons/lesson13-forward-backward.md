---
published: true
title: "Uncertainty in HMMs: Forward-Backward"
subtitle: "Lesson 13 · 2027-03-03 · Week 8"
nocite: |
  @durbin1998
---

## Core question

Given an observed sequence, how much support does each position have for each hidden label?

```{.python .execute}
#| label: setup
#| code-fold: true
import itertools
import math
```

## Genome labels and uncertainty

A genome annotation marks intervals as belonging to one kind of region or another. Near a boundary, several placements can explain the same bases almost equally well. A single decoded path gives coordinates but hides that uncertainty. We want a probability for each label at each position, using the entire observed sequence.

Keep the B/I hidden Markov model from Lessons 11 and 12. The observations are $S$ for C/G and $W$ for A/T; the hidden labels are background $B$ and island-like $I$. This reduced alphabet describes GC enrichment, not the CG-pair signal itself. All calculations condition on the model and a fixed sequence length, with no end-state probability.

$$
p_B=p_I=0.5,\qquad
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

::: {#def-posterior-state}
## Posterior state probability
For observations $x=x_1\ldots x_n$ with $P(x)>0$, define
$$
\gamma_k(i)=P(\pi_i=k\mid x).
$$
The output is one probability distribution over states for every position.
:::

For $\texttt{SSWS}$, Viterbi returns $\texttt{IIBB}$ with joint probability $0.0129024$. That establishes the best *complete path*. To assess the final label, we must add the probabilities of all paths ending in $I$ and compare that sum with all paths ending in $B$.

## The forward algorithm

::: {#def-forward-variable}
## Forward probability
The forward probability
$$
f_k(i)=P(x_1\ldots x_i,\pi_i=k)
$$
sums over every hidden prefix ending at state $k$.
:::

It includes the emission of $x_i$. The recurrence has the Viterbi table's shape, with a sum replacing its maximum:
$$
f_k(1)=p_ke_k(x_1),\qquad
f_k(i)=e_k(x_i)\sum_{\ell\in Q}f_\ell(i-1)a_{\ell k}.
$$
Finally, $P(x)=\sum_k f_k(n)$, because every full path ends in exactly one state.

::: {#exm-forward-ssws}
## Summing paths for SSWS
For the second symbol,
$$
f_B(2)=0.3[0.15(0.8)+0.4(0.4)]=0.084.
$$
Viterbi kept only the larger predecessor contribution and obtained $0.048$. Forward retains both contributions.

| State | $S$, $i=1$ | $S$, $i=2$ | $W$, $i=3$ | $S$, $i=4$ |
|---|---:|---:|---:|---:|
| $B$ | 0.15 | 0.084 | 0.10752 | 0.0293184 |
| $I$ | 0.40 | 0.216 | 0.02928 | 0.0312576 |

Thus $P(\texttt{SSWS})=0.060576$. The best path accounts for $0.0129024/0.060576\approx21.3\%$ of the conditional path probability.
:::

Small total sequence probability is normal: a probability model spreads its mass across many possible strings. The useful comparison here is between path probabilities for the *same* observed string.

## The backward algorithm

Forward probabilities use the observations up through position $i$. The observations to its right also inform the label, because neighboring hidden states tend to persist.

::: {#def-backward-variable}
## Backward probability
The backward probability
$$
b_k(i)=P(x_{i+1}\ldots x_n\mid\pi_i=k)
$$
is the probability of the remaining observations, given the state at $i$.
:::

This quantity excludes $x_i$, which the forward factor already contains. At $i=n$ there are no remaining observations, so $b_k(n)=1$. For earlier positions,
$$
b_k(i)=\sum_{\ell\in Q}a_{k\ell}e_\ell(x_{i+1})b_\ell(i+1).
$$
The next emission belongs to the *destination* state $\ell$, and the loop runs from right to left.

::: {#exm-backward-ssws}
## Using the suffix
At position 3, only the final $S$ remains:
$$
b_B(3)=0.8(0.3)+0.2(0.8)=0.40,\qquad
b_I(3)=0.4(0.3)+0.6(0.8)=0.60.
$$

| State | $i=1$ | $i=2$ | $i=3$ | $i=4$ |
|---|---:|---:|---:|---:|
| $B$ | 0.08896 | 0.248 | 0.40 | 1 |
| $I$ | 0.11808 | 0.184 | 0.60 | 1 |

The suffix is more likely starting from $I$ at position 3, although the symbol at position 3 is $W$. That symbol is absent from the backward probability; its evidence enters through the forward factor.
:::

## Combining the two directions

::: {#prp-forward-backward-posterior}
## Forward-backward identity
For every position $i$,
$$
P(x,\pi_i=k)=f_k(i)b_k(i),\qquad
P(x)=\sum_k f_k(i)b_k(i).
$$
Consequently, $\gamma_k(i)=f_k(i)b_k(i)/P(x)$.
:::

::: {.proof}
Sketch. Split the observations just after $i$. Conditional on $\pi_i$, the future is independent of the observed prefix. The probability of the prefix and state is $f_k(i)$; the conditional probability of the suffix is $b_k(i)$. Multiplying gives their joint probability. Summing over the mutually exclusive states gives $P(x)$, and division gives the conditional probability.
:::

::: {#exm-posterior-ssws}
## Two answers at the last position
At position 2, $\gamma_I(2)=0.216(0.184)/0.060576\approx0.656$.

| Position | Symbol | $P(I\mid x)$ | Most probable label |
|---|---|---:|---|
| 1 | S | 0.780 | I |
| 2 | S | 0.656 | I |
| 3 | W | 0.290 | B |
| 4 | S | 0.516 | I |

Choosing the largest posterior at each position gives $\texttt{IIBI}$. Viterbi gave $\texttt{IIBB}$. At the last position, all paths ending in $I$ together beat all paths ending in $B$, even though the strongest individual path ends in $B$.
:::

A 0.516 posterior is weak support for that last island label. Reporting the probability makes this visible. It is also conditional on our chosen emissions and transition probabilities; it does not measure uncertainty about whether those assumptions are right.

::: {#def-posterior-decoding}
## Posterior decoding
Posterior decoding chooses $\widehat\pi_i=\arg\max_k\gamma_k(i)$ separately at each position.
:::

For any proposed labels $z_1,\ldots,z_n$, the expected number correct, conditional on $x$, is $\sum_i\gamma_{z_i}(i)$. Choosing a largest probability in each column maximizes that sum. Viterbi instead maximizes the probability that the entire path is correct. Posterior decoding can violate path constraints when transitions are forbidden.

### Expected island positions and boundary support

Posterior probabilities can answer quantitative questions without committing to one path. Let $N_I$ be the number of positions labeled I in the hidden path. Write it as a sum of indicators, one for each position. Linearity of conditional expectation gives
$$
\mathbb E[N_I\mid x]=\sum_{i=1}^n\gamma_I(i).
$$
For SSWS, the displayed posteriors give approximately $0.780+0.656+0.290+0.516=2.242$ island positions. A noninteger expectation is natural: it averages the integer counts across possible paths. Viterbi reports two island positions and posterior decoding reports three; neither count is required to equal this average.

Uncertainty about a boundary involves two neighboring labels. For example, the probability of an I-to-B transition between positions 2 and 3 is
$$
P(\pi_2=I,\pi_3=B\mid x)
=\frac{f_I(2)a_{IB}e_B(W)b_B(3)}{P(x)}
=\frac{0.216(0.4)(0.7)(0.40)}{0.060576}
\approx0.399.
$$
The factorization joins a prefix ending in I, the specified transition and next emission, and the remaining suffix. It counts every complete path having that boundary. Multiplying the separate marginal probabilities $\gamma_I(2)\gamma_B(3)$ would give about 0.466 and would incorrectly assume those labels are independent conditional on the observations.

An annotation can therefore report a likely region while assigning only moderate support to its exact endpoint. State posteriors describe positions; transition posteriors describe boundaries. Both sum over alternative paths rather than treating one decoded path as known.

## Implementation and checks

::: {.algorithm}
**Algorithm: forward-backward inference.**

1. Fill the forward table left to right and sum its final column.
2. Fill the backward table right to left, starting with ones.
3. At every position, multiply corresponding forward and backward entries and divide by the sequence probability.
4. Check that posterior columns sum to one and that every forward-backward cut gives the same sequence probability.
:::

The code uses ordinary probabilities so each line matches the hand tables. It is for short validation examples. Long sequences require the log-space modification below.

```{.python .execute}
#| label: code-forward-backward
def forward_backward(
    observations: str, states: list[str],
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    emission: dict[str, dict[str, float]],
) -> tuple[list[dict[str, float]], list[dict[str, float]], float]:
    """Return forward and backward tables and the short-sequence likelihood."""
    if not observations:
        return [], [], 1.0
    forward = [{k: initial[k] * emission[k][observations[0]] for k in states}]
    for symbol in observations[1:]:
        column = {}
        for current in states:
            total = 0.0
            for previous in states:
                total += forward[-1][previous] * transition[previous][current]
            column[current] = emission[current][symbol] * total
        forward.append(column)

    backward = [{k: 1.0 for k in states} for _ in observations]
    for i in range(len(observations) - 2, -1, -1):
        for current in states:
            total = 0.0
            for following in states:
                total += (transition[current][following]
                          * emission[following][observations[i + 1]]
                          * backward[i + 1][following])
            backward[i][current] = total
    probability = sum(forward[-1].values())
    if probability == 0:
        raise ValueError("Impossible observations or numerical underflow.")
    return forward, backward, probability


states = ["B", "I"]
initial = {"B": 0.5, "I": 0.5}
transition = {"B": {"B": 0.8, "I": 0.2}, "I": {"B": 0.4, "I": 0.6}}
emission = {"B": {"S": 0.3, "W": 0.7}, "I": {"S": 0.8, "W": 0.2}}
forward, backward, probability = forward_backward(
    "SSWS", states, initial, transition, emission
)
posterior = []
for left, right in zip(forward, backward):
    cut_probability = sum(left[k] * right[k] for k in states)
    assert math.isclose(cut_probability, probability)
    column = {k: left[k] * right[k] / probability for k in states}
    assert math.isclose(sum(column.values()), 1.0)
    posterior.append(column)
probability, [round(column["I"], 3) for column in posterior]
```

For an independent check, enumerate all 16 paths and use the joint factorization directly.

```{.python .execute}
#| label: code-forward-enumeration-check
total = 0.0
island_mass = [0.0] * 4
for path in itertools.product(states, repeat=4):
    mass = initial[path[0]] * emission[path[0]]["S"]
    for i, symbol in enumerate("SWS", start=1):
        mass *= transition[path[i - 1]][path[i]] * emission[path[i]][symbol]
    total += mass
    for i, state in enumerate(path):
        if state == "I":
            island_mass[i] += mass

assert math.isclose(total, probability)
for i in range(4):
    assert math.isclose(island_mass[i] / total, posterior[i]["I"])
print("Enumeration agrees with the likelihood and all four posteriors.")
```

Each table has $nK$ entries and each entry sums $K$ terms. Time is $\mathcal O(nK^2)$ and storage is $\mathcal O(nK)$. A likelihood alone needs only a rolling forward column and $\mathcal O(K)$ storage.

### Stable sums in log space

Let $F_k(i)=\ln f_k(i)$. Forward now requires
$$
F_k(i)=\ln e_k(x_i)+
\operatorname{LSE}_{\ell}\bigl(F_\ell(i-1)+\ln a_{\ell k}\bigr).
$$
For finite maximum $m=\max_j z_j$,
$$
\operatorname{LSE}(z)=\ln\sum_j e^{z_j}
=m+\ln\sum_j e^{z_j-m}.
$$
Subtracting $m$ makes the largest exponent zero. If all entries are $-\infty$, the result is $-\infty$. The same identity stabilizes backward sums. Log posteriors are $F_k(i)+\ln b_k(i)-\ln P(x)$.

For $z=(-1000,-1000)$, direct exponentiation underflows; the stable answer is $-1000+\ln2$. Adding the two logs would give $-2000$, the log of a product, which answers the wrong question.

:::: {.callout-note collapse="true" title="Going deeper: learning from uncertain labels"}
With labeled paths, transition estimates use observed counts. Baum-Welch replaces those counts with expected counts under the current model. For example,
$$
\xi_{k\ell}(i)=
\frac{f_k(i)a_{k\ell}e_\ell(x_{i+1})b_\ell(i+1)}{P(x)}
$$
is the posterior probability of a particular transition. Sum over $i$ to estimate how often that transition was used, then normalize each source-state row. Emission counts similarly weight each observation by $\gamma_k(i)$.

Repeat inference and re-estimation. Exact unregularized updates do not decrease training likelihood, but can settle at a local optimum; the optimization proof is beyond this lesson. Better training likelihood alone does not establish better genome annotation. Initial parameters, state definitions, and held-out evaluation still matter.
::::

## Limitations

1. **Posterior labels need not form a legal path.** For a constrained annotation, choose a decoding objective that also respects allowed transitions.
2. **Probability-space code underflows.** Use log-sum-exp or a consistent scaling scheme for long sequences; a zero computed likelihood need not mean a mathematically impossible sequence.
3. **Confidence is model-dependent.** An inadequate state space can give confident wrong labels.
4. **Training and evaluation are separate.** Re-estimating parameters on the evaluation sequence changes what a held-out check measures.

## Exercises

::: {#exr-forward-hand}
Compute forward and backward tables for $\texttt{SW}$. Check $P(x)$ at both cuts.
:::

::: {#exr-posterior-filter}
For $\texttt{SSWS}$, compute $P(\pi_2=I\mid x_1x_2)$ by normalizing the forward column. Compare it with 0.656 and explain what information changes the answer.
:::

::: {#exr-forward-bug}
A backward implementation uses $e_k(x_i)$ in place of $e_\ell(x_{i+1})$. Identify which observation it counts twice when combined with the forward table.
:::

::: {#exr-forward-logsum}
Evaluate $\operatorname{LSE}(-1000,-1001)$ without directly exponentiating either input. Explain why replacing the sum by a maximum changes the inference.
:::

::: {#refs}
**References**
:::
