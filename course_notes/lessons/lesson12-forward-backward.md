---
published: true
title: "Uncertainty in HMMs: Forward-Backward"
subtitle: "Lesson 12 · 2027-03-01 · Week 8"
nocite: |
  @durbin1998
---

## Core question

How much support does each base have for belonging to an island-like region?

```{.python .execute}
#| label: setup
#| code-fold: true
import math
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, "../style")
import coursefigs as cf
cf.use_style()
```

## Genome labels and uncertainty

Lesson 11 chose one complete annotation of the assembled DNA: for `CGAC`, the model's best path put the first two bases in an island-like region and the last two in background. That gives a proposed boundary after base 2. How strongly do the letters support it? Other paths can explain the same sequence, including a path that puts all four bases in the island. Selecting the best path alone does not tell us how much support these alternatives share.

The input is still the DNA sequence and the fixed region model. We now want a probability for each label at each position, using all the observed bases. Such probabilities can identify an uncertain boundary or a weakly supported interval before we treat the decoded map as evidence about a genomic region [@durbin1998].

Keep the B/I hidden Markov model from Lessons 10 and 11. The observations are $S$ for C/G and $W$ for A/T; the hidden labels are background $B$ and island-like $I$. This reduced alphabet describes GC enrichment, not the CG-pair signal itself. All calculations condition on the model and a fixed sequence length, with no end-state probability.

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

$Q=\{B,I\}$ is the state set. Write $p_k$ for the initial-state probability, $a_{k\ell}$ for the probability of moving from state $k$ to state $\ell$, and $e_k(b)$ for the probability that state $k$ emits symbol $b$. Rows of $A$ are source states; rows of $E$ are emitting states. The hidden states form a first-order Markov chain. Given a hidden path, emissions are independent and each emission depends only on its state. These are the assumptions that let us combine evidence on either side of a position.

::: {#def-posterior-state}
## Posterior state probability
For observations $x=x_1\ldots x_n$ with $P(x)>0$, define
$$
\gamma_k(i)=P(\pi_i=k\mid x).
$$
The output is one probability distribution over states for every position.
:::

The DNA `CGAC` gives $x=\texttt{SSWS}$. Viterbi returns $\texttt{IIBB}$ with joint probability $0.0129024$, establishing the best *complete path*. To assess the final label, we must add the probabilities of all paths ending in $I$ and compare that sum with all paths ending in $B$. This changes the operation we need: retain the combined probability of alternative explanations rather than just the largest one.

## The forward algorithm

::: {#def-forward-variable}
## Forward probability
The forward probability
$$
f_k(i)=P(x_1\ldots x_i,\pi_i=k)
$$
sums over every hidden prefix ending at state $k$.
:::

It includes the emission of $x_i$. To compute the entry for state $k$, partition the possible hidden prefixes by their state $\ell$ at position $i-1$. Those alternatives are mutually exclusive. Their masses therefore add, even though many different paths can belong to each alternative.

A prefix ending in $\ell$ has total mass $f_\ell(i-1)$. The transition to $k$ multiplies it by $a_{\ell k}$, and the destination state emits $x_i$ with probability $e_k(x_i)$. That emission factor is common to every predecessor. Thus
$$
f_k(1)=p_ke_k(x_1),\qquad
f_k(i)=e_k(x_i)\sum_{\ell\in Q}f_\ell(i-1)a_{\ell k}.
$$
For $i=1$, there is no predecessor: the initial state emits the first symbol. Finally, $P(x)=\sum_k f_k(n)$, because every full path ends in exactly one state. Viterbi keeps a single best path in each entry, while forward retains the combined mass of all prefixes; changing the operation changes what the entry means [@durbin1998].

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

Could we get a state probability simply by normalizing a forward column? Yes, but it answers a question about only the prefix:
$$
\frac{f_I(2)}{f_B(2)+f_I(2)}
=\frac{0.216}{0.084+0.216}=0.72
=P(\pi_2=I\mid SS).
$$
This is the support for I before seeing the remaining WS. Those later bases also inform position 2 because region labels tend to persist. To use them, we need the probability of that suffix under each possible state at the current position.

::: {#def-backward-variable}
## Backward probability
The backward probability
$$
b_k(i)=P(x_{i+1}\ldots x_n\mid\pi_i=k)
$$
is the probability of the remaining observations, given the state at $i$.
:::

This quantity excludes $x_i$, which the forward factor already contains. At $i=n$ there are no remaining observations, so $b_k(n)=1$ under our fixed-length model.

For an earlier position, condition on the next hidden state $\pi_{i+1}=\ell$. Its probability given the current state is $a_{k\ell}$. Once $\ell$ is fixed, the next emission has probability $e_\ell(x_{i+1})$, and the later observations have probability $b_\ell(i+1)$. The model makes that emission independent of the later observations conditional on the next state. It also makes that remaining suffix independent of the earlier state $k$ once $\ell$ is known. Therefore
$$
b_k(i)=\sum_{\ell\in Q}a_{k\ell}e_\ell(x_{i+1})b_\ell(i+1).
$$
The terms sum over mutually exclusive choices of the next state. The next emission belongs to the *destination* state $\ell$; using $e_k(x_i)$ would include the current observation again. The loop runs from right to left so each required suffix probability has already been computed [@durbin1998].

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

At position 2, the suffix contains both W and S. Starting from B, it can first stay in B or switch to I:
$$
b_B(2)=\underbrace{0.8(0.7)(0.40)}_{\text{next state B}}+
\underbrace{0.2(0.2)(0.60)}_{\text{next state I}}=0.248.
$$
Each product has a transition, the W emission from the next state, and that state's probability of the final S. The factors $0.40$ and $0.60$ already sum over both choices for the final hidden state; the recurrence does not select a single suffix path.

The suffix is more likely starting from $I$ at position 3, although the symbol at position 3 is $W$. That symbol is absent from the backward probability; its evidence enters through the forward factor.
:::

## Combining the two directions

Fixing the state at position $i$ separates the calculation into a past and a future. Every prefix ending in $k$ can be joined to every suffix starting from $k$. The forward value adds the prefix probabilities; the backward value adds the conditional suffix probabilities. Multiplying the two sums accounts for every compatible combination once. Because forward includes $x_i$ and backward starts at $x_{i+1}$, no observation is counted twice.

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

@fig-forward-backward-path-groups shows the comparison responsible for the disagreement. The largest individual path belongs to the B-ending group, but the I-ending paths have slightly more probability when added together. Summing several explanations answers a different question from selecting one explanation.

```{.python .execute}
#| label: fig-forward-backward-path-groups
#| echo: false
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "Conditional probability mass of SSWS paths grouped by their final label. The red segment is the best complete path IIBB, while all I-ending paths together have slightly greater mass. A best path ending in B can coexist with a weak posterior preference for I at the final base."
#| fig-alt: "Two horizontal probability bars: B-ending paths total 0.484, including a red IIBB segment of 0.213; I-ending paths total 0.516. A vertical reference at 0.5 falls between the bar endpoints."
sequence_probability = 0.060576
best_path_mass = 0.0129024 / sequence_probability
end_B_mass = 0.0293184 / sequence_probability
end_I_mass = 0.0312576 / sequence_probability
fig, ax = plt.subplots()
ax.barh(1, best_path_mass, color=cf.EMPHASIS)
ax.barh(1, end_B_mass - best_path_mass, left=best_path_mass, color=cf.LIGHT)
ax.barh(0, end_I_mass, color="C0")
ax.text(best_path_mass / 2, 1, "IIBB", ha="center", va="center", color=cf.PAPER)
ax.text((best_path_mass + end_B_mass) / 2, 1, "other B paths",
        ha="center", va="center", fontsize=9)
ax.text(end_I_mass / 2, 0, "all I paths", ha="center", va="center", color=cf.PAPER)
for y, mass in [(1, end_B_mass), (0, end_I_mass)]:
    ax.text(mass + 0.025, y, f"{mass:.3f}", va="center", fontsize=9)
ax.axvline(0.5, color=cf.RULE, linewidth=0.8, linestyle=":")
ax.set(xlabel="Probability conditional on SSWS", xlim=(0, 1),
       yticks=[0, 1], yticklabels=["Ends in I", "Ends in B"])
plt.show()
```

::: {#def-posterior-decoding}
## Posterior decoding
Posterior decoding chooses $\widehat\pi_i=\arg\max_k\gamma_k(i)$ separately at each position.
:::

For any proposed labels $z_1,\ldots,z_n$, the expected number correct, conditional on $x$, is $\sum_i\gamma_{z_i}(i)$. Without constraints tying the reported labels together, choosing a largest probability in each column maximizes that expected number of correct labels. Viterbi instead maximizes the probability that the entire path is correct. Posterior decoding can violate path constraints when transitions are forbidden.

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

## Computing and comparing the probabilities

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

### A disagreement inside a region

Append one S to the running observation, giving SSWSS. The most probable complete path is now IIIII, with joint probability 0.005308416. The posterior at the interior W is only $\gamma_I(3)\approx0.323$, so positionwise decoding gives IIBII. A single path can win by staying in I throughout while the combined mass of paths using B at the middle position is larger. This is why the two decoding objectives can give different annotations even away from sequence boundaries.

```{.python .execute}
#| label: fig-forward-backward-evidence
#| fig-cap: "Island-state probabilities for SSWSS using the observed prefix or the entire sequence. Later S observations raise support for I at the interior W, but that support remains below one half."
#| fig-alt: "Two curves over S S W S S compare normalized forward values with full posteriors. Both dip at the central W; the full posterior there is about 0.323, above the forward-only value about 0.214."
observations = "SSWSS"
left, right, likelihood = forward_backward(
    observations, states, initial, transition, emission
)
filtered = []
smoothed = []
for prefix, suffix in zip(left, right):
    filtered.append(prefix["I"] / sum(prefix.values()))
    smoothed.append(prefix["I"] * suffix["I"] / likelihood)
positions = list(range(1, len(observations) + 1))
fig, ax = plt.subplots()
ax.plot(positions, filtered, "o--", label="Observed prefix only")
ax.plot(positions, smoothed, "s-", label="Entire sequence")
ax.axhline(0.5, color=cf.RULE, linewidth=0.8)
ax.set_xticks(positions, [f"{i}\n{symbol}" for i, symbol in zip(positions, observations)])
ax.set(xlabel="Position and observed symbol", ylabel="Probability of island state I",
       ylim=(0, 1))
ax.legend()
plt.show()
```

At the last position the two curves meet: the prefix is the whole sequence, and every backward value is one. At earlier positions their separation measures how the remaining observations change the state probability under this model.

With $K=|Q|$ states, each table has $nK$ entries and each entry sums $K$ terms. Time is $\mathcal O(nK^2)$ and storage is $\mathcal O(nK)$. A likelihood alone needs only a rolling forward column and $\mathcal O(K)$ storage.

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
Factor $e^m$ out of the sum and then take its logarithm to obtain the second expression. Subtracting $m$ makes the largest exponent zero. If all entries are $-\infty$, the result is $-\infty$. The same identity stabilizes backward sums. Log posteriors are $F_k(i)+\ln b_k(i)-\ln P(x)$.

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
