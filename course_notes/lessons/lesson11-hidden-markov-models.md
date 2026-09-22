---
published: true
title: "Hidden Markov Models"
subtitle: "Lesson 11 · 2027-02-22 · Week 7"
nocite: |
  @durbin1998, @compeau2015
---

## Core question

How can a sequence model represent a region label that we cannot observe directly?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import itertools
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## Observations and hidden states

Sliding a window along a genome gives overlapping decisions, and neighboring windows may disagree about the same bases. Choosing a window size also fixes how sharply a boundary can be located. Instead, we can model the sequence as alternating hidden regions and let the labels persist or change from one position to the next.

The input is the observed base sequence; a possible explanation is a sequence of region labels. Today we define a model that assigns probabilities to those explanations. Scoring one proposed path and summing over tiny examples will separate the questions that the next two lessons solve efficiently.

The classifier in Lesson 10 decides whether an entire window came from an island or background model. Along a chromosome, the region can change. We observe bases, but the region label at each position is unknown.

To keep the arithmetic short, collapse the DNA alphabet into $S$ for $C$ or $G$ and $W$ for $A$ or $T$. This collapse preserves GC enrichment but discards the distinction between CG and GC. The two-state example teaches hidden-region inference; it does not reproduce the dinucleotide classifier from Lesson 10.

::: {#def-hmm-data}
## Observed sequence and hidden path
An observed sequence is $x=x_1\ldots x_n$ over an alphabet $\Sigma$. A hidden path is $\pi=\pi_1\ldots\pi_n$ over a finite state set $Q$. State $\pi_i$ is the unobserved label associated with observation $x_i$.
:::

Our state set is $Q=\{B,I\}$: background and island. The observation $x=\texttt{SSW}$ might have hidden path $\pi=\texttt{IIB}$.

## The B/I hidden Markov model

::: {#def-hidden-markov-model}
## Hidden Markov model
A hidden Markov model (HMM) consists of initial probabilities $p_k=P(\pi_1=k)$, transition probabilities $a_{k\ell}=P(\pi_i=\ell\mid\pi_{i-1}=k)$, and emission probabilities $e_k(b)=P(x_i=b\mid\pi_i=k)$. The hidden states form a first-order Markov chain, and conditional on the entire hidden path, the observations are independent, and the distribution of $x_i$ depends only on $\pi_i$.
:::

We use these illustrative parameters through Lessons 11 to 13:
$$
p_B=p_I=0.5,
\qquad
A=
\begin{array}{c|cc}
 & B&I\\\hline
B&0.8&0.2\\
I&0.4&0.6
\end{array},
\qquad
E=
\begin{array}{c|cc}
 & S&W\\\hline
B&0.3&0.7\\
I&0.8&0.2
\end{array}.
$$

The emission table says that $S$ is more common in islands. The transition table says regions tend to persist, especially background. Every row is a probability distribution.

```{.python .execute}
#| label: fig-hidden-markov-model-diagram
#| fig-cap: "The B/I hidden Markov model used in Lessons 11 to 13. Transitions describe region changes; emission tables describe the observed symbol within each region."
#| fig-alt: "Two circular states labeled B and I with self-loops and arrows in both directions; beneath each state is a table of S and W emission probabilities."
#| fig-width: 5
#| fig-height: 3
transitions = {
    ("B", "B"): 0.8, ("B", "I"): 0.2,
    ("I", "B"): 0.4, ("I", "I"): 0.6,
}
emissions = {
    "B": {"S": 0.3, "W": 0.7},
    "I": {"S": 0.8, "W": 0.2},
}

fig, ax = plt.subplots()
cf.state_diagram(ax, ["B", "I"], transitions, emissions)
plt.show()
```

If a state remains with probability $r$, its expected run length is
$$
1+r+r^2+\cdots=\frac{1}{1-r}.
$$
The illustrative island state therefore has mean length $1/(1-0.6)=2.5$ symbols. Real CpG islands are much longer, so real base-level parameters would use an island self-loop close to one.


### State durations

Let $L$ be the length of a completed run in a state with self-transition probability $r<1$. A run of exactly $\ell$ positions requires $\ell-1$ self-transitions followed by one exit:
$$
P(L=\ell)=r^{\ell-1}(1-r),\qquad \ell=1,2,\ldots.
$$
For the island state, the probabilities of lengths 1, 2, and 3 are $0.4$, $0.24$, and $0.144$. The long tail makes the mean 2.5 even though length 1 is most probable. This is the geometric duration assumption. A run reaching the end of a recorded sequence may continue beyond the observation, so its observed length need not be its completed duration.

## Joint probability of sequence and path

The model assigns a joint probability to the observed symbols and the states that emitted them.

::: {#prp-hmm-joint-factorization}
## HMM joint probability
For any observation $x$ and path $\pi$ of length $n$,
$$
P(x,\pi)
=p_{\pi_1}e_{\pi_1}(x_1)
\prod_{i=2}^{n}
a_{\pi_{i-1}\pi_i}e_{\pi_i}(x_i).
$$
:::

::: {.proof}
Sketch. Apply the probability chain rule in the order $\pi_1,x_1,\pi_2,x_2,\ldots$. The Markov assumption reduces the probability of $\pi_i$ to a term conditioned on $\pi_{i-1}$. The emission assumption reduces the probability of $x_i$ to a term conditioned on $\pi_i$.
:::

::: {#exm-hmm-path-probability}
## Three candidate paths for SSW
For $x=\texttt{SSW}$ and $\pi=\texttt{IIB}$,
$$
P(x,\pi)=0.5(0.8)(0.6)(0.8)(0.4)(0.7)=0.05376.
$$
Two constant-state paths have probabilities
$$
P(\texttt{SSW},\texttt{BBB})
=0.5(0.3)(0.8)(0.3)(0.8)(0.7)=0.02016,
$$
$$
P(\texttt{SSW},\texttt{III})
=0.5(0.8)(0.6)(0.8)(0.6)(0.2)=0.02304.
$$
The mixed path wins among these three because it uses island emissions for the two $S$ symbols and a background emission for $W$, while paying for one switch.
:::

::: {.algorithm}
**Algorithm: score a specified HMM path.**

1. Multiply the initial probability and first emission probability.
2. At each later position, multiply by the transition into its state and the emission from that state.
3. Return the product.
:::

```{.python .execute}
#| label: code-hidden-markov-model-joint
def joint_probability(
    observations: str,
    path: str,
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    emission: dict[str, dict[str, float]],
) -> float:
    """Compute the joint probability of observations and a state path."""
    if len(observations) != len(path):
        raise ValueError("Observations and path must have equal length.")
    if not observations:
        return 1.0
    probability = initial[path[0]] * emission[path[0]][observations[0]]
    for i in range(1, len(observations)):
        probability *= transition[path[i - 1]][path[i]]
        probability *= emission[path[i]][observations[i]]
    return probability


initial = {"B": 0.5, "I": 0.5}
transition = {
    "B": {"B": 0.8, "I": 0.2},
    "I": {"B": 0.4, "I": 0.6},
}
emission = {
    "B": {"S": 0.3, "W": 0.7},
    "I": {"S": 0.8, "W": 0.2},
}

joint_probability("SSW", "IIB", initial, transition, emission)
```

## Summing over hidden paths

A path of length $n$ through $K$ states has $K^n$ possibilities. To obtain the probability of the observations, we sum over them:
$$
P(x)=\sum_{\pi\in Q^n}P(x,\pi).
$$
For $x=\texttt{SSW}$, all eight paths sum to $P(x)=0.1368$. The best path, $\texttt{IIB}$, contains only $0.05376/0.1368\approx39.3\%$ of that probability. A best explanation and the total evidence are different quantities.

::: {.algorithm}
**Algorithm: validate by enumerating paths.**

1. Generate all $K^n$ paths.
2. Compute each path joint probability.
3. Sum the probabilities for $P(x)$ or retain the largest for a best-path check.
:::

```{.python .execute}
#| label: code-hidden-markov-model-enumerate
def enumerate_paths(states: str, n: int) -> list[str]:
    """List state paths for a tiny validation example."""
    return ["".join(path) for path in itertools.product(states, repeat=n)]


paths = enumerate_paths("BI", 3)
probabilities = [
    joint_probability("SSW", path, initial, transition, emission)
    for path in paths
]
sum(probabilities), paths[probabilities.index(max(probabilities))]
```

Enumeration is useful as a test for short sequences, but its $\mathcal{O}(nK^n)$ cost is unusable on genomic data. Lesson 12 finds the most likely path in $\mathcal{O}(nK^2)$ time. Lesson 13 computes the sum and position-wise probabilities at the same cost.


::: {#prp-hmm-normalization}
## Normalization at a fixed length
If the initial probabilities and every transition and emission row sum to one, then
$$
\sum_{x\in\Sigma^n}\sum_{\pi\in Q^n}P(x,\pi)=1.
$$
:::

::: {.proof}
Sketch. Sum over the last observation first; its emission probabilities sum to one for every state. Then sum over the last state; the transition row sums to one. Repeat backward until only the initial-state sum remains, which is also one.
:::

We condition on a fixed sequence length $n$. The example has no end state, so its joint probability includes $n$ emissions and $n-1$ transitions. Scoring one specified path takes $\mathcal O(n)$ time. Extra factors for stopping would define a different model.

## Parameters and estimation

An HMM with $K$ states and $M$ symbols displays $K+K^2+KM$ probability entries. Normalization leaves
$$
(K-1)+K(K-1)+K(M-1)
$$
free parameters. The B/I model has six free parameters. With labeled paths, estimates come from normalized counts of first states, state transitions, and symbols emitted in each state. With only observed sequences, those state counts are unavailable; Lesson 13 sketches Baum-Welch estimation.

::: {#exm-hmm-supervised-counts}
## Estimating emissions from a labeled path
Suppose a short labeled training pair is
$$
x=\texttt{SSWWS},\qquad \pi=\texttt{IIBBB}.
$$
State $I$ emits two $S$ symbols and no $W$ symbols; state $B$ emits one $S$ and two $W$ symbols. Without pseudocounts, the fitted emissions are
$$
\widehat e_I(S)=1,\quad \widehat e_I(W)=0,
\qquad
\widehat e_B(S)=1/3,\quad \widehat e_B(W)=2/3.
$$
The path contains transitions $II,IB,BB,BB$, giving $\widehat a_{II}=\widehat a_{IB}=1/2$ and $\widehat a_{BB}=1$. These exact zeros and ones reflect a tiny sample, not biological certainty. Pseudocounts prevent these small-sample estimates from ruling out unobserved events, as in Lesson 9.
:::

The conditional-independence statement has a precise consequence. If the path is known, learning that $x_{i-1}=S$ does not change the emission distribution at position $i$ beyond what $\pi_i$ already tells us. Dependence can still appear after the path is hidden: two adjacent $S$ symbols are correlated because both are more likely to lie inside the same persistent $I$ run. HMMs create dependence among observations by mixing over correlated hidden states.

## Limitations

1. **Emissions are conditionally independent.** Once the state is known, the model forgets the previous symbol. A two-state, single-base HMM cannot directly express extra probability for $CG$; the state must also carry the previous base or emit pairs.
2. **State durations are geometric.** The highest-probability run length is always one, even when observed regions have a typical nontrivial length.
3. **The state set is chosen by the analyst.** Two states force every position into background or island and omit promoter structure.
4. **Parameters can fit the wrong population.** Estimates from one genome may fail on another. Evaluation should hold out long genomic regions, as in Lesson 10.

## Exercises

::: {#exr-hmm-joint}
Compute $P(\texttt{SWS},\texttt{IBI})$ under the B/I model. Identify every initial, transition, and emission factor.
:::

::: {#exr-hmm-total}
Enumerate the four paths for observation $\texttt{SW}$, compute their joint probabilities, and verify that their sum is a valid probability.
:::

::: {#exr-hmm-duration}
Choose $a_{II}$ so that the mean island length is 200 bases. State the corresponding probability of leaving $I$ at each base.
:::

::: {#exr-hmm-memory}
Describe an HMM state space that remembers both region ($B$ or $I$) and the previous DNA base. How many states does it have, and what dependence can it represent that the two-state model cannot?
:::

::: {#refs}
**References**
:::
