---
published: true
title: "Hidden Markov Models"
subtitle: "Lesson 10 · 2027-02-17 · Week 6"
nocite: |
  @durbin1998, @compeau2015
---

## Core question

Where do island-like regions begin and end along the DNA?

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

The window classifier from Lesson 9 tells us whether a chosen stretch of DNA looks island-like. To make a map along a chromosome, we also need to locate where each region begins and ends. Sliding a window along the assembled DNA gives overlapping decisions: a window crossing a boundary mixes island and background bases, and neighboring windows may disagree about the same positions. Choosing a window size also limits how sharply a boundary can be located [@durbin1998].

The input is the same kind of data, a DNA letter string with genomic coordinates. The desired output is a region label at each position; a consecutive run of island labels becomes a predicted interval. A single C or G could belong to either region, so a letter cannot reveal its label. Neighboring letters can supply evidence because a region extends across several bases. We model that persistence by allowing the unobserved label to stay the same or change between positions.

To keep the arithmetic short, collapse the DNA alphabet into $S$ for $C$ or $G$ and $W$ for $A$ or $T$. This collapse preserves GC enrichment but discards the distinction between CG and GC. The two-state example teaches hidden-region inference; it does not reproduce the dinucleotide classifier from Lesson 9.

::: {#def-hmm-data}
## Observed sequence and hidden path
An observed sequence is $x=x_1\ldots x_n$ over an alphabet $\Sigma$. A hidden path is $\pi=\pi_1\ldots\pi_n$ over a finite state set $Q$. State $\pi_i$ is the unobserved label associated with observation $x_i$.
:::

Our state set is $Q=\{B,I\}$: background and island-like sequence. For example, DNA `CGA` becomes the observation $x=\texttt{SSW}$. A proposed path $\pi=\texttt{IIB}$ assigns its first two bases to an island and its third to background; $\texttt{BBB}$ proposes no island. Both paths can produce the same letters. A *hidden Markov model* assigns probabilities to these competing explanations, using how often each region emits S or W and how often the region changes.

@fig-hidden-markov-model-observed-hidden keeps the recorded bases separate from the proposed region labels. Converting C/G to S is a known recoding step. Choosing between I and B is an inference: the same observed row is compatible with both proposed maps.

```{.python .execute}
#| label: fig-hidden-markov-model-observed-hidden
#| echo: false
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "Recorded DNA, its S/W recoding, and two proposed region maps at the same three positions. The red I cells mark a candidate island; those labels are an explanation of the letters, not an additional sequencing measurement."
#| fig-alt: "Four aligned rows: recorded DNA C G A, recoded observations S S W, candidate map I I B, and alternative map B B B. The first two I cells in the candidate map are red."
fig, ax = plt.subplots()
cf.sequence_tracks(
    ax,
    {"Recorded DNA": "CGA", "Recoded letters": "SSW",
     "Candidate map": "IIB", "Alternative map": "BBB"},
    highlighted_row="Candidate map", highlighted_positions=[0, 1],
)
plt.show()
```

## The B/I hidden Markov model

::: {#def-hidden-markov-model}
## Hidden Markov model
A hidden Markov model (HMM) consists of initial probabilities $p_k=P(\pi_1=k)$, transition probabilities $a_{k\ell}=P(\pi_i=\ell\mid\pi_{i-1}=k)$, and emission probabilities $e_k(b)=P(x_i=b\mid\pi_i=k)$. The hidden states form a first-order Markov chain, and conditional on the entire hidden path, the observations are independent, and the distribution of $x_i$ depends only on $\pi_i$.
:::

We use these illustrative parameters through Lessons 10 to 12:
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

The emission table says that $S$ is more common in islands. The transition table says regions tend to persist, especially background. Every row is a probability distribution. To generate a sequence, first draw a state from $p$ and emit a symbol from its emission row. Then draw the next state from the current state's transition row and emit again. Each position therefore involves two draws: a hidden region and an observed symbol.

```{.python .execute}
#| label: fig-hidden-markov-model-diagram
#| code-fold: true
#| code-summary: "Draw the hidden-region model"
#| fig-cap: "The B/I hidden Markov model used in Lessons 10 to 12. Transitions describe region changes; emission tables describe the observed symbol within each region."
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

Independent emissions may seem to contradict the neighboring patterns we want to model. The condition matters. If the path is known, learning that $x_{i-1}=S$ does not change the emission distribution at position $i$ beyond what $\pi_i$ already tells us. Dependence can still appear after the path is hidden: two adjacent $S$ symbols are correlated because both are more likely to lie inside the same persistent $I$ run.

### Why observations can be dependent

We can quantify the dependence between observations using the same B/I model. The first observation is S with probability
$$
P(X_1=S)=0.5(0.3)+0.5(0.8)=0.55.
$$
Before observing anything, the second hidden state is B with probability $0.5(0.8)+0.5(0.4)=0.6$, and I with probability 0.4. Thus $P(X_2=S)=0.6(0.3)+0.4(0.8)=0.50$.

Now sum the joint probability of SS over its four possible hidden paths:
$$
\begin{aligned}
P(SS)&=0.5(0.3)[0.8(0.3)+0.2(0.8)]\\
&\quad+0.5(0.8)[0.4(0.3)+0.6(0.8)]\\
&=0.06+0.24=0.30.
\end{aligned}
$$
Independence would instead give $0.55(0.50)=0.275$. Observing the first S raises the probability of a second S to $0.30/0.55\approx0.545$. The first observation supplies evidence for I; persistence of that hidden state carries some of the evidence to the next position.

Conditional independence concerns observations after the hidden path is fixed. Hiding a persistent path mixes different emission distributions, which can create dependence among the observed letters.

### How long does a region persist?

Let $L$ be the length of a completed run in a state with self-transition probability $r<1$. A run of exactly $\ell$ positions requires $\ell-1$ self-transitions followed by an exit:
$$
P(L=\ell)=r^{\ell-1}(1-r),\qquad \ell=1,2,\ldots.
$$
For the island state, run lengths 1, 2, and 3 have probabilities $0.4$, $0.24$, and $0.144$. Length one is most probable, but longer runs contribute to the average. The run reaches position $j$ with probability $P(L\ge j)=r^{j-1}$. Counting one unit for every position reached gives
$$
\mathbb E[L]=\sum_{j\ge1}P(L\ge j)
=1+r+r^2+\cdots=\frac{1}{1-r}.
$$
Thus the toy model has mean background length 5 and mean island length 2.5 symbols. These short runs keep hand calculations manageable; parameters for a base-level model of long regions would need self-transitions closer to one.

This is the *geometric duration* assumption. The exit probability stays $1-r$ even after a long run, so the model has no preferred minimum biological length. Increasing $r$ lengthens the mean but preserves that constant exit chance. A run cut off by the end of our recorded sequence may continue beyond it and is not a completed duration.

## Joint probability of sequence and path

We can now translate the two draws at each position into the probability of one proposed explanation. We condition on a fixed observed length $n$: the model has no extra factor for stopping. A length-$n$ path contains $n$ emissions and $n-1$ state transitions.

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
These factors are $p_I$, $e_I(S)$, $a_{II}$, $e_I(S)$, $a_{IB}$, and $e_B(W)$ in sequence order. The path starts in an island, stays there for the second symbol, then switches to background for W.
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

Scoring one supplied path takes only $\mathcal O(n)$ multiplications. The inference problem is that no path is supplied. With $K=|Q|$ states, every one of the $n$ positions has $K$ choices, giving $K^n$ candidate paths.

When we observe only $x$, every possible hidden path is an alternative way to produce it. These alternatives are mutually exclusive, so their joint probabilities add:
$$
P(x)=\sum_{\pi\in Q^n}P(x,\pi).
$$
For $x=\texttt{SSW}$, all eight paths sum to $P(x)=0.1368$. The best path, $\texttt{IIB}$, contains only $0.05376/0.1368\approx39.3\%$ of that probability. A best explanation and the total evidence are different quantities.

::: {.algorithm}
**Algorithm: sum over a tiny hidden-path space.**

1. Generate all $K^n$ paths.
2. Compute each path joint probability.
3. Sum the probabilities to obtain $P(x)$ and compare each path with this total.
:::

```{.python .execute}
#| label: code-hidden-markov-model-enumerate
def enumerate_paths(states: str, n: int) -> list[str]:
    """List all state paths of a specified length."""
    return ["".join(path) for path in itertools.product(states, repeat=n)]

paths = enumerate_paths("BI", 3)
probabilities = [
    joint_probability("SSW", path, initial, transition, emission)
    for path in paths
]
sum(probabilities), paths[probabilities.index(max(probabilities))]
```

The total converts path scores into probabilities after seeing SSW:
$$
P(\pi\mid\texttt{SSW})=\frac{P(\texttt{SSW},\pi)}{P(\texttt{SSW})}.
$$
The denominator is shared by every path, so the most probable joint path is also the most probable path conditional on the observations. But choosing it still discards the alternatives visible in the figure.

```{.python .execute}
#| label: fig-hidden-markov-model-path-mass
#| fig-cap: "Probabilities of all eight hidden paths after observing SSW, sorted from largest to smallest. IIB is the best explanation, but the remaining paths together hold most of the probability."
#| fig-alt: "Eight bars labeled by length-three B/I paths; IIB is the tallest and is highlighted, with height about 0.393."
#| code-fold: true
ranked = sorted(zip(paths, probabilities), key=lambda item: item[1], reverse=True)
path_names = [path for path, probability in ranked]
path_mass = [probability / sum(probabilities) for path, probability in ranked]
fig, ax = plt.subplots()
bars = ax.bar(path_names, path_mass)
bars[0].set_color(cf.EMPHASIS)
ax.set_xlabel("Hidden path")
ax.set_ylabel("Probability given SSW")
for tick in ax.get_xticklabels():
    tick.set_fontfamily(cf.MONO)
plt.show()
```

At 20 positions the two-state model already has $2^{20}=1{,}048{,}576$ paths. At 100 positions it has about $1.27\times10^{30}$. Scoring each separately costs $\mathcal{O}(nK^n)$ operations. The obstacle is the number of explanations, rather than the cost of scoring one. Dynamic programming will let us combine paths sharing a state at the same position instead of enumerating their complete histories.

## Parameters and estimation

If region labels are supplied for training, we can count first states, state transitions, and symbols emitted in each state, then normalize the appropriate counts. When only sequences are supplied, the state counts themselves are unknown. Estimation must then account for several possible paths, which motivates the inference methods that follow.

::: {#exm-hmm-supervised-counts}
## Estimating emissions from a labeled path
Suppose the assembled DNA fragment is `CGATC`, and an annotation marks its first two positions as island and the remaining three as background. Converting C/G to S and A/T to W gives the observations below; expanding the annotated interval into one label per base gives the path. This supplied annotation is training information, not a path inferred from these five bases:

$$
x=\texttt{SSWWS},\qquad \pi=\texttt{IIBBB}.
$$
State $I$ emits two $S$ symbols and no $W$ symbols; state $B$ emits one $S$ and two $W$ symbols. Without pseudocounts, the fitted emissions are
$$
\widehat e_I(S)=1,\quad \widehat e_I(W)=0,
\qquad
\widehat e_B(S)=1/3,\quad \widehat e_B(W)=2/3.
$$
The annotation source defines what the model learns to call an island. Labels based on a sequence-composition rule teach that rule; they do not by themselves establish a gene-regulatory function.

The path contains transitions $II,IB,BB,BB$, giving $\widehat a_{II}=\widehat a_{IB}=1/2$ and $\widehat a_{BB}=1$. These exact zeros and ones reflect a tiny sample, not biological certainty. Pseudocounts prevent these small-sample estimates from ruling out unobserved events, as in Lesson 8.
:::

The two-state S/W model represents GC enrichment and persistence. To retain the CpG transition signal as well, a state can include both the region label and the current DNA base. The eight states $B_A,B_C,B_G,B_T,I_A,I_C,I_G,I_T$ then distinguish, for example, a C-to-G transition inside an island from the same transition in background [@durbin1998]. The two-state simplification lets us study hidden paths without that additional bookkeeping.

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

::: {#exr-hmm-normalization}
Show that $\sum_x\sum_\pi P(x,\pi)=1$ at a fixed length $n$ by summing out the last observation and last hidden state, then repeating backward. Which row-normalization assumptions does each step use?
:::

::: {#exr-hmm-memory}
Describe an HMM state space that remembers both region ($B$ or $I$) and the previous DNA base. How many states does it have, and what dependence can it represent that the two-state model cannot?
:::

::: {#refs}
**References**
:::
