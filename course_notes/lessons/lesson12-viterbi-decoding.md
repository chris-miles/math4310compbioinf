---
published: true
title: "Viterbi Decoding"
subtitle: "Lesson 12 · 2027-03-01 · Week 8"
nocite: |
  @durbin1998, @compeau2015
---

## Core question

What is the most likely hidden path through an HMM?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import itertools
import math
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## The decoding problem

A hidden-region model assigns a probability to every possible annotation, but a genome browser needs a concrete set of intervals. We therefore want one complete path of labels through the sequence. Deciding each label from its base alone would ignore the model's preference for regions that continue across neighboring positions.

We retain the B/I model and observed S/W string from Lesson 11. Viterbi decoding chooses the most probable complete path, balancing emissions against the cost of switching states. The output is an annotation under this model; it will not tell us how confidently each individual base should receive its label.

Lesson 11 assigned a joint probability $P(x,\pi)$ to every observation and hidden path. Given only $x$, decoding chooses the path with the largest joint probability.

::: {#def-viterbi-decoding}
## Viterbi decoding
For an HMM and observed sequence $x=x_1\ldots x_n$ with $P(x)>0$, the *Viterbi path* is
$$
\pi^*=\underset{\pi\in Q^n}{\operatorname{argmax}}\ P(x,\pi).
$$
:::

Because $P(x)$ is fixed once $x$ is observed, this also maximizes $P(\pi\mid x)=P(x,\pi)/P(x)$. With $K$ states, direct enumeration evaluates $K^n$ paths. The B/I model already has more than a million paths at $n=20$.

## Best prefixes

Many complete paths share a prefix. Once two paths reach the same state at the same position, the better prefix stays better after either is extended by the same transition and emission. We therefore keep only one score per state and position.

::: {#def-viterbi-score}
## Viterbi score
For state $k$ at position $i$, define
$$
v_k(i)=
\max_{\pi_1,\ldots,\pi_{i-1}}
P(x_1\ldots x_i,\pi_1\ldots\pi_{i-1},\pi_i=k).
$$
It is the probability of the most likely partial path that emits $x_1\ldots x_i$ and ends in $k$.
:::

At the first position,
$$
v_k(1)=p_k e_k(x_1).
$$
For later positions, the previous state could be any $\ell\in Q$:
$$
v_k(i)=e_k(x_i)\max_{\ell\in Q}
\left[v_\ell(i-1)a_{\ell k}\right].
$$
A backpointer records which $\ell$ attained the maximum. The final maximum gives the best score, and following pointers backward reconstructs the path.

::: {.algorithm}
**Algorithm: Viterbi decoding.**

1. Initialize $v_k(1)=p_ke_k(x_1)$ for every state $k$.
2. For $i=2,\ldots,n$ and every state $k$, compute the recurrence and store the maximizing previous state.
3. Choose the state with largest $v_k(n)$.
4. Follow backpointers from position $n$ to position 1, then reverse the states.
:::

## A four-position trellis

Use the unchanged B/I model from Lesson 11 and observation $x=\texttt{SSWS}$.

::: {#exm-viterbi-table}
## Filling the Viterbi table
At position 1,
$$
v_B(1)=0.5(0.3)=0.15,\qquad
v_I(1)=0.5(0.8)=0.4.
$$
For the second $S$,
$$
v_B(2)=0.3\max\{0.15(0.8),0.4(0.4)\}=0.048,
$$
$$
v_I(2)=0.8\max\{0.15(0.2),0.4(0.6)\}=0.192.
$$
Continuing gives

| state | $S$, $i=1$ | $S$, $i=2$ | $W$, $i=3$ | $S$, $i=4$ |
|---|---:|---:|---:|---:|
| $B$ | 0.150 | 0.048 | 0.053760 | **0.0129024** |
| $I$ | 0.400 | 0.192 | 0.023040 | 0.0110592 |

The final maximum is in $B$. Its backpointers give
$$
B_4\leftarrow B_3\leftarrow I_2\leftarrow I_1,
$$
so $\pi^*=\texttt{IIBB}$. Direct multiplication gives $P(\texttt{SSWS},\texttt{IIBB})=0.0129024$.
:::

The last $S$ is assigned to background. Switching from $B$ back to $I$ costs a factor of $0.2$, and one island-favored emission does not compensate.

```{.python .execute}
#| label: fig-viterbi-trellis
#| fig-cap: "The Viterbi trellis for SSWS. Each circle holds the best path probability ending at that state and position; the red circles form the decoded path IIBB."
#| fig-alt: "A two-row, four-column trellis labeled B and I with observations S S W S above; a red path passes through I, I, B, B."
#| fig-width: 6
#| fig-height: 3
values = [
    [0.15, 0.048, 0.05376, 0.0129024],
    [0.4, 0.192, 0.02304, 0.0110592],
]
path = [(1, 0), (1, 1), (0, 2), (0, 3)]

fig, ax = plt.subplots()
backpointers = [
    ((1, 0), (0, 1)), ((1, 0), (1, 1)),
    ((1, 1), (0, 2)), ((1, 1), (1, 2)),
    ((0, 2), (0, 3)), ((1, 2), (1, 3)),
]
cf.trellis(
    ax, ["B", "I"], list("SSWS"), values,
    backpointers=backpointers, path=path,
)
plt.show()
```


::: {#exm-viterbi-greedy}
## Emissions alone choose a different path
Choosing the state with the larger emission probability at each position gives $\texttt{IIBI}$ for $\texttt{SSWS}$. Its joint probability is
$$
0.5(0.8)(0.6)(0.8)(0.4)(0.7)(0.2)(0.8)=0.0086016,
$$
which is smaller than $0.0129024$ for $\texttt{IIBB}$. The emission-only choice gains a factor $0.8/0.3$ at the final $S$, but loses a factor $0.2/0.8$ by switching instead of staying in background. The combined factor is $2/3$.
:::

Viterbi keeps a best prefix for each possible ending state because the next transition depends on that state. It cannot replace the whole column with a single winning prefix: a prefix that loses now may have a much better transition to a state needed later. The proof identifies exactly which histories can be merged safely.

### Evidence needed to introduce a region

Consider a run of $r\geq1$ observed S symbols between two positions whose hidden labels are fixed to B. Compare two candidate paths inside this interval: keeping all $r$ positions in B, or changing them all to I and returning to B afterward. Factors outside the interval, including the flanking emissions, cancel in the probability ratio.

The all-B candidate uses $r+1$ transitions of probability 0.8 and $r$ emissions of probability 0.3. The island candidate uses a B-to-I transition of probability 0.2, $r-1$ I self-transitions of probability 0.6, an I-to-B transition of probability 0.4, and $r$ emissions of probability 0.8. Consequently,
$$
\frac{P(\text{island candidate},x)}{P(\text{all-B candidate},x)}
=\frac{0.2(0.6)^{r-1}(0.4)}{(0.8)^{r+1}}
\left(\frac{0.8}{0.3}\right)^r
=\frac{2^{r-1}}{3}.
$$
For run lengths 1, 2, and 3, the ratios are $1/3$, $2/3$, and $4/3$. Three consecutive S emissions make this island candidate preferable; one or two do not compensate for introducing and ending the region.

This is a comparison of two paths with fixed flanking states, not a minimum-length rule for every Viterbi island. At a sequence boundary the initial probability replaces one transition, and uncertain flanking labels can change the comparison. The calculation shows explicitly how evidence from several emissions can pay for a state change. Viterbi performs all such compatible comparisons together.

## Log-space implementation

Long products underflow. For example, $0.25^{1000}$ is smaller than the representable positive range of a standard Python float. Define $V_k(i)=\log v_k(i)$. The recurrence becomes
$$
V_k(i)=\log e_k(x_i)+
\max_{\ell\in Q}\left[V_\ell(i-1)+\log a_{\ell k}\right].
$$
Products become sums, and the maximizing path is unchanged.

```{.python .execute}
#| label: code-viterbi-decoding
def viterbi(
    observations: str,
    states: list[str],
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    emission: dict[str, dict[str, float]],
) -> tuple[str, float]:
    """Return the Viterbi path and its natural-log joint probability."""
    if not observations:
        return "", 0.0

    def safe_log(probability: float) -> float:
        """Map an impossible event to negative infinity."""
        return -math.inf if probability == 0.0 else math.log(probability)

    scores = {
        state: safe_log(initial[state]) + safe_log(emission[state][observations[0]])
        for state in states
    }
    backpointers = []

    for symbol in observations[1:]:
        next_scores = {}
        pointers = {}
        for current in states:
            candidates = {
                previous: scores[previous]
                + safe_log(transition[previous][current])
                for previous in states
            }
            best_previous = max(candidates, key=candidates.get)
            next_scores[current] = (
                candidates[best_previous] + safe_log(emission[current][symbol])
            )
            pointers[current] = best_previous
        scores = next_scores
        backpointers.append(pointers)

    if all(score == -math.inf for score in scores.values()):
        raise ValueError("The model cannot emit this observation sequence.")

    final_state = max(scores, key=scores.get)
    path = [final_state]
    for pointers in reversed(backpointers):
        path.append(pointers[path[-1]])
    path.reverse()
    return "".join(path), scores[final_state]


states = ["B", "I"]
initial = {"B": 0.5, "I": 0.5}
transition = {
    "B": {"B": 0.8, "I": 0.2},
    "I": {"B": 0.4, "I": 0.6},
}
emission = {
    "B": {"S": 0.3, "W": 0.7},
    "I": {"S": 0.8, "W": 0.2},
}

path, log_probability = viterbi(
    "SSWS", states, initial, transition, emission
)
path, math.exp(log_probability)
```


The trellis above shows ordinary probabilities for the hand calculation. In code, a zero probability becomes $-\infty$: adding it makes that candidate impossible, and a finite candidate always beats it. If every final score is $-\infty$, the model gives the observed sequence probability zero; returning an arbitrary path would conceal this failure.

For example, change the emissions so $B$ can emit only $W$ and $I$ can emit only $S$, leaving the positive transitions unchanged. The observation $\texttt{SSWS}$ then forces $\texttt{IIBI}$. If we also forbid the transition $B\to I$, the same observation becomes impossible because its last symbol requires that transition. These two cases test different behavior: one must return the unique compatible path, and the other must report that no path exists.

## Correctness and cost

::: {#thm-viterbi-correctness}
## Viterbi recurrence
For every position $i$ and state $k$, the recurrence computes the largest joint probability among all partial paths ending at $k$. Therefore
$$
\max_k v_k(n)=\max_{\pi\in Q^n}P(x,\pi).
$$
Following the stored maximizing predecessors returns a path attaining this value.
:::

::: {.proof}
Sketch. The claim holds at $i=1$ by initialization. Assume it holds at $i-1$. Every path ending at $k$ at position $i$ has some previous state $\ell$, and its probability is its prefix probability times $a_{\ell k}e_k(x_i)$. For a fixed $\ell$, replacing the prefix by the best prefix ending at $\ell$ cannot decrease the probability. Maximizing over $\ell$ therefore considers the best possible final transition. This proves the claim by induction. Backpointers record the choices that attained those maxima.
:::

There are $nK$ table entries. Each entry considers $K$ previous states, so the time cost is $\mathcal{O}(nK^2)$. Scores for only two columns require $\mathcal{O}(K)$ space. Reconstructing the path needs $\mathcal{O}(nK)$ stored backpointers.

The trellis is a directed acyclic graph. Viterbi is a longest-path dynamic program after transition and emission probabilities are replaced by their logs. This is the same optimal-substructure argument used for alignment, with states replacing alignment moves.

### Validation by enumeration

For short observations, compare the algorithm against all paths. On $\texttt{SSWS}$ there are 16 paths. The three largest joint probabilities are
$$
P(x,\texttt{IIBB})=0.0129024,\quad
P(x,\texttt{IIII})=0.0110592,\quad
P(x,\texttt{IIBI})=0.0086016.
$$
The dynamic program agrees with brute force. Two additional checks catch boundary errors: a length-one sequence must return $\arg\max_k p_ke_k(x_1)$, and deterministic emissions must force the only compatible path.

### Ties and alternate optimal paths

The recurrence may have equal candidates. A single backpointer then returns one optimal path, but other optimal paths exist. To count or enumerate all optimal paths, store every predecessor attaining the maximum and trace the resulting backpointer graph. This can produce exponentially many paths even though their optimal score is computed in polynomial time.

::: {#exm-viterbi-tie}
## A model in which every path ties
Let both states have initial probability $1/2$, every transition probability $1/2$, and identical emission probabilities. For any fixed observation of length $n$, every state path has the same joint probability. Viterbi still returns one path because the implementation has a deterministic tie rule. The returned labels contain no evidence about the hidden state; the tie is a property of the model, not a computational failure.
:::

### From a path to intervals

A decoded path becomes a set of intervals by merging consecutive positions with the same state. For $\texttt{IIBB}$, positions 1 through 2 form one island interval and positions 3 through 4 form background. On real coordinates, off-by-one errors enter when converting one-based mathematical positions to zero-based, half-open file formats such as BED. A validation example should include a state change at the first or last possible boundary and check the exact reported coordinates. In zero-based half-open coordinates, the island interval in this example is $[0,2)$, which includes indices 0 and 1.

## Limitations

1. **The best path can carry little total probability.** For this observation, $P(x)=0.060576$, so the Viterbi path contains about $21.3\%$ of the probability mass. Lesson 13 sums across paths.
2. **A hard path hides uncertainty.** Positions with nearly tied state probabilities receive one label. Posterior decoding asks how likely each state is at each position.
3. **The answer is conditional on the model.** Incorrect emission, transition, or state definitions can produce a precisely computed but biologically poor path.
4. **Ties need a convention.** The code returns the first maximizing state. Tests should accept any optimal path unless a tie rule is part of the specification.

## Exercises

::: {#exr-viterbi-hand}
Fill both rows of the Viterbi table for $x=\texttt{SWS}$ and recover the path by traceback.
:::

::: {#exr-viterbi-log}
Recompute $v_B(2)$ and $v_I(2)$ in base-2 log space. Verify that the maximizing predecessors match @exm-viterbi-table.
:::

::: {#exr-viterbi-proof}
Explain why keeping only the largest score at each state is safe, but keeping only the largest score across all states at a position is not.
:::

::: {#exr-viterbi-sparse}
Suppose only $E$ of the $K^2$ state transitions have positive probability. State the time cost of Viterbi when the implementation visits only allowed transitions.
:::

::: {#refs}
**References**
:::
