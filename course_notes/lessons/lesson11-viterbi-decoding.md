---
published: true
title: "Viterbi Decoding"
subtitle: "Lesson 11 · 2027-02-22 · Week 7"
nocite: |
  @durbin1998, @compeau2015
---

## Core question

Which map of island-like regions best explains the observed DNA?

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

An assembled DNA sequence supplies letters and positions, but no island labels. Suppose the observed stretch is `CGAC`. Under the model from Lesson 10, one possible annotation puts its first two bases in an island and the last two in background; another puts all four in an island. Which complete annotation is better supported? We want one label per base, so consecutive equal labels can be combined into intervals on a genome map.

The hidden-region model assigns a probability to every such explanation. Deciding each label from its base alone would ignore the model's preference for regions that continue across neighboring positions. We must compare whole paths, including the changes of region that each explanation requires.

We retain the model from Lesson 10. Observations are $S$ for C or G and $W$ for A or T, so `CGAC` becomes $x=\texttt{SSWS}$. Hidden labels are background $B$ and island-like $I$, so $Q=\{B,I\}$ and the number of states is $K=2$. This simplified alphabet records GC enrichment, not the dinucleotide dependence used in the earlier CpG classifier. *Viterbi decoding* chooses the most probable complete path, balancing emissions against state transitions. Its output is an annotation under this model, rather than a confidence value for each base.

Recall the fixed illustrative parameters:
$$
p_B=p_I=0.5,\qquad
A=\begin{array}{c|cc}&B&I\\\hline B&0.8&0.2\\I&0.4&0.6\end{array},\qquad
E=\begin{array}{c|cc}&S&W\\\hline B&0.3&0.7\\I&0.8&0.2\end{array}.
$$
Here $p_k$ is the initial state probability, $a_{\ell k}$ is the probability of moving from $\ell$ to $k$, and $e_k(b)$ is the probability of observing $b$ in state $k$. Rows of $A$ are departure states; rows of $E$ are emitting states. The joint probability of observations and a proposed path is
$$
P(x,\pi)=p_{\pi_1}e_{\pi_1}(x_1)
\prod_{i=2}^n a_{\pi_{i-1}\pi_i}e_{\pi_i}(x_i).
$$
Given only $x$, decoding chooses the path with the largest such product.

::: {#def-viterbi-decoding}
## Viterbi decoding
For an HMM and observed sequence $x=x_1\ldots x_n$ with $P(x)>0$, the *Viterbi path* is
$$
\pi^*=\underset{\pi\in Q^n}{\operatorname{argmax}}\ P(x,\pi).
$$
:::

Because $P(x)$ is fixed once $x$ is observed, this also maximizes $P(\pi\mid x)=P(x,\pi)/P(x)$. With $K$ states, direct enumeration evaluates $K^n$ paths. The B/I model already has more than a million paths at $n=20$.

## Best prefixes

The alignment dynamic program kept one best score for each pair of consumed prefixes. Here the corresponding subproblem is a sequence prefix ending in a specified hidden state. Different histories can reach the same state at the same position. Once there, the HMM assigns either history the same probability factors for any specified continuation: future transitions depend on the current state, and emissions depend on the states that produce them. Multiplying by these same nonnegative factors preserves which prefix is better. We therefore keep only one best score per state and position.

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
For each previous state, multiply its best prefix probability by the transition into $k$. The current symbol is emitted from the destination state $k$, so $e_k(x_i)$ is common to every predecessor candidate and can sit outside the maximum. This step keeps the best history ending at $k$; it does not choose one best state for the whole column.

A backpointer records which $\ell$ attained the maximum. The final maximum gives the best score, and following pointers backward reconstructs the path.

::: {.algorithm}
**Algorithm: Viterbi decoding.**

1. Initialize $v_k(1)=p_ke_k(x_1)$ for every state $k$.
2. For $i=2,\ldots,n$ and every state $k$, compute the recurrence and store the maximizing previous state.
3. Choose the state with largest $v_k(n)$.
4. Follow backpointers from position $n$ to position 1, then reverse the states.
:::

## A four-position trellis

Use the unchanged B/I model from Lesson 10 and observation $x=\texttt{SSWS}$.

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
At the W observation, the best route into background is
$$
v_B(3)=0.7\max\{0.048(0.8),0.192(0.4)\}=0.05376.
$$
The winning predecessor is I: its stronger earlier prefix outweighs the smaller I-to-B transition probability. Continuing gives

| state | $S$, $i=1$ | $S$, $i=2$ | $W$, $i=3$ | $S$, $i=4$ |
|---|---:|---:|---:|---:|
| $B$ | 0.150 | 0.048 | 0.053760 | **0.0129024** |
| $I$ | 0.400 | 0.192 | 0.023040 | 0.0110592 |

For the last S, the two destination scores compare different best histories:
$$
v_B(4)=0.3\max\{0.05376(0.8),0.02304(0.4)\}=0.0129024,
$$
$$
v_I(4)=0.8\max\{0.05376(0.2),0.02304(0.6)\}=0.0110592.
$$
The best route into B stays in B; the best route into I stays in I, giving the competing path IIII. The final maximum is in B. Its backpointers give
$$
B_4\leftarrow B_3\leftarrow I_2\leftarrow I_1,
$$
so $\pi^*=\texttt{IIBB}$. Direct multiplication gives $P(\texttt{SSWS},\texttt{IIBB})=0.0129024$.
:::

The last S is assigned to background even though $e_I(S)>e_B(S)$. Comparing the complete histories favors IIBB over IIII. The next example shows why retaining both ending states is essential.

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
## Why we keep both states
For the shorter observation $\texttt{WSS}$, the two rows are

| State | $W$, $i=1$ | $S$, $i=2$ | $S$, $i=3$ |
|---|---:|---:|---:|
| $B$ | 0.35 | 0.084 | 0.02016 |
| $I$ | 0.10 | 0.056 | 0.02688 |

At position 2, the best path ending in I is BI, with probability $0.35(0.2)(0.8)=0.056$. It loses to BB at 0.084. But another S favors continuing in I: BII reaches $0.056(0.6)(0.8)=0.02688$, beating BBB at $0.084(0.8)(0.3)=0.02016$. Keeping only the column winner would discard BI before it could become the best complete path.
:::

```{.python .execute}
#| label: fig-viterbi-column-loser
#| fig-cap: "Viterbi scores for WSS. The optimal path BII passes through the smaller entry in column 2, so a best prefix must be retained for each ending state."
#| fig-alt: "A two-row three-column trellis for W S S. Red nodes follow B, I, I; the middle I score 0.056 is smaller than the B score 0.084."
#| fig-width: 6
#| fig-height: 3
fig, ax = plt.subplots()
cf.trellis(
    ax, ["B", "I"], list("WSS"),
    [[0.35, 0.084, 0.02016], [0.10, 0.056, 0.02688]],
    backpointers=[((0, 0), (0, 1)), ((0, 0), (1, 1)),
                  ((0, 1), (0, 2)), ((1, 1), (1, 2))],
    path=[(0, 0), (1, 1), (1, 2)], fmt="{:.5g}",
)
plt.show()
```

### How several bases support a region

The transition probabilities express a biological modeling choice: region labels tend to persist. How much sequence evidence pays for entering and leaving an island? Consider $r\geq1$ S symbols between two positions whose labels are fixed to B. Compare labeling the whole run B with labeling it I. The flanking emissions and all factors outside the interval cancel, leaving
$$
\frac{P(\text{island candidate},x)}{P(\text{all-B candidate},x)}
=\frac{0.2(0.6)^{r-1}(0.4)}{(0.8)^{r+1}}
\left(\frac{0.8}{0.3}\right)^r
=\frac{2^{r-1}}{3}.
$$
The numerator pays for entering I, staying there, and returning to B. Each S favors I through the emission ratio $0.8/0.3$. For run lengths 1, 2, and 3, the combined ratios are $1/3$, $2/3$, and $4/3$: the three-S run favors the island candidate. This comparison explains how several individually ambiguous bases can support a region together. It assumes fixed flanking labels; Viterbi also compares paths with different boundaries.

## Log-space implementation

Long products underflow. For example, $0.25^{1000}$ is smaller than the representable positive range of a standard Python float. Use natural logarithms and define $V_k(i)=\ln v_k(i)$. The recurrence becomes
$$
V_k(i)=\ln e_k(x_i)+
\max_{\ell\in Q}\left[V_\ell(i-1)+\ln a_{\ell k}\right].
$$
Products become sums, and the maximizing path is unchanged because the logarithm is increasing. Initialization is $V_k(1)=\ln p_k+\ln e_k(x_1)$. The code stores the current scores and the previous-state choices needed for traceback.

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
We can check these probabilities without using the recurrence.

::: {.algorithm}
**Algorithm: check decoding by exhaustive paths.**
Enumerate all length-$n$ state sequences. For each, multiply its initial, transition, and emission factors directly. Sort the resulting joint probabilities and compare their maximum with the dynamic-programming output.
:::

```{.python .execute}
#| label: code-viterbi-enumeration
observations = "SSWS"
ranked_paths = []
for labels in itertools.product(states, repeat=len(observations)):
    candidate = "".join(labels)
    probability = initial[labels[0]] * emission[labels[0]][observations[0]]
    for i in range(1, len(observations)):
        probability *= transition[labels[i - 1]][labels[i]]
        probability *= emission[labels[i]][observations[i]]
    ranked_paths.append((candidate, probability))
ranked_paths.sort(key=lambda item: item[1], reverse=True)
assert ranked_paths[0][0] == path
assert math.isclose(ranked_paths[0][1], math.exp(log_probability))
for candidate, probability in ranked_paths[:3]:
    print(candidate, round(probability, 10))
```

The dynamic program agrees with direct enumeration on this instance. Two additional checks catch boundary errors: a length-one sequence must return $\arg\max_k p_ke_k(x_1)$, and deterministic emissions must force the only compatible path.

A decoded path becomes intervals by merging consecutive equal labels. For IIBB, positions 1 through 2 form the predicted island-like region. This is a concrete annotation under the model, but it leaves a new question: how much support does each label have? The best path accounts for only $0.0129024/0.060576\approx21.3\%$ of the probability conditional on SSWS. The next lesson adds the probabilities of alternative paths to assess individual positions.

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
