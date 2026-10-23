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
import random
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## The decoding problem

An assembled DNA sequence supplies letters and positions, but no island labels. Suppose the observed stretch is `CGAC`. Under the model from Lesson 10, one possible annotation puts its first two bases in an island and the last two in background; another puts all four in an island. Which complete annotation is better supported? Deciding each label from its base alone would ignore the model's preference for regions that continue across neighboring positions, so we must compare whole paths, including the changes of region that each one requires. The winning path gives one label per base, and consecutive equal labels combine into intervals on a genome map.

We retain the model from Lesson 10: observations are $S$ for C or G and $W$ for A or T, so `CGAC` becomes $x=\texttt{SSWS}$, and hidden labels are background $B$ and island-like $I$, so $Q=\{B,I\}$ and $K=2$. *Viterbi decoding* chooses the most probable complete path, balancing emissions against state transitions, and computes it with the same kind of table-filling dynamic program that found optimal alignments. Its output is an annotation under this model, not a confidence value for each base.

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

Because $P(x)$ is fixed once $x$ is observed, this also maximizes $P(\pi\mid x)=P(x,\pi)/P(x)$. With $K$ states, direct enumeration evaluates $K^n$ paths, more than a million for the B/I model at $n=20$.

## Best prefixes

Decoding is an alignment problem in a new setting. Draw a grid, called a *trellis*, with one row per state and one column per position. A hidden path picks one node in each column, so it is a left-to-right route through the grid, just as a global alignment was a route through the alignment table $F(i,j)$. Taking logarithms turns the path probability into a sum of node and edge scores:
$$
\ln P(x,\pi)=
\underbrace{\ln p_{\pi_1}+\ln e_{\pi_1}(x_1)}_{\text{first node}}
+\sum_{i=2}^{n}\Bigl[
\underbrace{\ln a_{\pi_{i-1}\pi_i}}_{\text{edge into column } i}
+\underbrace{\ln e_{\pi_i}(x_i)}_{\text{node in column } i}
\Bigr].
$$
Read from left to right, the sum starts at a node in column 1 and then, at each later column, adds the score of the edge taken and the emission score of the node reached. Decoding therefore asks for the highest-scoring route through the trellis, as global alignment asked for the highest-scoring route through its table, and we solve it the same way: one best score per cell, filled column by column.

The alignment program kept one best score for each pair of consumed prefixes. Here the subproblem is a sequence prefix ending in a specified hidden state, which is one node of the trellis. Different histories can reach the same node. Once there, the HMM assigns both histories the same probability factors for any continuation, because future transitions depend only on the current state and each future emission depends only on its own state. Multiplying two prefix probabilities by the same nonnegative factors preserves which one is larger, so we keep only one best score per state and position.

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
The numerator pays for entering I, staying there, and returning to B. Each S favors I through the emission ratio $0.8/0.3$. For run lengths 1, 2, and 3, the combined ratios are $1/3$, $2/3$, and $4/3$: roughly three S symbols are needed to pay for entering an island and leaving it again. This comparison explains how several individually ambiguous bases can support a region together. It assumes fixed flanking labels; Viterbi also compares paths with different boundaries.

The same accounting runs in reverse inside an island. Put $r$ W symbols between two positions fixed to I, and compare keeping the whole run in I with dropping it into B:
$$
\frac{P(\text{all-I candidate},x)}{P(\text{gap candidate},x)}
=\frac{(0.6)^{r+1}}{0.4(0.8)^{r-1}(0.2)}
\left(\frac{0.2}{0.7}\right)^r .
$$
For one W the ratio is $(0.36/0.08)(2/7)\approx1.29$, so the island continues through a single W. For two W symbols it is $(0.216/0.064)(4/49)\approx0.28$, and the island breaks. With fixed flanks, then, an island under this model absorbs an isolated A or T but not two in a row.

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

In the language of the log-score display, the trellis is a directed acyclic graph and Viterbi finds its highest-scoring path; the proof above is the alignment table's optimal-substructure argument with states in place of alignment moves.

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

The dynamic program agrees with direct enumeration on this instance. Two additional checks catch boundary errors: a length-one sequence must return $\arg\max_k p_ke_k(x_1)$, and if each state deterministically emits its own symbol, the emissions must force the only compatible path.

## Decoding a simulated sequence

The four-position example shows the recurrence at work but cannot show whether decoding recovers regions; for that we need a sequence with known true labels. Lesson 10 simulated 300 positions from a sticky variant of the model, with $p$ and $E$ unchanged and transitions
$$
A_{\text{sticky}}=
\begin{array}{c|cc}
 & B&I\\\hline
B&0.98&0.02\\
I&0.05&0.95
\end{array}.
$$
The chunk below regenerates that sequence and decodes it with `viterbi` under the sticky transitions that generated it. For comparison, it also makes window calls in the spirit of Lesson 9: center a window of odd width $w$ at position $i$, add the emission log ratios of its symbols, and call position $i$ island when the total is positive:
$$
\Lambda_w(i)=\sum_{j=i-h}^{i+h}\log_2\frac{e_I(x_j)}{e_B(x_j)}>0,
\qquad h=\frac{w-1}{2}.
$$
Each S adds $\log_2(0.8/0.3)\approx1.42$ bits and each W adds $\log_2(0.2/0.7)\approx-1.81$ bits, so a window is called island when more than about 56% of its symbols are S; windows are truncated at the sequence ends.

```{.python .execute}
#| label: fig-viterbi-simulated
#| code-fold: true
#| code-summary: "Simulate, decode, and call windows"
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "One simulated sequence of 300 positions from the sticky variant, decoded three ways. Ticks mark observed S symbols; grey bars and shading show the true island runs. Viterbi (red) recovers the six longer islands with nearly correct boundaries, while short windows (w = 5) flicker and call false islands and long windows (w = 31) miss or blur the shorter islands."
#| fig-alt: "Five rows over positions 1 to 300: tick marks for S symbols, then grey bars for true islands, red bars for the Viterbi path, and blue bars for window calls of width 5 and 31. Light grey shading behind all rows marks the true islands."
def simulate_hmm(
    n: int,
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    emission: dict[str, dict[str, float]],
    rng: random.Random,
) -> tuple[str, str]:
    """Draw observations and a hidden path of length n (as in Lesson 10)."""

    def draw(distribution: dict[str, float]) -> str:
        """Draw one outcome with the probabilities stored as values."""
        return rng.choices(list(distribution), weights=list(distribution.values()))[0]

    state = draw(initial)
    path = [state]
    observations = [draw(emission[state])]
    for _ in range(1, n):
        state = draw(transition[state])
        path.append(state)
        observations.append(draw(emission[state]))
    return "".join(observations), "".join(path)


def window_calls(
    observations: str, width: int, emission: dict[str, dict[str, float]]
) -> str:
    """Label each position by the sign of its centered window's log ratio."""
    half = width // 2
    labels = []
    for i in range(len(observations)):
        window = observations[max(0, i - half):i + half + 1]
        score = 0.0
        for symbol in window:
            score += math.log2(emission["I"][symbol] / emission["B"][symbol])
        labels.append("I" if score > 0 else "B")
    return "".join(labels)


sticky_transition = {
    "B": {"B": 0.98, "I": 0.02},
    "I": {"B": 0.05, "I": 0.95},
}
x_sim, path_sim = simulate_hmm(
    300, initial, sticky_transition, emission, random.Random(6)
)
viterbi_sim, _ = viterbi(x_sim, states, initial, sticky_transition, emission)
rows = {
    "True states": path_sim,
    "Viterbi path": viterbi_sim,
    "Windows, w = 5": window_calls(x_sim, 5, emission),
    "Windows, w = 31": window_calls(x_sim, 31, emission),
}
colors = [cf.RULE, cf.EMPHASIS, cf.SECOND, cf.SECOND]

fig, ax = plt.subplots()
true_island = [i + 1 for i, state in enumerate(path_sim) if state == "I"]
ax.bar(true_island, 6.1, bottom=-0.5, width=1, color=cf.LIGHT, linewidth=0)
s_positions = [i + 1 for i, symbol in enumerate(x_sim) if symbol == "S"]
ax.vlines(s_positions, 4.6, 5.4, color=cf.INK, linewidth=0.6)
# One row per labeling; each island position is a unit-width bar.
for row, (labels, color) in enumerate(zip(rows.values(), colors)):
    island = [i + 1 for i, label in enumerate(labels) if label == "I"]
    ax.bar(island, 0.6, bottom=2.7 - row, width=1, color=color, linewidth=0)
ax.set_yticks([5, 3, 2, 1, 0], ["Observed S"] + list(rows))
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set(xlim=(0.5, 300.5), ylim=(-0.5, 5.6), xlabel="Position")
plt.show()
```

@fig-viterbi-simulated compares the three answers with the truth. Viterbi recovers the six longer islands, usually within a few positions of their true boundaries, and every island it calls overlaps a true one. It misses the two-position island at the start, and it extends the last island 11 positions into background that happened to emit S at 7 of them. Both behaviors follow from the price of switching: under the sticky transitions, the fixed-flank comparison for a run of $r$ S symbols becomes
$$
\frac{0.02\,(0.95)^{r-1}(0.05)}{(0.98)^{r+1}}\left(\frac{0.8}{0.3}\right)^r,
$$
which first exceeds one at $r=8$. A short burst of S therefore cannot pay to open an island, while a long S-rich run can.

The windows have no switching price, only a width. Width 5 reacts quickly but flickers, splitting one island in two and calling seven short false islands where a few S symbols cluster in background. Width 31 removes the flicker but dilutes shorter islands with surrounding background, missing those near positions 100 and 220 and nearly all of the one at 183 to 202. The HMM replaces the window width by transition probabilities, which state how long regions are expected to last, and lets each boundary fall wherever the whole sequence supports it best.

Merging consecutive equal labels turns a decoded path into intervals: for IIBB, positions 1 and 2 form the predicted island, and the Viterbi row of @fig-viterbi-simulated lists six intervals. A single path cannot say how firmly each label is supported. It draws a boundary as sharply as an island's interior, although letters near a boundary often fit several nearby endpoints almost equally well, and on SSWS the best path holds only $0.0129024/0.060576\approx21.3\%$ of the probability conditional on the observations. Measuring the support for each label requires adding the probabilities of the alternative paths, which is the computation of the next lesson.

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
