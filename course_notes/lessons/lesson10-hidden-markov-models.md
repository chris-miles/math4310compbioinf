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
import random
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## Observations and hidden states

The window classifier from Lesson 9 tells us whether a chosen stretch of DNA looks island-like. Because many CpG islands lie near promoters, locating these regions gives us candidate DNA intervals to investigate for roles in controlling when nearby genes are copied into RNA. To make a map along a chromosome, we also need to locate where each region begins and ends, and sliding a window along the DNA does this poorly. Suppose we call a ten-base window island-like when most of its bases are C or G. In `GCGCCATGCGTATAA`, the first ten bases contain eight C or G and pass, the last ten contain three and fail, and the five bases `ATGCG` shared by the two windows receive both labels. Changing the width trades one failure for another: a narrow window reacts to chance clusters of C and G, and a wide one blurs every boundary it straddles [@durbin1998; @compeau2015].

The input is the same kind of data, a DNA letter string with genomic coordinates. The output we want is a map: a label, island or background, at every position, with each maximal run of island labels reported as an interval, for example bases 1,204 to 1,890 of a chromosome. Real CpG islands are typically a few hundred to a few thousand bases long [@durbin1998], so such a map is a short list of intervals separated by long stretches of background. No single letter can place itself on this map, because a C or G occurs in both kinds of region. Compeau and Pevzner compare the situation to a dealer who now and then secretly swaps a fair coin for a loaded one: every flip is visible, the coin is not, and a long run of heads makes the loaded coin likely without proving it [@compeau2015]. Because either region can produce any base, no map is certain, and we instead grade competing maps by how well they explain the letters. Neighbors supply the evidence, since a region extends across many bases, and we model that persistence by letting the unobserved label either stay the same or change from one position to the next.

To keep the arithmetic short, collapse the DNA alphabet into $S$ for $C$ or $G$ and $W$ for $A$ or $T$. The collapse keeps GC enrichment but loses the distinction between CG and GC that the Lesson 9 classifier relied on.

::: {#def-hmm-data}
## Observed sequence and hidden path
An observed sequence is $x=x_1\ldots x_n$ over an alphabet $\Sigma$. A hidden path is $\pi=\pi_1\ldots\pi_n$ over a finite state set $Q$. State $\pi_i$ is the unobserved label associated with observation $x_i$.
:::

Our state set is $Q=\{B,I\}$: background and island-like sequence. For example, DNA `CGA` becomes the observation $x=\texttt{SSW}$. A proposed path $\pi=\texttt{IIB}$ assigns its first two bases to an island and its third to background; $\texttt{BBB}$ proposes no island. Both paths can produce the same letters. A *hidden Markov model* assigns probabilities to these competing explanations, using how often each region emits S or W and how often the region changes.

@fig-hidden-markov-model-observed-hidden keeps the recorded bases separate from the proposed region labels: converting C/G to S is a known recoding, whereas choosing between I and B is an inference.

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

The emission table says that $S$ is more common in islands, and the transition table says that regions tend to persist, background especially. Every row is a probability distribution. The model works like a machine that prints one symbol per position after two random choices: which region to occupy, drawn from $p$ at the first position and from the previous region's transition row afterward, and which symbol to print, drawn from the current region's emission row [@compeau2015]. Only the printed symbols reach us.

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

Independent emissions may seem to contradict the neighboring patterns we want to model. The condition matters: if the path is known, learning that $x_{i-1}=S$ does not change the emission distribution at position $i$ beyond what $\pi_i$ already tells us. Dependence can still appear after the path is hidden: two adjacent $S$ symbols are correlated because both are more likely to lie inside the same persistent $I$ run.

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

### Simulating sequences from the model

Running the machine on a computer shows what kind of data the model treats as typical. If simulated sequences look nothing like a genome, its inferences about a genome deserve suspicion.

::: {.algorithm}
**Algorithm: simulate an HMM.**

1. Draw $\pi_1$ with probabilities $p_k$, then draw $x_1$ with probabilities $e_{\pi_1}(b)$.
2. For $i=2,\ldots,n$, draw $\pi_i$ with probabilities $a_{\pi_{i-1}\ell}$, $\ell\in Q$, then draw $x_i$ with probabilities $e_{\pi_i}(b)$.
3. Return the observations $x$ and the hidden path $\pi$.
:::

With a mean island length of 2.5 symbols, 300 positions simulated from our parameters contain dozens of island runs, most too short to look like regions. To see the model's structure at a scale closer to a chromosome map, we also simulate a *sticky variant* that keeps $p$ and $E$ and raises both self-transitions:
$$
A_{\text{sticky}}=
\begin{array}{c|cc}
 & B&I\\\hline
B&0.98&0.02\\
I&0.05&0.95
\end{array},
$$
so the mean background and island lengths are $1/0.02=50$ and $1/0.05=20$. These are still far shorter than real islands but let several islands fit in 300 positions. Every hand calculation continues to use the original table $A$.

```{.python .execute}
#| label: code-hidden-markov-model-simulate
def simulate_hmm(
    n: int,
    initial: dict[str, float],
    transition: dict[str, dict[str, float]],
    emission: dict[str, dict[str, float]],
    rng: random.Random,
) -> tuple[str, str]:
    """Draw observations and a hidden path of length n from an HMM."""

    def draw(distribution: dict[str, float]) -> str:
        """Draw one outcome with the probabilities stored as values."""
        outcomes = list(distribution)
        weights = list(distribution.values())
        return rng.choices(outcomes, weights=weights)[0]

    state = draw(initial)
    path = [state]
    observations = [draw(emission[state])]
    for _ in range(1, n):
        # The next region depends only on the current region.
        state = draw(transition[state])
        path.append(state)
        observations.append(draw(emission[state]))
    return "".join(observations), "".join(path)


initial = {"B": 0.5, "I": 0.5}
transition = {
    "B": {"B": 0.8, "I": 0.2},
    "I": {"B": 0.4, "I": 0.6},
}
emission = {
    "B": {"S": 0.3, "W": 0.7},
    "I": {"S": 0.8, "W": 0.2},
}
sticky_transition = {
    "B": {"B": 0.98, "I": 0.02},
    "I": {"B": 0.05, "I": 0.95},
}
x_toy, path_toy = simulate_hmm(300, initial, transition, emission, random.Random(3))
x_sim, path_sim = simulate_hmm(
    300, initial, sticky_transition, emission, random.Random(6)
)

# How often does S follow S, compared with S following W?
after_S = []
after_W = []
for previous, current in zip(x_sim, x_sim[1:]):
    if previous == "S":
        after_S.append(current)
    else:
        after_W.append(current)
(
    round(x_sim.count("S") / len(x_sim), 3),
    round(after_S.count("S") / len(after_S), 3),
    round(after_W.count("S") / len(after_W), 3),
)
```

```{.python .execute}
#| label: fig-hidden-markov-model-simulated
#| code-fold: true
#| code-summary: "Draw the two simulated sequences"
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "Three hundred positions simulated from the B/I model (A) and from the sticky variant with the same emissions (B). Each black tick is an observed S; grey bars and shading mark the hidden island runs, which the simulation records but an observer would not see. S symbols crowd into islands without being confined to them, so a dense run of S is evidence for an island even though no single letter is."
#| fig-alt: "Two panels over positions 1 to 300, each with a row of vertical tick marks above a row of grey bars. In panel A the bars are many and narrow. In panel B about six wider bars sit under dense clusters of ticks, and ticks between them are sparse."
fig, axes = plt.subplots(2, 1, sharex=True)
panels = [
    ("A", "B/I model (mean island 2.5)", x_toy, path_toy),
    ("B", "Sticky variant (mean island 20)", x_sim, path_sim),
]
for ax, (letter, title, symbols, hidden) in zip(axes, panels):
    island = [i + 1 for i, state in enumerate(hidden) if state == "I"]
    # Light shading behind the ticks and a solid bar below them mark hidden islands.
    ax.bar(island, 1.35, width=1, align="center", color=cf.LIGHT, linewidth=0)
    ax.bar(island, 0.3, width=1, align="center", color=cf.RULE, linewidth=0)
    s_positions = [i + 1 for i, symbol in enumerate(symbols) if symbol == "S"]
    ax.vlines(s_positions, 0.45, 1.35, color=cf.INK, linewidth=0.6)
    ax.set(ylim=(0, 1.35))
    ax.set_yticks([0.15, 0.9], ["hidden I", "observed S"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title(f"{letter}   {title}", loc="left")
axes[1].set(xlim=(0.5, 300.5), xlabel="Position")
plt.show()
```

Panel B of @fig-hidden-markov-model-simulated shows what the model expects a chromosome map to look like. S ticks crowd into the shaded islands, as $e_I(S)=0.8$ requires, but background still emits S with probability 0.3 and islands emit W with probability 0.2, so stray S symbols appear between islands and gaps appear inside them. A boundary shows up as a change in tick density, not at any single letter. Panel A, with the original transitions, has many one- and two-position runs: the table that keeps hand calculations short is a poor description of regions hundreds of bases long.

Persistence also produces the dependence computed exactly above. In the sticky simulation, S makes up 45% of all symbols but 64% of the symbols that follow an S and only 30% of those that follow a W, although each symbol was drawn from its own state's emission row alone; neighboring positions usually share a hidden state. Recovering the shaded intervals from the ticks is the decoding problem, and because this sequence comes with its true answer, we will use it to check each decoding method.

## Joint probability of sequence and path

We can now translate the two choices at each position into the probability of one proposed explanation. We condition on a fixed observed length $n$: the model has no extra factor for stopping. A length-$n$ path contains $n$ emissions and $n-1$ state transitions, and for any observation $x$ and path $\pi$ of length $n$ the joint probability multiplies them:
$$
P(x,\pi)
=p_{\pi_1}e_{\pi_1}(x_1)
\prod_{i=2}^{n}
a_{\pi_{i-1}\pi_i}e_{\pi_i}(x_i).
$$ {#eq-hmm-joint}

::: {.proof}
Sketch. The machine first follows a path and then prints symbols along it, so
$$
P(x,\pi)=P(\pi)\,P(x\mid\pi).
$$
The Markov assumption makes each state depend only on the previous one, and the emission assumption makes each symbol depend only on its own state:
$$
P(\pi)=p_{\pi_1}\prod_{i=2}^{n}a_{\pi_{i-1}\pi_i},
\qquad
P(x\mid\pi)=\prod_{i=1}^{n}e_{\pi_i}(x_i).
$$
Multiplying the two products and grouping factors by position gives @eq-hmm-joint.
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

# initial, transition, and emission hold the B/I parameters set above.
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
The denominator is shared by every path, so the most probable joint path is also the most probable path conditional on the observations. But choosing it still discards the alternatives visible in @fig-hidden-markov-model-path-mass.

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
The path contains transitions $II,IB,BB,BB$, giving $\widehat a_{II}=\widehat a_{IB}=1/2$ and $\widehat a_{BB}=1$.
:::

The exact zeros and ones reflect a tiny sample, not biological certainty, and pseudocounts keep such estimates from ruling out unobserved events. The annotation source also defines what the model learns to call an island: labels drawn by a sequence-composition rule teach that rule and by themselves say nothing about gene regulation.

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
