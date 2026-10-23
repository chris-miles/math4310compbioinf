---
published: true
title: "Affine Gaps and Seed-and-Extend"
subtitle: "Lesson 5 · 2027-01-27 · Week 3"
nocite: |
  @durbin1998, @gotoh1982, @altschul1990, @compeau2015
---

## Core question

How should we score a missing stretch and search many sequences for a match?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import numpy as np

import coursefigs as cf

cf.use_style()
```

## Gap opening and extension

Suppose sequencing and gene prediction (identifying DNA regions likely to encode proteins) produce a candidate protein whose function is unknown. Its amino-acid string can be searched against a database of recorded proteins to find shared regions and candidate functions. We need a scoring rule for each proposed correspondence, and a search strategy to find useful candidates among many records.

First consider the scoring rule. A stretch of missing bases or amino acids may reflect one insertion or deletion event [@durbin1998]. The linear gap score from earlier lessons charges by the number of missing letters alone: one length-$\ell$ gap and $\ell$ separate one-letter gaps both score $\ell g$. It cannot express a preference for one longer event over several separate events. An *affine* gap score gives starting a gap a different cost from continuing it.

Once that score is defined, the database creates a second problem. We do not know which record contains a useful match, and filling a complete local-alignment table against every record repeats expensive work. A *seed* is a short matching word that selects a location for closer comparison. The gap model determines how we rank candidate alignments; the seed rule determines which candidates reach that comparison.

::: {#def-affine-gap}
## Affine gap penalty
For penalties $d_o>0$ and $d_e>0$, the cost of a gap of length $\ell\ge1$ is
$$
G(\ell)=d_o+(\ell-1)d_e.
$$
Here $d_o$ is the gap-open penalty and $d_e$ is the gap-extension penalty. Usually $d_o>d_e$.
:::

We subtract this positive cost from an alignment score, just as a negative gap score $g$ was added in earlier lessons. With $d_o=4$ and $d_e=1$, gaps of lengths 1, 2, and 3 cost 4, 5, and 6, while two separate one-letter gaps cost 8 because each represents its own event and pays its own opening.

The difference is visible even in the toy strings AAAA and AA. Both correspondences in @fig-affine-gap-runs have two matches and two gap columns. Linear scoring ties them; affine scoring charges 5 for the single two-letter gap and 8 for the two separate gaps. The repeated A bases leave the correspondence ambiguous, so this preference comes from the gap model rather than from additional letter evidence.

```{.python .execute}
#| label: fig-affine-gap-runs
#| code-fold: true
#| code-summary: "Compare one gap run with two gap runs"
#| fig-cap: "Two alignments of AAAA with AA use the same number of gap columns, in one run (left) or two runs (right). With opening penalty 4 and extension penalty 1, their gap costs are 5 and 8, expressing a preference for fewer separate events."
#| fig-alt: "AAAA over A--A on the left, and AAAA over -A-A on the right. Both pairs have two A matches, but their gaps form one run and two runs respectively."
#| fig-width: 6
#| fig-height: 1.8
fig, axes = plt.subplots(1, 2)
cf.alignment(axes[0], "AAAA", "A--A")
cf.alignment(axes[1], "AAAA", "-A-A")
plt.show()
```

The cost of the next column now depends on the preceding column. Appending `C/-` after a paired-letter column costs $d_o$, while appending it after another letter/gap column costs only $d_e$. Two prefix alignments can therefore have different continuation costs even when they end at the same $(i,j)$. Keeping only their overall maximum can discard a lower-scoring prefix whose next extension is cheaper. We retain a best score for each possible ending type [@durbin1998].

## Three alignment states

::: {#def-affine-states}
## Affine alignment states
For prefixes ending at $(i,j)$, let

- $M(i,j)$ be the best score ending with $x_i$ aligned to $y_j$;
- $I_x(i,j)$ be the best score ending with $x_i$ aligned to a gap;
- $I_y(i,j)$ be the best score ending with a gap aligned to $y_j$.
:::

The recurrences are
$$
\begin{aligned}
M(i,j)&=s(x_i,y_j)+\max\{M,I_x,I_y\}(i-1,j-1),\\
I_x(i,j)&=\max\{M(i-1,j)-d_o,\ I_x(i-1,j)-d_e\},\\
I_y(i,j)&=\max\{M(i,j-1)-d_o,\ I_y(i,j-1)-d_e\}.
\end{aligned}
$$ {#eq-affine-recurrence}

The first candidate in a gap state opens a gap; the second extends one. Entering $I_x$ consumes the next letter of $x$ but none of $y$, so its predecessors lie at $(i-1,j)$. Entering $I_y$ consumes only a letter of $y$, and entering $M$ consumes one from each. These index changes are the vertical, horizontal, and diagonal moves from Lesson 3, with the ending type now recorded as well.

This version forbids immediate switches between $I_x$ and $I_y$; other formulations allow them, and the convention determines which alignments are admissible. The name $M$ denotes any paired-letter column, match or mismatch.

```{.python .execute}
#| label: fig-affine-state-transitions
#| fig-cap: "Allowed ending-state transitions for affine alignment. The highlighted loop continues a gap for extension cost rather than paying another opening cost."
#| fig-alt: "Three states M, I sub x, and I sub y, each with a self-loop. M connects in both directions to each gap state; there is no edge between the gap states."
#| fig-width: 6
#| fig-height: 4
paired_state, x_gap_state, y_gap_state = "$M$", "$I_x$", "$I_y$"
state_edges = [
    (paired_state, paired_state),
    (paired_state, x_gap_state), (x_gap_state, paired_state),
    (paired_state, y_gap_state), (y_gap_state, paired_state),
    (x_gap_state, x_gap_state), (y_gap_state, y_gap_state),
]
fig, ax = plt.subplots()
cf.digraph(
    ax, state_edges,
    pos={paired_state: (0, 0), x_gap_state: (-1.5, 1.5),
         y_gap_state: (1.5, 1.5)},
    highlight_edges=[(x_gap_state, x_gap_state)],
)
plt.show()
```

For global alignment, set $M(0,0)=0$,
$$
I_x(i,0)=-d_o-(i-1)d_e,\qquad I_y(0,j)=-d_o-(j-1)d_e,
$$
for $i,j\geq1$, and set every other border state except $M(0,0)$ to $-\infty$. The final global score is $\max\{M(n,m),I_x(n,m),I_y(n,m)\}$.

::: {#exm-affine-gap}
## A worse prefix with a better continuation
Use $x=\texttt{ACG}$, $y=\texttt{C}$, match $+1$, mismatch $-1$, $d_o=4$, and $d_e=1$. At $(2,1)$, two endings compete:
$$
\begin{array}{c|c|c}
\text{alignment of AC and C} & \text{ending state} & \text{score}\\
\texttt{AC}/\texttt{-C} & M & -4+1=-3\\
\texttt{AC}/\texttt{C-} & I_x & -1-4=-5
\end{array}
$$
The paired-letter prefix wins by two points. But appending G/- opens a new gap after the first alignment and extends the existing gap after the second:
$$
I_x(3,1)=\max\{M(2,1)-4,\ I_x(2,1)-1\}
=\max\{-7,-6\}=-6.
$$
Discarding the lower prefix score $-5$ would lose this better continuation. Scoring the two complete alignments under both gap models shows that the choice of model changes the answer:
$$
\begin{array}{c|c|c}
\text{alignment} & \text{linear, } g=-4 & \text{affine, } d_o=4,\ d_e=1\\
\hline
\texttt{ACG}/\texttt{-C-} & 1-4-4=-7 & 1-4-4=-7\\
\texttt{ACG}/\texttt{C--} & -1-4-4=-9 & -1-5=-6
\end{array}
$$
Linear gaps prefer ACG/-C-, which keeps the C/C match and pays for two separate one-letter gaps. Affine gaps prefer ACG/C--, accepting a mismatch so that the two missing letters form one run and pay one opening.

The best affine score is $-6$; ACG/--C also attains it. That tie leaves the placement ambiguous, but both optima use one opening rather than two. The model encodes a preference for one insertion or deletion event; it does not establish which event occurred.
:::

::: {.callout-note title="Key idea"}
With affine gaps, the cost of the next column depends on how the alignment so far ends. Each cell therefore keeps one best score per ending state, because a prefix that scores lower now can still lead to the best complete alignment.
:::

::: {#prp-affine-correctness}
## Correctness of the three-state recurrence
For every state and cell, the recurrence in @eq-affine-recurrence gives the best score among alignments of the two prefixes that end in that state.
:::

::: {.proof}
Induct on $i+j$. An alignment ending in $M$ removes a paired-letter column and may come from any state. An alignment ending in $I_x$ either opens that gap after an $M$ column or extends an existing $I_x$ gap; the $I_y$ case is symmetric. Replacing the remaining prefix by a better prefix ending in the same state preserves the next opening or extension cost and would improve the whole alignment. Thus an optimal alignment uses an optimal same-state prefix. These cases are exhaustive under the stated transition convention, and appending the indicated column realizes each candidate.
:::

Each of the three tables has $(n+1)(m+1)$ entries. An interior entry compares at most three predecessors, so the time and memory costs remain $\mathcal O(nm)$. Keeping gap states avoids trying every possible gap length anew at each cell [@durbin1998; @compeau2015].

::: {.algorithm}
**Algorithm: score a global alignment with affine gaps.**

1. Create three tables for paired-letter, letter/gap, and gap/letter endings; initialize impossible states to $-\infty$.
2. Set the origin and gap-only borders using the stated global boundary rules.
3. In row order, fill all three states at each interior cell using @eq-affine-recurrence.
4. Return the largest state score at $(n,m)$.
:::

Here `gap_in_y` implements $I_x$, because a letter of $x$ faces a gap in $y$; `gap_in_x` implements $I_y$.

```{.python .execute}
#| label: code-affine-global-score
def affine_global_score(
    x: str, y: str, match: int = 1, mismatch: int = -1,
    gap_open: int = 4, gap_extend: int = 1,
) -> float:
    """Return a global score with no direct switches between opposite gaps."""
    n, m = len(x), len(y)
    paired = [[float("-inf")] * (m + 1) for _ in range(n + 1)]
    gap_in_y = [[float("-inf")] * (m + 1) for _ in range(n + 1)]
    gap_in_x = [[float("-inf")] * (m + 1) for _ in range(n + 1)]
    paired[0][0] = 0
    for i in range(1, n + 1):
        gap_in_y[i][0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        gap_in_x[0][j] = -gap_open - (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            substitution = match if x[i - 1] == y[j - 1] else mismatch
            paired[i][j] = substitution + max(
                paired[i - 1][j - 1],
                gap_in_y[i - 1][j - 1],
                gap_in_x[i - 1][j - 1],
            )
            gap_in_y[i][j] = max(
                paired[i - 1][j] - gap_open,
                gap_in_y[i - 1][j] - gap_extend,
            )
            gap_in_x[i][j] = max(
                paired[i][j - 1] - gap_open,
                gap_in_x[i][j - 1] - gap_extend,
            )
    return max(paired[n][m], gap_in_y[n][m], gap_in_x[n][m])


affine_global_score("ACG", "C"), affine_global_score("ACG", "C", gap_extend=4)
```

The outputs are `(-6, -7)`. The affine score agrees with @exm-affine-gap, and setting opening and extension costs both to 4 recovers the best linear score. Impossible border states remain $-\infty$ so that they can never win a maximum.

The function returns a score. Recovering the alignment would require traceback that records both a predecessor cell and its state. A pointer to $(i-1,j)$ alone cannot tell whether that step opened a gap from $M$ or extended it from $I_x$.

## Database search and seeds

If a query has length $n$ and a database contains total length $N$, exact dynamic programming takes $\mathcal O(nN)$ time. For one query this may already mean billions of cell updates. Database search tools first look for short promising matches and run expensive alignment only near those locations.

::: {#def-seed}
## Seed
A *seed* is a short pattern from the query used to find candidate positions in a target. An exact $k$-mer seed requires the same length-$k$ word in both sequences.
:::

::: {.algorithm}
**Algorithm: seed and extend.**

1. Build an index from target $k$-mers to their positions.
2. Look up each query $k$-mer in the index.
3. Extend each hit along its ungapped offset by adding pair scores in both directions. Keep the best endpoints; a common stopping rule quits after the running score falls a chosen amount below its best value.
4. For promising hits, allow gaps by filling a local-alignment table in a window or band around the hit.
5. Return the scored candidate alignments.
:::

Index lookup avoids considering most query-target position pairs. For query `ACGTAC` and target `AGGTTC`, the 3-mer call should return no hits. The 2-mer call should return one hit, `[(2, 2)]`, because GT begins at zero-based position 2 in both strings.

```{.python .execute}
#| label: code-exact-seeds
def exact_seed_hits(query: str, target: str, k: int) -> list[tuple[int, int]]:
    """Return positions of all shared exact length-k words."""
    if k < 1:
        raise ValueError("k must be positive.")
    index = {}
    for target_start in range(len(target) - k + 1):
        word = target[target_start:target_start + k]
        index.setdefault(word, []).append(target_start)
    hits = []
    for query_start in range(len(query) - k + 1):
        word = query[query_start:query_start + k]
        for target_start in index.get(word, []):
            hits.append((query_start, target_start))
    return hits

exact_seed_hits("ACGTAC", "AGGTTC", 3), exact_seed_hits("ACGTAC", "AGGTTC", 2)
```

The function returns pairs of zero-based start positions `(query_start, target_start)`, not alignments. The difference `target_start - query_start` identifies a candidate ungapped offset, along which extension is cheap because the base pairing is fixed. An insertion or deletion changes the offset, so gapped extension must be able to leave that diagonal; a window or band limits the cells examined, and an alignment that leaves it can still be missed.

Exact DNA words keep the search rule easy to inspect. The original BLAST algorithm also admits high-scoring protein words that need not match exactly, then extends its word hits in the same way [@altschul1990; @durbin1998].

## Seed sensitivity

The ungapped ACGTAC/AGGTTC comparison has four matches and two mismatches, scoring 2 under $+1/-1$. In @fig-affine-seed-interruptions, the mismatches split those matches into runs of lengths 1, 2, and 1, so no 3-mer is shared at this offset (and checking all the words shows none at any other). The 3-mer index finds no candidate while the 2-mer index finds GT: a promising final score cannot help if the seed stage never sends the location to extension.

```{.python .execute}
#| label: fig-affine-seed-interruptions
#| code-fold: true
#| code-summary: "Show how mismatches interrupt an exact seed"
#| fig-cap: "Four paired bases match in ACGTAC/AGGTTC, but the two mismatch dots split them into runs of lengths 1, 2, and 1. The exact 3-mer search misses this pair; a 2-mer search finds the central GT run."
#| fig-alt: "ACGTAC above AGGTTC, with match marks at positions one, three, four, and six and red mismatch dots at positions two and five."
#| fig-width: 5
#| fig-height: 1.7
fig, ax = plt.subplots()
cf.alignment(ax, "ACGTAC", "AGGTTC", labels=("query", "target"))
plt.show()
```

::: {#def-seed-sensitivity}
## Seed sensitivity
For a specified class of true alignments, seed sensitivity is the fraction that contain at least one seed accepted by the search rule.
:::

For a tiny sensitivity check, take ACGTAC as the query and suppose the true targets are ACGTTC, ACATAC, AGGTTC, and ATGTTC. The first two contain a shared 3-mer; the last two do not. This class has exact 3-mer sensitivity $2/4$, while all four contain a shared 2-mer. The fraction refers to this specified class, rather than a universal property of $k$.

Shorter seeds improve sensitivity, but they also produce more chance hits, and every hit costs an extension. We can count the chance hits under a simple background. If both sequences are independent, uniformly random DNA, two length-$k$ words agree with probability $4^{-k}$. A query of length $n$ and a target of length $N$, both at least $k$, offer $(n-k+1)(N-k+1)$ pairs of starting positions. Let $h$ count the exact seed hits. Adding the agreement probability over all position pairs gives
$$
\mathbb E[h]=(n-k+1)(N-k+1)\,4^{-k}.
$$
The first two factors count the comparisons the index makes implicitly; the last is the chance that any one of them succeeds. Overlapping words make hit events dependent, but an expected count still adds over events. Lengthening the seed by one letter multiplies the expected number of chance hits by $4^{-(k+1)}/4^{-k}=1/4$, while the number of position pairs barely changes. For a length-100 query against a length-1,000 target, $k=4$ gives $97\cdot997/256\approx378$ expected chance hits and $k=8$ gives $93\cdot993/65536\approx1.4$.

```{.python .execute}
#| label: fig-affine-seed-dotplot
#| code-fold: true
#| code-summary: "Plot seed hits for a query that shares one region with a random target"
#| fig-cap: "Exact seed hits for a length-100 query against a random length-1,000 target; the query is a copy of target positions 401 to 500 with about 20% of its bases substituted. Red hits lie on the true offset and grey hits arise by chance. With k = 4 the true diagonal is surrounded by about as many chance hits as the formula predicts; with k = 8 almost all chance hits disappear, and so do the true seeds in stretches where substitutions leave no 8-letter match."
#| fig-alt: "Two dot plots of query position against target position. The left panel, k equals 4, has several hundred grey points scattered everywhere and a red near-vertical line of points near target position 400 to 500. The right panel, k equals 8, has one grey point and a broken red line covering only the first 60 query positions."
#| fig-height: 3.5
rng = np.random.default_rng(7)
target = "".join(rng.choice(list("ACGT"), size=1000))
true_offset = 400
# The query copies 100 target bases, substituting each with probability 0.2.
query_letters = []
for base in target[true_offset:true_offset + 100]:
    if rng.random() < 0.8:
        query_letters.append(base)
    else:
        query_letters.append(rng.choice([b for b in "ACGT" if b != base]))
query = "".join(query_letters)

fig, axes = plt.subplots(1, 2, sharey=True)
for ax, k in zip(axes, (4, 8)):
    hits = exact_seed_hits(query, target, k)
    true_hits = [(q, t) for q, t in hits if t - q == true_offset]
    chance_hits = [(q, t) for q, t in hits if t - q != true_offset]
    ax.plot([t for q, t in chance_hits], [q for q, t in chance_hits], "o",
            markersize=2.5, color=cf.RULE, label="chance offset")
    ax.plot([t for q, t in true_hits], [q for q, t in true_hits], "o",
            markersize=2.5, color=cf.EMPHASIS, label="true offset")
    ax.set(title=f"$k={k}$: {len(chance_hits)} chance, {len(true_hits)} true",
           xlabel="Target position", xlim=(0, 1000), ylim=(0, 100))
axes[0].set(ylabel="Query position")
axes[0].invert_yaxis()
axes[1].legend(loc="lower right", markerscale=2)
plt.show()
```

@fig-affine-seed-dotplot runs `exact_seed_hits` on such a pair, with one planted region of 80% identity. At $k=4$ the 376 chance hits agree with the predicted 378, since off the true diagonal the two sequences are still unrelated, and each would send the search to an extension that finds nothing. At $k=8$ one chance hit remains, but the true diagonal has lost its seeds wherever substitutions fall closer together than eight letters. A region with slightly lower identity could lose all of them, and an alignment containing no accepted seed is never proposed; unlike full dynamic programming, the search can miss the best alignment. Spaced seeds, which require matches only at selected positions within a longer window, are one way to keep the low chance-hit rate of a long seed while tolerating some substitutions.

:::: {.callout-note collapse="true" title="Going deeper: a bounded-mismatch guarantee"}

The expected-hit calculation addresses workload. For sensitivity, a different question is whether mismatches can destroy every accepted word.

::: {#prp-seed-guarantee}
## An exact seed survives bounded mismatches
An ungapped alignment of length $L$ with $t<L$ mismatches contains a matching run of length at least
$$
\left\lceil\frac{L-t}{t+1}\right\rceil.
$$
It therefore contains an exact $k$-mer seed for every $k$ no larger than that bound.
:::

::: {.proof}
The $t$ mismatches divide the $L-t$ matching positions into at most $t+1$ runs, including possibly empty runs at the ends. If every run were shorter than their average $(L-t)/(t+1)$, their total would be less than $L-t$. At least one run reaches that average; its integer length reaches the ceiling.
:::

For $L=20$ and $t=3$, some run has at least five matches, guaranteeing a 5-mer seed. An 8-mer seed has no such guarantee. Gaps need different accounting because they interrupt the correspondence between positions.

::::

For fixed $k$, building the index takes $\mathcal O(N)$ expected time with a hash table, and looking up the query's words and listing $h$ hits takes $\mathcal O(n+h)$ before any extension. Repeated words can make $h$ as large as $\mathcal O(nN)$, so the index speeds up typical searches but not every input.

## Limitations

1. **Affine gaps still use a simple duration model.** Every extension has the same cost, regardless of gap length or position.
2. **Seeds can miss real alignments.** Distributed mismatches or gaps can destroy every accepted word.
3. **Repeats flood the candidate list.** Common seeds occur at many positions and erase the speed advantage unless the method filters or downweights them.
4. **Statistical significance depends on the search.** The same raw score can be less surprising in a larger database or in repetitive sequence.

## Exercises

::: {#exr-affine-cost}
With $d_o=5$ and $d_e=1$, compare the cost of one length-4 gap with four separate length-1 gaps. State the event interpretation.
:::

::: {#exr-affine-cell}
For the ACG/C example, compute $M(3,1)$ from its predecessor states. Explain why it ties $I_x(3,1)$ even though the two alignments pair the target C with different query letters.
:::

::: {#exr-affine-linear}
Set $d_o=d_e$. Explain why the score of every gap becomes linear in its length, even though the algorithm still uses three tables. Does changing the penalties also remove the restriction on immediate switches between opposite gap states?
:::

::: {#exr-seed-counterexample}
Construct two length-8 strings with a positive ungapped score under match $+1$, mismatch $-1$, but no shared exact 3-mer. Mark the mismatches.
:::

::: {#exr-seed-choice}
A target database is highly repetitive. Predict how decreasing $k$ affects the number of seed hits, runtime, and sensitivity.
:::

::: {#refs}
**References**
:::
