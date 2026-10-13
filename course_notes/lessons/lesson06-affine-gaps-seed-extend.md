---
published: true
title: "Affine Gaps and Seed-and-Extend"
subtitle: "Lesson 6 · 2027-02-01 · Week 4"
nocite: |
  @durbin1998, @gotoh1982, @altschul1990, @compeau2015
---

## Core question

How do we model gaps as events, then search a database without filling every table?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt

import coursefigs as cf

cf.use_style()
```

## Gap opening and extension

Suppose sequencing and gene prediction (identifying DNA regions likely to encode proteins) produce a candidate protein whose function is unknown. A database search compares its amino-acid sequence with known records and returns regions worth investigating. Unlike the two-sequence problems so far, we do not know in advance which record contains the useful match. Filling a complete local-alignment table against every record repeats expensive work; seeds first select promising locations.

Within a candidate alignment, a stretch of missing residues may reflect one insertion or deletion event. An affine gap score distinguishes that event from several isolated gaps. These are separate choices: the gap model ranks the alignments we examine, while the seed rule determines which locations we examine at all.

A linear gap score $\ell g$ charges the same amount for one length-$\ell$ gap as for $\ell$ separate one-letter gaps. A single insertion or deletion often affects a run of adjacent bases or amino acids. We therefore distinguish starting a gap from continuing it.

::: {#def-affine-gap}
## Affine gap penalty
For penalties $d_o>0$ and $d_e>0$, the cost of a gap of length $\ell\ge1$ is
$$
G(\ell)=d_o+(\ell-1)d_e.
$$
Here $d_o$ is the gap-open penalty and $d_e$ is the gap-extension penalty. Usually $d_o>d_e$.
:::

We subtract this positive cost from an alignment score, just as a negative gap score $g$ was added in earlier lessons. With $d_o=4$ and $d_e=1$, gaps of lengths 1, 2, and 3 cost 4, 5, and 6. Two separate one-letter gaps cost 8. The second opening represents a second event.

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

This version forbids immediate switches between $I_x$ and $I_y$ in either direction. Other formulations allow them; the transition convention determines which alignments are admissible. The name $M$ denotes a paired-letter column, which can be either a match or a mismatch.

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
Discarding the lower prefix score $-5$ would lose this better continuation. The two resulting alignments are ACG/-C- (score $-7$) and ACG/C-- (score $-6$). A linear cost of 4 per gap character instead ranks them $-7$ and $-9$, favoring the isolated gaps. The affine model accepts a mismatch to group the missing bases into one run.

The best affine score is $-6$; ACG/--C also attains it. That tie leaves the placement ambiguous, but both optima use one opening rather than two. The model encodes a preference for one insertion or deletion event; it does not establish which event occurred.
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

The outputs are `(-6, -7)`. The affine score agrees with @exm-affine-gap; setting opening and extension costs both to 4 gives the best linear score. Impossible border states remain $-\infty$ so they cannot win a maximum. The same initialization gives score 0 for two empty strings and $-d_o-(\ell-1)d_e$ for a length-$\ell$ string against an empty one.

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
5. Report high-scoring alignments and their search-calibrated significance.
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

The function returns pairs of zero-based start positions `(query_start, target_start)`, not alignments. The difference `target_start - query_start` identifies a candidate ungapped offset. Extending along it is cheap because the base pairing is fixed. An insertion or deletion changes the offset, so a later gapped extension must be able to move away from that diagonal. A window or band limits the cells examined; an alignment leaving that region can still be missed.

This example uses DNA exact words to make the search rule inspectable. The original BLAST algorithm also admits high-scoring protein words that need not match exactly, then extends its word hits [@altschul1990; @durbin1998]. Exact seeding is a simpler member of the same search architecture.

## Seed sensitivity

The ungapped ACGTAC/AGGTTC comparison has four matches and two mismatches, scoring 2 under $+1/-1$. Yet the 3-mer index finds no candidate, while the 2-mer index finds GT. A promising final score cannot help if the seed stage never sends the location to extension.

::: {#def-seed-sensitivity}
## Seed sensitivity
For a specified class of true alignments, seed sensitivity is the fraction that contain at least one seed accepted by the search rule.
:::

For a tiny sensitivity check, take ACGTAC as the query and suppose the true targets are ACGTTC, ACATAC, AGGTTC, and ATGTTC. The first two contain a shared 3-mer; the last two do not. This class has exact 3-mer sensitivity $2/4$, while all four contain a shared 2-mer. The fraction refers to this specified class, rather than a universal property of $k$.

Longer exact seeds generate fewer random hits and speed up extension, but they miss more alignments containing substitutions or gaps. Shorter seeds improve sensitivity and produce more candidates. Spaced seeds require matches at selected positions rather than every position in a contiguous word; they can detect some alignments that contiguous seeds miss.

The index finds every accepted shared word when it retains every occurrence. It cannot propose an alignment containing no accepted seed. Unlike full dynamic programming, this search rule can therefore miss the unrestricted best alignment.

### Seed length and candidate counts

Under an iid uniform DNA background, two independent length-$k$ words agree with probability $4^{-k}$. A query of length $n$ and target of length $N$, both at least $k$, offer $(n-k+1)(N-k+1)$ pairs of starting positions. If $h$ counts exact seed hits, adding the probabilities of these position pairs gives
$$
\mathbb E[h]=(n-k+1)(N-k+1)4^{-k}.
$$
Overlapping words make hit events dependent, but expectation still adds. Independence of the two sequences and the uniform background are the assumptions behind each event's probability.

For a length-100 query and a length-10,000 target, $k=4$ gives an expected $97(9997)/256\approx3788$ random hits. At $k=8$, the expectation is $93(9993)/65536\approx14.2$. These are counts of candidate position pairs, not counts of distinct matching records or final reported alignments.

The random-hit calculation addresses workload. For sensitivity, a different question is whether mismatches can destroy every accepted word.

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

## Database size and E-values

A raw local-alignment score is less surprising in a large search than in a small one. BLAST reports an E-value: the expected number of chance hits scoring at least this well in a search of the stated size under its background model. An E-value of $0.01$ means one such chance hit per 100 comparable searches on average. It is not the probability that the reported sequences are unrelated.

For fixed $k$, building the index costs $\mathcal O(N)$ expected time with a hash table. Looking up query words and listing $h$ hits costs $\mathcal O(n+h)$ before extension. Repeated words can make $h$ as large as $\mathcal O(nN)$, so indexing does not guarantee a speedup on every input.

Database composition, low-complexity sequence, score parameters, and database size all affect significance calibration. The extension score and the seed rule answer different questions: a sensitive score cannot recover a target that the seed stage never visits.

## Limitations

1. **Affine gaps still use a simple duration model.** Every extension has the same cost, regardless of gap length or position.
2. **State conventions differ.** Allowing direct transitions between opposite gap states changes the set of represented alignments and can change scores under unusual parameter choices.
3. **Seeds can miss real alignments.** Distributed mismatches or gaps can destroy every accepted word.
4. **Repeats flood the candidate list.** Common seeds occur at many positions and erase the speed advantage unless the method filters or downweights them.
5. **E-values depend on the search.** Changing the database size or background composition changes significance even when the raw alignment score stays fixed.

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
