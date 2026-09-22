---
published: true
title: "Affine Gaps and Seed-and-Extend"
subtitle: "Lesson 6 · Week 3"
---

## Core question

How do we model gaps as events, then search a database without filling every table?

## Gap opening and extension

A linear gap score $kg$ charges the same amount for one length-$k$ gap as for $k$ separate one-letter gaps. A single insertion or deletion often affects a run of adjacent bases or amino acids. We therefore distinguish starting a gap from continuing it.

::: {#def-affine-gap}
## Affine gap penalty
For penalties $d>0$ and $e>0$, the cost of a gap of length $k\ge1$ is
$$
G(k)=d+(k-1)e.
$$
Here $d$ is the gap-open penalty and $e$ is the gap-extension penalty. Usually $d>e$.
:::

With $d=2$ and $e=1$, gaps of lengths 1, 2, and 3 cost 2, 3, and 4. Two separate one-letter gaps cost 4. The model treats opening a second gap as a second event.

The recurrence must remember whether its last column is already inside a gap. One table entry can no longer summarize all relevant history.

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
I_x(i,j)&=\max\{M(i-1,j)-d,\ I_x(i-1,j)-e\},\\
I_y(i,j)&=\max\{M(i,j-1)-d,\ I_y(i,j-1)-e\}.
\end{aligned}
$$ {#eq-affine-recurrence}

The first candidate in a gap state opens a gap; the second extends one. This version forbids immediate switches between $I_x$ and $I_y$ in either direction. Other formulations allow them; the scoring convention determines which alignments are admissible.

For global alignment, set $M(0,0)=0$,
$$
I_x(i,0)=-d-(i-1)e,\qquad I_y(0,j)=-d-(j-1)e,
$$
for $i,j\geq1$, and set every other border state except $M(0,0)$ to $-\infty$. The final global score is $\max\{M(n,m),I_x(n,m),I_y(n,m)\}$.

::: {#exm-affine-gap}
## Opening and extending
Take $x=\texttt{ACGT}$, $y=\texttt{AT}$, match $+1$, mismatch $-1$, $d=2$, and $e=1$. After matching A/A,
$$
I_x(2,1)=M(1,1)-2=-1.
$$
At the next letter,
$$
I_x(3,1)=\max\{M(2,1)-2, I_x(2,1)-1\}=-2.
$$
The extension candidate wins. The alignment
$$
\begin{array}{c}\texttt{ACGT}\\\texttt{A--T}\end{array}
$$
scores $1-3+1=-1$. A linear penalty of 2 per gap character gives $1-4+1=-2$.
:::

::: {#prp-affine-correctness}
## Correctness of the three-state recurrence
For every state and cell, the recurrence in @eq-affine-recurrence gives the best score among alignments of the two prefixes that end in that state.
:::

::: {.proof}
Induct on $i+j$. An alignment ending in $M$ removes a paired-letter column and may come from any state. An alignment ending in $I_x$ either opens that gap after an $M$ column or extends an existing $I_x$ gap; the $I_y$ case is symmetric. These cases are exhaustive under the stated transition convention, and appending the indicated column realizes each candidate.
:::

Three tables multiply the constant amount of work, but the asymptotic cost remains $\mathcal O(nm)$ time and memory.

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
3. Extend each hit in both directions, first cheaply and then with local alignment for promising hits.
4. Report high-scoring alignments and their search-calibrated significance.
:::

Index lookup avoids considering most query-target position pairs.

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

The original BLAST algorithm follows this architecture, using short word hits followed by local extension [@altschul1990].

## Seed sensitivity

Let the query be $\texttt{ACGTAC}$ and the target be $\texttt{AGGTTC}$. Their ungapped alignment has four matches and two mismatches:
$$
\begin{array}{c}
\texttt{ACGTAC}\\
\texttt{AGGTTC}
\end{array}
$$
With match $+1$ and mismatch $-1$, its score is 2.

The query 3-mers are ACG, CGT, GTA, TAC. The target 3-mers are AGG, GGT, GTT, TTC. There is no shared 3-mer, so an exact 3-mer seeding rule never sends this pair to extension. With $k=2$, both contain GT and the candidate is examined.

::: {#def-seed-sensitivity}
## Seed sensitivity
For a specified class of true alignments, seed sensitivity is the fraction that contain at least one seed accepted by the search rule.
:::

Longer exact seeds generate fewer random hits and speed up extension, but they miss more alignments containing substitutions or gaps. Shorter seeds improve sensitivity and produce more candidates. Spaced seeds require matches at selected positions rather than every position in a contiguous word; they can detect some alignments that contiguous seeds miss.

::: {#prp-seed-guarantee}
## What exact seeding guarantees
An exact $k$-mer index discovers every occurrence of an exact shared $k$-mer, provided indexing and lookup enumerate every occurrence. It gives no candidate location for an alignment that contains no accepted seed.
:::

::: {.proof}
Every exact shared $k$-mer occurs as a query word and as an indexed target word, so lookup returns that pair of positions. What happens after the hit depends on the extension rule. If no accepted seed exists, lookup cannot propose that location.
:::

Unlike dynamic programming, seed and extend is a heuristic for the unrestricted best-alignment problem. Its limitation is explicit and testable: construct a true alignment with no accepted seed.

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
With $d=5$ and $e=1$, compare the cost of one length-4 gap with four separate length-1 gaps. State the event interpretation.
:::

::: {#exr-affine-cell}
For the ACGT/AT example, compute $M(4,2)$ from its three predecessor-state values and verify the score of the displayed alignment.
:::

::: {#exr-affine-linear}
Set $d=e$. Explain why the score of every gap becomes linear in its length, even though the algorithm still uses three tables.
:::

::: {#exr-seed-counterexample}
Construct two length-8 strings with a positive ungapped score under match $+1$, mismatch $-1$, but no shared exact 3-mer. Mark the mismatches.
:::

::: {#exr-seed-choice}
A target database is highly repetitive. Predict how decreasing $k$ affects the number of seed hits, runtime, and sensitivity.
:::

## Further reading

Durbin et al. derive linear and affine gap models and their dynamic programs in Chapter 2 [@durbin1998]. Gotoh gave the standard quadratic-time affine-gap algorithm [@gotoh1982]. Altschul et al. describe BLAST's local-search strategy and statistical reporting [@altschul1990]. Compeau and Pevzner motivate affine penalties as a model of gap events [@compeau2015].
