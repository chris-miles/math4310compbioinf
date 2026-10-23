---
published: true
title: "Eulerian Paths, Repeats, and Assembly Errors"
subtitle: "Lesson 7 · 2027-02-03 · Week 4"
nocite: |
  @compeau2015, @pevzner2001
---

## Core question

Do a genome's short words determine its full sequence?

```{.python .execute}
#| label: setup
#| code-fold: true
import sys
sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## Eulerian paths

Assembly starts with overlapping reads whose positions in the genome are unknown. Lesson 6 represented their short words in a de Bruijn graph: each word supplies an edge, and matching suffixes and prefixes allow edges to follow one another. For the ideal 3-mer multiset of `TACGACG`, those edges must be ordered to use TAC once, ACG twice, CGA once, and GAC once. Following that order spells a candidate DNA sequence.

Our first task is to find such an ordering without guessing at branches. Throughout, the multiset is the ideal genomic one, in which each genome position contributes exactly one $k$-mer occurrence, so a complete reconstruction uses each edge exactly once. This is an *Eulerian path* problem: simple degree conditions decide whether a solution exists, and an algorithm constructs one in linear time.

Finding one complete path answers the computational question but leaves a biological one: do the words force that sequence, or could another path spell a different genome? Repeats can leave the data compatible with several genomes, and with real reads, errors and coverage gaps add or remove edges. An assembler then reports contigs, stretches that every compatible genome shares.

For example, `ACGTCGACG` and `ACGACGTCG` contain exactly the same 3-mers with the same counts, so they have the same de Bruijn graph. @fig-euler-repeat-graph draws that graph twice, each time with a different loop from CG taken first, and spells the resulting genome beneath it. The two loops append TCG and ACG, and taking them in the other order swaps those letters in the genome. A route through this graph is a valid reconstruction, but the words alone cannot say which of the two genomes supplied them.

```{.python .execute}
#| label: fig-euler-repeat-graph
#| fig-cap: "One graph, two genomes. Each panel highlights the loop from CG that a complete traversal takes first and, beneath it, the string the traversal spells, with the three letters appended by the highlighted loop in red. Swapping the order of the two loops swaps TCG and ACG in the spelled genome while using the same edges."
#| fig-alt: "Two copies of the same directed multigraph, each above a nine-letter string. The upper graph highlights CG to GT to TC to CG above ACGTCGACG, with letters 4 to 6 (TCG) on red; the lower highlights CG to GA to AC to CG above ACGACGTCG, with letters 4 to 6 (ACG) on red."
#| fig-width: 6
#| fig-height: 6
#| code-fold: true
#| code-summary: "Compare two reconstructions from the same short words"
repeat_edges = [("AC", "CG"), ("AC", "CG"), ("CG", "GT"),
                ("GT", "TC"), ("TC", "CG"), ("CG", "GA"), ("GA", "AC")]
positions = {"AC": (0, 0), "CG": (1.6, 0), "GT": (3, .9),
             "TC": (4.2, .2), "GA": (2.5, -1.1)}
first_loops = [
    [("CG", "GT"), ("GT", "TC"), ("TC", "CG")],
    [("CG", "GA"), ("GA", "AC"), ("AC", "CG", 0)],
]
spellings = ["ACGTCGACG", "ACGACGTCG"]
titles = ["TCG loop first", "ACG loop first"]
fig, axes = plt.subplots(4, 1, gridspec_kw={"height_ratios": [3, 0.7, 3, 0.7]})
for panel in range(2):
    graph_ax, string_ax = axes[2 * panel], axes[2 * panel + 1]
    cf.digraph(graph_ax, repeat_edges, pos=positions,
               highlight_edges=first_loops[panel])
    graph_ax.set_title(titles[panel])
    # positions 4-6 of the genome (Python indices 3-5) come from the first loop
    cf.sequence_tracks(string_ax, {"spells": spellings[panel]},
                       highlighted_row="spells", highlighted_positions=[3, 4, 5],
                       show_positions=False)
plt.show()
```

Recall that $\operatorname{Comp}_k(x)$ lists the length-$k$ substrings of $x$, including repeated occurrences. In $D(K)$, each word in the multiset $K$ supplies an edge from its length-$(k-1)$ prefix to its length-$(k-1)$ suffix, identical words supply parallel edges, and following an edge appends the word's last base to the sequence spelled so far.

::: {#def-eulerian-trail}
## Eulerian path and cycle
An *Eulerian path* uses every edge exactly once. An *Eulerian cycle* is an Eulerian path whose first and last nodes agree.
:::

Nodes may repeat; the restriction falls on edges, because they represent the $k$-mer occurrences that a reconstruction must consume.

::: {#prp-assembly-eulerian}
## Reconstruction is an Eulerian path problem
A string $x$ with $\operatorname{Comp}_k(x)=K$ determines an Eulerian path in $D(K)$. Every Eulerian path in $D(K)$ spells a string with composition $K$.
:::

::: {.proof}
Sketch. Consecutive genomic $k$-mers overlap in $k-1$ symbols, so their edges meet in order. Conversely, traversing an edge appends its final symbol. Using every edge occurrence once gives exactly $K$.
:::

The proposition reduces exact string reconstruction to a graph traversal. We still need to know when such a traversal exists, how to find it, and whether a different traversal would spell a different string.

## Degree conditions

Write $d^+(v)$ and $d^-(v)$ for out-degree and in-degree, counting parallel edges.

::: {#thm-directed-euler-path}
## Eulerian path criterion
Suppose all nonzero-degree nodes of a finite directed multigraph are connected when directions are ignored. It has an Eulerian path exactly when either every node is balanced, $d^+(v)=d^-(v)$, in which case every Eulerian path is a cycle; or there are two nodes $s\neq t$ with
$$
d^+(v)-d^-(v)=\begin{cases}
\phantom{-}1, & v=s \text{ (the start)},\\
-1, & v=t \text{ (the end)},\\
\phantom{-}0, & \text{every other node } v.
\end{cases}
$$
:::

::: {.proof}
Sketch. Every visit to an internal node pairs an entering edge with a leaving edge. Only the start and end can have an unpaired edge, proving necessity.

For sufficiency, first make the graph balanced by adding $t\to s$ in the open-path case. Follow unused outgoing edges from a node $v$. At any other node, the walk has used one more incoming than outgoing edge, so balance guarantees an unused exit. It can get stuck only back at $v$, producing a closed walk.

Removing that closed walk leaves balanced unused degrees. If an unused edge already touches the closed walk, we have a place to continue. Otherwise, weak connectivity gives an undirected route from the walk to a remaining edge. The first edge leaving the visited nodes must be unused, so it supplies that contact. Balance of the unused edges then guarantees an unused outgoing edge at some visited node. Start there, make another closed walk, and insert it at that node. Each insertion uses at least one new edge, so the finite process consumes every edge. In the open-path case, cut the completed cycle at the added edge and delete it to obtain a path from $s$ to $t$.
:::

The proof explains why local degree counts suffice once connectivity is checked: every unfinished excursion has somewhere to rejoin the growing walk.

::: {#exm-euler-degrees}
## Degree check for `TACGACG`

| Node | $d^-$ | $d^+$ | difference |
|---|---:|---:|---:|
| `TA` | 0 | 1 | 1 |
| `AC` | 2 | 2 | 0 |
| `CG` | 2 | 1 | -1 |
| `GA` | 1 | 1 | 0 |

The path must start at `TA` and end at `CG`.
:::

### Interpreting degree differences

A genome walk contributes one entry and one exit each time it passes through an internal $(k-1)$-mer, and these cancel in the degree difference. Only the first node keeps an unpaired exit and the last an unpaired entry; for a circular genome, every node balances.

This accounting gives a fast rejection test. If two nodes each have one extra outgoing edge, no ordering of the given edges can produce one linear string. The failure could come from a missing $k$-mer, an extra erroneous one, or data from more than one genome, and the degrees alone do not say which.

Degrees also do not test connectivity. Two separate cycles, one on nodes $A,B$ and one on $C,D$, balance every node, yet a walk stays trapped in the cycle where it begins. The weak-connectivity condition in @thm-directed-euler-path rules out this case.

## Constructing a path

Hierholzer's algorithm follows the proof. It postpones adding a node to the answer until no unused outgoing edge remains.

::: {.algorithm}
**Algorithm: construct an Eulerian path.**

1. Start at the node with degree difference 1, or any node with an outgoing edge for a cycle.
2. If the current node has an unused edge, remove it and push its endpoint on a stack.
3. Otherwise, pop the node into the output.
4. Reverse the output and verify that it contains $|E|+1$ nodes.
:::

```{.python .execute}
#| label: code-euler-path
def eulerian_path(edges: list[tuple[str, str]]) -> list[str]:
    """Return an Eulerian node path, or raise ValueError."""
    if not edges:
        return []
    adjacency = {}
    indegree = {}
    outdegree = {}
    for source, target in edges:
        adjacency.setdefault(source, []).append(target)
        adjacency.setdefault(target, [])
        outdegree[source] = outdegree.get(source, 0) + 1
        outdegree.setdefault(target, 0)
        indegree[target] = indegree.get(target, 0) + 1
        indegree.setdefault(source, 0)

    starts = [v for v in adjacency if outdegree[v] - indegree[v] == 1]
    ends = [v for v in adjacency if indegree[v] - outdegree[v] == 1]
    bad = [v for v in adjacency if abs(outdegree[v] - indegree[v]) > 1]
    if bad or (len(starts), len(ends)) not in [(0, 0), (1, 1)]:
        raise ValueError("Degree conditions fail.")

    start = starts[0] if starts else edges[0][0]
    stack, reverse_path = [start], []
    while stack:
        node = stack[-1]
        if adjacency[node]:
            stack.append(adjacency[node].pop())
        else:
            reverse_path.append(stack.pop())
    path = reverse_path[::-1]
    if len(path) != len(edges) + 1:
        raise ValueError("Edges are disconnected.")
    return path

def spell_node_path(path: list[str]) -> str:
    """Spell a string from overlapping node labels."""
    return "" if not path else path[0] + "".join(v[-1] for v in path[1:])

edges = [("TA", "AC"), ("AC", "CG"), ("CG", "GA"),
         ("GA", "AC"), ("AC", "CG")]
path = eulerian_path(edges)
path, spell_node_path(path)
```

Each traversal pushes one endpoint, and each stack entry is popped once. The adjacency lists therefore give $\mathcal O(|V|+|E|)$ time and space. The output is reversed because a node is finalized only after all its outgoing excursions have been completed.

### Stack order and path order

Take edges $A\to B$, $B\to A$, and $A\to C$. The degree conditions require starting at $A$ and ending at $C$. Suppose the algorithm first takes $A\to C$, reaching a dead end while two edges remain. A greedy walk that commits this order as its final answer would fail.

The stack algorithm instead postpones the final ordering:

| Action | Stack | Output so far, before reversal |
|---|---|---|
| start | A | empty |
| take $A\to C$ | A, C | empty |
| pop C | A | C |
| take $A\to B\to A$ | A, B, A | C |
| pop until empty | empty | C, A, B, A |

Reversing the output gives $A,B,A,C$, which uses every edge once. The earlier dead end becomes the end of the final path. Popping a vertex records where its completed portion belongs after the remaining excursions have been inserted.

## Repeats and ambiguity

For `ACGTCGACG` at $k=3$,
$$K=\{\!\{\texttt{ACG},\texttt{ACG},\texttt{CGT},\texttt{GTC},\texttt{TCG},\texttt{CGA},\texttt{GAC}\}\!\}.$$
The node `CG` occurs three times in the genome. Between those visits lie two closed excursions: `CG → GT → TC → CG`, which appends `TCG`, and `CG → GA → AC → CG`, which appends `ACG`. The graph records both excursions but not their order. After the first edge `AC → CG`, taking the TCG excursion first spells `ACGTCGACG`, and taking the ACG excursion first spells `ACGACGTCG`, using the same edges (@fig-euler-repeat-graph).

::: {#exm-repeat-ambiguity}
## Two strings with one composition
Both `ACGTCGACG` and `ACGACGTCG` have composition $K$. Data consisting only of $K$ cannot favor either reconstruction.
:::

The swap works for any node that a complete walk visits three times. Write such a walk as
$$
P\;E_1\;E_2\;Q,
$$
where $P$ is the part up to the first visit of the node $v$, $E_1$ and $E_2$ are the closed excursions from $v$ back to $v$, and $Q$ is the rest. Because each excursion begins and ends at $v$, $P\,E_2\,E_1\,Q$ is also a walk using exactly the same edges. By @prp-assembly-eulerian it spells a string with the same composition, which differs from the first whenever $E_1$ and $E_2$ spell different letters. A $(k-1)$-mer occurring three times in the genome supplies such a node: in `ACGTCGACG`, CG occurs at positions 2, 5, and 8.

A $(k-1)$-mer that occurs only twice gives one excursion and nothing to reorder, but two such repeats can still cause ambiguity when they interleave. If the walk visits $u$, then $w$, then $u$, then $w$, write it as $P\,A\,B\,C\,Q$, where $A$ and $C$ run from $u$ to $w$ and $B$ runs from $w$ back to $u$. Then $P\,C\,B\,A\,Q$ uses the same edges, and it spells a different string whenever $A$ and $C$ spell different letters. `TACGACG` shows why that last condition matters. Its repeated 2-mers interleave, with AC at positions 2 and 5 and CG at 3 and 6, but both segments from AC to CG are the single word ACG, so the swap returns the same string. A triple repeat or two interleaved repeats whose swapped segments differ is therefore enough to make the genome impossible to recover from its $k$-mers.

### Existence and uniqueness

The degree theorem says whether at least one traversal exists; assembly also asks whether all traversals spell the same string. A branching node alone does not settle this, because one choice may get stuck before using every edge, or two choices may spell the same letters, as in `TACGACG`.

Once two strings produce the same $K$, every algorithm receiving only $K$ sees identical data, so no tie-breaking rule or cleverer traversal can recover which string generated it. Additional information must cross the repeat: longer k-mers, paired reads with a known separation, or a long read spanning both flanks.

At $k=4$, the two strings in @exm-repeat-ambiguity no longer have the same composition. For `ACGTCGACG`, the 4-mers are `ACGT, CGTC, GTCG, TCGA, CGAC, GACG`. For `ACGACGTCG`, they are `ACGA, CGAC, GACG, ACGT, CGTC, GTCG`. The first contains `TCGA` and the second `ACGA`, so one extra base of context separates the candidates. A longer repeat would need a correspondingly larger $k$.

## Errors, tips, and bubbles

A wrong measured base changes several overlapping k-mers at once. The resulting edges may look like a second route through a region, so an assembler must decide whether that route is a measurement artifact or a sequence actually present in the sample.

::: {#def-tip-bubble}
## Tips and bubbles
A *tip* is a short dead-end path. A *bubble* consists of two distinct paths that diverge from one node and later rejoin, with no shared internal nodes.
:::

::: {#exm-euler-error-bubble}
## One substitution creates a bubble
Suppose some reads report `ACGTC` and one reports `ACTTC`. The third base differs. At $k=3$, their k-mer lists and node paths are

| Read | 3-mers | Node path |
|---|---|---|
| `ACGTC` | `ACG, CGT, GTC` | `AC → CG → GT → TC` |
| `ACTTC` | `ACT, CTT, TTC` | `AC → CT → TT → TC` |

All three words in the second read contain the changed base, so the two paths share only their ends, `AC` and `TC`: the unchanged two-base flanks are enough for them to rejoin. If the changed base were at a read's end, the erroneous path could instead stop before reaching a shared flank.
:::

```{.python .execute}
#| label: fig-euler-error-bubble
#| fig-cap: "The reads ACGTC and ACTTC give two paths from AC to TC at k = 3; the ACTTC route is highlighted. A bubble locates conflicting sequence evidence but does not identify its cause."
#| fig-alt: "Two directed paths from AC at the left to TC at the right. The upper path passes through CG and GT; the lower highlighted path passes through CT and TT."
#| fig-width: 6
#| fig-height: 3
first_edges = [("AC", "CG"), ("CG", "GT"), ("GT", "TC")]
alternative_edges = [("AC", "CT"), ("CT", "TT"), ("TT", "TC")]
bubble_positions = {
    "AC": (0, 0), "CG": (1, 0.65), "GT": (2, 0.65),
    "TC": (3, 0), "CT": (1, -0.65), "TT": (2, -0.65),
}
fig, ax = plt.subplots()
cf.digraph(
    ax, first_edges + alternative_edges,
    pos=bubble_positions, highlight_edges=alternative_edges,
)
plt.show()
```

Suppose twenty independently sampled reads spanning these flanks support `ACGTC`, and one supports `ACTTC`. Under a model with rare, independent sequencing errors, this imbalance supports keeping the first route. But the same graph could come from two biological versions of the region, or from similar repeat copies sampled at different depths. Read counts alone cannot distinguish those explanations. Base quality and reads spanning additional flanking sequence can supply evidence beyond this local graph [@compeau2015].

Bubble detection answers only the structural question of where paths disagree. Choosing a route requires assumptions about the sample and the measurement process, and when the evidence stays ambiguous, ending the contiguous sequences at the bubble avoids an unsupported join.

### From a graph to contigs

After filtering and graph simplification, real data still contain branches supported by substantial evidence. To identify routes rather than count observations, aggregate repeated observations of the same $k$-mer into one edge and keep the count as an attribute of that edge. In this graph, a maximal path whose internal nodes each have exactly one incoming and one outgoing edge involves no choice, so every genome compatible with the graph contains the sequence it spells. That sequence is a *contig*.

Return to the repeat graph of @fig-euler-repeat-graph. With one edge per distinct 3-mer, AC has one incoming and one outgoing edge, and only CG branches, with two edges in and two out. The non-branching paths therefore run from CG back to CG:
$$
\texttt{CG}\to\texttt{GT}\to\texttt{TC}\to\texttt{CG}\ \text{spells}\ \texttt{CGTCG},
\qquad
\texttt{CG}\to\texttt{GA}\to\texttt{AC}\to\texttt{CG}\ \text{spells}\ \texttt{CGACG}.
$$
Both candidate genomes contain both contigs, which overlap only in the repeated CG:

```text
ACGTCGACG        ACGACGTCG
 CGTCG            CGACG
    CGACG            CGTCG
```

An assembler reports these two contigs and leaves their order open, which is exactly what the 3-mers determine. In @fig-euler-error-bubble, the two outgoing neighbors of `AC` are genuinely different, so aggregating counts leaves the bubble intact, and contigs stop at its branch points. A contig ends where a route encounters branching, a coverage gap, or a chromosome end. An isolated component that forms a directed cycle needs separate treatment: choose a starting node and report the circular sequence rather than wait for a branch [@compeau2015].

A contig is not necessarily a chromosome, and an assembly with fewer contigs is not automatically more accurate. Joining across an unresolved repeat can produce a longer but false sequence.

A simplified short-read workflow therefore has several distinct steps:

1. Count read-derived k-mers and remove counts poorly supported under an error model.
2. Build a de Bruijn graph with one edge per retained distinct k-mer and keep observation counts as edge attributes.
3. Simplify tips and bubbles when the evidence supports one path.
4. Emit unbranched contigs and use longer-range data to order or connect them.

The degree theorem and traversal algorithm solve the ideal composition problem exactly. The contig workflow answers a different biological question: which stretches can the measured reads support without choosing among unresolved joins?

## Exercises

::: {#exr-euler-check}
For $A\to B$, $B\to C$, $C\to A$, and $A\to C$, compute degrees and identify required endpoints.
:::

::: {#exr-euler-connectivity}
Give a balanced directed graph with no Eulerian cycle using all edges. Which theorem hypothesis fails?
:::

::: {#exr-euler-ambiguity}
Verify the two strings in @exm-repeat-ambiguity have the same 3-mer multiset. Compare their 4-mer multisets.
:::

::: {#exr-euler-error}
Change the last base of `ACGTC` to make `ACGTA`. Draw its 3-mer path together with that of `ACGTC`. Which edges survive, and why does this example produce two terminal alternatives rather than the rejoining bubble in @exm-euler-error-bubble?
:::

::: {#exr-euler-code}
Run `eulerian_path` on two disconnected cycles. Explain the final length check.
:::

::: {#refs}
**References**
:::
