---
published: true
title: "Eulerian Paths, Repeats, and Errors"
subtitle: "Lesson 8 · Week 4"
---

## Core question

When does a de Bruijn graph determine a genome, and what do repeats and errors do to that answer?

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

::: {#def-eulerian-trail}
## Eulerian path and cycle
An *Eulerian path* uses every edge exactly once. An *Eulerian cycle* is an Eulerian path whose first and last nodes agree.
:::

Vertices may repeat. The restriction is on edges, which represent the $k$-mer occurrences that reconstruction must consume.

::: {#prp-assembly-eulerian}
## Reconstruction is an Eulerian path problem
A string $x$ with $\operatorname{Comp}_k(x)=K$ determines an Eulerian path in $D(K)$. Every Eulerian path in $D(K)$ spells a string with composition $K$.
:::

::: {.proof}
Sketch. Consecutive genomic $k$-mers overlap in $k-1$ symbols, so their edges meet in order. Conversely, traversing an edge appends its final symbol. Using every edge occurrence once gives exactly $K$.
:::

This equivalence assumes an ideal spectrum. It does not guarantee existence for noisy data or uniqueness.

## Degree conditions

Write $d^+(v)$ and $d^-(v)$ for out-degree and in-degree, counting parallel edges.

::: {#thm-directed-euler-path}
## Eulerian path criterion
Suppose all nonzero-degree nodes of a finite directed multigraph are connected when directions are ignored. It has an Eulerian path exactly when either:

1. every node is balanced, giving an Eulerian cycle; or
2. one start $s$ has $d^+(s)-d^-(s)=1$, one end $t$ has $d^-(t)-d^+(t)=1$, and all other nodes are balanced.
:::

::: {.proof}
Sketch. Every visit to an internal node pairs one entering edge with one leaving edge, proving necessity. For sufficiency, add $t\to s$ if needed, follow unused edges until a cycle closes, then splice in closed walks from nodes incident to unused edges. Remove the added edge.
:::

Connectivity matters. Two disjoint cycles are balanced but cannot share one walk.

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

A genome walk contributes one entry and one exit whenever it passes through an internal $(k-1)$-mer. These contributions cancel in the degree difference. The first node has one unpaired exit, and the last has one unpaired entry. For a circular genome, every exit is paired with an entry and all nodes balance.

This accounting argument is a fast rejection test. If two nodes have one extra outgoing edge, no ordering of the given edge multiset can produce one linear string. It does not tell us whether the failure came from a missing k-mer, an extra erroneous k-mer, or data from more than one genome. Those are different biological explanations for the same graph obstruction.

Degrees also do not test connectivity. Imagine one cycle on nodes $A,B$ and another on $C,D$. All four nodes balance, but no edge joins the cycles. Any walk remains trapped in the component where it begins. The weak-connectivity condition in @thm-directed-euler-path rules out this case.

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


::: {#exm-hierholzer-stack}
## Hierholzer's algorithm on the shared graph
Use the edge order
$$
\texttt{TA→AC},\quad \texttt{AC→CG},\quad \texttt{CG→GA},\quad
\texttt{GA→AC},\quad \texttt{AC→CG}.
$$
Starting at `TA`, the stack grows through
$$
[\texttt{TA}], [\texttt{TA,AC}], [\texttt{TA,AC,CG}],
[\texttt{TA,AC,CG,GA}], [\texttt{TA,AC,CG,GA,AC,CG}].
$$
The last `CG` has no unused exit, so it moves to the output. The same is then true of `AC`, `GA`, and the preceding nodes. The output is built backward:
$$
[\texttt{CG,AC,GA,CG,AC,TA}].
$$
Reversing it gives the Eulerian path and spells `TACGACG`.
:::

Why append nodes only when stuck? A greedy walk can close a cycle while unused edges remain elsewhere. The delayed output lets the algorithm insert that cycle at the correct occurrence of its starting node. The stack is the data structure that performs the splicing described in the proof.

For an Eulerian graph, every edge is removed exactly once, and every node placed on the stack is popped exactly once. This accounts for the linear running time directly. The final length check in the code verifies that the traversal used all components rather than merely completing one balanced component.

Each traversal pushes one endpoint, so total time and space are $\mathcal O(|V|+|E|)$. Choosing a different outgoing edge may return a different valid assembly.

## Repeats and ambiguity

For `ACGTCGACG` at $k=3$,
$$K=\{\!\{\texttt{ACG},\texttt{ACG},\texttt{CGT},\texttt{GTC},\texttt{TCG},\texttt{CGA},\texttt{GAC}\}\!\}.$$
The graph branches at `CG`.

```{.python .execute}
#| label: fig-euler-repeat-graph
#| fig-cap: "A repeat creates a branch at `CG`. Either outgoing edge can be taken first in a complete Eulerian path."
#| fig-alt: "A directed multigraph with two edges from AC to CG and branches from CG toward GT and GA."
#| fig-width: 6
#| fig-height: 3.6
repeat_edges = [("AC", "CG"), ("AC", "CG"), ("CG", "GT"),
                ("GT", "TC"), ("TC", "CG"), ("CG", "GA"), ("GA", "AC")]
positions = {"AC": (0, 0), "CG": (1.6, 0), "GT": (3, .9),
             "TC": (4.2, .2), "GA": (2.5, -1.1)}
fig, ax = plt.subplots()
cf.digraph(ax, repeat_edges, pos=positions, highlight_nodes=["CG"])
plt.show()
```

::: {#exm-repeat-ambiguity}
## Two strings with one composition
Both `ACGTCGACG` and `ACGACGTCG` have composition $K$. Data consisting only of $K$ cannot favor either reconstruction.
:::

A repeat causes ambiguity when the available context cannot pair its incoming flanks with the correct outgoing flanks. Larger $k$ supplies more context, but demands longer error-free overlaps and fragments more readily under missing coverage.


### Existence and uniqueness

The degree theorem answers whether at least one traversal exists. Assembly also asks whether all traversals spell the same string. A node with two outgoing edges does not automatically imply two assemblies: one choice may become stuck before using all edges, or two edge choices may carry identical labels and spell the same symbols. Conversely, the graph in @exm-repeat-ambiguity has two complete traversals that spell different strings.

We can prove non-identifiability by exhibiting two strings with the same input. Once both strings produce $K$, every algorithm receiving only $K$ sees identical data in the two cases. No tie-breaking rule, faster traversal, or clever implementation can recover which string generated it. Additional information must cross the repeat: longer k-mers, paired reads with a known separation, or a long read spanning both flanks.

At $k=4$, the two strings in @exm-repeat-ambiguity no longer have the same composition. For `ACGTCGACG`, the 4-mers are `ACGT, CGTC, GTCG, TCGA, CGAC, GACG`. For `ACGACGTCG`, they are `ACGA, CGAC, GACG, ACGT, CGTC, GTCG`. The first contains `TCGA`; the second contains `ACGA`. Here one extra base of context separates the candidates. This calculation is specific to the example, not a theorem that larger k always gives a unique assembly.

## Errors, tips, and bubbles

Changing read `ACGTC` to `ACTTC` replaces 3-mers
$$\{\!\{\texttt{ACG},\texttt{CGT},\texttt{GTC}\}\!\}$$
with
$$\{\!\{\texttt{ACT},\texttt{CTT},\texttt{TTC}\}\!\}.$$

::: {#def-tip-bubble}
## Tips and bubbles
A *tip* is a short dead-end path. A *bubble* is a pair of paths that diverge from one node and later rejoin.
:::

An end error often makes a tip; an internal substitution can make a low-count bubble. Biological variation can also make a bubble, so lower count alone does not prove error. Assemblers combine depth, path length, base quality, and long-range support. They simplify supported artifacts and emit maximal unbranched sequences called *contigs* when branches remain [@pevzner2001].


### From a graph to contigs

After filtering and graph simplification, real data still contain branches supported by substantial evidence. An assembler can report every maximal path whose internal nodes have in-degree and out-degree one. Such a path has no local choice, so it spells a contig. The contig ends when the graph reaches a branch, a coverage gap, or a chromosome end.

This output separates what the reads determine from what they leave ambiguous. A contig is not necessarily a chromosome, and an assembly with fewer contigs is not automatically more accurate. Joining across an unresolved repeat can produce a longer but false sequence. Paired reads and long reads help because they connect graph regions farther apart than one k-mer overlap.

A simplified short-read workflow therefore has several distinct steps:

1. Count read-derived k-mers and remove counts poorly supported under an error model.
2. Build a de Bruijn graph while preserving multiplicity.
3. Simplify tips and bubbles when the evidence supports one path.
4. Emit unbranched contigs and use longer-range data to order or connect them.

The traversal guarantee applies to the supplied graph. Filtering and simplification make additional claims about sequencing noise and biological variation.

## Limitations

1. **Balance gives existence, not uniqueness.** The repeat graph has two reconstructions.
2. **Long repeats exceed local context.** Increasing $k$ helps only when reads support those longer exact overlaps.
3. **Errors break the exact model.** Cleaning changes the graph and therefore requires an error model.
4. **Coverage does not equal genomic multiplicity.** Sampling variation can blur one-copy and repeated regions.
5. **The genome may be unidentifiable.** Contigs report the unbranched stretches when data cannot resolve branches.

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
Draw the 3-mer paths for `ACGTC` and `ACTTC` together. Mark divergence and rejoining.
:::

::: {#exr-euler-code}
Run `eulerian_path` on two disconnected cycles. Explain the final length check.
:::

## Further reading

Compeau and Pevzner develop Eulerian assembly and repeat ambiguity [@compeau2015]. Pevzner, Tang, and Waterman give the graph basis for assembly with sequencing errors [@pevzner2001].
