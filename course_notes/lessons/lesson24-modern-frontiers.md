---
published: false
title: "Modern Frontiers in Bioinformatics"
subtitle: "Lesson 24 · 2027-04-21 · Week 15"
nocite: |
  @compeau2015, @goh2026expertise, @hickey2024pangenome
---

## Core question

Which genome differences can longer DNA reads and references representing several genomes help us recover?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. A sequence-reconstruction comparison connects familiar methods to newer data. An optional protein-family extension introduces alternative representations. It does not survey every current tool.
:::

## Recovering sequence that a short observation cannot place

Suppose several genome regions contain the same repeated DNA sequence but have different sequences on either side. A read lying entirely inside the repeat supplies letters without identifying which copy produced them. Lessons 6 and 7 turned such fragments into an assembly graph and showed why several reconstructions can fit the same observations. The biological question is whether we can recover the arrangement of these regions, including differences among genomes. A more elaborate algorithm cannot infer a missing connection from the repeated letters alone.

A *long read* still records a base sequence inferred from an instrument's signal, but it can include more surrounding DNA in one observation. A read that crosses a whole repeat and reaches distinct sequence on both sides provides a connection that shorter reads may miss [@compeau2015]. Consider the two reconstructions from Lesson 7:

| candidate | reconstructed sequence |
|---|---|
| 1 | `ACGTCGACG` |
| 2 | `ACGACGTCG` |

Both produce the same multiset of short 3-mers: `ACG` twice, and `CGT`, `GTC`, `TCG`, `CGA`, `GAC` once each. Those data cannot choose between the candidates. Now suppose an additional read reports `GTCGA`. It contains `GT` before the repeated `CG` and `A` after it. Candidate 1 contains that complete five-base string; candidate 2 does not. Assuming an error-free read from the same genome, in the displayed orientation, the added observation selects candidate 1 among these two possibilities. Replacing the 3-mer traversal algorithm while retaining only the original data would not supply this connection. Here five bases suffice because the example is tiny; what matters is spanning the ambiguity, rather than a universal read-length cutoff.

```{.python .execute}
#| label: setup
#| code-fold: true
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

```{.python .execute}
#| label: fig-frontiers-spanning-read
#| code-fold: true
#| fig-width: 6
#| fig-height: 2.2
#| fig-cap: "Two reconstructions with identical 3-mer multisets; the red bracket locates the added read GTCGA in candidate 1. The longer observation connects flanks that the short-word data leave unpaired."
#| fig-alt: "Two rows of DNA letters, ACGTCGACG and ACGACGTCG; a bracket under GTCGA appears only in the first row."
candidates = ["ACGTCGACG", "ACGACGTCG"]
spanning_read = "GTCGA"
fig, axes = plt.subplots(2, 1)
for i, (ax, sequence) in enumerate(zip(axes, candidates), start=1):
    cf.sequence_windows(ax, sequence, word=spanning_read, label=f"Candidate {i}")
    ax.set_ylim(-1.2, 0.6)
plt.show()
```

The reference representation can also change. A *pangenome* combines sequence from several genomes; a pangenome graph represents shared stretches and alternative sequence paths [@hickey2024pangenome]. Instead of asking where a read fits one linear reference, we can ask which path its letters support. Today extends familiar sequence and graph reasoning to these richer observations and representations, with an optional protein-family comparison. Each extension still needs an explicit input, an output we can inspect, and evidence connecting that output to a biological claim.

For a small graph example, suppose two genomes contain `ACGTT` and `ACATT`. They share a prefix `AC` and suffix `TT`, with alternatives `G` and `A` between them. Concatenating the node strings along the two paths recovers the two sequences. An error-free read `CATT` follows the lower alternative, whereas `CGTT` follows the upper one. A reference containing only `ACGTT` would describe the first read using a mismatch; the graph can represent its alternative sequence explicitly. This graph depicts known sequence alternatives, rather than the unknown ordering in the assembly example above.

```{.python .execute}
#| label: fig-frontiers-variation-paths
#| code-fold: true
#| fig-width: 6
#| fig-height: 2.6
#| fig-cap: "A toy variation graph represents ACGTT and ACATT as paths sharing the AC and TT nodes; the path compatible with read CATT is highlighted. Alternative sequence can be part of the reference instead of appearing only as a mismatch to one chosen genome."
#| fig-alt: "The AC node splits toward G above and A below, and both alternatives join TT; the lower AC-to-A-to-TT path is red."
variation_edges = [("AC", "G"), ("G", "TT"), ("AC", "A"), ("A", "TT")]
variation_positions = {"AC": (0, 0), "G": (1.5, 0.7),
                       "A": (1.5, -0.7), "TT": (3, 0)}
fig, ax = plt.subplots()
cf.digraph(ax, variation_edges, pos=variation_positions,
           highlight_edges=[("AC", "A"), ("A", "TT")], node_size=750)
plt.show()
```

## Long reads and genome graphs

- Revisit assembly when two copies of a repeated region have different unique flanks. Short reads contained entirely inside the repeat cannot identify the copy that produced them.
- Use a small repeat example to compare a collection of short reads with a read spanning a unique flank, the repeat, and the next unique flank. Show which edge ordering the spanning observation supports and which ambiguity remains if no read spans the whole repeat.
- Compare the information gained from a longer observation with changes in sequencing error assumptions. More context and more accurate base calls address different obstacles.
- Introduce a pangenome graph as a representation of shared sequence and alternatives across genomes. Trace two supplied paths through a small variation graph; compare the information represented with one linear reference, without deriving a new graph algorithm.

## Optional extension: protein families and learned representations

- Pose a sequence-family question: can a new protein share conserved features with a family even when no single pairwise alignment gives strong evidence?
- Start from a supplied small multiple alignment. Estimate a separate residue distribution at each conserved position, adding pseudocounts so unseen residues are not assigned probability zero. Compare a new same-length sequence $s=s_1\cdots s_L$ with a positive background residue distribution $b$ using $\sum_{i=1}^{L}\log(p_i(s_i)/b(s_i))$, where $p_i$ is the residue distribution at position $i$. This log-odds score extends the comparisons from Lessons 2 and 9 to position-specific family evidence.
- Explain why an ordinary chain through those positions cannot handle an inserted or missing residue. Introduce profile-HMM match, insertion, and silent deletion states as an extension connecting alignment to the familiar hidden-state model [@compeau2015]. Keep training and a full decoding implementation outside the lesson.
- Contrast this explicit model with a supplied learned representation of protein sequences. Identify its training observations and prediction target; ask what a neighboring vector or predicted label supports about family membership or function. Before developing this optional illustration, inspect a suitable published primary source and select a validated example; no learned-representation example is required for the core lesson.

## Worked comparisons and biological interpretation

- Carry the repeat example through two reconstructions and show the observation responsible for a resolved placement. State which conclusions rely on the repeat being spanned and on the error assumptions.
- If the optional protein example is developed, compare a family-model score and a representation-based result with independent validation evidence. Sequence similarity or representation proximity can suggest an experiment; it does not alone establish a protein's activity in a cell.

## Validation and limitations

- Distinguish additional information in the data from a more flexible model. An algorithm's optimum is a guarantee under its objective; biological interpretation also needs evidence that the model and observations address the intended question.
- Discuss training overlap and distribution shifts for the supplied representation example. State which validation observations were independent of model development.
- Connect these limits to earlier lessons: repeat ambiguity, uncertain hidden labels, batch confounding, and the distinction between observed and biological states each change what an analysis can establish.

## Exercises

- Choose one method discussed in class. State its literal input, output, one assumption, and one independently checkable biological claim.
- Identify an assumption shared with a method from an earlier lesson and explain how the newer method changes it.
- Optional exercise or Flex activity: audit a prepared AI-produced expression analysis. Check biological-replicate identities, recompute one statistic, inspect treatment-by-batch confounding and multiple-testing decisions, and revise the biological conclusion to match the evidence [@goh2026expertise].

::: {#refs}
**References**
:::
