---
published: false
title: "Comparing Biological Conditions"
subtitle: "Lesson 17 · 2027-03-24 · Week 11"
nocite: |
  @holmeshuber2019
---

## Core question

Does a treatment change a gene's measured abundance beyond the variation among independent cultures?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The outline develops one-feature inference; Lesson 18 begins with the resulting p-values and addresses testing many features together.
:::

## Comparing treated and control cultures

Suppose we grow three independent cell cultures under control conditions and three under a treatment, then measure RNA abundance for each culture. As in Lessons 15 and 16, sequencing produces a gene-count table; a fixed normalization and transformation give the expression matrix we will compare. A PCA plot can show a broad difference among samples. Here the biological question is narrower: how much does one specified gene's measured abundance differ between conditions, and could variation among cultures account for that difference?

For a toy example, take the following measurements in fixed normalized units. Each row represents a separate culture, with its condition recorded alongside the measurement.

| culture | condition | gene abundance |
|---|---|---:|
| c1 | control | 2 |
| c2 | control | 3 |
| c3 | control | 4 |
| t1 | treated | 5 |
| t2 | treated | 6 |
| t3 | treated | 7 |

The treated mean exceeds the control mean by $3$ units. That is an observed effect, but its evidential strength depends on the experimental design. Independent cultures are *biological replicates*: they let us observe variation across separate experimental units. Measuring the same culture repeatedly gives *technical replicates*, which help assess the assay but do not supply more independent cultures [@holmeshuber2019]. We will compare the observed effect with the effects obtained by relabeling whole cultures under a model with no condition difference.

Consider just the three control cultures. Assaying c1, c2, and c3 once each gives the values $(2,3,4)$ from three independent cultures. Suppose instead that we assay only c1 three times and obtain $(1.9,2.0,2.1)$. Both designs produce three recorded values, but the second observes assay variation around one culture's abundance. It cannot tell us how much independently grown cultures vary. The drawing follows the physical source of each measurement so that the number of records cannot be mistaken for the number of biological replicates.

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

```{.python .execute}
#| label: fig-conditions-culture-replication
#| code-fold: true
#| fig-width: 6
#| fig-height: 3.4
#| fig-cap: "Three control cultures assayed once each, compared with one culture assayed three times; labels give illustrative normalized abundances. Three recorded values provide different biological evidence when they come from different numbers of independent cultures."
#| fig-alt: "One panel shows separate cultures c1, c2, c3 leading to measurements two, three, four; another shows only c1 leading to three measurements 1.9, 2.0, 2.1."
fig, ax = plt.subplots()
cf.replicate_sampling(
    ax, biological_values=[2, 3, 4], technical_values=[1.9, 2.0, 2.1],
    culture_labels=["c1", "c2", "c3"],
)
plt.show()
```

## One gene and its condition labels

Let $X\in\mathbb R^{n\times p}$ be the supplied expression matrix, with cultures as rows and genes as columns. For gene $j$, retain $x_i=X_{ij}$ and the condition label $g_i\in\{C,T\}$. The inputs to the comparison are these measurements, labels, and the design metadata; the outputs are a mean difference and a p-value describing its position in a specified null reference distribution. Culture identity, pairing, and processing batch matter because they determine which relabelings that reference can justify. The measurement scale remains fixed throughout, so the numerical contrast keeps a consistent biological interpretation.

## The question and effect size

- Estimate the population mean contrast $\Delta_j=\mu_{T,j}-\mu_{C,j}$ with $\widehat\Delta_j=\bar x_{T,j}-\bar x_{C,j}$, where $\mu_{T,j}$ and $\mu_{C,j}$ are the unknown condition means and bars denote averages over the observed biological replicates.
- Interpret the sign and magnitude on the supplied measurement scale. A difference of means on a log scale differs from a difference of mean raw abundances; do not silently call either a fold change.
- Separate biological importance from evidence against a null model. A small p-value alone does not state how large the contrast is.

## The null model

- For the unpaired example, assume that under the null, biological-replicate measurements have the same distribution in both conditions and are independent. This makes the labels exchangeable: conditional on the pooled observations and group sizes, every allowed assignment is equally likely.
- Equality of means alone does not justify unrestricted permutations when the distributions differ. The equal-distribution and independence assumptions are essential for this unpaired comparison.

## The technique

- Use the two-sided statistic $S=|\widehat\Delta_j|$. Relabel complete biological replicates while keeping the group sizes fixed, recompute $S$, and compare the observed value with this null reference.
- For the complete allowed assignment set $\mathcal A$, calculate
  $$
  p=\frac{\#\{a\in\mathcal A:S(a)\ge S_{\mathrm{obs}}\}}{|\mathcal A|}.
  $$
  Include the observed assignment and ties. This is a tail probability under the null assignment model, not the probability that the null is true.
- Enumerate the complete tiny example exactly and interpret the fraction of allowed assignments as the null reference.

## Worked example

- Suppose three control cultures give $(2,3,4)$ and three treated cultures give $(5,6,7)$ in fixed normalized units. The observed means are $3$ and $6$, so $\widehat\Delta=3$ and $S_{\mathrm{obs}}=3$.
- Choose the three treated labels among the six observations. There are $\binom63=20$ equally weighted assignments under the null. The assignment with treated values $(2,3,4)$ gives contrast $-3$; all other assignments except the observed one have smaller absolute contrast. Thus $p=2/20=0.1$. Because every assignment has a complementary assignment with the opposite signed contrast, $0.1$ is the smallest attainable two-sided p-value in this six-sample design.
- Plot the measurements by culture and the 20 null statistics. The observed contrast is large on this scale, but six biological replicates give a coarse reference distribution. Failure to reject at $0.05$ does not establish equal expression.

## Why it works, and what it costs

- Explain calibration by ranking the observed statistic among exchangeable assignments. With ties, the upper-tail rule is conservative: under the null, $P(p\le\alpha)\le\alpha$. Exact permutation p-values are discrete and need not be uniform.
- Count $\binom n{n_T}$ assignments for $n_T$ treated samples. Each direct mean comparison costs $O(n)$.

## Where the model breaks

- Duplicate the assay rows for each culture and explain why treating them as independent biological replicates changes the reference without adding independent evidence.
- Return to the treatment-by-batch example from Lesson 16. Confounding and dependence can invalidate unrestricted relabeling even when the code runs.
- Count-specific models and a full RNA-sequencing analysis are outside this lesson. The common lesson is to connect an effect, independent experimental units, and a defensible null comparison.

## Exercises

- Enumerate the 20 assignments and verify the exact p-value. Explain what changing to a prespecified one-sided statistic would change.
- Identify the independent experimental unit in a study with four organisms and three repeated measurements per organism.
- Explain why duplicating technical measurements does not provide independent evidence for the population mean contrast.

## Optional extensions

- For larger examples, sample $B$ assignments uniformly and independently and use $(1+b)/(B+1)$, where $b$ counts sampled statistics at least as large as the observed statistic. Discuss simulation variability and the smallest attainable reported p-value.
- Distinguish the equal-distribution observational null from a randomization test of no treatment effect justified by randomized treatment assignment. Use only the assignments allowed by the experimental design.
- For paired samples, consider swaps within pairs; for blocked experiments, consider permutations within suitable blocks. Optional exercise: describe the allowed assignments for a matched comparison and explain why unrestricted relabeling changes the reference.

## Bridge to multiple testing

This lesson produces an effect estimate and one calibrated p-value for a specified feature. Lesson 18 starts with a vector of such p-values and asks how to choose discoveries when thousands of comparisons are considered together. It develops false discovery proportions and Benjamini-Hochberg rather than repeating the single-feature test.

::: {#refs}
**References**
:::
