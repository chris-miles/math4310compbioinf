---
published: false
title: "Multiple Testing and False Discovery Rates"
subtitle: "Lesson 18 · 2027-03-29 · Week 12"
nocite: |
  @holmeshuber2019
---

## Core question

Which measured gene changes should we follow up when thousands of tests can produce chance discoveries?

Lesson 17 compared one gene between treated and control cultures. An expression experiment measures many genes in those same samples, so we can repeat the question for every column of the matrix. We want a list of genes to investigate as possible treatment responses. Even when a gene has no condition-associated change, sampling variation can give it a small p-value. Repeating the comparison thousands of times creates many opportunities for such chance results [@holmeshuber2019].

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## From individual comparisons to a discovery list

Suppose an experiment supplies a gene name, an estimated expression difference, and a calibrated p-value for each of $m=20000$ genes. If every null hypothesis is true and we flag each $p_j\le0.05$, the expected number of falsely flagged genes is at most $20000(0.05)=1000$. Equality holds for continuous uniform null p-values. This is an average over repeated experiments, not an upper bound on the number in one list; adding more measurements has also added more tests.

We therefore need a rule for the whole list. A *false discovery* is a gene we flag even though its tested null hypothesis is true. We cannot identify those genes from their p-values alone, but under stated testing assumptions we can control the expected fraction of false discoveries among the flagged genes. This is the purpose of the false discovery rate, or FDR. Effect estimates remain essential: a list with controlled statistical error can still include changes too small to matter for the biological question.

## The input table and desired output

Write $p_j$ for the p-value and $\widehat\Delta_j$ for the estimated condition difference of gene $j$, where $j=1,\ldots,m$. The input is the table $(j,\widehat\Delta_j,p_j)$; the desired output is a subset of gene indices to flag at a chosen FDR target $q$. We assume that the individual tests are calibrated under their null hypotheses. A multiple-testing rule cannot repair a test invalidated by confounded batches or dependent biological replicates. The illustrative values below and the later simulation teach the list-selection rule; they are supplied inputs, not repeated copies of Lesson 17's six-culture experiment.

## The question

- Choose a rejection set with a controlled fraction of false discoveries.

## The model

- Valid null p-values satisfy $P(p_j\le a)\le a$ for $0\le a\le1$; continuous exact null p-values are uniform, while permutation p-values can be discrete and conservative. State the independence or suitable positive-dependence conditions used for BH.

## The technique

- Define $\mathrm{FDP}=V/\max(R,1)$, where $V$ is the number of false rejections and $R$ the total number of rejections; $\mathrm{FDR}=E[\mathrm{FDP}]$. Distinguish a realized fraction from its expectation.
- FWER versus FDR; Benjamini--Hochberg step-up; the sorted-p-value picture; a differential-expression mini-case with effect size versus p-value.

## Worked example

- Use the illustrative list $(0.001,0.006,0.021,0.055,0.060,0.120,0.400,0.800)$ and apply BH by hand at $q=0.1$. These values are supplied for learning the algorithm; they do not come from the six-sample permutation example.
- Sort as $p_{(1)}\le\cdots\le p_{(m)}$ and find the largest rank $k$ with $p_{(k)}\le kq/m$. For this list, rank 4 fails because $0.055>4(0.1)/8=0.05$, but rank 5 passes because $0.060\le5(0.1)/8=0.0625$, so $k=5$. Reject all ranks through $k$, even when an earlier rank fails its own inequality.
- In a pedagogical model-based simulation, generate independent continuous $\mathrm{Uniform}(0,1)$ p-values for labeled true nulls, and a stated synthetic distribution concentrated near zero for labeled alternatives. Known labels let us count false rejections and compare realized FDP across repeated runs.

The sorted-p-value picture shows why the rule uses the largest qualifying rank. In this eight-gene example, rank $5$ satisfies its threshold, so the discovery list contains the first five genes. The cutoff applies to a list of hypotheses; it does not say that each selected gene has a $10\%$ probability of being a false discovery.

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

```{.python .execute}
#| label: fig-fdr-ranked-cutoff
#| code-fold: true
#| fig-cap: "The supplied sorted p-values and Benjamini-Hochberg thresholds at q = 0.1. The largest qualifying rank selects the first five genes for follow-up."
#| fig-alt: "Eight p-values plotted against rank, with a rising threshold line and a red ring around rank five; ranks six through eight lie above their thresholds."
p_values = np.array([0.001, 0.006, 0.021, 0.055, 0.060, 0.120, 0.400, 0.800])
ranks = np.arange(1, len(p_values) + 1)
thresholds = ranks * 0.1 / len(p_values)
fig, ax = plt.subplots()
ax.plot(ranks, p_values, "o", label="Sorted p-values")
ax.plot(ranks, thresholds, "--", label="BH threshold")
ax.plot(5, p_values[4], "o", markerfacecolor="none", markeredgewidth=1.5,
        markeredgecolor=cf.EMPHASIS, markersize=9)
ax.set(xlabel="Sorted rank", ylabel="p-value (log scale)", xticks=ranks,
       yscale="log")
ax.legend()
plt.show()
```

## Why it works, and what it costs

- BH controls FDR at $q$ under independence; proof sketch of the load-bearing inequality.

## Where the model breaks

- Dependence between genes.
- Significant is not large. For a volcano-plot mini-case, supply explicitly illustrative $(\widehat\Delta,p)$ pairs such as $(0.2,0.001)$, $(3,0.12)$ and $(-2,0.006)$, where $\widehat\Delta$ is the effect estimate on a fixed measurement scale. These pairs are distinct from both the p-value-only simulation and the six-culture experiment.
- Coarse permutation p-values limit resolution, but do not universally prevent BH discoveries: thresholds depend on rank. If all $m$ p-values equal $0.1$ and $q=0.1$, rank $m$ satisfies the cutoff and BH rejects all; the minimum need not satisfy $q/m$.
- If every p-value comes from the six-culture two-sided design and is at least $0.1$, BH makes no discoveries at $q<0.1$; at $q=0.1$, a rejection requires all $m$ values to equal $0.1$.
- Real gene-level p-values require an appropriate calibrated test and adequate biological replication. The illustrative values and labeled simulation teach multiple testing without claiming to provide a complete expression-analysis workflow.

## Exercises

- Apply BH to a given list.
- Simulate the null and count false discoveries with and without BH.

::: {#refs}
**References**
:::
