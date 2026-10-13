---
published: true
title: "PCA in Practice"
subtitle: "Lesson 16 · 2027-03-22 · Week 11"
nocite: |
  @holmeshuber2019
---

## Core question

Does a separation in a PCA plot reflect biology or how the samples were measured?

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

## Which matrix describes the experiment?

Suppose treated samples lie on one side of a PCA plot and controls on the other. That could reflect a shared RNA response to treatment. It could also reflect the experiment: one set of samples may have more sequenced reads, or may have been processed differently. PCA finds the largest variation in the numbers it receives; deciding what produced that variation requires the count table and the sample records [@holmeshuber2019].

In an RNA-sequencing experiment, RNA is sampled, copied into DNA, sequenced in fragments, and the reads are matched to genes. The count $C_{ij}$ records reads assigned to gene $j$ in sample $i$. We retain samples as rows and genes as columns. Lesson 15 held a transformed matrix fixed to study the PCA optimization; here we ask which representation makes the samples comparable and how to interpret the resulting coordinates.

The aim is to investigate whether treated samples have a coordinated expression pattern different from controls. Sample identity, condition, and processing batch belong in a separate *metadata* table, meaning records that describe each sample alongside its gene measurements. A *batch* is a group processed together, perhaps on the same day or instrument. A shared processing difference can affect many genes and create large variance without representing the treatment response. These records help us relate PCA directions to possible biological or technical causes [@holmeshuber2019].

Write $C\in\mathbb Z_{\geq0}^{n\times p}$ for counts from $n\geq2$ samples and $p$ genes. For a synthetic example, keep metadata separate from the numeric counts:

| sample | condition | sequencing run |
|---|---|---|
| s1 | control | 1 |
| s2 | control | 2 |
| s3 | treated | 1 |
| s4 | treated | 2 |

| sample | gene 1 | gene 2 | gene 3 | total assigned reads |
|---|---:|---:|---:|---:|
| s1 | 10 | 20 | 70 | 100 |
| s2 | 20 | 40 | 140 | 200 |
| s3 | 20 | 20 | 60 | 100 |
| s4 | 40 | 40 | 120 | 200 |

Within each condition, run 2 has exactly twice as many assigned reads in this synthetic example. A PCA plot can separate these samples sharply even though that difference reflects sequencing depth. To investigate treatment-associated expression, we need a representation that compares gene abundances at comparable depths.

## Sequencing depth, relative abundance, and logarithms

Let $L_i=\sum_j C_{ij}>0$ be the total assigned count for sample $i$. If a sampled read belongs to gene $j$ with probability $\rho_{ij}$, where $\sum_j\rho_{ij}=1$, its mean count at fixed depth is
$$
\mathbb E[C_{ij}\mid L_i]=L_i\rho_{ij}.
$$
Here $\rho_{ij}$ is the fraction of sampled reads expected to belong to that gene. Doubling the number of reads doubles the expected count without changing that fraction [@holmeshuber2019]. A simple normalization is counts per million:
$$
R_{ij}=10^6\frac{C_{ij}}{L_i}.
$$
Each row sums to one million. In the table, s1 and s2 both become $(100000,200000,700000)$, while s3 and s4 become $(200000,200000,600000)$. The representation removes their depth difference and keeps their different composition.

This representation measures composition: the share of reads assigned to each gene. If every RNA abundance doubles, those shares remain the same. If one highly abundant gene increases, other shares can decrease even when their absolute abundances are unchanged. Interpreting relative changes as absolute changes therefore requires additional information. Other normalization methods use assumptions about a stable reference set or the majority of genes [@holmeshuber2019]. Compute row totals before filtering genes, so removing a column does not silently redefine the reference total.

Counts also span a wide range. On the original scale, changing a value from 1000 to 2000 contributes a much larger squared difference than changing it from 10 to 20. For this illustration, choose
$$
X_{ij}=\log_2(1+R_{ij}).
$$
For large $R_{ij}$, doubling it adds approximately one to $X_{ij}$, so multiplicative differences become approximately additive. A twofold change in an abundant gene and a twofold change in a less abundant gene then make comparable contributions. The offset of one count per million makes zeros finite; near zero, it also strongly affects the difference being measured. This illustrative transform does not make all measurements equally reliable.

::: {.algorithm}
**Algorithm: prepare a count matrix for this illustrative PCA.**

1. Match sample rows to metadata; verify the count units and orientation.
2. Compute positive sample totals on the full count table and divide each row by its total, multiplying by $10^6$.
3. Apply $\log_2(1+R)$, then retain a stated set of genes.
4. Subtract each retained gene's mean across samples. Optionally divide by its standard deviation if unit-variance weighting is the intended objective.

For each representation, fit PCA using the Lesson 15 algorithm: center the columns, compute the SVD, and return scores $T=ZV$, loadings $V$, and component variances $d_k^2/(n-1)$. Record which matrix was supplied to each fit.
:::

```{.python .execute}
#| label: code-pca-count-preparation
def pca_summary(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return scores, loadings, and component variances for centered columns."""
    Z = X - X.mean(axis=0)
    _, singular_values, vt = np.linalg.svd(Z, full_matrices=False)
    loadings = vt.T
    for k in range(loadings.shape[1]):
        largest = np.argmax(np.abs(loadings[:, k]))
        if loadings[largest, k] < 0:
            loadings[:, k] *= -1
    scores = Z @ loadings
    variances = singular_values ** 2 / (X.shape[0] - 1)
    return scores, loadings, variances


counts = np.array([
    [10, 20, 70],
    [20, 40, 140],
    [20, 20, 60],
    [40, 40, 120],
], dtype=float)
totals = counts.sum(axis=1)
relative_counts = 1_000_000 * counts / totals[:, None]
X = np.log2(1 + relative_counts)
Z = X - X.mean(axis=0)
print("log-CPM values:\n", np.round(X, 3))
raw_scores, _, raw_variances = pca_summary(counts)
normalized_scores, _, normalized_variances = pca_summary(X)
print("raw explained fractions:", np.round(raw_variances / raw_variances.sum(), 3))
print("log-CPM fractions:",
      np.round(normalized_variances / normalized_variances.sum(), 3))
```

Raw PC1 explains 93.3% of the spread and separates the two sequencing depths. After normalization, the two control rows coincide, as do the two treated rows; log-CPM PCA has one nonzero component, accounting for 100% of its spread. That remaining contrast is a difference in composition between the conditions. The exact coincidences arise from the table's noiseless ratios; real replicates vary. These PCA percentages describe different input matrices, so a larger percentage does not validate a representation.

```{.python .execute}
#| label: fig-pca-depth-comparison
#| echo: false
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "The same count table summarized by PCA before and after depth normalization and a logarithm. Blue circles denote controls and orange triangles treated samples. Raw PC1 separates sequencing depths; the normalized representation retains the condition-associated composition difference."
#| fig-alt: "Left panel shows four samples in two columns: s1 and s3 at low raw PC1, s2 and s4 at high raw PC1. Right panel shows coincident controls at negative PC1 and coincident treated samples at positive PC1, with PC2 zero."
fig, axes = plt.subplots(1, 2)
for i in range(4):
    treated = i >= 2
    axes[0].scatter(raw_scores[i, 0], raw_scores[i, 1],
                    color="C1" if treated else "C0",
                    marker="^" if treated else "o", s=30)
    axes[0].annotate(f"s{i + 1}: depth {int(totals[i])}",
                     raw_scores[i, :2], xytext=(0, 7 if treated else -15),
                     textcoords="offset points", ha="center", fontsize=9)
for i, label in [(0, "control: s1, s2"), (2, "treated: s3, s4")]:
    axes[1].scatter(normalized_scores[i, 0], 0,
                    color="C1" if i == 2 else "C0",
                    marker="^" if i == 2 else "o", s=30)
    axes[1].annotate(label, (normalized_scores[i, 0], 0),
                     xytext=(0, 9), textcoords="offset points",
                     ha="center", fontsize=9)
axes[0].set(xlabel="Raw PC1 (93.3%)", ylabel="Raw PC2 (6.7%)")
axes[1].set(xlabel="log-CPM PC1 (100%)", ylabel="log-CPM PC2 (0%)")
axes[0].set_title("Raw counts")
axes[1].set_title("Normalized and transformed")
axes[0].set(xlim=(-80, 90), ylim=(-85, 85))
axes[1].set(xlim=(-1.15, 1.15), ylim=(-1.15, 1.15), yticks=[-1, 0, 1])
for ax, label in zip(axes, ["A", "B"]):
    ax.set_aspect("equal", adjustable="box")
    ax.text(-0.12, 1.04, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold")
plt.show()
```

## Scaling changes the optimization problem

For $n\geq2$ samples and $p$ retained genes, recall the centered matrix $Z_{ij}=X_{ij}-\bar x_j$ and covariance $S=Z^\top Z/(n-1)$. Unscaled PCA chooses a unit loading $v\in\mathbb R^p$ to maximize $\|Zv\|^2/(n-1)$. Genes with larger numerical variance can dominate this objective.

For each gene with nonzero sample standard deviation,
$$
\sigma_j=\sqrt{\frac{1}{n-1}\sum_i Z_{ij}^2},
\qquad W_{ij}=\frac{Z_{ij}}{\sigma_j}.
$$
Each column of $W$ has variance one. Constant genes have $\sigma_j=0$ and must be omitted before division. Standardization makes each retained gene contribute equal marginal variance; it does not make their measurements equally reliable.

The change can be read directly from a sample distance. After standardization,
$$
\|w_i-w_k\|^2=\sum_j\frac{(X_{ij}-X_{kj})^2}{\sigma_j^2}.
$$
Each difference is measured in units of that gene's standard deviation. A change of one standard deviation contributes one to the squared distance, whatever the original units. Thus scaling changes the weighting of gene differences and the directions PCA prefers. It can expose a coordinated pattern or amplify unreliable low-count measurements. Filtering should have a stated measurement-based purpose, such as excluding genes with too little coverage; choosing genes because they already separate the plotted conditions makes that separation partly a selection artifact.

### Revisit the Lesson 15 matrix

Use the same supplied transformed values:
$$
X=\begin{pmatrix}4&7&5\\5&6&3\\7&10&3\\8&9&5\end{pmatrix},
\qquad
Z=\begin{pmatrix}-2&-1&1\\-1&-2&-1\\1&2&-1\\2&1&1\end{pmatrix}.
$$
Its covariance is
$$
S=\frac13\begin{pmatrix}10&8&0\\8&10&0\\0&0&4\end{pmatrix}.
$$
The first direction combines genes 1 and 2 equally, with variance 6 and explained fraction $3/4$. Suppose we multiply gene 3 by ten. Its variance becomes $400/3$; it now supplies PC1 and explains $400/(400+20)=20/21$ of the total. We have changed no samples or correlations, but changed which variation the objective favors.

After standardizing all three genes, the covariance is instead
$$
S_W=\begin{pmatrix}1&0.8&0\\0.8&1&0\\0&0&1\end{pmatrix}.
$$
Its eigenvalues are $1.8,1,0.2$, so PC1 again combines genes 1 and 2 but now explains $1.8/3=60\%$. Multiplying gene 3 by ten before standardizing gives the same $W$, because its standard deviation also increases tenfold. Standardization removes that change of units, at the cost of replacing the original variance weighting.

```{.python .execute}
#| label: code-pca-scaling-check
example_X = np.array([[4, 7, 5], [5, 6, 3], [7, 10, 3], [8, 9, 5]], dtype=float)
example_Z = example_X - example_X.mean(axis=0)
rescaled_Z = example_Z * np.array([1, 1, 10])
W = example_Z / example_Z.std(axis=0, ddof=1)
rescaled_W = rescaled_Z / rescaled_Z.std(axis=0, ddof=1)
assert np.allclose(W, rescaled_W)
assert np.allclose(W.T @ W / 3, [[1, 0.8, 0], [0.8, 1, 0], [0, 0, 1]])
for representation in (example_Z, rescaled_Z, W):
    squared_singular = np.linalg.svd(representation, compute_uv=False) ** 2
    print(np.round(squared_singular / squared_singular.sum(), 3))
```

The fractions are $(0.75,0.167,0.083)$, $(0.952,0.043,0.005)$, and $(0.6,0.333,0.067)$, agreeing with the hand calculation. A gene's large contribution can reflect the chosen units. Whether to preserve its original variance or give it equal weight is a scientific choice about which differences to summarize.

## A dominant component can be a processing batch

To isolate a second issue, use a new synthetic matrix on an already transformed scale. Each batch includes one control and one treated sample:

| sample | batch sign $b_i$ | condition sign $\tau_i$ | gene 1 | gene 2 | gene 3 |
|---|---:|---:|---:|---:|---:|
| s1 | −1 | −1 (control) | 7 | 7 | 9 |
| s2 | −1 | +1 (treated) | 7 | 7 | 11 |
| s3 | +1 | −1 (control) | 13 | 13 | 9 |
| s4 | +1 | +1 (treated) | 13 | 13 | 11 |

We planted an additive model
$$
x_i=\mu+b_i a+\tau_i c,
\qquad \mu=(10,10,10),\quad a=(3,3,0),\quad c=(0,0,1).
$$
Here $a$ is the per-gene batch shift and $c$ the per-gene condition shift; $x_i$ and these profiles are row vectors. They are orthogonal, $ac^\top=0$. The sample contrasts also satisfy $\sum_i b_i=\sum_i \tau_i=0$ and $\sum_i b_i\tau_i=0$. Thus centering gives $z_i=b_i a+\tau_i c$, and the two effects do not mix in the cross-product matrix:
$$
Z^\top Z=4a^\top a+4c^\top c
=\begin{pmatrix}36&36&0\\36&36&0\\0&0&4\end{pmatrix}.
$$
The batch direction $(1,1,0)/\sqrt2$ has variance $72/3=24$. The condition direction $(0,0,1)$ has variance $4/3$. Consequently PC1 explains $72/76=18/19\approx94.7\%$ and separates batches, while PC2 explains $1/19\approx5.3\%$ and separates conditions. A rule retaining 90% of variance with one component would discard the planted biological contrast.

Use the same PCA function to check the calculation. Its sign convention affects the displayed orientation, while preserving reconstructions and separations.

```{.python .execute}
#| label: code-pca-batch
batch_X = np.array([[7, 7, 9], [7, 7, 11], [13, 13, 9], [13, 13, 11]], dtype=float)
batch = np.array([-1, -1, 1, 1])
condition = np.array([-1, 1, -1, 1])
scores_before, loadings_before, variance_before = pca_summary(batch_X)
assert np.allclose(variance_before, [24, 4 / 3, 0])
assert np.allclose(scores_before[:, 0], 3 * np.sqrt(2) * batch)
assert np.allclose(scores_before[:, 1], condition)
np.round(scores_before[:, :2], 3)
```

The scores agree with the planted contrasts: the largest spread is technical by construction. In real data we do not know the planted effects. Plotting measured batch and condition metadata helps investigate what a component is associated with; loadings identify the genes involved, but neither association alone proves a cause [@holmeshuber2019].

## What batch adjustment assumes

In this balanced example, a batch's mean is $\mu+b_i a$ because its two condition signs cancel. Subtracting that mean and restoring the overall mean gives
$$
x_i^{\mathrm{adj}}=x_i-\bar x_{\mathrm{batch}(i)}+\bar x
=\mu+\tau_i c.
$$
It removes the known additive batch shift while preserving the condition contrast. This calculation relies on the balanced design and additive model; it is not evidence that subtracting batch means is safe in any experiment.

::: {.algorithm}
**Algorithm: center known batches in the balanced example.** Compute the overall gene means. For each batch, subtract its gene means from its sample rows and add the overall means. Fit PCA to the adjusted matrix and check which planted contrast remains.
:::

```{.python .execute}
#| label: code-pca-batch-adjustment
adjusted_X = batch_X.copy()
grand_mean = batch_X.mean(axis=0)
for group in np.unique(batch):
    selected = batch == group
    adjusted_X[selected] = batch_X[selected] - batch_X[selected].mean(axis=0) + grand_mean
scores_after, loadings_after, variance_after = pca_summary(adjusted_X)
expected_adjusted = np.column_stack([np.full(4, 10), np.full(4, 10), 10 + condition])
assert np.allclose(adjusted_X, expected_adjusted)
assert np.allclose(scores_after[:, 0], condition)
assert np.allclose(variance_after, [4 / 3, 0, 0])
adjusted_X
```

```{.python .execute}
#| label: fig-pca-batch-comparison
#| fig-width: 6
#| fig-height: 2.5
#| fig-cap: "PCA before and after batch-mean adjustment of the balanced additive example, using equal data-unit scales. Blue denotes control and orange treated; circles denote batch -1 and squares batch +1 before adjustment. PCA is refitted afterward. The larger technical contrast disappears while the smaller condition contrast remains."
#| fig-alt: "Before adjustment, four samples occupy two batch-separated columns with controls below treated samples. After adjustment, controls s1 and s3 coincide at minus1 and treated s2 and s4 coincide at plus1, all with PC2 zero."
fig, axes = plt.subplots(1, 2)
condition_colors = {-1: "C0", 1: "C1"}
batch_markers = {-1: "o", 1: "s"}
for i in range(4):
    axes[0].scatter(scores_before[i, 0], scores_before[i, 1],
                    color=condition_colors[condition[i]], marker=batch_markers[batch[i]], s=30)
    axes[0].annotate(f"s{i + 1}", scores_before[i, :2],
                     xytext=(4, 5), textcoords="offset points")
for sign, label in [(-1, "control\ns1, s3"), (1, "treated\ns2, s4")]:
    axes[1].scatter(sign, 0, color=condition_colors[sign], s=30)
    axes[1].annotate(label, (sign, 0), xytext=(8 * sign, 8),
                     textcoords="offset points", ha="right" if sign < 0 else "left",
                     fontsize=9)
axes[0].set(xlabel="PC1 (94.7%)", ylabel="PC2 (5.3%)")
axes[1].set(xlabel="PC1 (100%)", ylabel="PC2 (0%)")
for ax, label in zip(axes, ["A  Before", "B  Adjusted"]):
    ax.set(xlim=(-6, 6), ylim=(-2, 2))
    ax.set_aspect("equal", adjustable="box")
    ax.text(0, 1.04, label, transform=ax.transAxes, fontsize=11, fontweight="bold")
plt.show()
```

The adjusted example has rank one, so PC2 scores are all zero and its loading direction is arbitrary. The disappearance of batch separation is expected here; in real data it does not by itself prove that biology was preserved.

## When batch and condition coincide

Suppose all controls were processed in one batch and all treated samples in another. Then $b_i=\tau_i$, and the model becomes
$$
x_i=\mu+\tau_i(a+c).
$$
In the balanced experiment, comparing control and treatment within the same batch held processing fixed. That comparison is now missing. A treated sample differs from a control in both condition and batch, so the data reveal only their combined profile $a+c$.

For example, if $a+c=(3,3,1)$, both of these explanations produce the same measurements and the same PCA plot:

| explanation | batch profile $a$ | condition profile $c$ |
|---|---|---|
| entirely technical | $(3,3,1)$ | $(0,0,0)$ |
| entirely biological | $(0,0,0)$ | $(3,3,1)$ |

::: {#thm-pca-confounding}
## Confounded effects cannot be separated from these measurements
In the model $x_i=\mu+b_i a+\tau_i c$, if $b_i=\tau_i$ for every sample, the batch and condition profiles are unidentifiable: distinct pairs $(a,c)$ give the same data.
:::

::: {.proof}
For any row vector $h$, set $a'=a+h$ and $c'=c-h$. Then
$$
\mu+b_i a'+\tau_i c'
=\mu+\tau_i(a+h+c-h)
=\mu+\tau_i(a+c)=x_i.
$$
Thus even exact measurements cannot determine the two profiles separately.
:::

Subtracting batch means in this design subtracts $\mu+\tau_i(a+c)$, leaving only the restored grand mean. It removes the entire condition-associated contrast, including any true biological contribution. Leaving the contrast untouched can retain a technical contribution. Neither choice resolves the missing comparison.

Recording batch labels is essential, but it cannot repair this lack of independent comparisons. Processing both conditions within each batch, with biological replication, provides information that the fully confounded design lacks. Adjustment in less balanced data requires a model and enough design information to separate effects; it cannot manufacture that information [@holmeshuber2019].

## Reading a PCA result as evidence

The reported component variances describe how the chosen matrix's spread is distributed. A sharp drop can motivate a low-dimensional summary, but a variance cutoff is not a biological criterion: the batch example's condition contrast occupies the smaller component.

Sample labels and metadata help ask whether an unusual score follows treatment, batch, sequencing depth, or another measured feature. Loadings show which genes contribute to that direction. Refit under scientifically plausible preprocessing choices to see which relationships persist; changing an analysis until the desired groups separate is not validation [@holmeshuber2019].

The two central decisions precede the plot: which representation should define sample differences, and which experimental comparisons separate biology from processing. PCA summarizes the resulting variation through scores and loadings. The next lesson asks a different question, comparing biological conditions for one gene, with an effect size and a null model for sampling variation.

## Exercises

::: {#exr-pca-composition}
Two samples have absolute RNA abundances $(10,10,10)$ and $(10,10,30)$. Suppose reads sample these proportions. Calculate the expected relative abundances. Which genes look reduced after total-count normalization even though their absolute abundance did not change? What extra measurement could distinguish a global abundance change from a depth change?
:::

::: {#exr-pca-scaling}
Reproduce the three explained-fraction vectors in the scaling example. Derive the standardized sample distance formula and explain what is gained and lost by giving every retained gene unit variance. What must happen to a constant gene?
:::

::: {#exr-pca-balanced-batch}
In the planted batch example, verify $\sum_i b_i\tau_i=0$ and derive both score vectors by multiplying $Z$ by their loading directions. Change the batch profile from $(3,3,0)$ to $(1/2,1/2,0)$. Which contrast becomes PC1? Does this change which effect is biologically meaningful?
:::

::: {#exr-pca-confounding}
Give two different pairs of batch and condition profiles that produce the same fully confounded data. Explain why overlapping batches in an adjusted PCA plot does not prove that an adjustment preserved the condition effect. Propose a design change that supplies the missing comparison.
:::

::: {#refs}
**References**
:::
