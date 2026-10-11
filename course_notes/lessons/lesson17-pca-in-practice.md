---
published: true
title: "PCA in Practice"
subtitle: "Lesson 17 · 2027-03-24 · Week 11"
nocite: |
  @holmeshuber2019
---

## Core question

What can a PCA plot tell us, and what can it hide?

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

Lesson 16 found directions of maximum variance in a supplied numerical matrix. In an RNA-sequencing experiment, the available data begin earlier in the pipeline: RNA is sampled, copied into DNA, sequenced in fragments, and the reads are matched to genes. The count $C_{ij}$ records reads assigned to gene $j$ in sample $i$. It measures RNA through sampling and processing, rather than directly counting proteins or measuring a gene's function [@holmeshuber2019]. We retain samples as rows and genes as columns.

Suppose we want to investigate whether treated samples have a coordinated expression pattern different from controls. Sample identity, condition, and processing batch belong in a separate metadata table. A *batch* is a group processed together, perhaps on the same day or instrument. A shared processing difference can affect many genes and create large variance without representing the treatment response.

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

Within each condition, run 2 has exactly twice as many assigned reads in this noiseless example. A PCA of raw counts can use that larger spread to separate runs. We first need to state what differences the representation should retain.

## Sequencing depth, relative abundance, and logarithms

Let $L_i=\sum_j C_{ij}>0$ be the total assigned count for sample $i$. If a sampled read belongs to gene $j$ with probability $p_{ij}$, where $\sum_jp_{ij}=1$, its mean count at fixed depth is
$$
\mathbb E[C_{ij}\mid L_i]=L_ip_{ij}.
$$
Thus more sequencing can raise all counts even when the relative abundances are unchanged [@holmeshuber2019]. A simple normalization is counts per million:
$$
R_{ij}=10^6\frac{C_{ij}}{L_i}.
$$
Each row sums to one million. In the table, s1 and s2 both become $(100000,200000,700000)$, while s3 and s4 become $(200000,200000,600000)$. The representation removes their depth difference and keeps their different composition.

This is a modeling choice, not a universal correction. The row totals must be computed on the chosen count universe before subsequent gene filtering. Interpreting relative changes as changes in absolute RNA abundance requires further information: if every RNA abundance doubles, the proportions remain the same. If one highly abundant gene increases, other genes' proportions can decrease even when their absolute abundances are unchanged. More robust normalization methods use assumptions about a stable reference set or the majority of genes [@holmeshuber2019]. We use total-count normalization here to make the assumption visible.

Counts also span a wide range. On the original scale, changing a value from 1000 to 2000 contributes a much larger squared difference than changing it from 10 to 20. For this illustration, choose
$$
X_{ij}=\log_2(1+R_{ij}).
$$
For large $R_{ij}$, doubling it adds approximately one to $X_{ij}$, so multiplicative differences become approximately additive. The positive offset makes zeros finite. Here the offset is one *count per million*, not one raw read; changing units before adding an offset changes the transformation. At small $R_{ij}$ it is not a log fold change, and it does not guarantee that measurement variance is constant. Count-based variance-stabilizing transformations address that separate issue more directly [@holmeshuber2019].

::: {.algorithm}
**Algorithm: prepare a count matrix for this illustrative PCA.**

1. Match sample rows to metadata; verify the count units and orientation.
2. Compute positive sample totals on the full count table and divide each row by its total, multiplying by $10^6$.
3. Apply $\log_2(1+R)$, then retain a stated set of genes.
4. Subtract each retained gene's mean across samples. Optionally divide by its standard deviation if unit-variance weighting is the intended objective.
:::

```{.python .execute}
#| label: code-pca-count-preparation
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
assert np.allclose(relative_counts.sum(axis=1), 1_000_000)
assert np.allclose(X[0], X[1]) and np.allclose(X[2], X[3])
assert np.allclose(Z.mean(axis=0), 0)
np.round(X, 3)
```

The control rows now coincide, as do the treated rows. That is expected from our exact ratios, not evidence that real replicates have no variation. Nor does their separation test a treatment effect: this tiny example has no sampling noise. The next decisions concern which genes and which units of variation PCA should emphasize.

## Scaling changes the optimization problem

For $n\geq2$ samples and $p$ retained genes, recall the centered matrix $Z_{ij}=X_{ij}-\bar x_j$ and covariance $S=Z^\top Z/(n-1)$. Unscaled PCA chooses a unit loading $v\in\mathbb R^p$ to maximize $\|Zv\|^2/(n-1)$. Genes with larger numerical variance can dominate this objective.

For each gene with nonzero sample standard deviation,
$$
\sigma_j=\sqrt{\frac{1}{n-1}\sum_i Z_{ij}^2},
\qquad W_{ij}=\frac{Z_{ij}}{\sigma_j}.
$$
Each column of $W$ has variance one. Constant genes have $\sigma_j=0$ and must be omitted before division. Standardization makes each retained gene contribute equal marginal variance; it does not make their measurements equally reliable.

::: {#prp-pca-scaling-objective}
## The geometry selected by standardization
Let $G=\operatorname{diag}(\sigma_1,\ldots,\sigma_p)$, with positive diagonal entries. PCA on $W=ZG^{-1}$ maximizes
$$
\frac{\|Za\|^2}{n-1}
\quad\text{subject to}\quad a^\top G^2a=1.
$$
The squared distance between two samples becomes
$$
\|w_i-w_k\|^2=\sum_j\frac{(X_{ij}-X_{kj})^2}{\sigma_j^2}.
$$
:::

::: {.proof}
The standardized objective uses a unit $v$ and scores $Wv=ZG^{-1}v$. Set $a=G^{-1}v$. Then $v=Ga$, so $v^\top v=1$ becomes $a^\top G^2a=1$. The distance formula follows by subtracting the corresponding scaled coordinates.
:::

Thus scaling changes the weighting of gene differences rather than simply relabeling axes. A small difference in a low-variance gene can receive more weight than a larger difference in a high-variance gene. This may expose a coordinated pattern, but it can also amplify unreliable low-count measurements. Filtering should have a stated measurement-based purpose, such as excluding genes with too little coverage; choosing genes because they already separate the plotted conditions makes that separation partly a selection artifact.

### Revisit the Lesson 16 matrix

Use the same supplied transformed values:
$$
X=\begin{pmatrix}4&6&5\\5&7&3\\7&9&3\\8&10&5\end{pmatrix},
\qquad
Z=\begin{pmatrix}-2&-2&1\\-1&-1&-1\\1&1&-1\\2&2&1\end{pmatrix}.
$$
Its covariance is
$$
S=\frac13\begin{pmatrix}10&10&0\\10&10&0\\0&0&4\end{pmatrix}.
$$
The first direction combines genes 1 and 2 equally, with variance $20/3$ and explained fraction $5/6$. Suppose we multiply gene 3 by ten. Its variance becomes $400/3$; it now supplies PC1 and explains $400/(400+20)=20/21$ of the total. We have changed no samples or correlations, but changed which variation the objective favors.

After standardizing all three genes, the covariance is instead
$$
S_W=\begin{pmatrix}1&1&0\\1&1&0\\0&0&1\end{pmatrix}.
$$
Its eigenvalues are 2,1,0, so PC1 again combines genes 1 and 2 but now explains $2/3$. Multiplying gene 3 by ten before standardizing gives the same $W$, because its standard deviation also increases tenfold. Standardization removes that change of units, at the cost of replacing the original variance weighting.

```{.python .execute}
#| label: code-pca-scaling-check
example_X = np.array([[4, 6, 5], [5, 7, 3], [7, 9, 3], [8, 10, 5]], dtype=float)
example_Z = example_X - example_X.mean(axis=0)
rescaled_Z = example_Z * np.array([1, 1, 10])
W = example_Z / example_Z.std(axis=0, ddof=1)
rescaled_W = rescaled_Z / rescaled_Z.std(axis=0, ddof=1)
assert np.allclose(W, rescaled_W)
assert np.allclose(W.T @ W / 3, [[1, 1, 0], [1, 1, 0], [0, 0, 1]])
for representation in (example_Z, rescaled_Z, W):
    squared_singular = np.linalg.svd(representation, compute_uv=False) ** 2
    print(np.round(squared_singular / squared_singular.sum(), 3))
```

Expect fractions $(0.833,0.167,0)$, $(0.952,0.048,0)$, and $(0.667,0.333,0)$ from the hand calculation. These percentages describe three different representations. Their values do not rank the representations by biological validity.

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

We compute scores with the same SVD convention as Lesson 16: $Z=UDV^\top$, loadings are columns of $V$, scores are $T=ZV=UD$, and variances are $\lambda_k=d_k^2/(n-1)$. The helper below also fixes signs for display; a sign flip alone changes no reconstruction or separation.

::: {.algorithm}
**Algorithm: obtain comparable PCA summaries.** Center gene columns, compute an SVD, read loadings from $V$ and scores from $ZV$, and report score variances and their fractions of the total. Record the representation supplied to each fit; comparing different preprocessing choices compares different objectives.
:::

```{.python .execute}
#| label: code-pca-batch
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
#| fig-cap: "PCA of the planted batch-by-condition matrix before and after balanced additive batch adjustment, using equal data-unit scales. Blue denotes control and orange treated; circles denote batch −1 and squares batch +1 on the left. PCA is refitted after adjustment: adjusted PC1 is the original condition direction. Coincident samples are labeled together. Removing the larger batch contrast preserves the smaller planted condition contrast only under this example's design and model."
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
    ax.text(0, 1.04, label, transform=ax.transAxes, fontweight="bold")
plt.show()
```

The adjusted example has rank one, so PC2 scores are all zero and its loading direction is arbitrary. The disappearance of batch separation is expected here; in real data it does not by itself prove that biology was preserved.

### When batch and condition coincide

Suppose all controls were processed in one batch and all treated samples in another. Then $b_i=\tau_i$, and the model becomes
$$
x_i=\mu+\tau_i(a+c).
$$
Only the combined profile $a+c$ is observable. For any row vector $q$, replacing $a$ by $a+q$ and $c$ by $c-q$ gives exactly the same measurements. The separate technical and biological effects are therefore *unidentifiable*: these data cannot distinguish them, even without noise. Subtracting batch means can remove the entire condition contrast, while leaving it untouched can retain a technical artifact.

Recording batch labels is essential, but it cannot repair this lack of independent comparisons. Processing both conditions within each batch, with biological replication, provides information that the fully confounded design lacks. Adjustment in less balanced data requires a model and enough design information to separate effects; it cannot manufacture that information [@holmeshuber2019].

## Reading a PCA result as evidence

A *scree plot* graphs component variances or explained fractions against component number. It describes how the chosen matrix's spread is distributed. A sharp drop can motivate a low-dimensional summary, but a variance cutoff is not a biological criterion: our example's treatment signal was in the smaller component.

Sample labels and metadata help ask whether an unusual score follows treatment, batch, sequencing depth, or another measured feature. Loadings show which genes contribute to that direction. Refit under scientifically plausible preprocessing choices to see which relationships persist; changing an analysis until the desired groups separate is not validation. Near-equal eigenvalues can make individual axes unstable under small changes, even when their joint subspace is stable [@holmeshuber2019].

For projecting a future sample, retain the fitted gene set, transformation, means, scales, and loadings. Subtract the training gene means and use the training scales before multiplying by $V$. Recomputing them with the new sample changes the coordinate system. A count-based normalization must likewise use a stated rule for that sample.

Normalization, logarithms, centering, and scaling each take $\mathcal O(np)$ work on an $n\times p$ matrix. A dense thin SVD costs roughly $\mathcal O(np\min(n,p))$; it is usually the larger computational step. The harder inferential issue is deciding which representation and experimental comparisons justify a biological interpretation. The next lesson tests gene-level differences and addresses how many simultaneous tests can produce misleading discoveries.

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

::: {#exr-pca-future-projection}
Fit the first two loading directions of the Lesson 16 matrix and project a new sample $(6,8,5)$ using the original column means. Compare with the result of refitting after adding it as a fifth sample. Why are these two different tasks?
:::

::: {#refs}
**References**
:::
