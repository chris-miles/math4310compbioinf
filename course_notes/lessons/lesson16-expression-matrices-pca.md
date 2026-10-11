---
published: true
title: "Gene Expression and PCA"
subtitle: "Lesson 16 · 2027-03-22 · Week 11"
nocite: |
  @holmeshuber2019
---

## Core question

How can we summarize coordinated differences among samples using many gene measurements, and what does PCA preserve?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: NumPy and figure style"
import sys

sys.path.insert(0, "../style")
import numpy as np
import matplotlib.pyplot as plt
import coursefigs as cf
cf.use_style()
```

## From RNA reads to a measurement matrix

Cells copy portions of DNA into RNA molecules called transcripts. Their abundances differ across cell types and conditions. In a typical RNA-sequencing workflow, RNA is collected from each sample and copied into complementary DNA (cDNA); fragments of this DNA are sequenced. Matching those reads to genes gives a count table. A row such as `(120, 40, 0)` for three genes records how many reads were assigned to each gene in that sample. These counts measure sampled RNA through an experimental pipeline; they are not direct counts of proteins or a complete measure of gene activity [@holmeshuber2019, secs. 7.2, 8.2].

We use **samples as rows and genes as columns**. Many RNA-seq files use the opposite orientation, so check the row and column names before transposing. The biological question is whether samples differ through coordinated changes in several genes, rather than unrelated changes in each gene. A low-dimensional summary may expose an unusual sample or a pattern associated with a condition. It cannot establish the cause of that pattern.

Raw counts depend on sequencing depth and have strongly unequal scales. Today the input is an explicitly supplied, normalized and transformed numerical matrix $X\in\mathbb R^{n\times p}$, with $n\ge2$ samples and $p$ genes. Lesson 17 examines the preprocessing choices. We hold them fixed here to understand the mathematical summary. The entries need not be integers or remain nonnegative after centering.

PCA looks for large variation in this chosen representation. Useful biological structure might dominate that variation, but technical differences can dominate it too. "Large variance" is an optimization criterion, not a synonym for biological importance.

## Centering and a direction in gene space

Different genes have different baseline levels. To summarize differences among samples, remove each gene's mean:
$$
\bar x_j=\frac1n\sum_{i=1}^n X_{ij},
\qquad Z_{ij}=X_{ij}-\bar x_j.
$$
Thus each column of $Z$ sums to zero. A row $z_i$ describes how sample i differs from the mean expression profile. We center genes across samples, rather than subtracting a separate mean from each sample.

Choose a unit vector $v\in\mathbb R^p$, so $\|v\|_2^2=\sum_jv_j^2=1$. Its entries weight the genes. The projected coordinate of sample i is
$$
t_i=z_i v=\sum_{j=1}^p Z_{ij}v_j,
\qquad t=Zv\in\mathbb R^n.
$$
The gene weights $v$ are *loadings*; the sample coordinates $t$ are *scores*. For example, positive weights on two genes combine their departures from their means into one score. A loading has one entry per gene, while a score vector has one entry per sample [@holmeshuber2019, secs. 7.4, 7.6.4].

::: {#def-pca-first-direction}
## First principal component
The first PCA loading direction solves
$$
v_1\in\operatorname*{argmax}_{\|v\|_2=1}
\frac{\|Zv\|_2^2}{n-1}.
$$
The corresponding scores are $t_1=Zv_1$. The objective is their sample variance: because the scores have mean zero, it is their sum of squared values divided by $n-1$.
:::

The unit-length constraint is essential: multiplying $v$ by ten would otherwise multiply the objective by one hundred without choosing a better direction. The denominator $n-1$ is the usual sample-variance convention; replacing it by n would rescale all variances without changing the maximizing directions.

### Why large spread also gives a good reconstruction

A sample's projection back into gene coordinates is $t_i v^\top$. Its residual $z_i-t_i v^\top$ is perpendicular to $v$ because its inner product with v is zero. Pythagoras therefore gives
$$
\|Z\|_F^2=\|Zv\|_2^2+\|Z-(Zv)v^\top\|_F^2,
$$
where $\|M\|_F^2=\sum_{i,j}M_{ij}^2$ is the squared Frobenius norm. The total on the left is fixed. Maximizing the retained spread is therefore equivalent to minimizing squared reconstruction error for a line through the centered origin [@holmeshuber2019, sec. 7.3.2]. This is a geometric guarantee about the input matrix, not a guarantee that omitted variation is unimportant biology.

## Why an eigenvector solves the problem

Define the sample covariance matrix
$$
S=\frac{Z^\top Z}{n-1}\in\mathbb R^{p\times p}.
$$
Its entry $S_{jk}$ averages products of the centered measurements for genes j and k. Positive values indicate that their departures tend to have the same sign; negative values indicate opposite signs. The diagonal $S_{jj}$ is the variance of gene j. No probability distribution for the measurements is needed to compute this descriptive matrix.

The PCA objective is $v^\top Sv$. Recall that an eigenvector satisfies $Sw=\lambda w$: multiplication by S leaves its direction unchanged. A real symmetric matrix such as $S$ has an orthonormal basis of eigenvectors $w_1,\ldots,w_p$. Here the eigenvalues are nonnegative, since $v^\top Sv=\|Zv\|^2/(n-1)\ge0$. Order them as $\lambda_1\ge\cdots\ge\lambda_p\ge0$.

::: {#thm-pca-variance}
## Maximum-variance direction
A unit eigenvector of S with largest eigenvalue maximizes the PCA objective, and its score variance is $\lambda_1$.
:::

::: {.proof}
Write a unit v in the orthonormal eigenbasis as $v=\sum_j\alpha_jw_j$, with $\sum_j\alpha_j^2=1$. Then
$$
v^\top Sv=\sum_j\lambda_j\alpha_j^2
\le\lambda_1\sum_j\alpha_j^2=\lambda_1.
$$
Taking $v=w_1$ attains the bound.
:::

Further components maximize the same objective subject to being perpendicular to earlier loading directions. Choosing successive eigenvectors gives variances $\lambda_2,\lambda_3,\ldots$. Their scores have zero pairwise sample covariance, because
$$
\frac{t_j^\top t_k}{n-1}=v_j^\top Sv_k=0\quad(j\ne k).
$$
This is an algebraic property of the summary; it does not say the underlying genes or biological mechanisms are independent.

## A four-sample example by hand

Suppose these numbers have already been put on a comparable transformed expression scale. They are synthetic measurements, not raw read counts:

| sample | gene 1 | gene 2 | gene 3 |
|---|---:|---:|---:|
| s1 | 4 | 6 | 5 |
| s2 | 5 | 7 | 3 |
| s3 | 7 | 9 | 3 |
| s4 | 8 | 10 | 5 |

The column means are $(6,8,4)$, so
$$
Z=\begin{pmatrix}
-2&-2&1\\
-1&-1&-1\\
1&1&-1\\
2&2&1
\end{pmatrix},
\qquad
S=\frac13\begin{pmatrix}10&10&0\\10&10&0\\0&0&4\end{pmatrix}.
$$
Genes 1 and 2 have identical centered profiles, while gene 3 is perpendicular to that shared profile. Check that
$$
v_1=\frac1{\sqrt2}(1,1,0)^\top,\quad \lambda_1=\frac{20}{3};
\qquad
v_2=(0,0,1)^\top,\quad \lambda_2=\frac43.
$$
The remaining direction $(1,-1,0)^\top/\sqrt2$ has eigenvalue zero: there is no variation in the difference between the first two centered genes.

The score vectors are
$$
t_1=(-2\sqrt2,-\sqrt2,\sqrt2,2\sqrt2)^\top,
\qquad t_2=(1,-1,-1,1)^\top.
$$
Thus PC1 orders the samples by the coordinated shift in genes 1 and 2. PC2 retains the gene-3 contrast between the outer and inner samples in this ordering. One coordinate cannot preserve both patterns.

The first component's fraction of explained variance is
$$
\frac{\lambda_1}{\sum_j\lambda_j}
=\frac{20}{24}=\frac56\approx83.3\%.
$$
"Explained" means retained squared spread of these transformed measurements. It does not mean 83.3% of a treatment effect or of all biological information. To relate the ordering to a condition, we would need sample metadata and appropriate validation.

## SVD computes scores and loadings together

The singular value decomposition (SVD) writes an $n\times p$ matrix as
$$
Z=UDV^\top=\sum_{j=1}^{q}d_j u_jv_j^\top,
\qquad q=\min(n,p).
$$
In this reduced form, $U$ is $n\times q$, $V$ is $p\times q$, and $D$ is a $q\times q$ diagonal matrix. The columns of $U$ and $V$ are orthonormal, and the singular values satisfy $d_1\ge\cdots\ge d_q\ge0$. Each term is a rank-one matrix: one sample pattern multiplied by one gene-weight pattern. Rank counts independent directions; centering makes the sum of the row vectors zero, so at most $\min(n-1,p)$ singular values are positive [@holmeshuber2019, secs. 7.6.1–7.6.4].

Multiplying out connects the SVD to PCA:
$$
S=\frac1{n-1}VD^2V^\top,
\qquad ZV=UD.
$$
The columns of $V$ are loading directions; the score matrix is $T=ZV=UD$, with $T_{ij}$ the coordinate of sample i on component j. For each nonzero component, $\lambda_j=d_j^2/(n-1)$. The explained-variance fraction is therefore $d_j^2/\sum_kd_k^2$.

For our example, $d_1=\sqrt{20}$, $d_2=2$, and $d_3=0$. Multiplying any loading and its scores by −1 changes neither the projected data nor the variance. Repeated eigenvalues permit additional rotations within their shared subspace; individual directions need not be unique. We fix signs only to make displayed outputs reproducible.

::: {.algorithm}
**Algorithm: PCA of supplied transformed data.**

1. Subtract the column means to form $Z$.
2. Compute its SVD and order components by descending singular value.
3. Read loadings from $V$ and scores from $ZV$.
4. Report $d_j^2/(n-1)$ and the retained fraction of total squared singular values.
5. For $r$ components, reconstruct $Z_r=T_rV_r^\top$; add back the column means if reconstructing $X$.
:::

```{.python .execute}
#| label: code-pca-hand-example
X = np.array([[4, 6, 5], [5, 7, 3], [7, 9, 3], [8, 10, 5]], dtype=float)
means = X.mean(axis=0)
Z = X - means
U, singular_values, Vt = np.linalg.svd(Z, full_matrices=False)
loadings = Vt.T.copy()
for component in range(loadings.shape[1]):
    pivot = np.argmax(np.abs(loadings[:, component]))
    if loadings[pivot, component] < 0:
        loadings[:, component] *= -1
        U[:, component] *= -1
scores = Z @ loadings
variances = singular_values**2 / (X.shape[0] - 1)
total_spread = np.sum(singular_values**2)
if total_spread == 0:
    raise ValueError("Explained-variance fractions are undefined for constant samples.")
explained = singular_values**2 / total_spread
Z_one = scores[:, :1] @ loadings[:, :1].T
print("variances:", np.round(variances, 3))
print("PC1 loadings:", np.round(loadings[:, 0], 3))
print("PC1 scores:", np.round(scores[:, 0], 3))
print("explained:", np.round(explained, 3))
print("rank-one squared error:", np.sum((Z - Z_one)**2))
```

Expect variances $(6.667,1.333,0)$, first loadings $(0.707,0.707,0)$, first scores $(-2.828,-1.414,1.414,2.828)$, and explained fractions $(0.833,0.167,0)$. The rank-one reconstruction keeps the first two columns exactly and replaces the third centered column by zero. Its squared error is $1+1+1+1=4$, independently of the software.

```{.python .execute}
#| label: fig-pca-score-projection
#| fig-cap: "Sample coordinates for the two nonzero components and their projection onto PC1. Retaining 83.3% of the spread still discards the entire gene-3 pattern."
#| fig-alt: "Four sample points at PC1 scores minus2.828, minus1.414,1.414,2.828 and PC2 scores1,minus1,minus1,1. Vertical segments connect them to their rank-one projections at PC2 zero."
#| fig-width: 6
#| fig-height: 3.5
fig, ax = plt.subplots()
ax.axhline(0, color=cf.RULE, linewidth=1)
for i, sample in enumerate(["s1", "s2", "s3", "s4"]):
    ax.plot([scores[i, 0], scores[i, 0]], [0, scores[i, 1]],
            color=cf.RULE, linewidth=1)
    ax.annotate(sample, scores[i, :2], xytext=(5, 5), textcoords="offset points")
ax.scatter(scores[:, 0], scores[:, 1], s=25, label="sample")
ax.scatter(scores[:, 0], np.zeros(4), s=25, marker="x",
           color=cf.EMPHASIS, label="PC1 reconstruction")
ax.set_xlabel("PC1 score (83.3% of variance)")
ax.set_ylabel("PC2 score (16.7% of variance)")
ax.set_aspect("equal", adjustable="box")
ax.margins(x=0.15, y=0.35)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.35), ncol=2)
plt.show()
```

## Best low-rank reconstruction

For $0\le r\le q$ retained directions define $Z_r=\sum_{j=1}^r d_j u_jv_j^\top$. This solves a second explicit optimization problem:
$$
\min_{\operatorname{rank}(A)\le r}\|Z-A\|_F^2
=\|Z-Z_r\|_F^2
=\sum_{j>r}d_j^2.
$$
This is the squared-error form of the Eckart–Young result. In our example, r=1 leaves error 4 and r=2 leaves zero. More generally, the discarded squared singular values quantify precisely what the compression loses.

### Optional: the general reconstruction argument

::: {.proof}
Sketch. The rows of any rank-at-most-$r$ approximation $A$ lie in some subspace $W$ of dimension at most $r$. Projecting each row of $Z$ orthogonally onto $W$ is at least as close as its proposed row in $A$. Choose an orthonormal basis $w_1,\ldots,w_r$ for $W$, padding $W$ if necessary. Its retained squared spread is $\sum_k\|Zw_k\|^2$. Extend $d_j=0$ for $j>q$. In a full orthonormal eigenbasis $v_1,\ldots,v_p$ of $Z^\top Z$, the retained spread equals $\sum_j d_j^2\beta_j$, where $\beta_j=\sum_{k=1}^r(v_j^\top w_k)^2$. Each $\beta_j$ is the squared length of the projection of a unit eigenvector onto $W$, so $0\le\beta_j\le1$; summing over the full eigenbasis gives $\sum_j\beta_j=r$. The largest possible sum puts weight one on the $r$ largest squared singular values. The leading right singular vectors attain that bound. The projection identity then gives the displayed minimum and residual.
:::

If $p>q$, the omitted eigenvalues in this argument are zero. Take $0\le r\le q$; retaining more than the number of positive singular values gives no further improvement. Adding the means back yields $X_r=Z_r+\mathbf1\bar x^\top$. The rank guarantee concerns the centered variation $Z$, not necessarily this matrix with its baseline restored.

A dense reduced SVD costs roughly $\mathcal O(np\min(n,p))$ operations. When genes greatly outnumber samples, explicitly forming the $p\times p$ covariance matrix wastes storage; SVD can work directly with $Z$. The formulas identify the quantities we need without requiring that large matrix to be built.

## What the summary can and cannot establish

PCA supplies a reproducible description of the chosen input: gene weights, sample coordinates, and squared-error guarantees. A high loading identifies a gene contributing to that direction on the chosen scale. It does not establish that the gene causes a sample difference.

Without centering, the squared spread includes the common mean profile, so the first direction can describe abundance shared by every sample rather than differences among samples. Rescaling genes also changes the objective: a highly variable gene may dominate unless the analysis deliberately changes its weight. Lesson 17 examines these choices, sequencing depth, and sample metadata.

An outlier or a technical batch can account for large variance. Conversely, a biologically important change in a few genes may explain little total variance. PCA discards neither "noise" nor "signal" by definition; it discards the directions with less squared spread. The useful inference comes from relating scores and loadings to the experimental design and checking plausible alternative explanations.

## Exercises

::: {#exr-pca-hand-check}
Verify the displayed covariance matrix and both nonzero eigenpairs. Compute the variance of $t_1$ directly from its four entries and check the eigenvalue formula.
:::

::: {#exr-pca-projection}
Write the rank-one centered reconstruction explicitly, then add the means back. Which sample differences disappear? Explain why its squared error is 4.
:::

::: {#exr-pca-centering}
Add 100 to gene 3 in every row of $X$. Show that centered PCA is unchanged. Explain why an uncentered SVD need not be unchanged.
:::

::: {#exr-pca-score-covariance}
Using orthonormal loading eigenvectors, prove that different score columns have sample covariance zero. Why is that insufficient to establish independence of biological mechanisms?
:::

::: {#refs}
**References**
:::
