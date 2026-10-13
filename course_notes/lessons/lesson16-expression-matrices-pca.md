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

We use **samples as rows and genes as columns**. Many RNA-seq files use the opposite orientation, so check the row and column names before transposing. The biological question is whether samples differ through coordinated changes in several genes. For example, a change in cell proliferation can affect many genes together [@holmeshuber2019, sec. 7.2]. Looking at one gene at a time can miss that shared pattern. We want a few coordinates per sample that retain as much of the measured variation as possible.

Raw counts depend on sequencing depth and have strongly unequal scales. Today the input is an explicitly supplied, normalized and transformed numerical matrix $X\in\mathbb R^{n\times p}$, with $n\ge2$ samples and $p$ genes. Lesson 17 examines the preprocessing choices. We hold them fixed here to understand the mathematical summary. The entries need not be integers or remain nonnegative after centering.

Principal component analysis (PCA) constructs those coordinates by combining genes. We first need to specify which combinations count as good summaries.

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

The unit-length constraint is essential: multiplying $v$ by ten would otherwise multiply the objective by one hundred without choosing a better direction. It lets us compare directions on equal terms.

### Why large spread also gives a good reconstruction

A score tells us where a sample lies along a line, but discards its position perpendicular to that line. To reconstruct the sample from its score, multiply the coordinate $t_i$ by the direction: $\hat z_i=t_i v^\top$. The residual $z_i-\hat z_i$ is perpendicular to $v$, since $(z_i-t_i v^\top)v=t_i-t_i=0$. Each sample therefore forms a right triangle with the centered origin and its projection. Adding Pythagoras over the samples gives
$$
\|Z\|_F^2=\|Zv\|_2^2+\|Z-(Zv)v^\top\|_F^2,
$$
where $\|M\|_F^2=\sum_{i,j}M_{ij}^2$ is the squared Frobenius norm. Rotating the line redistributes a fixed total: a direction with more spread along it has less squared distance left perpendicular to it. Thus the line that maximizes score variance also minimizes reconstruction error [@holmeshuber2019, sec. 7.3.2]. PCA's two descriptions, a direction of large spread and a line that fits the samples well, solve the same problem.

## Why an eigenvector solves the problem

Define the sample covariance matrix
$$
S=\frac{Z^\top Z}{n-1}\in\mathbb R^{p\times p}.
$$
Its entry $S_{jk}$ averages products of the centered measurements for genes j and k. Positive values indicate that their departures tend to have the same sign; negative values indicate opposite signs. The diagonal $S_{jj}$ is the variance of gene j. No probability distribution for the measurements is needed to compute this descriptive matrix.

The PCA objective is $v^\top Sv$. Recall that an eigenvector satisfies $Sv_k=\lambda_kv_k$: multiplication by $S$ leaves its direction unchanged. A real symmetric matrix such as $S$ has an orthonormal basis of eigenvectors $v_1,\ldots,v_p$. Its eigenvalues are nonnegative, since $v^\top Sv=\|Zv\|^2/(n-1)\ge0$. Order them as $\lambda_1\ge\cdots\ge\lambda_p\ge0$.

Along eigenvector $v_k$, the score variance is exactly $\lambda_k$. A general unit direction combines these perpendicular directions. Its variance is a weighted average of their eigenvalues, with weights equal to its squared coefficients. Because those weights sum to one, mixing in a direction with smaller variance cannot improve on the largest eigenvalue. The calculation below makes that argument precise.

::: {#thm-pca-variance}
## Maximum-variance direction
A unit eigenvector of S with largest eigenvalue maximizes the PCA objective, and its score variance is $\lambda_1$.
:::

::: {.proof}
Write a unit $v$ in the orthonormal eigenbasis as $v=\sum_k\alpha_kv_k$, with $\sum_k\alpha_k^2=1$. Orthogonality removes the cross terms, giving
$$
v^\top Sv=\sum_k\lambda_k\alpha_k^2
\le\lambda_1\sum_k\alpha_k^2=\lambda_1.
$$
Taking $v=v_1$ puts all the weight on $\lambda_1$ and attains the bound.
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
| s1 | 4 | 7 | 5 |
| s2 | 5 | 6 | 3 |
| s3 | 7 | 10 | 3 |
| s4 | 8 | 9 | 5 |

The column means are $(6,8,4)$, so
$$
Z=\begin{pmatrix}
-2&-1&1\\
-1&-2&-1\\
1&2&-1\\
2&1&1
\end{pmatrix},
\qquad
S=\frac13\begin{pmatrix}10&8&0\\8&10&0\\0&0&4\end{pmatrix}.
$$
Genes 1 and 2 tend to increase together, but are not identical. Each alone has variance $10/3$. Combining them equally gives a larger variance:
$$
S\frac{(1,1,0)^\top}{\sqrt2}
=6\frac{(1,1,0)^\top}{\sqrt2}.
$$
The joint direction therefore retains variance 6, compared with $10/3$ for either gene alone. The complete eigenpairs, in decreasing order, are
$$
v_1=\frac1{\sqrt2}(1,1,0)^\top,\quad \lambda_1=6;
\qquad
v_2=(0,0,1)^\top,\quad \lambda_2=\frac43.
$$
The third direction is $v_3=(1,-1,0)^\top/\sqrt2$, with $\lambda_3=2/3$. It records the smaller disagreements between genes 1 and 2.

The score vectors are
$$
t_1=(-3,-3,3,3)^\top/\sqrt2,
\qquad t_2=(1,-1,-1,1)^\top.
$$
For instance, sample s1 has score $(-2-1)/\sqrt2=-3/\sqrt2$ on PC1. Its reconstruction is $(-3/2,-3/2,0)$, so its residual is $(-1/2,1/2,1)$, with squared length $3/2$. PC1 separates s1 and s2 from s3 and s4 through the coordinated shift in genes 1 and 2. It loses both the gene-3 contrast and the smaller gene-1 versus gene-2 differences.

The first component's fraction of explained variance is
$$
\frac{\lambda_1}{\sum_j\lambda_j}
=\frac{6}{6+4/3+2/3}=\frac34=75\%.
$$
"Explained" means retained squared spread of these transformed measurements. To relate the separation to a condition, we would need sample metadata.

@fig-pca-score-projection shows the same samples using only genes 1 and 2. Projecting onto the first gene discards all variation in the second. Rotating toward their shared direction keeps more of both genes: the retained sum of squares rises from 10 to 18, and the residual falls from 10 to 2. Their total remains 20.

```{.python .execute}
#| label: fig-pca-score-projection
#| echo: false
#| fig-cap: "The first two centered genes projected onto a single-gene axis and onto their joint PCA direction. Red crosses mark reconstructed points and grey segments mark residuals. A shared expression direction retains more spread with less reconstruction error."
#| fig-alt: "Two panels show points at minus2,minus1; minus1,minus2; 1,2; and 2,1. Left projections lie on the horizontal axis; right projections lie on the diagonal, with shorter perpendicular residuals."
#| fig-width: 6
#| fig-height: 3.3
two_genes = np.array([[-2, -1], [-1, -2], [1, 2], [2, 1]], dtype=float)
fig, axes = plt.subplots(1, 2)
cf.pca_projection(axes[0], two_genes, [1, 0], ["s1", "s2", "s3", "s4"])
cf.pca_projection(axes[1], two_genes, [1, 1], ["s1", "s2", "s3", "s4"])
axes[0].set_title("Gene 1: retained 10, residual 10")
axes[1].set_title("PC1: retained 18, residual 2")
for ax, label in zip(axes, ["A", "B"]):
    ax.text(-0.12, 1.04, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold")
axes[1].set_ylabel("")
plt.show()
```

## SVD computes scores and loadings together

The singular value decomposition (SVD) writes an $n\times p$ matrix as
$$
Z=UDV^\top=\sum_{j=1}^{q}d_j u_jv_j^\top,
\qquad q=\min(n,p).
$$
In this reduced form, $U$ is $n\times q$, $V$ is $p\times q$, and $D$ is a $q\times q$ diagonal matrix. The columns $u_k$ of $U$ are orthonormal, as are the columns $v_k$ of $V$, and the singular values satisfy $d_1\ge\cdots\ge d_q\ge0$. Each term is a rank-one pattern: entry $(i,j)$ is $d_k u_{ik}v_{jk}$. Samples differ in how strongly they express that pattern, while the same gene weights are used for every sample [@holmeshuber2019, secs. 7.6.1, 7.6.3].

Rank counts independent directions. Centering makes the $n$ sample row vectors add to the zero vector (each column sums to zero), so at most $\min(n-1,p)$ singular values are positive. Even a matrix with thousands of genes cannot have more than $n-1$ independent centered sample contrasts.

Multiplying out connects the SVD to PCA:
$$
S=\frac1{n-1}VD^2V^\top,
\qquad ZV=UD.
$$
The columns of $V$ are loading directions; the score matrix is $T=ZV=UD$, with $T_{ij}$ the coordinate of sample i on component j. For each nonzero component, $\lambda_j=d_j^2/(n-1)$. The explained-variance fraction is therefore $d_j^2/\sum_kd_k^2$.

The columns of $U$ are normalized sample-score patterns: $Zv_k=d_ku_k$. For our example, $d_1=\sqrt{18}$, $d_2=2$, $d_3=\sqrt2$, and $u_1=(-1,-1,1,1)^\top/2$. The first rank-one term is therefore
$$
d_1u_1v_1^\top=
\begin{pmatrix}
-3/2&-3/2&0\\
-3/2&-3/2&0\\
3/2&3/2&0\\
3/2&3/2&0
\end{pmatrix}.
$$
This piece describes two low-score samples and two high-score samples sharing the same gene profile. Adding the other two pieces recovers their remaining differences. Multiplying both $u_k$ and $v_k$ by $-1$ leaves the piece unchanged; we fix signs only to make displayed outputs reproducible.

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
X = np.array([[4, 7, 5], [5, 6, 3], [7, 10, 3], [8, 9, 5]], dtype=float)
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
explained = singular_values**2 / total_spread
Z_one = scores[:, :1] @ loadings[:, :1].T
print("variances:", np.round(variances, 3))
print("PC1 loadings:", np.round(loadings[:, 0], 3))
print("PC1 scores:", np.round(scores[:, 0], 3))
print("explained:", np.round(explained, 3))
print("rank-one squared error:", np.sum((Z - Z_one)**2))
```

The output gives variances $(6,1.333,0.667)$, first loadings $(0.707,0.707,0)$, first scores $(-2.121,-2.121,2.121,2.121)$, and explained fractions $(0.75,0.167,0.083)$. Each reconstructed row has squared residual $3/2$, so the total error is 6, agreeing with the hand calculation.

## Best low-rank reconstruction

For $0\le r\le q$ retained directions define $Z_r=\sum_{j=1}^r d_j u_jv_j^\top$. This solves a second explicit optimization problem:
$$
\min_{\operatorname{rank}(A)\le r}\|Z-A\|_F^2
=\|Z-Z_r\|_F^2
=\sum_{j>r}d_j^2.
$$
This is the squared-error form of the Eckart-Young result. In our example, $r=1$ leaves error $4+2=6$, $r=2$ leaves error 2, and $r=3$ leaves zero. Each added component restores a pattern perpendicular to the ones already retained, reducing the error by its squared singular value.

::: {.callout-note collapse="true" title="Going deeper: why no other rank-r matrix fits better"}

::: {.proof}
Sketch. The rows of any rank-at-most-$r$ approximation $A$ lie in some subspace $W$ of dimension at most $r$. Projecting each row of $Z$ orthogonally onto $W$ is at least as close as its proposed row in $A$. Choose an orthonormal basis $w_1,\ldots,w_r$ for $W$, padding $W$ if necessary. Its retained squared spread is $\sum_k\|Zw_k\|^2$. Extend $d_j=0$ for $j>q$. In a full orthonormal eigenbasis $v_1,\ldots,v_p$ of $Z^\top Z$, the retained spread equals $\sum_j d_j^2\beta_j$, where $\beta_j=\sum_{k=1}^r(v_j^\top w_k)^2$. Each $\beta_j$ is the squared length of the projection of a unit eigenvector onto $W$, so $0\le\beta_j\le1$; summing over the full eigenbasis gives $\sum_j\beta_j=r$. The largest possible sum puts weight one on the $r$ largest squared singular values. The leading right singular vectors attain that bound. The projection identity then gives the displayed minimum and residual.
:::

If $p>q$, the omitted eigenvalues in this argument are zero. Take $0\le r\le q$; retaining more than the number of positive singular values gives no further improvement. Adding the means back yields $X_r=Z_r+\mathbf1\bar x^\top$. The rank guarantee concerns the centered variation $Z$, not necessarily this matrix with its baseline restored.

:::

A dense reduced SVD costs roughly $\mathcal O(np\min(n,p))$ operations. When genes greatly outnumber samples, explicitly forming the $p\times p$ covariance matrix wastes storage; SVD can work directly with $Z$. The formulas identify the quantities we need without requiring that large matrix to be built.

## What the summary can and cannot establish

Centering defines differences from a common baseline. Orthogonal projection divides their squared spread into retained variation and residual error. The eigenvector argument selects the best directions, and SVD computes their loadings, scores, and reconstruction errors together.

Those guarantees concern the matrix supplied to PCA. Without centering, the squared spread includes the common mean profile. Rescaling genes changes which differences count most. Lesson 17 examines these choices, sequencing depth, and sample metadata.

PCA retains directions with large squared spread; a technical batch or an outlier can supply that spread. A high loading identifies a gene contributing to a direction on the chosen scale, without establishing a causal role. The biological inference comes from relating these gene patterns and sample coordinates to the experiment.

## Exercises

::: {#exr-pca-hand-check}
Verify the displayed covariance matrix and all three eigenpairs. Compute the variance of $t_1$ directly from its four entries and check the eigenvalue formula.
:::

::: {#exr-pca-projection}
Add the means back to the displayed rank-one reconstruction. Which sample differences disappear? Explain why its squared error is 6, and check that retaining PC2 reduces the error to 2.
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
