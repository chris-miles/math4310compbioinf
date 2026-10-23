---
published: false
title: "Multiple Testing and False Discovery Rates"
subtitle: "Lesson 18 · 2027-03-29 · Week 12"
nocite: |
  @holmeshuber2019
---

## Core question

Which genes should we follow up when we test thousands of them at once and some small p-values must arise by chance?

```{.python .execute}
#| label: setup
#| include: false
import itertools
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## From one comparison to a screen

Lesson 17 asked whether a drug changes the expression of one stress-response gene in a human cell line. Comparing treated and control cultures produced two numbers for that gene: an effect size $\widehat\Delta_j$, the difference in mean expression, and a permutation p-value $p_j$ measuring how unusual that difference is among relabelings of the cultures. The sequencing run behind that comparison measured far more than one gene. RNA sequencing counts reads from every expressed gene in each culture, so the expression matrix $X$ has a column for every measured gene, roughly 20,000 protein-coding genes for human cells, and nothing stops us from repeating the comparison on every column.

The result is a table with one row per gene. On the $\log_2(1+\text{CPM})$ scale of Lesson 16, $\widehat\Delta_j=1$ means the treated cultures have about twice the typical abundance of the controls. Here are the first rows of a simulated screen that we will use throughout:

```{.python .execute}
#| label: code-fdr-simulate-screen
#| code-fold: true
#| code-summary: "Simulate the screen and compute permutation p-values"
def permutation_p_values(X: np.ndarray, n_control: int) -> tuple[np.ndarray, np.ndarray]:
    """Effect sizes and exact two-sided permutation p-values for every gene column."""
    n = X.shape[0]
    relabelings = []  # weights giving each relabeling's difference in means
    for treated in itertools.combinations(range(n), n - n_control):
        weights = np.full(n, -1 / n_control)
        weights[list(treated)] = 1 / (n - n_control)
        relabelings.append(weights)
    observed = X[n_control:].mean(axis=0) - X[:n_control].mean(axis=0)
    relabeled = np.abs(np.array(relabelings) @ X)  # relabelings x genes
    # Fraction of relabelings at least as large as observed (observed included).
    return observed, (relabeled >= np.abs(observed) - 1e-12).mean(axis=0)


def simulate_screen(rng: np.random.Generator, m: int, m_changed: int,
                    n_per_group: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Log2 expression, control rows first; returns X and the true-null indicator."""
    noise_sd = rng.uniform(0.15, 0.8, size=m)  # genes vary differently across cultures
    delta = np.zeros(m)
    changed = rng.permutation(m)[:m_changed]
    delta[changed] = rng.choice([-1, 1], size=m_changed) * rng.uniform(0.3, 2.0, size=m_changed)
    control = rng.normal(0.0, noise_sd, size=(n_per_group, m))
    treated = rng.normal(delta, noise_sd, size=(n_per_group, m))
    return np.vstack([control, treated]), delta == 0


rng = np.random.default_rng(18)
X, is_null = simulate_screen(rng, m=5000, m_changed=500)
delta_hat, p_values = permutation_p_values(X, n_control=6)
print("gene       delta_hat   p_value")
for j in range(5):
    print(f"gene{j:04d}  {delta_hat[j]:9.2f}  {p_values[j]:8.4f}")
```

The biological goal is a *discovery list*: genes that respond to the drug and deserve follow-up, such as a targeted measurement in new cultures or a knockdown experiment. Each follow-up costs weeks of work, so the list should not be dominated by genes whose small p-values arose by chance. With $m$ genes tested, the input is the table $(j,\widehat\Delta_j,p_j)$ for $j=1,\ldots,m$ and the output is a set of selected genes, chosen by a rule that tells us something about the errors the set contains.

## What a p-value does when nothing changed

A p-value is useful because we know how it behaves when the gene's null hypothesis is true, that is, when the drug does not change the gene's measured expression. Lesson 17 computed $p_j$ as the fraction of relabelings whose relabeled difference is at least as large in absolute value as the observed one, and @prp-permutation-calibration showed that such a p-value is *calibrated*: whenever the null hypothesis for gene $j$ holds,
$$
P(p_j\le a)\le a\qquad\text{for every } 0\le a\le 1 ,
$$
where the probability is over repeated experiments with fresh cultures. When equality holds for every $a$, the null p-value is *uniform* on $[0,1]$, as likely to fall in $[0,0.05]$ as in $[0.50,0.55]$. Permutation p-values take a limited set of values, so for them the inequality can be strict.

Suppose we flag every gene with $p_j\le t$ for some threshold $t$. Let $\mathcal H_0$ be the set of genes whose null hypothesis is true and $m_0=|\mathcal H_0|$. We do not know which genes these are, and that is the whole difficulty. The number of falsely flagged genes is
$$
V(t)=\sum_{j\in\mathcal H_0}\mathbf 1[p_j\le t],
$$
where the indicator $\mathbf 1[\cdot]$ is 1 when its condition holds and 0 otherwise. Recall from Lesson 2 that an expectation is a long-run average over repeated experiments. The average of a sum is the sum of the averages, and the average of an indicator is the probability of its event, so
$$
E[V(t)]=\sum_{j\in\mathcal H_0}P(p_j\le t)\le m_0\,t .
$$ {#eq-fdr-expected-false}

Each null gene has at most a $t$ chance of slipping under the threshold, and there are $m_0$ of them. If all $m=20{,}000$ genes are null and $t=0.05$, we expect $20{,}000\times0.05=1{,}000$ flagged genes, every one a false lead, from an experiment in which the drug did nothing. Each test behaved exactly as designed; the false leads come from asking the same question 20,000 times.

### A simulated screen with known answers

Real data never reveal which genes are null, so we simulate a screen where we know. Six control and six treated cultures are measured for $m=5000$ genes. For 500 of them the drug shifts mean $\log_2$ expression by between $0.3$ and $2$ in either direction; the other 4500 are null. Each gene gets the exact two-sided permutation p-value of Lesson 17, computed over all $\binom{12}{6}=924$ relabelings.

```{.python .execute}
#| label: code-fdr-flag-005
flagged = p_values <= 0.05
print("p <= 0.05:", flagged.sum(), "genes, of which null:", (flagged & is_null).sum())
print("p <= 0.05 / m (Bonferroni):", (p_values <= 0.05 / len(p_values)).sum(), "genes")
```

Flagging every gene with $p_j\le0.05$ selects 634 genes, and 225 of them are null, matching the $4500\times0.05=225$ predicted by @eq-fdr-expected-false. About 35% of this list would be wasted follow-up.

@fig-fdr-screen shows two standard pictures of such a screen. The *volcano plot* places each gene by its effect size and by $-\log_{10}p_j$, so genes with strong evidence rise toward the top. The *p-value histogram* is the more useful diagnostic for multiple testing [@holmeshuber2019]. Because null p-values are uniform, the null genes spread evenly across the bars and form a flat floor about $m_0$ times the bar width high, while genes that truly changed pile up near zero. The first bar mixes real responses with its share of the floor, and the floor's height says roughly how many small p-values to distrust.

```{.python .execute}
#| label: fig-fdr-screen
#| code-fold: true
#| code-summary: "Draw the volcano plot and p-value histogram"
#| fig-width: 6.4
#| fig-height: 3
#| fig-cap: "The simulated screen of 5000 genes, colored by the truth that real data hide. (A) Volcano plot, changed genes as orange triangles; no permutation p-value with six cultures per group falls below 2/924, which produces the flat top row. (B) Stacked p-value histogram, bar width 0.05, with the expected null count per bar dashed. Null genes form a flat floor and changed genes a spike near zero, so the first bar contains both."
#| fig-alt: "Left: scatter of effect size against minus log10 p-value; grey null points cluster near zero effect, orange triangles spread outward and upward to a flat ceiling. Right: p-value histogram with a flat grey floor, a tall orange spike in the first bar, and a dashed line at the floor height."
null_color, changed_color = cf.RULE, cf.SERIES[1]
fig, (ax_volcano, ax_hist) = plt.subplots(1, 2)
ax_volcano.scatter(delta_hat[is_null], -np.log10(p_values[is_null]), s=4, color=null_color)
ax_volcano.scatter(delta_hat[~is_null], -np.log10(p_values[~is_null]), s=6,
                   color=changed_color, marker="^")
ax_volcano.set(xlabel=r"effect size $\widehat\Delta_j$ (log$_2$ units)",
               ylabel=r"$-\log_{10} p_j$")
ax_hist.hist([p_values[is_null], p_values[~is_null]], bins=np.arange(0, 1.0001, 0.05),
             stacked=True, color=[null_color, changed_color], label=["null", "changed"])
ax_hist.axhline(is_null.sum() * 0.05, color=cf.INK, linestyle="--", linewidth=1)
ax_hist.set(xlabel=r"p-value $p_j$", ylabel="number of genes", xlim=(0, 1))
ax_hist.legend(loc="upper right")
ax_volcano.set_title("A", loc="left", fontweight="bold", fontsize=11)
ax_hist.set_title("B", loc="left", fontweight="bold", fontsize=11)
plt.show()
```

## Controlling the chance of any false discovery

A first repair is a threshold strict enough that, with high probability, the list contains no null gene at all. Write $V$ and $R$ for the numbers of false and total selections. The *family-wise error rate* (FWER) is $P(V\ge1)$, and dividing the threshold by the number of tests controls it.

::: {#prp-bonferroni}
## Bonferroni threshold
If every null p-value is calibrated and we select the genes with $p_j\le\alpha/m$, then $P(V\ge1)\le\alpha$, whatever the dependence between genes.
:::

::: {.proof}
A count of at least one is at least one, so $V\ge\mathbf 1[V\ge1]$. Averaging both sides over repeated experiments gives $P(V\ge1)\le E[V]$, and @eq-fdr-expected-false with $t=\alpha/m$ gives $E[V]\le m_0\alpha/m\le\alpha$.
:::

For 20,000 genes at $\alpha=0.05$ the threshold is $2.5\times10^{-6}$. With six cultures per condition the smallest possible permutation p-value is $2/924\approx0.0022$, nearly a thousand times larger, so the Bonferroni list is empty however strongly a gene responds. In our screen the threshold is $10^{-5}$, and Bonferroni selects nothing although 500 genes changed, many of them by a factor of two or more.

FWER suits settings where one false claim is costly, such as declaring a forensic DNA match. A discovery screen feeds further experiments that retest each selected gene, so a few false leads are tolerable if they are a small share of the list. We should control a fraction instead.

## False discovery proportion and rate

::: {#def-false-discovery-rate}
## False discovery proportion and false discovery rate
A *false discovery* is a selected gene whose null hypothesis is true. If a selection rule selects $R$ genes, $V$ of them false discoveries, the *false discovery proportion* is
$$
\mathrm{FDP}=\frac{V}{\max(R,1)},
$$
and the *false discovery rate* is its expectation, $\mathrm{FDR}=E[\mathrm{FDP}]$.
:::

The $\max(R,1)$ only handles an empty list, whose FDP is defined to be $0$. FDP is the fraction of wasted follow-ups in one particular list; it belongs to one experiment and cannot be computed from real data because $V$ is hidden. FDR is the average of FDP over repetitions of the experiment, and a selection rule can guarantee something about it. Thresholding our screen at $0.05$ gave FDP $=225/634\approx0.35$, far too high for a follow-up list.

## From the histogram to Benjamini–Hochberg

Although we cannot compute FDP, @eq-fdr-expected-false lets us estimate it. Select every gene with $p_j\le t$, and let $R(t)$ be the observed number selected. On average at most $m_0t\le mt$ of them are null, which suggests the estimate
$$
\widehat{\mathrm{FDP}}(t)=\frac{m\,t}{R(t)} .
$$
In the histogram this is the floor's share of everything to the left of $t$, with the floor drawn as if every gene were null. Using $m$ for the unknown $m_0$ overstates the floor and so errs toward caution.

We want the largest threshold, and hence the longest list, whose estimated false fraction is at most a target $q$. Because $R(t)$ only changes when $t$ passes an observed p-value, it is enough to try thresholds equal to the observed p-values. Sort them as $p_{(1)}\le p_{(2)}\le\cdots\le p_{(m)}$. At $t=p_{(i)}$ exactly $i$ genes are selected (assuming no ties), and $\widehat{\mathrm{FDP}}\le q$ becomes
$$
\frac{m\,p_{(i)}}{i}\le q
\qquad\Longleftrightarrow\qquad
p_{(i)}\le\frac{i\,q}{m}.
$$ {#eq-fdr-bh-line}

The right side of @eq-fdr-bh-line is a line through the origin with slope $q/m$ in the rank $i$. Choosing the largest rank on or below it is the procedure of @benjamini1995.

::: {.algorithm}
**Algorithm: Benjamini–Hochberg selection.**
Given p-values $p_1,\ldots,p_m$ and a target $q$, sort them so that $p_{(1)}\le\cdots\le p_{(m)}$. Find the largest rank $k$ with $p_{(k)}\le kq/m$. If no rank qualifies, select nothing; otherwise select the genes with the $k$ smallest p-values.
:::

Why the largest rank, rather than stopping at the first failure? The estimate $mt/R(t)$ does not increase steadily with $t$. A gene whose p-value sits just above its line may be followed closely by several more, and each extra selection divides the same estimated floor among more genes. BH judges the list as a whole, so a gene that fails its own inequality is still selected when a later rank passes. Because it searches for the last crossing, BH is called a *step-up* procedure.

```{.python .execute}
#| label: code-fdr-bh
def benjamini_hochberg(p_values: np.ndarray, q: float) -> np.ndarray:
    """Boolean array marking the genes that Benjamini-Hochberg selects at target q."""
    m = len(p_values)
    order = np.argsort(p_values)  # gene indices from smallest p-value to largest
    k = 0
    for i in range(1, m + 1):  # i is the rank, starting at 1 as in the text
        if p_values[order[i - 1]] <= i * q / m:
            k = i  # keep the largest rank that passes
    selected = np.zeros(m, dtype=bool)
    selected[order[:k]] = True
    return selected
```

Sorting dominates the cost, $\mathcal O(m\log m)$ comparisons, which is negligible beside computing the p-values.

## Benjamini–Hochberg on ten genes

::: {#exm-bh-ten-genes}
## Selecting from ten genes at $q=0.05$
Suppose a screen of $m=10$ genes returned the values below, sorted by p-value. The line $iq/m$ rises by $0.005$ per rank.

| rank $i$ | gene | $\widehat\Delta$ | $p_{(i)}$ | $iq/m$ | pass? |
|---:|---|---:|---:|---:|---|
| 1 | g07 | 0.15 | 0.0004 | 0.005 | yes |
| 2 | g02 | $-1.30$ | 0.0042 | 0.010 | yes |
| 3 | g09 | 0.90 | 0.0180 | 0.015 | no |
| 4 | g04 | 1.10 | 0.0190 | 0.020 | yes |
| 5 | g01 | $-0.60$ | 0.0320 | 0.025 | no |
| 6 | g10 | 0.40 | 0.0410 | 0.030 | no |
| 7 | g05 | 2.10 | 0.0470 | 0.035 | no |
| 8 | g03 | 0.20 | 0.2200 | 0.040 | no |
| 9 | g08 | $-0.10$ | 0.5300 | 0.045 | no |
| 10 | g06 | 0.05 | 0.8600 | 0.050 | no |

Compare ranks 3 and 4. Stopping after three genes uses $t=0.018$, so the estimated number of null genes on the list is $mt=10\times0.018=0.18$ and the estimated false fraction is $0.18/3=0.060$, above the target. A fourth gene raises the threshold only to $0.019$, and $0.19/4=0.0475\le0.05$. Rank 4 passes and no later rank does, so $k=4$ and BH selects g07, g02, g09, and g04, including g09, which failed its own inequality.

The same table gives three different lists: thresholding each gene at $0.05$ selects seven genes, Bonferroni at $0.05/10=0.005$ selects two, and BH at $q=0.05$ selects four.
:::

The function agrees with the hand calculation.

```{.python .execute}
#| label: code-fdr-ten-genes
genes = np.array(["g07", "g02", "g09", "g04", "g01", "g10", "g05", "g03", "g08", "g06"])
toy_p = np.array([0.0004, 0.0042, 0.0180, 0.0190, 0.0320,
                  0.0410, 0.0470, 0.2200, 0.5300, 0.8600])
genes[benjamini_hochberg(toy_p, 0.05)]
```

@fig-fdr-ranked-cutoff draws the rules as lines over the sorted p-values. Bonferroni is a flat line at the height of the first BH threshold and the unadjusted rule a flat line at $0.05$. BH's line rises with rank, so it becomes more lenient as more genes lie below it. Panel B shows the step-up rule at work on the simulated screen.

```{.python .execute}
#| label: fig-fdr-ranked-cutoff
#| code-fold: true
#| code-summary: "Draw sorted p-values against the thresholds"
#| fig-width: 6.4
#| fig-height: 3
#| fig-cap: "Sorted p-values against rank. (A) The ten-gene table with the BH line iq/m at q = 0.05 (solid) and the Bonferroni and unadjusted thresholds (dashed); filled points are selected. Rank 3 lies above its line, but rank 4, ringed in red, is the last rank on or below it. (B) The 450 smallest p-values of the simulated screen and the BH line at q = 0.1. The smallest attainable value, 2/924, lies above the line for about the first hundred ranks, yet those genes are selected because later ranks pass; the red line marks the cutoff."
#| fig-alt: "Left: ten points rising with rank, a sloped solid line and two horizontal dashed lines; the first four points are filled and the fourth is circled in red. Right: a step-shaped curve of sorted p-values crossing a straight sloped line, with a red vertical line at the last crossing."
fig, (ax_toy, ax_screen) = plt.subplots(1, 2)
ranks = np.arange(1, 11)
toy_selected = benjamini_hochberg(toy_p, 0.05)
ax_toy.plot(ranks, ranks * 0.05 / 10, color=cf.SERIES[0], label="BH: $iq/m$")
ax_toy.hlines([0.05, 0.005], 0.5, 10.5, color=cf.RULE, linestyle="--", linewidth=1)
ax_toy.text(10.7, 0.05, "unadjusted", va="center", fontsize=8)
ax_toy.text(10.7, 0.005, "Bonferroni", va="center", fontsize=8)
ax_toy.plot(ranks[toy_selected], toy_p[toy_selected], "o", color=cf.INK)
ax_toy.plot(ranks[~toy_selected], toy_p[~toy_selected], "o", markerfacecolor="none",
            markeredgecolor=cf.INK, markeredgewidth=1)
ax_toy.plot(4, toy_p[3], "o", markersize=12, markerfacecolor="none",
            markeredgecolor=cf.EMPHASIS, markeredgewidth=1.5)
ax_toy.set(xlabel="rank $i$", ylabel="p-value", xticks=ranks, ylim=(-0.003, 0.1),
           xlim=(0.5, 13))
ax_toy.legend(loc="upper left")

screen_ranks = np.arange(1, 451)
ax_screen.step(screen_ranks, np.sort(p_values)[:450], where="post", color=cf.INK,
               label="sorted $p_{(i)}$")
ax_screen.plot(screen_ranks, screen_ranks * 0.1 / len(p_values), color=cf.SERIES[0],
               label="BH: $iq/m$")
ax_screen.axvline(benjamini_hochberg(p_values, 0.1).sum(), color=cf.EMPHASIS)
ax_screen.set(xlabel="rank $i$", ylabel="p-value", ylim=(0, 0.012))
ax_screen.legend(loc="upper left")
ax_toy.set_title("A", loc="left", fontweight="bold", fontsize=11)
ax_screen.set_title("B", loc="left", fontweight="bold", fontsize=11)
plt.show()
```

```{.python .execute}
#| label: code-fdr-screen-bh
selected = benjamini_hochberg(p_values, 0.1)
print("BH at q = 0.1:", selected.sum(), "genes, of which null:", (selected & is_null).sum())
```

On the screen, BH at $q=0.1$ selects 328 genes, 27 of them null, an FDP of $27/328\approx0.08$. Compared with the unadjusted list, it keeps 301 of the 409 real responses while cutting the false leads from 225 to 27.

## Why Benjamini–Hochberg controls the FDR

BH was derived from an estimate, so we still need to show that it keeps its promise.

::: {#thm-benjamini-hochberg}
## Benjamini–Hochberg controls the FDR
Suppose the null p-values are calibrated and independent of each other and of the non-null p-values. Then BH at target $q$ satisfies
$$
\mathrm{FDR}\le\frac{m_0}{m}\,q\le q .
$$
:::

The argument rests on a cancellation. When BH makes $R$ selections, it selects exactly the genes with $p_j\le Rq/m$: rank $R$ passes, and if a gene beyond rank $R$ also had $p_j\le Rq/m$, rank $R+1$ would pass too. A null gene $j$ therefore contributes $\mathbf 1[p_j\le Rq/m]/R$ to the FDP. If $R$ were a fixed number $k$, that contribution would average to $(kq/m)/k=q/m$. A longer list makes each null gene more likely to be selected, in exact proportion to how much it dilutes each false selection.

::: {.proof}
Sketch. Write $\mathrm{FDP}=\sum_{j\in\mathcal H_0}\mathbf 1[p_j\le Rq/m]/\max(R,1)$. Since $R$ depends on $p_j$, it is not fixed. Let $R_j\ge1$ be the number of genes BH would select if $p_j$ were replaced by $0$, which depends only on the other p-values. If gene $j$ is selected, lowering $p_j$ to $0$ leaves the same genes in the first $R$ ranks and the same values above them, so $R_j=R$; conversely, if $p_j\le R_jq/m$, restoring $p_j$ from $0$ keeps the first $R_j$ ranks passing and cannot make a higher rank pass, so gene $j$ is selected with $R=R_j$. Gene $j$'s term therefore equals $\mathbf 1[p_j\le R_jq/m]/R_j$. Hold the other p-values fixed, so that $R_j=k$. By independence and calibration, $P(p_j\le kq/m)\le kq/m$, so the term averages to at most $q/m$ for every $k$. Summing over the $m_0$ null genes gives $\mathrm{FDR}\le m_0q/m$.
:::

The theorem bounds an average over experiments. The next chunk runs the simulated screen 200 times with $m=2000$ genes, 200 of them changed, and records the realized FDP of each list.

```{.python .execute}
#| label: code-fdr-repeat
rng = np.random.default_rng(1995)
fdp_bh, fdp_unadjusted = [], []
for run in range(200):
    X_run, null_run = simulate_screen(rng, m=2000, m_changed=200)
    _, p_run = permutation_p_values(X_run, n_control=6)
    for selected, record in [(benjamini_hochberg(p_run, 0.1), fdp_bh),
                             (p_run <= 0.05, fdp_unadjusted)]:
        record.append((selected & null_run).sum() / max(selected.sum(), 1))
print(f"BH: mean FDP {np.mean(fdp_bh):.3f}, range {min(fdp_bh):.3f} to {max(fdp_bh):.3f}")
print(f"p <= 0.05: mean FDP {np.mean(fdp_unadjusted):.3f}")
```

```{.python .execute}
#| label: fig-fdr-repeated-fdp
#| code-fold: true
#| code-summary: "Draw the realized FDP across repeated experiments"
#| fig-cap: "Realized false discovery proportion in 200 simulated experiments, each with 2000 genes of which 200 changed. Lists from BH at q = 0.1 (blue) vary from run to run but average below q (red line); thresholding each gene at 0.05 (orange) gives lists about a third false. FDR control describes the average over experiments, not any single list."
#| fig-alt: "Two histograms of false discovery proportion. The blue BH histogram spans about 0.02 to 0.16 and is centered just below a red vertical line at 0.1; the orange unadjusted histogram is centered near 0.35."
fig, ax = plt.subplots()
bins = np.arange(0, 0.5001, 0.0125)
ax.hist(fdp_bh, bins=bins, color=cf.SERIES[0], label="BH, $q=0.1$")
ax.hist(fdp_unadjusted, bins=bins, color=cf.SERIES[1], label=r"$p_j\leq 0.05$")
ax.axvline(0.1, color=cf.EMPHASIS)
ax.set(xlabel="realized FDP of one experiment's list", ylabel="number of experiments")
ax.legend(loc="upper right")
plt.show()
```

The BH lists average an FDP of 0.084, a little below the bound $m_0q/m=0.09$, partly because permutation p-values are slightly conservative. Individual runs range from 0.025 to 0.16. A lab that runs the experiment once gets one draw from the blue histogram, not its mean.

::: {.callout-note title="Key idea"}
BH selects the longest list whose estimated false fraction $m\,p_{(k)}/k$ is at most $q$. It guarantees that false discoveries make up at most a fraction $q$ of the list on average over repeated experiments.
:::

## What a controlled list does and does not say

The FDR describes a list, not a gene. Selecting four genes at $q=0.05$ in @exm-bh-ten-genes does not mean that $0.2$ of them are false, nor that each has a 5% chance of being false: the evidence against the null is far stronger for g07 than for g09, which entered only because g04 followed it. An unselected gene is not shown to be unchanged either. Gene g05 missed the list, yet its estimate $\widehat\Delta=2.10$ is the largest in the table, and more cultures might have put it on the list.

A small p-value also does not mean a large change. Gene g07 has the strongest evidence but an estimated change of $0.15$ $\log_2$ units, about 11%; its cultures vary so little that a small shift is easy to detect. Such genes sit high and near the center of the volcano plot in @fig-fdr-screen. Follow-up lists are therefore usually filtered on both axes, requiring evidence that the change is real and a change large enough to matter for the biological question.

Real screens need care with the theorem's assumptions. BH takes each p-value at face value. If cultures were sequenced in batches confounded with the drug, or technical replicates were relabeled as if they were separate cultures (pseudoreplication), null p-values come out too small across the board, the histogram's floor rises toward zero, and the guarantee is lost; no list-level procedure repairs a miscalibrated test. Independence is the second assumption. Co-regulated genes rise and fall together, so their p-values are dependent. BH still controls the FDR under a form of positive dependence that covers many such cases, and a more conservative variant works under any dependence [@benjamini2001], but strong dependence makes the realized FDP vary more from run to run.

Finally, the design must be able to produce small p-values. With three cultures per condition, the design of Lesson 17, every permutation p-value is at least $0.1$. BH at $q=0.05$ then selects nothing for any number of genes, because no p-value can meet even the largest threshold $mq/m=q$. This is one reason expression studies use more replicates, or model-based tests that share information across genes.

## Exercises

::: {#exr-fdr-by-hand}
Apply BH at $q=0.1$ to the p-values $0.001$, $0.006$, $0.021$, $0.055$, $0.060$, $0.120$, $0.400$, $0.800$. Identify the selected genes, point out any gene selected despite failing its own inequality, and compare with Bonferroni at $0.1$.
:::

::: {#exr-fdr-equal-p}
Suppose all $m$ p-values equal the same value $p_0$. Show that BH at target $q$ selects every gene if $p_0\le q$ and none otherwise. Compare with Bonferroni, and check that both reduce to the single test $p\le q$ when $m=1$.
:::

::: {#exr-fdr-all-null}
When every null hypothesis is true, explain why the FDP of any list is either $0$ or $1$, and deduce that the FDR then equals the family-wise error rate $P(V\ge1)$.
:::

::: {#exr-fdr-interpret}
A collaborator reports: "BH at $q=0.05$ gave 200 genes, so exactly 10 of them are false, and the 4800 genes not on the list are unaffected by the treatment." Explain what is wrong with each part of the claim and state what the BH result does support.
:::

::: {#refs}
**References**
:::
