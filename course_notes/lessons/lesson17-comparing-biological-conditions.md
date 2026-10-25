---
published: true
title: "Comparing Biological Conditions"
subtitle: "Lesson 17 · 2027-03-24 · Week 11"
nocite: |
  @holmeshuber2019
---

## Core question

Does a treatment change one gene's expression by more than independent cultures already differ from one another?

```{.python .execute}
#| label: setup
#| include: false
import sys
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## One gene in treated and control cultures

Suppose a laboratory wants to know whether a drug switches on a particular stress-response gene in a human cell line. They grow six cultures separately from the same frozen stock, add the drug to three, and leave three as controls. RNA from each culture is then sequenced, and the reads are counted, normalized, and log-transformed as in Lesson 16. The result is an expression matrix $X$ with one row per culture and one column per gene, on the $\log_2(1+\text{CPM})$ scale, together with a metadata table recording each culture's condition.

Lesson 16 used PCA to summarize all genes at once. Here the drug is thought to act on one named gene, so we keep that gene's column and ask about it alone. In this example the column holds the following toy values.

| culture | condition | expression |
|---|---|---:|
| c1 | control | 2 |
| c2 | control | 3 |
| c3 | control | 4 |
| t1 | treated | 5 |
| t2 | treated | 6 |
| t3 | treated | 7 |

```{.python .execute}
#| label: code-conditions-data
culture_names = ["c1", "c2", "c3", "t1", "t2", "t3"]
culture_values = [2, 3, 4, 5, 6, 7]
treated_cultures = {3, 4, 5}  # list positions of t1, t2, t3
```

Panel A of @fig-conditions-one-gene plots these six numbers, one point per culture. Every treated culture lies above every control, and the treated average exceeds the control average by 3 units. Yet c1 and c3 were treated identically and still differ by 2, and six new cultures would give six different numbers. We want to know whether a 3-unit difference is larger than culture-to-culture variation alone tends to produce, and to answer with a number we can compute by hand.

```{.python .execute}
#| label: fig-conditions-one-gene
#| code-fold: true
#| code-summary: "Plot the six cultures"
#| fig-width: 6
#| fig-height: 3.2
#| fig-cap: "One gene's expression in six cultures; blue circles are labeled control, orange triangles treated, and grey bars mark the labeled means. (A) The experiment's labels: the treated mean is 3 units higher, while the controls alone span 2 units. (B, C) Two relabelings of the same six values used later in the permutation test; moving one culture across the line shrinks the difference."
#| fig-alt: "Three panels of the same six values 2 to 7 named c1 to c3 and t1 to t3. Panel A labels t1, t2, t3 treated with difference 3; panel B labels c3, t1, t3 treated with difference 5/3; panel C labels c1, c2, c3 treated with difference minus 3."
labelings = [
    ("A  observed", {3, 4, 5}, "3"),
    ("B  c3, t1, t3 treated", {2, 3, 5}, "5/3"),
    ("C  reversed", {0, 1, 2}, "-3"),
]
fig, axes = plt.subplots(1, 3, sharey=True)
for ax, (title, treated_set, difference) in zip(axes, labelings):
    groups = {0: [], 1: []}
    for i, (name, value) in enumerate(zip(culture_names, culture_values)):
        position = 1 if i in treated_set else 0
        groups[position].append(value)
        ax.plot(position, value, "^" if position else "o",
                color="C1" if position else "C0", markersize=7)
        ax.text(position + 0.1, value, name, va="center", fontsize=9)
    for position, group in groups.items():
        mean = sum(group) / len(group)
        ax.plot([position - 0.35, position - 0.1], [mean, mean], color=cf.RULE)
    ax.set(xticks=[0, 1], xticklabels=["control", "treated"], xlim=(-0.5, 1.5),
           ylim=(1, 8), title=title + "\n" + rf"$\widehat\Delta(a)={difference}$")
axes[0].set_ylabel(r"Expression, $\log_2(1+\mathrm{CPM})$")
plt.show()
```

## Biological and technical replicates

Why grow three control cultures rather than measure one control culture three times? The question concerns cells of this line in general, not one particular dish. Each culture has its own small differences in cell density, growth phase, and handling, which change its RNA content, and only several independently grown cultures per condition show how large those differences are [@holmeshuber2019].

::: {#def-replicates}
## Biological and technical replicates
An *experimental unit* is the smallest entity to which a condition is independently applied; here it is one culture. *Biological replicates* are distinct experimental units in the same condition. *Technical replicates* are repeated measurements of one unit, such as sequencing one culture's RNA several times.
:::

@fig-conditions-culture-replication contrasts the two designs for the controls. Assaying c1, c2, and c3 once each gives $(2,3,4)$; assaying only c1 three times might give $(1.9,2.0,2.1)$. Both designs produce three rows of data, but the tight spread of the second reflects only the assay's repeatability and says nothing about how far a second culture would fall from c1. Technical replicates help assess an instrument; they cannot show whether a drug effect extends beyond the dishes we grew. For comparing conditions, the count that matters is the number of experimental units: $n_C=3$ control and $n_T=3$ treated cultures, with $n=n_C+n_T=6$.

```{.python .execute}
#| label: fig-conditions-culture-replication
#| code-fold: true
#| code-summary: "Draw the two replication designs"
#| fig-width: 6
#| fig-height: 3.4
#| fig-cap: "Three control cultures assayed once each, compared with one culture assayed three times; boxes give illustrative expression values. Three recorded values carry different biological evidence when they come from different numbers of independent cultures."
#| fig-alt: "One panel shows separate cultures c1, c2, c3 leading to measurements two, three, four; another shows only c1 leading to three measurements 1.9, 2.0, 2.1."
fig, ax = plt.subplots()
cf.replicate_sampling(
    ax, biological_values=[2, 3, 4], technical_values=[1.9, 2.0, 2.1],
    culture_labels=["c1", "c2", "c3"],
)
plt.show()
```

## Effect size and variability

Fix the gene and write $x_i$ for culture $i$'s entry in its column, with condition label $g_i\in\{C,T\}$. The condition means are
$$
\bar x_C=\frac1{n_C}\sum_{i:\,g_i=C}x_i,
\qquad
\bar x_T=\frac1{n_T}\sum_{i:\,g_i=T}x_i,
$$
and the observed *effect size* is their difference,
$$
\widehat\Delta=\bar x_T-\bar x_C .
$$
For our cultures, $\bar x_C=(2+3+4)/3=3$ and $\bar x_T=(5+6+7)/3=6$, so $\widehat\Delta=3$.

$\widehat\Delta$ estimates an unknown quantity. Imagine every culture this lab could grow and assay under each condition, and let $\mu_C$ and $\mu_T$ be the average measurements over the possible control and treated cultures. The biological question concerns $\Delta=\mu_T-\mu_C$. Our six cultures are a small sample from those populations, and new cultures would give a different $\widehat\Delta$.

The sign says which way the gene moved and the magnitude says by how much, on the log scale. A difference of 3 in $\log_2$ units corresponds to a factor of $2^3=8$: averaging the $\log_2$ values within each condition and exponentiating the difference gives an eightfold ratio of typical $1+\text{CPM}$ values, which is close to an eightfold increase for a gene expressed well above one count per million. A difference of 3 on the raw CPM scale would be negligible for most genes, so an effect is interpretable only together with its scale.

How far would $\widehat\Delta$ move between repetitions? That depends on how much cultures within one condition vary. Measure the spread in condition $C$ by the sample variance
$$
s_C^2=\frac{1}{n_C-1}\sum_{i:\,g_i=C}(x_i-\bar x_C)^2 ,
$$
and define $s_T^2$ in the same way. Its square root, the *standard deviation* $s_C$, is the typical distance of a culture from its condition mean; dividing by $n_C-1$ rather than $n_C$ is the convention when the mean is estimated from the same values. For the controls, $s_C^2=\bigl((2-3)^2+0^2+(4-3)^2\bigr)/2=1$, and the treated cultures also give $s_T^2=1$, so the effect is three standard deviations. Had the controls been $(0,3,6)$, with $s_C=3$, the same $\widehat\Delta$ would be far less convincing. We need a comparison that weighs the effect against the spread, using only the six numbers and the way the experiment assigned its labels.

## A null hypothesis and relabeled cultures

Suppose the lab chose the three treated cultures at random, say by drawing names from a hat, and suppose the drug has no effect on this gene. Each culture would then have given the same measurement whether or not it received the drug, since its RNA content was set by its own growth history. The six numbers $(2,3,4,5,6,7)$ would be fixed; only their labels would be random.

::: {#def-permutation-null}
## Null hypothesis and relabelings
The *null hypothesis* $H_0$ states that the treatment does not change this gene's measured expression in any culture. A *relabeling* $a$ is a choice of $n_T$ of the $n$ cultures to call treated. The set $\mathcal A$ contains every relabeling with the experiment's group sizes, so $|\mathcal A|=\binom{n}{n_T}$.
:::

Under $H_0$ with randomized treatment, the observed labeling is a uniformly random element of $\mathcal A$; here each of the $\binom63=20$ relabelings has probability $1/20$. Recomputing the effect size $\widehat\Delta(a)$ for each relabeling $a$ gives the 20 differences the hat alone could have produced if the drug did nothing. If the observed $\widehat\Delta$ lies far out among them, either the drug changed the gene or an unusual draw occurred.

Some conditions, such as tumor versus normal tissue, cannot be randomized. The null hypothesis then says that all six values were generated the same way, with no condition difference. In that case which three values carry the "treated" label tells us nothing about the values themselves, so every split is equally plausible; the labels are called *exchangeable*. The computation is unchanged, but its justification moves from the design to that assumption.

To measure how far out a difference lies, use the two-sided statistic $S(a)=|\widehat\Delta(a)|$, since a drug could lower expression as well as raise it.

::: {#def-permutation-p-value}
## Permutation p-value
**Input:** one gene's measurements $x_1,\ldots,x_n$ and the observed labeling $a_{\mathrm{obs}}\in\mathcal A$.

**Output:** the exact permutation p-value
$$
p=\frac{\#\{a\in\mathcal A:\ S(a)\ge S_{\mathrm{obs}}\}}{|\mathcal A|},
\qquad S_{\mathrm{obs}}=S(a_{\mathrm{obs}}).
$$
:::

In words, $p$ is the fraction of allowed relabelings whose absolute difference is at least the observed one. The observed labeling belongs to $\mathcal A$, so $p\ge 1/|\mathcal A|$, and ties count toward the numerator. Because each relabeling has probability $1/|\mathcal A|$ under $H_0$, this fraction is a probability: $p$ is the chance that the hat alone produces a difference at least as large as the observed one.

::: {.algorithm}
**Algorithm: exact permutation test for one gene.**
Given measurements $x_1,\ldots,x_n$ and the observed set of treated cultures, compute $S_{\mathrm{obs}}=|\bar x_T-\bar x_C|$. For each of the $\binom{n}{n_T}$ subsets $a$ of size $n_T$, compute $\widehat\Delta(a)$ as the mean over $a$ minus the mean over the remaining cultures, and add one to a counter when $|\widehat\Delta(a)|\ge S_{\mathrm{obs}}$. Return the counter divided by $\binom{n}{n_T}$.
:::

```{.python .execute}
#| label: code-conditions-permutation
def mean_difference(values: list[float], treated: set[int]) -> float:
    """Treated mean minus control mean under one labeling of the cultures."""
    treated_values = [values[i] for i in range(len(values)) if i in treated]
    control_values = [values[i] for i in range(len(values)) if i not in treated]
    return (sum(treated_values) / len(treated_values)
            - sum(control_values) / len(control_values))


def permutation_p_value(values: list[float],
                        treated: set[int]) -> tuple[float, list[float]]:
    """Exact two-sided permutation p-value and the relabeled differences."""
    observed = abs(mean_difference(values, treated))
    relabeled_differences = []
    at_least_as_large = 0
    for assignment in combinations(range(len(values)), len(treated)):
        difference = mean_difference(values, set(assignment))
        relabeled_differences.append(difference)
        # the small tolerance counts exact ties despite floating-point rounding
        if abs(difference) >= observed - 1e-9:
            at_least_as_large += 1
    return at_least_as_large / len(relabeled_differences), relabeled_differences


p_value, relabeled_differences = permutation_p_value(culture_values, treated_cultures)
print("observed difference:", mean_difference(culture_values, treated_cultures))
print("relabelings:", len(relabeled_differences), "  p-value:", p_value)
```

## Worked example: all twenty relabelings

::: {#exm-conditions-permutation}
## Twenty relabelings of six cultures
The six values sum to 27. If a relabeling's treated cultures have values summing to $t$, its control cultures sum to $27-t$, so
$$
\widehat\Delta(a)=\frac{t}{3}-\frac{27-t}{3}=\frac{2t-27}{3}.
$$
For example, calling c3, t1, and t3 treated gives $t=4+5+7=16$ and $\widehat\Delta(a)=5/3$, so moving one culture across the line already pulls the difference well below 3 (@fig-conditions-one-gene, panel B). The condition $|\widehat\Delta(a)|\ge3$ holds exactly when $t\ge18$ or $t\le9$. The largest possible treated sum is $5+6+7=18$, reached only by the observed labeling, and the smallest is $2+3+4=9$, reached only by the reversed labeling that calls c1, c2, and c3 treated. Two of the twenty relabelings qualify, so
$$
p=\frac{2}{20}=0.1 .
$$
:::

The code reports the same 20 relabelings and $p=0.1$. In @fig-conditions-permutation-null, most relabelings mix high and low values and give differences near zero; only the observed labeling and its reversal reach the tails.

```{.python .execute}
#| label: fig-conditions-permutation-null
#| code-fold: true
#| code-summary: "Plot the twenty relabeled differences"
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "The twenty values of the relabeled difference for the six-culture experiment, one dot per relabeling. Red dots have absolute difference at least the observed 3; they are the observed labeling and its reversal, giving p = 2/20."
#| fig-alt: "Stacked dots at ten values from minus 3 to 3 in steps of two thirds, with stack heights 1, 1, 2, 3, 3, 3, 3, 2, 1, 1. The single dots at minus 3 and plus 3 are red and labeled."
stack_heights = {}
fig, ax = plt.subplots()
for difference in sorted(relabeled_differences):
    key = round(difference, 6)
    stack_heights[key] = stack_heights.get(key, 0) + 1
    extreme = abs(difference) >= 3 - 1e-9
    ax.plot(difference, stack_heights[key], "o", markersize=9,
            color=cf.EMPHASIS if extreme else "C0")
ax.annotate("observed:\nt1, t2, t3 treated", (3, 1), xytext=(2.15, 2.6),
            ha="center", fontsize=9, arrowprops={"arrowstyle": "-", "color": cf.RULE})
ax.annotate("reversed:\nc1, c2, c3 treated", (-3, 1), xytext=(-2.15, 2.6),
            ha="center", fontsize=9, arrowprops={"arrowstyle": "-", "color": cf.RULE})
ax.set(xlabel=r"Relabeled difference $\widehat\Delta(a)$",
       ylabel="Relabelings", xticks=range(-3, 4), yticks=[1, 2, 3],
       xlim=(-3.5, 3.5), ylim=(0.4, 3.6))
plt.show()
```

The observed labeling is as extreme as six numbers allow, yet $p=0.1$. With equal group sizes, reversing every label is also an allowed relabeling and turns $\widehat\Delta$ into $-\widehat\Delta$, so the two-sided count always includes the observed labeling and its reversal. With $n_C=n_T$ cultures per condition, the smallest attainable p-value is therefore
$$
p_{\min}=\frac{2}{\binom{2n_C}{n_C}},
$$
which is $2/20=0.1$ for $n_C=3$, $2/70\approx0.029$ for $n_C=4$, and $2/252\approx0.008$ for $n_C=5$. With three cultures per condition, no data can give this test a p-value below 0.1. That limit belongs to the experimental design and can be computed before any RNA is sequenced.

::: {.callout-note title="Key idea"}
A permutation p-value ranks the observed difference among the differences produced by relabeling whole experimental units. The number of independent units, not the number of measurements, determines how finely that ranking can resolve evidence.
:::

## Why the permutation p-value is calibrated

A p-value is useful only if small values are rare when nothing is happening; a test that reported $p\le0.05$ in a third of experiments with an inactive drug would flag inactive drugs constantly. For the permutation p-value, this property follows from a counting argument.

::: {#thm-permutation-calibration}
## Permutation p-values are calibrated
Suppose $H_0$ holds and the observed labeling was drawn uniformly at random from $\mathcal A$. Then for every $\alpha\in[0,1]$,
$$
P(p\le\alpha)\le\alpha .
$$
:::

::: {.proof}
Under $H_0$ the measurements do not depend on which labeling was drawn, so the statistics $S(a)$, $a\in\mathcal A$, form one fixed list. The p-value that labeling $a$ would report is
$$
p(a)=\frac{\#\{a'\in\mathcal A:\ S(a')\ge S(a)\}}{|\mathcal A|},
$$
and the reported p-value is $p(a_{\mathrm{obs}})$. Collect the labelings that would report $p\le\alpha$ in $\mathcal R=\{a\in\mathcal A: p(a)\le\alpha\}$; if $\mathcal R$ is empty there is nothing to prove. Otherwise let $b\in\mathcal R$ have the smallest statistic. Every member of $\mathcal R$ has a statistic at least $S(b)$, so
$$
|\mathcal R|\le\#\{a':\ S(a')\ge S(b)\}=p(b)\,|\mathcal A|\le\alpha\,|\mathcal A| .
$$
Since $a_{\mathrm{obs}}$ is uniform on $\mathcal A$,
$$
P(p\le\alpha)=\frac{|\mathcal R|}{|\mathcal A|}\le\alpha.
$$
:::

A p-value with this property is called *calibrated*. The randomness here is in the draw of labels, with the measurements held fixed: with an inactive drug, a rule that reports a change whenever $p\le0.05$ does so in at most 5% of experiments. The bound can be loose, since in the six-culture design no labeling gives $p\le0.05$ at all.

In the example, $p=0.1$ answers one question: how often the hat alone produces a difference this large. The estimated change, roughly eightfold, is a separate number, and $p=0.1$ says only that three cultures per condition cannot separate that change from a lucky draw at the 0.05 level. It is not the probability that the drug has no effect, which would require a model of how often drugs affect genes.

Exact enumeration costs $\mathcal O(n)$ work for each of the $\binom{n}{n_T}$ relabelings, about $1.4\times10^{11}$ relabelings for twenty cultures per condition, so large designs sample relabelings at random instead.

## Variability and the number of cultures

The six-culture experiment cannot show how replicate spread enters the test, because it already sits at the floor $p_{\min}=0.1$. Consider instead two experiments with five cultures per condition, so that $p_{\min}=2/252$. Both have $\bar x_C=4$, $\bar x_T=5$, and $\widehat\Delta=1$, but the within-condition standard deviation is about 0.3 in the first and about 1.1 in the second.

```{.python .execute}
#| label: fig-conditions-spread
#| code-fold: true
#| code-summary: "Compare two spreads with the same difference in means"
#| fig-width: 6
#| fig-height: 4.6
#| fig-cap: "Two experiments with five cultures per condition and the same difference in means, 1 (top), and the 252 relabeled differences for each (bottom); red lines mark the observed difference and its negative. Larger culture-to-culture spread widens the relabeling distribution, so the same difference becomes far less unusual."
#| fig-alt: "Top left: tightly clustered control values near 4 and treated values near 5 with no overlap. Top right: widely spread overlapping values. Bottom left: a narrow histogram with red lines at plus and minus 1 in its extreme tails, p = 0.008. Bottom right: a wide histogram with the red lines inside its bulk, p = 0.206."
experiments = {
    "Low spread": ([3.6, 3.8, 4.0, 4.2, 4.4], [4.6, 4.9, 5.0, 5.1, 5.4]),
    "High spread": ([2.6, 3.4, 3.9, 4.6, 5.5], [3.3, 4.4, 5.3, 5.8, 6.2]),
}
fig, axes = plt.subplots(2, 2)
for column, (title, (control, treated)) in enumerate(experiments.items()):
    values = control + treated
    p, differences = permutation_p_value(values, set(range(5, 10)))
    top, bottom = axes[0, column], axes[1, column]
    top.plot([0] * 5, control, "o", color="C0", alpha=0.8)
    top.plot([1] * 5, treated, "^", color="C1", alpha=0.8)
    top.set(xticks=[0, 1], xticklabels=["control", "treated"], xlim=(-0.6, 1.6),
            ylim=(2, 7), title=title)
    bottom.hist(differences, bins=np.arange(-2.3, 2.4, 0.2), color="C0",
                edgecolor="none", alpha=0.9)
    for observed in (-1, 1):
        bottom.axvline(observed, color=cf.EMPHASIS, linestyle="--")
    bottom.set(xlim=(-2.3, 2.3), ylim=(0, 60), title=f"p = {p:.3f}",
               xlabel=r"Relabeled difference $\widehat\Delta(a)$")
for observed in (-1, 1):
    axes[1, 0].annotate("1 relabeling", (observed, 1), xytext=(observed * 1.15, 14),
                        ha="left" if observed > 0 else "right", fontsize=8,
                        arrowprops={"arrowstyle": "-", "color": cf.RULE})
axes[0, 0].set_ylabel(r"Expression, $\log_2(1+\mathrm{CPM})$")
axes[1, 0].set_ylabel("Relabelings")
plt.show()
```

@fig-conditions-spread shows the permutation test making the comparison that the effect size alone could not. When cultures are tightly clustered, any relabeling that moves a culture across the line pulls the two means together, so the observed difference of 1 is matched only by itself and its reversal: $p=2/252\approx0.008$. When cultures overlap, swapping a high control with a low treated culture barely changes the means, many relabelings reach $\pm1$, and $p=52/252\approx0.21$. The relabeling distribution widens with the within-condition spread, and the observed difference is judged against that width. More cultures help in two ways: they lower $p_{\min}$, and they make each condition mean less sensitive to one unusual culture.

## When relabeling is not justified

Which relabelings are allowed is where the test can go wrong. Suppose each of the six cultures had been sequenced three times, with technical deviations of $-0.1$, $0$, and $+0.1$, and an analyst treated the eighteen assays as eighteen replicates.

```{.python .execute}
#| label: code-conditions-pseudoreplication
assay_offsets = [-0.1, 0.0, 0.1]  # repeatability of the assay on one culture
assay_values = []
for value in culture_values:
    for offset in assay_offsets:
        assay_values.append(value + offset)
treated_assays = set(range(9, 18))  # the nine assays of t1, t2, t3
assay_p, assay_differences = permutation_p_value(assay_values, treated_assays)
print("observed difference:", round(mean_difference(assay_values, treated_assays), 3))
print("relabelings:", len(assay_differences), "  p-value:", round(assay_p, 6))
```

The effect estimate is still 3, because nothing about the cultures changed, yet the p-value falls from $0.1$ to $2/48620\approx0.00004$. The new relabelings split one culture's three assays between conditions, which the experiment could never have done, since the drug was applied to whole cultures. Assays of one culture differ by about $0.1$, far less than cultures differ from one another, so the relabeling distribution is too narrow and the observed difference looks far more unusual than it is [@holmeshuber2019]. This mistake is called *pseudoreplication*. The remedy is to relabel exactly the units to which the condition was assigned, averaging any technical replicates within each unit first.

The same principle governs processing batches. If all treated cultures were sequenced on one day and all controls on another, the test still returns $p=0.1$, but every relabeling moves treatment and processing day together, so the data cannot separate a drug effect from a day effect (@thm-pca-confounding). Paired designs restrict the relabelings in yet another way (@exr-conditions-paired).

For one gene we now have an effect estimate and a calibrated p-value; the expression matrix has thousands of columns, and an experiment usually asks which of them changed. Lesson 18 starts from the resulting list of p-values and asks how to choose genes for follow-up when many comparisons are made at once.

## Exercises

::: {#exr-conditions-enumerate}
List all twenty relabelings of the six-culture example with their treated sums $t$ and differences $\widehat\Delta(a)$, and check the stack heights in @fig-conditions-permutation-null. What would the one-sided p-value be for a statistic $\widehat\Delta(a)$ chosen in advance to detect an increase?
:::

::: {#exr-conditions-spread}
Replace the cultures by controls $(1,3,5)$ and treated $(4,6,8)$. Show that $\widehat\Delta$ is unchanged, compute both sample variances, and find the exact two-sided permutation p-value by hand. Explain why it differs from $0.1$.
:::

::: {#exr-conditions-unit}
A study feeds four mice a test diet and four a control diet, then measures a gene in three tissue sections from each mouse. Identify the experimental unit, the number of allowed relabelings, and the smallest attainable two-sided p-value.
:::

::: {#exr-conditions-paired}
Cells from each of three donors are split into a control and a treated dish, giving (control, treated) measurements $(2,4)$, $(5,6)$, and $(3,5)$. Using the mean within-donor difference as the statistic and swapping labels only within donors, list the allowed relabelings and compute the exact two-sided p-value.
:::

::: {#exr-conditions-interpret}
A collaborator summarizes the six-culture experiment as "$p=0.1$, so there is a 10% chance the drug works, and the gene does not respond." Rewrite the summary correctly, reporting both the effect size and the p-value.
:::

::: {#refs}
**References**
:::
