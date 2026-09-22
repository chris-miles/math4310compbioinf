---
title: "Computing and Math Reference"
---

## Running the examples

The notes use Python 3.11 or later, with NumPy, SciPy, pandas, Matplotlib, and NetworkX. To install these packages, run:

```sh
python -m pip install numpy scipy pandas matplotlib networkx
```

Copy a lesson's imports and code into a notebook or script, and run the examples in order. To reproduce the figures, download the [course repository](https://github.com/chris-miles/math4310compbioinf) and run from `course_notes/lessons/` so the plotting helpers can be found.

You may use R or another language for assignments. Follow any instructions about implementing an algorithm yourself.

## Computational work and AI

Chris will rarely read your code line by line. Focus your write-up on the reasoning, checks, actual outputs, and interpretation that establish your answer. See the [tools and AI policy](../course-info.html#tools-data-and-ai-use) for expectations about assistance and disclosure.

## Notation and indexing

Mathematical sequence positions start at 1; Python indices start at 0. Thus $x_i$ corresponds to `x[i - 1]`. A Python slice `x[start:stop]` includes `start` and excludes `stop`.

The [notation page](notation.html) collects symbols used in the lessons. Mathematical explanations appear alongside the methods that need them.

## Big-O notation

Big-O describes how computational work grows with input size, ignoring constant factors. It is an eventual upper bound, not a runtime in seconds.

Formally, $T(n)\in\mathcal O(g(n))$ if constants $C>0$ and $n_0$ exist such that
$$
0\le T(n)\le Cg(n)\qquad\text{for all }n\ge n_0.
$$
For example, $3n+7\le 10n$ for $n\ge1$, so $3n+7\in\mathcal O(n)$. Comparing two length-$n$ sequences position by position takes linear work; filling their alignment table takes quadratic work.

```{.python .execute}
#| label: fig-reference-growth
#| code-fold: true
#| code-summary: "Plot growth rates"
#| fig-cap: "Doubling sequence length doubles linear work but quadruples quadratic work."
#| fig-alt: "Curves n and n squared for input sizes 1 to 20; the quadratic curve grows increasingly steep while the linear curve has constant slope."
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("../style/figstyle.mplstyle")
n = np.arange(1, 21)
fig, ax = plt.subplots()
ax.plot(n, n, label=r"Linear: $n$", linestyle="-")
ax.plot(n, n**2, label=r"Quadratic: $n^2$", linestyle="--")
ax.set_xlabel("Input size n")
ax.set_ylabel("Work (arbitrary units)")
ax.set_xticks([1, 5, 10, 15, 20])
ax.legend()
plt.show()
```
