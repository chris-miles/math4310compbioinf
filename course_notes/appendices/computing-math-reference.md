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

Big-O gives an eventual upper bound on computational work as input size grows. For a nonnegative work function $T$, we write $T(n)\in\mathcal O(g(n))$ if constants $C>0$ and $n_0$ exist such that
$$
0\le T(n)\le Cg(n)\qquad\text{for all }n\ge n_0.
$$
The constant $C$ allows a fixed multiplicative factor; the threshold $n_0$ allows exceptions at small input sizes. For example, $3n+7\le 10n$ for $n\ge1$, so $3n+7\in\mathcal O(n)$. Comparing two length-$n$ sequences position by position takes linear work; filling their alignment table takes quadratic work.

### Common growth rates

Common rates, ordered by increasing growth for sufficiently large $n$, are
$$
1,\quad \log n,\quad n,\quad n\log n,\quad n^2,\quad n^3,\quad 2^n.
$$
These are constant, logarithmic, linear, linearithmic, quadratic, cubic, and exponential growth. Any fixed logarithm base greater than 1 gives the same Big-O class: changing the base multiplies the logarithm by a constant. The plot uses base 2.

```{.python .execute}
#| label: fig-reference-growth
#| code-fold: true
#| code-summary: "Plot common growth rates"
#| fig-cap: "Six growth rates on a logarithmic vertical axis; each major tick multiplies work by the same factor. The exponential curve starts below the cubic curve but eventually overtakes it, so their order at a small input size does not determine their growth order."
#| fig-alt: "Six curves for input sizes 2 to 20 on a logarithmic work axis: log base 2 of n, n, n log base 2 of n, n squared, n cubed, and 2 to the n. The exponential curve is below the cubic curve at n equals 2 but above it by n equals 10."
import matplotlib.pyplot as plt
import numpy as np

plt.style.use("../style/figstyle.mplstyle")
n = np.arange(2, 21)
curves = [
    (np.log2(n), r"Logarithmic: $\log_2 n$", "o", "-"),
    (n, r"Linear: $n$", "s", "--"),
    (n * np.log2(n), r"Linearithmic: $n\log_2 n$", "^", "-."),
    (n**2, r"Quadratic: $n^2$", "D", ":"),
    (n**3, r"Cubic: $n^3$", "v", (0, (5, 1))),
    (2.0**n, r"Exponential: $2^n$", "P", (0, (3, 1, 1, 1))),
]
fig, ax = plt.subplots()
for work, label, marker, linestyle in curves:
    ax.plot(n, work, label=label, marker=marker, markevery=3,
            linestyle=linestyle)
ax.set_yscale("log")
ax.set_xlabel("Input size n")
ax.set_ylabel("Work (log scale)")
ax.set_xticks([2, 5, 10, 15, 20])
ax.legend(loc="upper left", ncol=2, fontsize=8)
plt.show()
```

### Small inputs can reverse the comparison

Suppose two algorithms do $T_A(n)=20n$ and $T_B(n)=n^2$ operations on the same task. At $n=5$, the linear algorithm does 100 operations and the quadratic algorithm only 25. They tie at $n=20$. For every $n>20$, the linear algorithm does less work.

```{.python .execute}
#| label: fig-reference-crossover
#| code-fold: true
#| code-summary: "Plot a change in the cheaper algorithm"
#| fig-height: 3
#| fig-cap: "The work curves 20n and n squared cross at n = 20. A larger constant can make a linear algorithm more expensive on small inputs even though it scales better."
#| fig-alt: "Linear-axis plot from n equals 1 to 40. The line 20n lies above n squared until their intersection at n equals 20 and work equals 400, then lies below the quadratic curve."
n = np.arange(1, 41)
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
fig, ax = plt.subplots()
# Keep linear and quadratic colors and markers from the growth-rate plot.
ax.plot(n, 20 * n, label=r"Linear: $20n$", color=colors[1],
        marker="s", markevery=5, linestyle="--")
ax.plot(n, n**2, label=r"Quadratic: $n^2$", color=colors[3],
        marker="D", markevery=5, linestyle=":")
ax.annotate("Equal work at n = 20", xy=(20, 400), xytext=(4, 950),
            arrowprops={"arrowstyle": "->", "color": "black"}, fontsize=9)
ax.set_xlabel("Input size n")
ax.set_ylabel("Operations")
ax.set_xticks([1, 10, 20, 30, 40])
ax.legend(loc="upper left")
plt.show()
```

The ratio $T_B(n)/T_A(n)=n/20$ grows without bound. No fixed multiplicative factor can make linear growth keep up with quadratic growth for all sufficiently large $n$. A plot over a finite range illustrates this comparison; the ratio establishes the eventual behavior.

Big-O is an upper bound, so $20n$ belongs to both $\mathcal O(n)$ and $\mathcal O(n^2)$. The linear bound is more informative. Constants still matter when choosing an algorithm for a particular input size, and operation counts alone do not specify elapsed time in seconds.
