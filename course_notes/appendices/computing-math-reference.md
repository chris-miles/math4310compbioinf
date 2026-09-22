---
title: "Computing and Math Reference"
---

Use this page when a lesson's code or notation is unfamiliar. The lessons introduce the ideas they need; these reminders are here to help you work through an example or check your own result.

## Languages and the lecture code

The notes use Python 3.11 or later. You are welcome to use R, Julia, MATLAB, or another language for your work, subject to any algorithm-specific instructions in the assignment. You do not need to translate your solution into Python.

If you want to copy and run the lecture examples, their main packages are NumPy, SciPy, pandas, and Matplotlib. Some lessons use NetworkX or a domain library. Install only what the example you are running needs. From a terminal, a typical setup is:

```sh
python -m pip install numpy scipy pandas matplotlib networkx
```

These are the notes' tools, not a required framework for submissions. You do not need Pandoc or the website builder to run an individual algorithm. Copy the setup imports and the relevant function into a script or notebook. Run a page's examples in order when later chunks reuse earlier variables.

Figures use the helper in [coursefigs.py](https://github.com/chris-miles/math4310compbioinf/blob/main/course_notes/style/coursefigs.py). Its relative import paths assume the working directory is the lesson folder. The mathematical functions generally do not need the figure helper; keep the plotting setup only when you are reproducing a figure.

When an assignment asks you to implement a method, a library call that already performs that method may bypass the task. Libraries for arrays, plotting, or loading files are usually fine. Follow the particular problem's instructions.

## What to show with a computational answer

Chris will rarely read your code line by line. The main evidence is your mathematical reasoning, validation, actual outputs, and interpretation. Submit runnable code when requested so that the results can be reproduced, and be ready to explain the relevant part of it.

A useful write-up states the input and expected answer for a small check, shows the actual output, and explains why they should agree. A known answer might come from a hand calculation, exhaustive enumeration on a tiny input, or a second formulation that does not reuse the same algorithm.

A passing test establishes only what it checks. An alignment can have the correct reported score but the wrong printed letters. Check both that removing gaps recovers the inputs and that independently scoring the printed columns recovers the reported score.

See the [course policy on tools, data, and AI](../course-info.html#tools-data-and-ai-use). When AI assistance is permitted, you may use it to draft or debug code. Disclose material assistance, check its output, and explain your conclusions in terms of the calculation rather than the tool's authority.

## Strings, indices, and intervals

The mathematics numbers sequence positions from 1: $x_1,\ldots,x_n$. Python indexes the corresponding string from 0, so mathematical $x_i$ is Python's <code>x[i - 1]</code>. Python's slice <code>x[start:stop]</code> includes <code>start</code> and excludes <code>stop</code>. R's usual indexing starts at 1; translate the mathematical ranges rather than copying Python bounds unchanged.

A length-$n$ string contains $\max(n-k+1,0)$ overlapping length-$k$ words for positive integer $k$. Repeated occurrences count separately. Different reads must be processed separately: concatenating them creates words crossing a boundary that was never observed.

## A small verification example

::: {.algorithm}
**Algorithm: count overlapping words.**

For every legal starting position, extract the length-$k$ word and add one to its dictionary count. Reject a nonpositive value of $k$.
:::

```{.python .execute}
#| label: code-reference-kmers
def count_kmers(sequence: str, k: int) -> dict[str, int]:
    """Count overlapping length-k substrings."""
    if k < 1:
        raise ValueError("k must be positive.")
    counts = {}
    for start in range(len(sequence) - k + 1):
        word = sequence[start:start + k]
        counts[word] = counts.get(word, 0) + 1
    return counts


observed = count_kmers("ACGACG", 2)
expected = {"AC": 2, "CG": 2, "GA": 1}
assert observed == expected
assert sum(observed.values()) == 5
observed
```

The explicit words are AC, CG, GA, AC, CG. The dictionary comparison checks both labels and multiplicities. The total-count check alone would not catch a program that assigned all five occurrences to AC.

## Asymptotic notation

For nonnegative functions $f,g$, we write $f(n)\in\mathcal O(g(n))$ if there are constants $C>0,n_0$ such that $f(n)\le Cg(n)$ for all $n\ge n_0$. This is an upper bound on growth, not an exact runtime prediction.

One pass over a length-$n$ string with constant work per position takes $\mathcal O(n)$ time. A table with $n$ rows and $m$ columns and constant work per cell takes $\mathcal O(nm)$ time. With two independently varying lengths, keep both variables.

For k-mer code, “linear in sequence length” usually treats $k$ as fixed. Creating and hashing an ordinary length-$k$ Python substring also takes work proportional to $k$. State which quantities are fixed when comparing costs.

## Counting

If a construction has $a$ choices followed by $b$ choices for each first choice, it has $ab$ possibilities. If outcomes fall into disjoint groups of sizes $a,b$, it has $a+b$ possibilities. The disjointness is why a recurrence can count paths by their final move.

The binomial coefficient $\binom nr=n!/[r!(n-r)!]$ counts ways to choose $r$ positions from $n$. A grid path containing $r$ horizontal and $s$ vertical moves has $\binom{r+s}{r}$ possible orders. Adding diagonal moves changes the counting problem.

## Probability

For $P(B)>0$, conditional probability is $P(A\mid B)=P(A\cap B)/P(B)$. Rearranging gives the product rule. Events are independent when $P(A\cap B)=P(A)P(B)$.

A random variable is a numerical or symbolic outcome; a realization is the particular value observed. In the sequence lessons, $X_i$ is a random base and $x_i$ is the letter in the file.

For fixed model parameters $\theta$, $P_\theta(x)$ assigns probabilities to possible data. For fixed observed data $x$, the same expression is the likelihood $L(\theta;x)$ as a function of $\theta$. Likelihoods do not have to sum to one over parameter choices.

Bayes' rule combines likelihood and a prior:
$$
P(M\mid x)=\frac{P(x\mid M)P(M)}{P(x)}.
$$
For two model classes, posterior odds equal likelihood ratio times prior odds. A large likelihood ratio can still favor a rare class only weakly after its prior is included.

## Logarithms and numerical checks

Logs turn products into sums: $\log(ab)=\log a+\log b$. They do not turn sums into sums of logs. Base 2 gives bits; natural logs are common in code. Keep the base consistent when interpreting a score.

Use approximate comparisons for floating-point results, such as Python's <code>math.isclose</code>, rather than expecting decimal calculations to match exactly. A probability table should have nonnegative entries and rows or columns summing to one according to its definition. An HMM forward column, for example, is a joint probability over a prefix and a state; it need not sum to one.

## Linear algebra

Vectors $u,v$ are orthogonal when $u^\mathsf Tv=0$. An orthonormal set also has unit-length vectors. Projection onto a unit vector $u$ is $(u^\mathsf Tx)u$.

A singular value decomposition writes $X=U\Sigma V^\mathsf T$, with orthonormal singular vectors and nonnegative singular values in decreasing order. PCA applies this decomposition to a centered data matrix. Check whether rows denote samples or genes before interpreting the vectors or multiplying matrices.

## Differential equations

A fixed point of $\dot x=f(x)$ satisfies $f(x^*)=0$. For a differentiable scalar equation, $f'(x^*)<0$ implies local asymptotic stability and $f'(x^*)>0$ implies instability. If the derivative is zero, this test is inconclusive; examine the sign of $f$ nearby.

A numerical trajectory is one solution from one starting point, approximated at a chosen step size and tolerance. Check steady states algebraically and compare numerical results at a smaller step size before interpreting a plot.

The [notation page](notation.html) collects the symbols used across lessons.
