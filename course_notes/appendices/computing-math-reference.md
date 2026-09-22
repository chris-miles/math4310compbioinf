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
