---
search: false
toc: false
number-sections: false
---

# Using These Notes {.unnumbered}

The numbered lessons are the main reading sequence.

Assignment pages are in the Assignments page in the top navigation. The editable [assignment source files](https://github.com/chris-miles/math4310compbioinf/tree/main/assignments) are on GitHub, so you can copy LaTeX, code, or starter text directly from the `.md` files when that is more convenient than copying from the rendered page.

Small datasets for examples and homework are in the GitHub [data folder](https://github.com/chris-miles/math4310compbioinf/tree/main/data). Assignment pages will link to the specific files they use.

The [Computing and Math Reference](appendices/computing-math-reference.html) collects setup information and reminders. Python is the language used by the examples; assignments also welcome R and other languages unless a problem says otherwise.

Code chunks are part of the notes. They are small mathematical examples, not finished software. Use the copy button to move an example into a notebook or script, change one input, and check whether the output changes in the way the definition predicts.

Try the chunk below. First predict the output by hand. Then expand the code, read the loop, and change one of the two strings.

```{.python .execute}
#| label: code-using-notes-hamming
#| code-fold: true
#| code-summary: "Show code"
def hamming_distance(x: str, y: str) -> int:
    """Count positions where two equal-length strings differ."""
    if len(x) != len(y):
        raise ValueError("Hamming distance needs equal-length strings.")

    distance = 0
    for a, b in zip(x, y):
        if a != b:
            distance += 1
    return distance


hamming_distance("ACGTT", "ACGCT")
```

The function counts substitutions at fixed positions. It gives a quick hand-checkable example before we allow gaps, shifted matches, and dynamic-programming tables.

When a lesson gives an algorithm, keep three questions in view:

1. What does each table, graph, matrix, or probability mean?
2. Why does the procedure return the claimed object?
3. What small example would catch a wrong implementation?

Assignments will ask for verification checks for the same reason. A program that runs is not yet evidence that the method was applied correctly.
