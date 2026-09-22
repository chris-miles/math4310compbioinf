---
published: true
title: "Problem Set 1: Comparing Sequences and Checking Results"
subtitle: "Lessons 1 and 2 · About 60 minutes"
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Python, R, Julia, MATLAB, or another language is fine. You may copy and adapt code from the notes and use AI to help write or debug it. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies.

Chris will rarely read your code line by line. Most of the evidence for your work should be in the write-up: a calculation with a known answer, a test that would catch a specific mistake, and an explanation of what an output means. Keep those explanations short and specific. For the hand calculations, make a first attempt without tools as exam practice; then use tools to check and correct it.

## 1. Two ways to count a change (25 points, about 12 minutes)

Let $x=\texttt{ACGTA}$ and $y=\texttt{CGTAA}$.

a. Compute the Hamming distance by marking disagreeing positions.

b. Give an edit sequence that transforms $x$ into $y$ using two edits. Prove that one edit cannot suffice. A program reporting a distance is not a substitute for this lower-bound argument.

c. Suppose these strings are short fragments from corresponding genomic regions. Explain in two or three sentences why a large Hamming distance need not mean many mutation events. Also state why a minimum edit distance cannot establish the actual evolutionary history.

**Evidence to submit:** marked strings, the edit sequence, and the lower-bound argument.

## 2. Matches as evidence (25 points, about 15 minutes)

Use the synthetic pair model from Lesson 2. Background bases are independent and uniform, so every ordered pair has probability $1/16$. Under the related model, each matching pair has probability $1/8$ and each mismatching pair has probability $1/24$.

Consider these two fixed ungapped alignments:

| Candidate | Top | Bottom |
|---|---|---|
| A | ACGT | ACGA |
| B | ACGT | AGGA |

a. Compute both log-odds scores in bits. Leaving a logarithm unevaluated is acceptable. For candidate A, also compute the likelihood ratio by multiplying the four pairwise ratios and verify that it equals $2^{\text{score}}$.

b. A program instead scores a match as $+1$ and a mismatch as $-1$. Does it rank these candidates in the same order? Does its score have the same numerical interpretation in bits under the supplied model?

c. A classmate says, “The higher score proves that candidate A is related.” Explain what the comparison actually establishes, and name one piece of information it does not supply.

**Evidence to submit:** both scores, one independent ratio calculation, and the interpretation. A screenshot of a calculator alone is not enough.

## 3. Repair a k-mer counter and test the repair (30 points, about 23 minutes)

This Python code has an intentional counting error. It should count every overlapping length-$k$ word, including repeated occurrences.

```python
def count_kmers(sequence: str, k: int) -> dict[str, int]:
    """Count overlapping length-k words; this version contains a bug."""
    if k < 1:
        raise ValueError("k must be positive.")
    counts = {}
    for start in range(len(sequence) - k):
        word = sequence[start:start + k]
        counts[word] = counts.get(word, 0) + 1
    return counts
```

a. Before running it, predict its output for $\texttt{AAAA}$ with $k=2$. Write the correct count separately. Identify the missing starting position and repair the loop. You may translate the code into your preferred language.

b. Run your repaired function on the four cases below. Report the input, expected output, actual output, and why the expectation is known.

- $\texttt{AAAA}$ with $k=2$.
- A DNA string of your choice of length 6 to 8, with $k=3$; list its overlapping words by hand.
- That same string with $k$ equal to its length.
- That same string with $k$ one greater than its length.

c. For valid $k$, check that the sum of all reported counts equals $\max(n-k+1,0)$. Explain why this check would detect the supplied bug, but would *not* detect a program that assigned the right total count to the wrong words.

**Evidence to submit:** your corrected function, the four output comparisons, and a two- or three-sentence explanation connecting the loop bounds to the legal starting positions. If a test fails, preserve the failed output and explain the correction.

## 4. What the counts cannot tell us (20 points, about 10 minutes)

A researcher compares $\texttt{AACAA}$ and $\texttt{ACAAA}$ using their 2-mer count dictionaries.

a. List both dictionaries by hand, then check them with your repaired function. A claim is made: “Identical counts mean identical sequences.” Use these outputs to assess the claim.

b. Repeat the comparison at $k=3$. Identify a particular word that distinguishes the strings. Explain why this result does not prove that increasing $k$ will always solve an assembly problem with short or erroneous reads.

**Evidence to submit:** the dictionaries and a concise interpretation of what was preserved and what was lost.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State “No outside assistance” if applicable.

Credit follows the evidence: correct calculations, a justified lower bound, tests with independently known answers, and conclusions limited to what the data support. A program that runs without errors is only the beginning of a check.
