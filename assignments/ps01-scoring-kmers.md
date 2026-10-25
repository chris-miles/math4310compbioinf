---
published: true
title: "Problem Set 1: Comparing Sequences and Checking Results"
subtitle: "Lessons 1 and 2 · About 60 minutes"
nocite: |
  @compeau2015
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Python, R, Julia, MATLAB, or another language is fine. You may copy and adapt code from the notes and use AI to help write or debug it. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies.

Chris will rarely read your code line by line. Most of the evidence for your work should be in the write-up: a calculation with a known answer, a test that would catch a specific mistake, and an explanation of what an output means. Keep those explanations short and specific. For the hand calculations in Problems 1, 2, and 4a, make a first attempt without tools as exam practice; then use tools to check and correct it.

## 1. Two ways to count a change (20 points, about 10 minutes)

Let $x=\texttt{ACGTA}$ and $y=\texttt{CGTAA}$.

a. Compute the Hamming distance by marking disagreeing positions.

b. Give an edit sequence that transforms $x$ into $y$ using two edits. Prove that one edit cannot suffice. A program reporting a distance is not a substitute for this lower-bound argument.

c. Suppose these strings are short fragments from corresponding genomic regions. Explain in two or three sentences why a large Hamming distance need not mean many mutation events. Also state why a minimum edit distance cannot establish the actual evolutionary history.

**Evidence to submit:** marked strings, the edit sequence, and the lower-bound argument.

## 2. Matches as evidence (25 points, about 13 minutes)

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

## 3. Repair a k-mer counter and test the repair (30 points, about 20 minutes)

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

c. For valid $k$, check that the sum of all reported counts equals $\max(n-k+1,0)$. Explain why, for $1\le k\le n$, this check would detect the supplied bug, but would *not* detect a program that assigned the right total count to the wrong words.

**Evidence to submit:** your corrected function, the four output comparisons, and a two- or three-sentence explanation connecting the loop bounds to the legal starting positions. If a test fails, preserve the failed output and explain the correction.

## 4. What word counts keep and what they hide (25 points, about 14 minutes)

a. By hand, list the 2-mer dictionaries of $\texttt{AACAA}$ and $\texttt{ACAAA}$, then check them with your repaired function. Use these outputs to assess the claim “identical counts mean identical sequences.” Name a 3-mer that occurs in one string but not the other.

b. Now count words in a real genome, the [phiX174 genome](../data/phix174-NC_001422.1.fasta) (5,386 bases, NCBI accession NC_001422.1), treated as linear as [the data notes](https://github.com/chris-miles/math4310compbioinf/blob/main/data/README.md) describe. The code below reads it into one string and shows how to generate a seeded random genome of the same length for part d.

```python
import random

def read_fasta(path: str) -> str:
    """Return the sequence of a one-record FASTA file as one string."""
    with open(path) as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">"))

phix = read_fasta("data/phix174-NC_001422.1.fasta")   # path from the repository root

rng = random.Random(4310)                               # choose and report your own seed
random_genome = "".join(rng.choice("ACGT") for _ in range(len(phix)))
```

Before counting, compute the expected number of occurrences of one particular 9-mer, such as $\texttt{ACGTTGCAA}$, in a genome of this length whose letters are independent and uniform. Then write down a guess for the largest count that *any* 9-mer reaches in such a random genome.

c. Count the 9-mers of phiX174 with your repaired function. Report the total count, checked as in Problem 3c, the number of distinct 9-mers, and the most frequent 9-mer with its count.

d. Generate 20 random genomes of length 5,386 with independent uniform letters, using a fixed seed that you report. Record the largest 9-mer count in each genome. Compare these 20 values with your guess from part b and with phiX174. Explain why the largest count in a random genome is far above the expected count from part b. Then say in one or two sentences what the phiX174 count does show and what it does not.

Adapted from @compeau2015, Ch. 1.

**Evidence to submit:** both dictionaries, the expected count and your written guess, the three phiX174 outputs, the 20 maxima with their seed, and the interpretation.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State “No outside assistance” if applicable.

Credit follows the evidence: correct calculations, a justified lower bound, tests with independently known answers, predictions compared honestly with outputs, and conclusions limited to what the data support. A program that runs without errors is only the beginning of a check.

::: {#refs}
**References**
:::
