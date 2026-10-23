---
published: false
title: "Problem Set 4: Assembling a Genome from Short Words"
subtitle: "Lessons 6 and 7 · About 85 minutes"
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Python, R, Julia, MATLAB, or another language is fine. You may copy and adapt code from the notes and use AI to help write or debug it. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies. Graph and assembly libraries would skip the steps being assessed, so build the graphs and paths with your own code or the code from Lessons 6 and 7. Dictionaries, sorting, and file reading are fine.

Chris will rarely read your code line by line, so put the evidence in the write-up: a calculation with a known answer, a check that would catch a specific mistake, and a short explanation of what an output means. Make a first attempt at Problems 1 and 2 without tools as exam practice; then use tools to check and correct it.

Problems 3 and 4 use the [phiX174 reference genome](../data/phix174-NC_001422.1.fasta) (5,386 bases, NCBI accession NC_001422.1). The genome is circular, but as in Lesson 6 we treat the FASTA sequence as a linear string and omit words that cross its end/start boundary. All reads in this set are synthetic, error-free unless stated otherwise, and already oriented to the strand given in the file. Report genome positions starting from 1.

## 1. Build a graph and check every edge (25 points, about 20 minutes)

A synthetic 11-base genome segment has the ideal 3-mer composition below. The words are listed alphabetically, not in genome order.

$$
K=\{\!\{\texttt{ACG},\texttt{ACT},\texttt{CAC},\texttt{CAC},\texttt{CCA},\texttt{CGC},\texttt{CTA},\texttt{GCC},\texttt{TAA}\}\!\}
$$

a. Before drawing anything, predict from the list alone the number of edges in $D(K)$, the number of nodes, and the length of every string with composition $K$. Give a one-sentence reason for each.

b. Draw $D(K)$, keeping parallel edges. Tabulate $d^-(v)$ and $d^+(v)$ for every node, and use the degree criterion from Lesson 7 to identify where an Eulerian path must start and end.

c. Find an Eulerian path by hand and spell its string. Then do an edge-use check: a tally showing that your path uses each word exactly as many times as it occurs in $K$.

d. Node `AC` has two different outgoing edges. Prove that only one string has composition $K$ by explaining what happens to a path that takes the other edge out of `AC` first.

**Evidence to submit:** your predictions, the drawing, the degree table, the spelled string with its tally, and the uniqueness argument.

## 2. Same words, different genomes (20 points, about 15 minutes)

a. Construct two different DNA strings, each of length at most 10, with the same 3-mer composition. The two strings must begin with the same 2-mer and end with the same 2-mer. Do not use `ACGTCGACG` and `ACGACGTCG` from Lesson 7, or that pair with letters renamed. Certify the pair by printing both compositions in sorted order, and name the node at which the two strings take their routes in a different order.

b. Find the smallest $k$ at which your two compositions differ. Name one $k$-mer that occurs in only one of the strings.

c. Two classmates make claims about strings $x$ and $y$ and an integer $k\geq2$.

- Claim 1: if $\operatorname{Comp}_k(x)=\operatorname{Comp}_k(y)$, then $x$ and $y$ have the same length.
- Claim 2: if $\operatorname{Comp}_k(x)=\operatorname{Comp}_k(y)$, then $\operatorname{Comp}_{k-1}(x)=\operatorname{Comp}_{k-1}(y)$.

For each claim, either prove it or give a counterexample. If a claim is false, add a condition under which it holds and prove that version. For Claim 2, explain in a sentence why a circular genome such as phiX174, whose FASTA file can begin at any position on the circle, is the natural source of trouble.

Adapted from @compeau2015.

**Evidence to submit:** the certified pair, the separating $k$ with its distinguishing word, and the two short proofs or counterexamples.

## 3. A repeat in a real genome (30 points, about 30 minutes)

Lesson 7 showed that a repeated $(k-1)$-mer can allow two strings with the same composition. Here we ask whether this happens in phiX174 at a particular $k$, and whether the Lesson 7 code notices.

a. Make a prediction before computing. Suppose the genome were a string of $n$ independent, uniformly random bases. Two fixed starting positions then carry the same $k$-mer with probability $4^{-k}$, so the expected number of position pairs sharing a $k$-mer is about
$$
\binom{n-k+1}{2}\,4^{-k}.
$$
Use this estimate to predict the length of the longest repeated substring of phiX174. Show the arithmetic for two neighboring values of $k$.

b. For $k=10,\ldots,14$, count the distinct $k$-mers that occur at least twice in phiX174. Report the longest repeated words and the start positions of each copy. Compare the result with your prediction in one or two sentences.

c. For $k=13$, build the ideal edge list with one prefix-to-suffix edge per genome position, in genome order. Run Lesson 7's `eulerian_path` and `spell_node_path` on it, or a translation that, like the lesson code, removes the most recently added unused edge from a node's adjacency list. Repeat for $k=14$. For each $k$, report whether the spelled string equals the reference. If it differs, report its length, the first position where it differs, and whether its $k$-mer composition equals the reference's. Does the code have a bug?

d. Explain the $k=13$ result using the repeat positions from part (b). Which two stretches of the genome can trade places, and why does the order in which the repeat copies occur along the genome allow the trade? How many different strings have the reference's $13$-mer composition? Justify your count. (Compare `TACGACG` in Lesson 7, where a repeated node does not create a second string.) If your translation returned the reference at $k=13$, use this explanation to construct the other string and certify it as in part (c).

e. In two or three sentences: would this ambiguity remain if we had error-free 150-base reads covering phiX174? Describe a genome for which the analogous ambiguity would remain even with such reads.

**Evidence to submit:** the prediction with its arithmetic, the table for $k=10,\ldots,14$, the repeat positions, the outputs and composition checks at $k=13$ and $14$, and the explanation of which stretches trade places.

## 4. Reading a k-mer count table (25 points, about 20 minutes)

Real reads contain errors and observe each genomic word many times. The table below comes from synthetic reads sampled from phiX174 with the generator shown after it: 1,616 reads of length $L=100$, start positions uniform over the linear sequence, and each base independently replaced, with probability $\varepsilon=0.005$, by one of the other three bases chosen uniformly. This gives nominal coverage $C=30$. All 21-mers inside each read were counted; no words cross read boundaries. The last column is available only because the reads are synthetic; an assembler would not know it.

| Count in reads | Distinct 21-mers | Of these, in the reference |
|---|---:|---:|
| 1 | 12,245 | 9 |
| 2 | 247 | 7 |
| 3–5 | 13 | 13 |
| 6–10 | 77 | 77 |
| 11–15 | 338 | 338 |
| 16–20 | 1,676 | 1,676 |
| 21–25 | 2,106 | 2,106 |
| 26–30 | 986 | 986 |
| 31–37 | 152 | 152 |
| Total | 17,840 | 5,364 |

: 21-mer counts in synthetic reads from phiX174 (generator seed 4310). The reference contains 5,366 distinct 21-mers; 2 never appear in any read. The largest count is 37. {#tbl-ps4-kmer-counts}

```python
# Generator used for Problem Set 4, Problem 4 (you do not need to rerun it)
import random

def simulate_reads(genome: str, n_reads: int, length: int,
                   error_rate: float, rng: random.Random) -> list[str]:
    """Sample same-strand reads with independent substitution errors."""
    reads = []
    for _ in range(n_reads):
        start = rng.randrange(len(genome) - length + 1)
        bases = list(genome[start:start + length])
        for i, base in enumerate(bases):
            if rng.random() < error_rate:
                bases[i] = rng.choice([b for b in "ACGT" if b != base])
        reads.append("".join(bases))
    return reads

# reference holds the phiX174 sequence as one string, read as in Lesson 6
reads = simulate_reads(reference, 1616, 100, 0.005, random.Random(4310))
```

a. Make two predictions before using the table. First, a reference 21-mer away from the genome ends should be counted about $C\,\frac{L-k+1}{L}(1-\varepsilon)^k$ times; explain each factor and evaluate it. Second, estimate the number of erroneous 21-mer observations: multiply the expected number of substitutions by the average number of 21-mers that contain a given read position, which you should derive. Compare both predictions with the table. What count would you expect for a 21-mer that occurs twice in the genome, and what does the table suggest about repeated 21-mers in phiX174?

b. An assembler keeps only 21-mers with count at least $t$. For $t=2$ and $t=3$, how many of the 5,366 reference 21-mers are absent from the filtered set, and how many erroneous 21-mers remain? Which threshold would you choose, and why?

c. All 31 reference 21-mers with count at most 5 start at positions 1–18 or 5,353–5,366. Explain why, and compute the expected count of the 21-mer starting at position 1. What would a count threshold do to the ends of the assembly? Why would reads from the actual circular molecule not show this pattern?

d. A classmate passes all 129,280 observations, repeats included, to `eulerian_path`, receives "Degree conditions fail," and concludes that sequencing errors caused the failure. Give a set of two or three error-free reads for which the same call fails, and explain what the Lesson 7 workflow does with repeated observations instead.

Adapted from @compeau2015.

**Evidence to submit:** both predictions with their comparison, the threshold counts and your choice, the end-of-genome explanation, and the error-free counterexample.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: correct hand constructions, certificates that a reader can check, predictions compared honestly with outputs, and conclusions limited to what the data support.

::: {#refs}
**References**
:::
