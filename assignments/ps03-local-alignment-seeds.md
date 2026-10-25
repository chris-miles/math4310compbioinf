---
published: true
title: "Problem Set 3: Local Alignment, Gap Costs, and Seeds"
subtitle: "Lessons 4 and 5 · About 80 minutes"
nocite: |
  @compeau2015
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Any language is fine, and you may adapt the Lesson 4 and 5 code and use AI to help write or debug it. Disclose that help briefly and say how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies. Unless a problem says otherwise, score a match $+1$, a mismatch $-1$, and each gap column $-2$. Make a first attempt at Problems 1, 2, and 3a without tools.

## 1. Finding a shared core and its coordinates (20 points, about 18 minutes)

Let $x=\texttt{CCGATG}$ and $y=\texttt{GATGAA}$.

a. Fill the local-alignment table $H(i,j)$ by hand. Most entries are zero; say in one sentence why.

b. Report the optimal local score, the local alignment, and its coordinates in the Lesson 4 convention: the cell $(i_0,j_0)$ where traceback stops and the cell $(i_1,j_1)$ where it starts. Write out which letters of $x$ and of $y$ the alignment uses.

c. Compute the global score of the same pair. In two or three sentences, explain what the local result tells a biologist that the global result hides, using the positions of the shared region in each string.

**Evidence to submit:** the filled table with the path marked, the score, alignment, and coordinates, the global score, and the explanation.

## 2. One long gap or several short ones (15 points, about 12 minutes)

A read differs from its reference by a deletion of three adjacent bases. Two candidate alignments have the same matched letters. Alignment A places the three missing bases in one gap of length 3. Alignment B splits them into three gaps of length 1 at different places.

a. Under the linear gap score of Lesson 3, with $-2$ per gap column, compute the gap cost of each alignment. Can linear scoring prefer either one?

b. Under affine costs $G(\ell)=d_o+(\ell-1)d_e$ from Lesson 5, write each alignment's gap cost in terms of $d_o$ and $d_e$, and find exactly when A is preferred. Evaluate both costs for $(d_o,d_e)=(4,1)$.

c. Explain in two or three sentences why the preference you found matches a biological view of how insertions and deletions arise, and why the Lesson 3 recurrence, which keeps one score per cell and looks back only at the three neighboring cells, cannot apply affine costs correctly.

**Evidence to submit:** the cost expressions, the condition for preferring A, the numbers at $(4,1)$, and the explanation.

## 3. Seeds that survive, and seeds that do not (20 points, about 15 minutes)

Lesson 5 shows that an ungapped alignment of length $L$ with $t$ mismatches contains an exact match of length at least $\lceil (L-t)/(t+1)\rceil$.

a. For $L=24$ and $t=3$, compute this bound. Then construct a reference segment of length 24 and a read that differs from it at exactly 3 positions so that the read and segment share **no** exact 7-letter word at the same offset. Your certificate is the list of mismatch positions and the lengths of the matching runs between them.

b. Run `exact_seed_hits` from Lesson 5, or your own version, on your pair for $k=5$, $6$, and $7$. Report the number of hits on the true offset, where the two starting positions are equal. Explain why no placement of 3 mismatches can remove every 6-letter seed.

Adapted from @compeau2015, Ch. 9.

**Evidence to submit:** the constructed pair with its certificate, the seed counts for three values of $k$, and the explanation.

## 4. Finding a gene in a related genome (30 points, about 22 minutes)

In Problem Set 2, NCBI's annotation told us where gene D sits in both phages. A database search has no such help, so here we locate a gene with seeds and local alignment, and check the annotation only at the end. Gene F of phiX174 encodes the major capsid protein; its coding sequence occupies positions 1001–2284 of the [phiX174 genome](../data/phix174-NC_001422.1.fasta). Use that coding sequence as the query, with $n=1284$, and the whole [G4 genome](../data/g4-NC_001420.2.fasta) as the target, with $N=5577$. Search both files as given, as linear strings: do not wrap around the circular genomes or search reverse complements. The scaffold below loads the sequences and cuts a G4 segment for part c.

```python
# --- supplied for Problem Set 3, Problem 4 ---
def read_fasta(path: str) -> str:
    """Return the sequence of a one-record FASTA file as one string."""
    with open(path) as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">"))

phix = read_fasta("data/phix174-NC_001422.1.fasta")
g4 = read_fasta("data/g4-NC_001420.2.fasta")
gene_f = phix[1000:2284]            # 1-based positions 1001-2284

def g4_segment(d: int, flank: int = 50) -> tuple[int, str]:
    """Return the 1-based G4 position of the segment's first letter, and the segment,
    covering G4 positions d+1-flank through d+len(gene_f)+flank."""
    first = max(1, d + 1 - flank)
    last = min(len(g4), d + len(gene_f) + flank)
    return first, g4[first - 1:last]

# Letters segment[j0:j1] of a local alignment are G4 positions first + j0 through first + j1 - 1.
```

a. Before running anything, compute the expected number of chance seed hits $\mathbb E[h]$ from Lesson 5 for $k=8$, $11$, and $14$, assuming both sequences have independent uniform letters.

b. Find the exact seed hits for $k=8$, $11$, $14$, and $20$. For each $k$, report the total number of hits, the offset `target_start - query_start` with the most hits (or "none" if there are no hits), and the number of hits at that offset. Compare the hits at all other offsets with your predictions from part a, and give one reason a real genome can disagree with the uniform model.

c. Let $d$ be the top offset found at $k=8$, $11$, or $14$. A homolog of gene F along that offset would occupy G4 positions $d+1$ through $d+1284$. Extend the hits by running local alignment of gene F against `g4_segment(d)`, which adds 50 bases on each side of that range. Report the score, which letters of gene F and which G4 positions the alignment uses, and the number of gap columns. More than one endpoint can attain the best score; say which one your code reports and why. NCBI annotates the G4 gene F coding sequence at positions 2600–3883. Compare your coordinates with the annotation and explain any difference at the ends.

d. A colleague proposes $k=20$, because a 150-letter read searched against $3\times10^9$ letters would then expect fewer than one chance hit. Compute that expectation. Then use your counts from part b to explain what $k=20$ would cost when the homolog is as divergent as G4's gene F.

**Evidence to submit:** the three expected counts written before counting, a table of seed results for four values of $k$, the local alignment's score and coordinates with your tie policy and the comparison to the annotation, and the calculation and explanation for $k=20$.

*Optional extension, not graded.* Local alignment may leave letters at either end of the query unaligned. Suppose every letter of gene F must be aligned, while the alignment may still begin and end anywhere in G4. Describe what changes in the local recurrence, in its first row and first column, and in the cell where traceback starts. This extension is adapted from @compeau2015, Ch. 5.

## 5. Predict, then test: when local alignment stops being local (15 points, about 12 minutes)

Local alignment finds the edges of a shared region only if unrelated columns lower the score on average.

a. For uniform DNA, compute the expected score of one unrelated column, $\mathbb E_R[s]=\sum_{a,b}q_aq_b\,s(a,b)$, under three schemes: match $+1$ and mismatch $-1$; match $+3$ and mismatch $-1$; match $+5$ and mismatch $-1$ (gaps $-2$ throughout).

b. **Before computing anything else**, predict how long the optimal local alignment of two *unrelated* random DNA strings of length 200 will be under the first and third schemes. Write your prediction down.

c. Generate two random strings of length 200 with any seed you choose, run local alignment under both schemes, and report the score and the length of the aligned region in each string. Compare with your prediction and explain the difference using part a.

**Evidence to submit:** the three expected scores, your written prediction, the actual outputs with the seed you used, and the explanation.

## Submission check

Include enough output to reproduce each conclusion without opening your code. End with two or three sentences naming any code, tools, AI assistance, or people you used, and one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: correct tables and coordinates, certificates that settle a construction, predictions made before computing and compared honestly with outputs, and conclusions limited to what the scores support.

::: {#refs}
**References**
:::
