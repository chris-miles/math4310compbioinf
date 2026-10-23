---
published: false
title: "Problem Set 3: Local Alignment, Gap Costs, and Seeds"
subtitle: "Lessons 4 and 5 · About 80 minutes"
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Any language is fine, and you may adapt the Lesson 4 and 5 code and use AI to help write or debug it. Disclose that help briefly and say how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies. Unless a problem says otherwise, score a match $+1$, a mismatch $-1$, and each gap column $-2$. Make a first attempt at the hand calculations without tools.

## 1. Finding a shared core and its coordinates (25 points, about 20 minutes)

Let $x=\texttt{CCGATG}$ and $y=\texttt{GATGAA}$.

a. Fill the local-alignment table $H(i,j)$ by hand. Most entries are zero; say in one sentence why.

b. Report the optimal local score, the local alignment, and its coordinates in the Lesson 4 convention: the cell $(i_0,j_0)$ where traceback stops and the cell $(i_1,j_1)$ where it starts. Write out which letters of $x$ and of $y$ the alignment uses.

c. Compute the global score of the same pair. In two or three sentences, explain what the local result tells a biologist that the global result hides, using the positions of the shared region in each string.

**Evidence to submit:** the filled table with the path marked, the score, alignment, and coordinates, the global score, and the explanation.

## 2. One long gap or several short ones (20 points, about 15 minutes)

A read differs from its reference by a deletion of three adjacent bases. Two candidate alignments have the same matched letters. Alignment A places the three missing bases in one gap of length 3. Alignment B splits them into three gaps of length 1 at different places.

a. Under the linear gap score of Lesson 3, with $-2$ per gap column, compute the gap cost of each alignment. Can linear scoring prefer either one?

b. Under affine costs $G(\ell)=d_o+(\ell-1)d_e$ from Lesson 5, write each alignment's gap cost in terms of $d_o$ and $d_e$, and find exactly when A is preferred. Evaluate both costs for $(d_o,d_e)=(4,1)$.

c. Explain in two or three sentences why the preference you found matches a biological view of how insertions and deletions arise, and why a single score table that remembers only the best score per cell cannot apply affine costs correctly.

**Evidence to submit:** the cost expressions, the condition for preferring A, the numbers at $(4,1)$, and the explanation.

## 3. Seeds that survive, and seeds that do not (35 points, about 30 minutes)

Lesson 5 shows that an ungapped alignment of length $L$ with $t$ mismatches contains an exact match of length at least $\lceil (L-t)/(t+1)\rceil$.

a. For $L=24$ and $t=3$, compute this bound. Then construct a reference segment of length 24 and a read that differs from it at exactly 3 positions so that the read and segment share **no** exact 7-letter word at the same offset. Your certificate is the list of mismatch positions and the lengths of the matching runs between them.

b. Write a function that lists every exact $k$-mer seed hit between a read and a reference, as pairs of starting positions. Run it on your pair for $k=5$, $6$, and $7$ and report the number of hits on the true diagonal (equal starting positions). Explain why no placement of 3 mismatches can remove every 6-letter seed.

c. Suppose a read of length 150 is searched against a genome of $3\times10^9$ letters with uniform, independent letters. Using $\mathbb E[h]=(n-k+1)(N-k+1)4^{-k}$, compute the expected number of chance seed hits for $k=6$, $11$, $16$, and $20$.

d. A colleague chooses $k=20$ because it gives fewer than one chance hit. Using parts a to c, explain what that choice costs for reads like yours, and what a mapping program must do to keep both errors (missed true placements and wasted extensions) small. Two or three sentences.

**Evidence to submit:** the constructed pair with its certificate, your seed function's output for three values of $k$, the four expected counts, and the explanation.

## 4. Predict, then test: when local alignment stops being local (20 points, about 15 minutes)

Local alignment finds the edges of a shared region only if unrelated columns lower the score on average.

a. For uniform DNA, compute the expected score of one unrelated column, $\mathbb E_R[s]=\sum_{a,b}q_aq_b\,s(a,b)$, under three schemes: match $+1$ and mismatch $-1$; match $+3$ and mismatch $-1$; match $+5$ and mismatch $-1$ (gaps $-2$ throughout).

b. **Before computing anything else**, predict how long the optimal local alignment of two *unrelated* random DNA strings of length 200 will be under the first and third schemes. Write your prediction down.

c. Generate two random strings of length 200 with any seed you choose, run local alignment under both schemes, and report the score and the length of the aligned region in each string. Compare with your prediction and explain the difference using part a.

**Evidence to submit:** the three expected scores, your written prediction, the actual outputs with the seed you used, and the explanation.

## Submission check

Include enough output to reproduce each conclusion without opening your code. End with two or three sentences naming any code, tools, AI assistance, or people you used, and one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: correct tables and coordinates, certificates that settle a construction, predictions made before computing, and conclusions limited to what the scores support.
