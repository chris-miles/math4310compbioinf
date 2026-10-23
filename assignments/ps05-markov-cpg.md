---
published: false
title: "Problem Set 5: Markov Chains and CpG Scores"
subtitle: "Lessons 8 and 9 · About 60 minutes"
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas. Midterm 1 does not cover Lessons 8 and 9, so it is fine to start this set after the exam.

Python, R, Julia, MATLAB, or another language is fine. You may copy and adapt code from the notes and use AI to help write or debug it. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies.

Chris will rarely read your code line by line, so put the evidence in the write-up: a calculation with a known answer, a check that would catch a specific mistake, and a short explanation of what an output means. Make a first attempt at Problems 1 and 3 without tools as exam practice; then use tools to check and correct it. Use base-2 logarithms throughout, so scores are in bits.

## 1. Probabilities that add to one (20 points, about 12 minutes)

a. Let a first-order Markov chain on $\Sigma=\{A,C,G,T\}$ have initial probabilities $q_b$ and a transition matrix $A=(a_{ab})$ whose rows each sum to one. Prove that for every $n\geq1$,
$$
\sum_{x\in\Sigma^n}P(x)=\sum_{x\in\Sigma^n}q_{x_1}\prod_{i=2}^{n}a_{x_{i-1}x_i}=1.
$$
Identify the step where the row sums are used.

b. A classmate estimates the transition matrix by dividing each pair count $n_{ab}$ by the total number of adjacent pairs in the training sequence, $N=\sum_{c,d}n_{cd}$, rather than by the row total. Call the result $M$, so $M_{ab}=n_{ab}/N$. With the uniform initial distribution $q_b=1/4$, show that the total probability assigned to all strings of length 2 is exactly $1/4$, whatever the training data (assume the training sequence has at least two bases, so $N>0$). Which step of your proof in part (a) fails?

c. Show that $M_{ab}=f_a\,\widehat a_{ab}$, where $f_a$ is the fraction of adjacent pairs that begin with $a$ and $\widehat a_{ab}$ is the row-normalized estimate from Lesson 8 (assume every base occurs at least once before another base, so each row total is positive). In two sentences, explain what extra factor this places in each window likelihood, and how that factor could distort a comparison between an island model and a background model fitted this way.

Adapted from @durbin1998.

**Evidence to submit:** the proof with the row-sum step marked, the length-2 calculation, and the factorization with its interpretation.

## 2. Does a chain predict phiX174 better than independent bases? (30 points, about 20 minutes)

Use the [phiX174 reference genome](../data/phix174-NC_001422.1.fasta) (5,386 bases, NCBI accession NC_001422.1), treated as a linear string. The first half, positions 1–2,693, is the training sequence; the second half, positions 2,694–5,386, is held out. Do not count the adjacent pair that spans the split. Fit an iid model and a first-order chain on the training half, each with one pseudocount per cell as in Lesson 8, and use the fitted base frequencies $\widehat q$ as the initial distribution for both models.

a. Make two predictions before computing. First, Lesson 9 attributes CG depletion to methylation of C in CG pairs of vertebrate DNA. phiX174 is a virus that infects the bacterium *Escherichia coli*. Predict whether the ratio $R(CG)=f(CG)/(f(C)f(G))$ from Lesson 8 is well below 1 for phiX174, and state what your prediction assumes. Second, predict whether the chain's held-out log-likelihood will exceed the iid model's by more than 0.1 bit per base.

b. Report the row of the fitted transition matrix for $C$, and check that every row sums to one. Report the number of departures from $C$ in the training half, and explain why it differs from the number of C bases in that half.

c. Report the log-likelihood per base of both models on the training half and on the held-out half. For the held-out half, also report the total difference in bits and the likelihood ratio it implies. Explain how a small difference per base and a large total can both describe the same comparison.

d. Compute $R(ab)$ for all 16 pairs over the whole genome. Report $R(CG)$, the most depleted pair, and $R(CG)$ computed separately within each half. In three or four sentences, assess both predictions from part (a): does the chain's advantage come mainly from CG, and how much weight should one genome's $R(CG)$ carry?

**Evidence to submit:** both predictions, the $C$ row with the departure count, the four log-likelihoods per base and the held-out total, the 16 ratios, and your assessment.

## 3. Prior odds and the shortest convincing window (25 points, about 15 minutes)

Use the synthetic island and background chains from Lesson 9 and their pair-score table $\beta(a,b)=\log_2\bigl(a^+_{ab}/a^-_{ab}\bigr)$. Both models use the same initial distribution, so a window's score $S(x)$ is the sum of its pair contributions.

a. A window scores $S=3$ bits. Classmate A says its island probability is $2^3/(1+2^3)\approx0.89$. Classmate B says that if 1% of windows are islands, the probability is about 0.075. Identify the assumption behind A's number and check B's.

b. With prior island probability 0.01 and equal costs for the two kinds of error, find the Bayes threshold in bits from Lesson 9.

c. Construct the shortest window whose score exceeds that threshold. Certify its score by listing its pair contributions, and prove that no shorter window exceeds the threshold.

d. Construct a window of length 6 containing at least four C or G bases whose score is negative. Which entries of the $\beta$ table make this possible, and what does your window say about calling a GC-rich stretch an island under these models?

**Evidence to submit:** the two posterior probabilities, the threshold, the certified shortest window with the bound on shorter ones, and the GC-rich negative window.

## 4. Reading a threshold table (25 points, about 13 minutes)

To evaluate the window classifier where the true labels are known, we simulated a synthetic 3,000-base sequence. Bases 726–1,125 and 2,026–2,325 were generated from Lesson 9's island chain, and all other bases from the background chain. The sequence was cut into consecutive, nonoverlapping windows (1–20, 21–40, and so on), and each window was scored with Lesson 9's `markov_log_ratio`. A window is *pure* if it lies entirely inside one region and *straddling* otherwise. A window is called an island when $S>c$.

```python
# Generator used for Problem Set 5, Problem 4 (you do not need to rerun it)
# transition_plus, transition_minus, and markov_log_ratio come from Lesson 9.
import random

def simulate_regions(n_bases: int, islands: list[tuple[int, int]],
                     rng: random.Random) -> tuple[str, list[bool]]:
    """Simulate a sequence whose transition matrix follows known island labels."""
    labels = [any(lo <= p <= hi for lo, hi in islands) for p in range(1, n_bases + 1)]
    sequence = rng.choices("ACGT", k=1)
    for position in range(2, n_bases + 1):
        matrix = transition_plus if labels[position - 1] else transition_minus
        row = matrix[sequence[-1]]
        sequence += rng.choices("ACGT", weights=[row[b] for b in "ACGT"], k=1)
    return "".join(sequence), labels

sequence, labels = simulate_regions(3000, [(726, 1125), (2026, 2325)], random.Random(4310))
```

With 20-base windows there are 33 pure island windows and 113 pure background windows:

| Threshold $c$ | TP | FN | FP | TN |
|---:|---:|---:|---:|---:|
| −2 | 31 | 2 | 21 | 92 |
| 0 | 27 | 6 | 10 | 103 |
| 2 | 27 | 6 | 2 | 111 |
| 4 | 23 | 10 | 0 | 113 |

: Pure 20-base windows classified at four thresholds. {#tbl-ps5-threshold-counts}

The four straddling 20-base windows are 721–740 ($S=8.5$), 1,121–1,140 ($S=-4.2$), 2,021–2,040 ($S=11.2$), and 2,321–2,340 ($S=-3.6$). With 50-base windows instead, all 12 pure island and 44 pure background windows are classified correctly at $c=0$. The four straddling 50-base windows are 701–750 ($S=8.0$), 1,101–1,150 ($S=3.1$), 2,001–2,050 ($S=12.0$), and 2,301–2,350 ($S=3.1$).

a. For the 20-base windows, compute sensitivity, specificity, and the total number of misclassified pure windows at each threshold. Treat the fraction of pure windows that are islands as the prior and compute the Bayes threshold. Which row is closest to it, and does that row have the fewest errors?

b. A classmate reads the $c=0$ row and says, "Specificity is 0.91, so 91% of the windows called islands really are islands." Assess the claim with a calculation from the table.

c. Another classmate concludes that 50-base windows are better in every respect. Count the background bases labeled island inside straddling windows at $c=0$ for each width, and explain the tradeoff between window width and locating an island's boundaries.

d. In one or two sentences: these sequences were simulated from the same two chains used to score them. Why would labeled genomic data probably give worse numbers?

**Evidence to submit:** the completed table of sensitivity, specificity, and errors, the Bayes threshold, the calculation for the classmate's claim, the base counts for the two widths, and your explanation.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: proofs that name the assumption they use, predictions compared honestly with outputs, certified constructions, and conclusions limited to what the models and data support.

::: {#refs}
**References**
:::
