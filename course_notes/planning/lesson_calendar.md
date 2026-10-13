# Lesson Calendar

Spring 2027, Monday and Wednesday, 80 minutes per meeting. The public schedule in `site/schedule.json` is the date source of truth.

The teaching term is January 11 through April 27. January 18 and February 15 are Monday holidays; spring break removes March 8 and 10. This leaves **27 meetings: 24 lessons, two midterms, and one Flex day on April 26**. Projects remain undecided; no presentation days are reserved.

Calendar checked against the [University calendar](https://registrar.utah.edu/academic-calendars/pdf-academic-calendars/main_and_online_2026_2027.pdf).

| Week | Date | Meeting |
|---:|---|---|
| 1 | 2027-01-11 | Lesson 1: Sequences, Distances, and Alignments |
| 1 | 2027-01-13 | Lesson 2: Log-Odds Scoring |
| 2 | 2027-01-20 | Lesson 3: Global Alignment: Dynamic Programming |
| 3 | 2027-01-25 | Lesson 4: Local Alignment |
| 3 | 2027-01-27 | Lesson 5: Affine Gaps and Seed-and-Extend |
| 4 | 2027-02-01 | Lesson 6: Reads, k-mers, and de Bruijn Graphs |
| 4 | 2027-02-03 | Lesson 7: Eulerian Paths, Repeats, and Assembly Errors |
| 5 | 2027-02-08 | Lesson 8: Probability and Markov Chains for DNA |
| 5 | 2027-02-10 | Lesson 9: CpG Islands and Likelihood Ratios |
| 6 | 2027-02-17 | Lesson 10: Hidden Markov Models |
| 7 | 2027-02-22 | Lesson 11: Viterbi Decoding |
| 7 | 2027-02-24 | Midterm 1: Lessons 1–7 |
| 8 | 2027-03-01 | Lesson 12: Uncertainty in HMMs: Forward-Backward |
| 8 | 2027-03-03 | Lesson 13: Evolutionary Distances and UPGMA |
| 9 | March 6 to 14 | Spring break; no class |
| 10 | 2027-03-15 | Lesson 14: Parsimony and Fitch’s Algorithm |
| 10 | 2027-03-17 | Lesson 15: Gene Expression and PCA |
| 11 | 2027-03-22 | Lesson 16: PCA in Practice |
| 11 | 2027-03-24 | Lesson 17: Comparing Biological Conditions |
| 12 | 2027-03-29 | Lesson 18: Multiple Testing and False Discovery Rates |
| 12 | 2027-03-31 | Midterm 2: Lessons 8–14 |
| 13 | 2027-04-05 | Lesson 19: Distances, Clustering, and Neighborhood Graphs |
| 13 | 2027-04-07 | Lesson 20: Single-Cell Data and Embeddings |
| 14 | 2027-04-12 | Lesson 21: ODE Models of Gene Regulation |
| 14 | 2027-04-14 | Lesson 22: Biological Switches and Cell States |
| 15 | 2027-04-19 | Lesson 23: Biological Dynamics and Noisy Measurements |
| 15 | 2027-04-21 | Lesson 24: Modern Frontiers in Bioinformatics |
| 16 | 2027-04-26 | Flex |

## Assessment timing

Midterm 1 is February 24 and covers Lessons 1–7; assembly homework is due February 17. Midterm 2 is March 31 and covers Lessons 8–14; phylogenetics homework is due March 24. Each exam follows its last relevant homework deadline by one week. Later lessons proceed during the study interval without entering exam coverage.

Ten problem sets are released after their required material is taught. The week containing an assignment entry in the site schedule is its release week. Canvas is the submission authority.

| Set | Release | Due | Coverage |
|---:|---|---|---|
| 1 | 2027-01-13 | 2027-01-25 | Scoring, distances, k-mers (1–2) |
| 2 | 2027-01-20 | 2027-02-03 | Global alignment, including traceback and correctness (3) |
| 3 | 2027-02-03 | 2027-02-10 | Local alignment, gap models, seeds (4–5) |
| 4 | 2027-02-08 | 2027-02-17 | Assembly (6–7) |
| 5 | 2027-02-22 | 2027-03-03 | Markov chains and CpG windows (8–9) |
| 6 | 2027-03-03 | 2027-03-17 | Hidden regions, decoding, uncertainty (10–12) |
| 7 | 2027-03-15 | 2027-03-24 | Evolutionary trees (13–14) |
| 8 | 2027-03-29 | 2027-04-07 | PCA, biological comparisons, multiple testing (15–18) |
| 9 | 2027-04-07 | 2027-04-14 | Clustering and embeddings (19–20) |
| 10 | 2027-04-14 (supplement April 19) | 2027-04-21 | Gene regulation, switches, and noisy measurements (21–23) |

PS10 should remain compact: release the regulation and stability problems April 14, then add a short supplied-data interpretation after Lesson 23 on April 19. If the set already fills its workload budget, keep that exercise as ungraded class practice. Its April 21 deadline leaves only two days for that final part; avoid requiring a new stochastic simulation implementation. April 26 remains Flex for review, completion, or an activity chosen later. Project format, checkpoints, deadlines, and assessment weights remain undecided.

## Boundaries of the new lessons

Lesson 17 uses gene-expression measurements to compare one feature across biological conditions: effect size, independent biological replication, variability, a null hypothesis, and a hand-checkable label-permutation test. It does not teach an assay pipeline or simultaneous testing. Lesson 18 starts from many such comparisons and studies false discoveries, FDP, FDR, and Benjamini–Hochberg.

Lessons 21–23 form a three-meeting dynamics block. Lesson 21 develops production/degradation ODEs; Lesson 22 connects nonlinear feedback and stability to biological switches and cell states; Lesson 23 separates process noise from measurement noise and connects hidden biological dynamics to observations through a simple state-space formulation. Intrinsic and extrinsic biological variability are briefly distinguished from the process/measurement distinction. Kalman filtering is outside scope. Lesson 24 retains Modern Frontiers.
