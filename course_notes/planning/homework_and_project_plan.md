# Homework and Project Plan

The date source of truth is [the lesson calendar](lesson_calendar.md), mirrored in site/schedule.json. The class meets Monday and Wednesday; holidays reduce the available meetings. There are ten problem sets, two in-class midterms, and one Flex day on April 26; project format remains undecided.

Schedule each midterm one week after the last homework deadline for its covered material. Later lessons can proceed during that interval without entering the exam coverage. Midterm 1 covers Lessons 1–7; Midterm 2 covers Lessons 8–14.

## Evidence and workload

Problem Set 1 has four problems and targets about one hour. Later sets can ask for more sustained reasoning, but keep the central implementation small. Avoid pairing a full new implementation with several unrelated derivations.

Each set should combine a small hand calculation or mathematical argument, a computational task when useful, and an interpretation tied to a concrete output. Students submit independently justified checks and actual outputs, not simply a program that ran. A counterexample or deliberately broken implementation can be useful when the student must explain what it diagnoses.

Students may use Python, R, Julia, MATLAB, or another language. The notes provide Python examples for students who want a starting point. An assignment should say when using a ready-made implementation would bypass the algorithm being assessed. Array libraries, plotting, and file parsing are support tools.

AI assistance is allowed as described in the course policy. Students disclose material assistance, verify outputs, and explain the mathematics. Chris will rarely read code line by line; the write-up should make the result assessable through calculations, checks, and interpretation. Code remains part of reproducibility. Do not claim that a take-home task is LLM-proof, and do not introduce a graded oral check or quiz without an explicit assessment decision.

## The ten sets

| Set | Lessons | Focus | Evidence |
|---:|---|---|---|
| 1 | 1–2 | Distances, scores, k-mers | Hand lower bound; likelihood-ratio check; repair a counting bug; interpret lost order |
| 2 | 3 | Global alignment | Tiny score table; traceback; independently rescore the printed alignment |
| 3 | 4–5 | Shared regions and sequence search | Local coordinates; gap-open/extension calculation; seed failure example |
| 4 | 6–7 | Assembly | Graph construction; edge-use check; two candidate genomes with identical input |
| 5 | 8–9 | Markov models and CpG windows | Row-normalized counts; held-out likelihood; prior/threshold interpretation |
| 6 | 10–12 | Hidden regions and uncertainty | Joint path probability; Viterbi versus forward; posterior check by enumeration |
| 7 | 13–14 | Evolutionary trees | Small distance reconstruction; character calculation; clock assumption |
| 8 | 15–18 | PCA, comparing conditions, multiple testing | Centering/projection check; one-feature effect and permutation calculation; small BH table |
| 9 | 19–20 | Clustering and embeddings | Compare distance choices; perturb an input; separate a plot from a biological claim |
| 10 | 21–23 | Gene regulation and noisy observations | Check a steady state; interpret stability; short supplied-data noise comparison |

PS6 should not require students to implement Viterbi, forward, backward, and model training from scratch in one set. Supply working code for some components and assess a small change plus an independent check. Baum-Welch remains optional enrichment.

PS2 supplies a working traceback function. Required work combines a hand traceback, compact score-table implementation, and independent rescoring/input-reconstruction checks; implementing traceback from scratch and comparing scoring systems are enrichment. This keeps the assignment within Lesson 3's main teaching scope.

PS5 and PS8 are released shortly before midterms covering earlier material. Keep both sets compact and allow students to begin the new work after the exam; the release date does not imply immediate completion. PS8's four lessons should produce one small connected analysis rather than four separate implementations.

Lesson 17 develops a single-feature comparison, biological replication, effect size, and a permutation null. Lesson 18 handles many comparisons, FDP/FDR, and BH. PS8 should reinforce that boundary without adding a full assay workflow.

For PS8, use library SVD rather than asking for an SVD implementation. Keep the multiple-testing part to a hand-checkable table or a short supplied simulation. Neither set should become several problem sets packed under one number.

## First-set deliverable

The student-facing first set is assignments/ps01-scoring-kmers.md. Its point allocation totals 100. The questions use inline toy strings so setup and data loading do not consume the one-hour budget.

The small public reference in data/phix174-NC_001422.1.fasta is available for later word-counting or assembly exercises. Use the provenance and circular-boundary conventions in data/README.md. Label reads sampled from that reference as synthetic; they are not experimental measurements.

## Optional project ideas

Projects remain undecided. The following ideas do not establish project requirements, teams, deadlines, or grade weights.

A feasible project uses one course method, a small public or synthetic dataset, an independent validation case, and a reproducible explanation of what the result does and does not establish.

A useful starter gives students the biological question, input format, baseline task, and an attainable validation check. Alignment, assembly, hidden-region models, evolutionary trees, and expression analysis are suitable tracks. Gene-regulation projects need more scaffolding because that material arrives late.

If a project is adopted, release curated choices early enough for students to explore the data and method. No project checkpoints, presentation dates, final submission deadline, or team requirement is currently scheduled.

If adopted, a project repository should contain runnable code, appropriate data provenance, a small check with a known answer, actual results, and a README understandable to a peer. Keep restricted data, student records, instructor solutions, and private reference extracts out of the public course repository.

## Scope controls

Provide file parsing or plotting support when those steps are not the learning objective. Prefer tiny examples for correctness and a small real dataset for interpretation. Do not require full raw-read processing, large-genome assembly, deep-learning training, or building a web application.

Assessment weights and collaboration requirements belong in the course assessment policy. This document does not change those weights or impose a new group-work requirement.

PS10 core problems are released April 14 and due April 21. A brief Lesson 23 supplied-data interpretation can be added April 19; if the set already fills its workload budget, keep that exercise as ungraded class practice. Avoid a new stochastic simulator with a two-day deadline.
