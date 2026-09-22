# Homework and Project Plan

The date source of truth is [the lesson calendar](lesson_calendar.md), mirrored in site/schedule.json. The class meets Monday and Wednesday; holidays reduce the available meetings. There are ten problem sets, two in-class midterms, and two project presentation days.

## Evidence and workload

Problem Set 1 has four problems and targets about one hour. Later sets can ask for more sustained reasoning, but keep the central implementation small. Avoid pairing a full new implementation with several unrelated derivations.

Each set should combine a small hand calculation or mathematical argument, a computational task when useful, and an interpretation tied to a concrete output. Students submit independently justified checks and actual outputs, not simply a program that ran. A counterexample or deliberately broken implementation can be useful when the student must explain what it diagnoses.

Students may use Python, R, Julia, MATLAB, or another language. The notes provide Python examples for students who want a starting point. An assignment should say when using a ready-made implementation would bypass the algorithm being assessed. Array libraries, plotting, and file parsing are support tools.

AI assistance is allowed as described in the course policy. Students disclose material assistance, verify outputs, and explain the mathematics. Chris will rarely read code line by line; the write-up should make the result assessable through calculations, checks, and interpretation. Code remains part of reproducibility. Do not claim that a take-home task is LLM-proof, and do not introduce a graded oral check or quiz without an explicit assessment decision.

## The ten sets

| Set | Lessons | Focus | Evidence |
|---:|---|---|---|
| 1 | 1–2 | Distances, scores, k-mers | Hand lower bound; likelihood-ratio check; repair a counting bug; interpret lost order |
| 2 | 3–4 | Global alignment | Tiny score table; traceback; independently rescore the printed alignment |
| 3 | 5–6 | Shared regions and sequence search | Local coordinates; gap-open/extension calculation; seed failure example |
| 4 | 7–8 | Assembly | Graph construction; edge-use check; two candidate genomes with identical input |
| 5 | 9–10 | Markov models and CpG windows | Row-normalized counts; held-out likelihood; prior/threshold interpretation |
| 6 | 11–13 | Hidden regions and uncertainty | Joint path probability; Viterbi versus forward; posterior check by enumeration |
| 7 | 14–15 | Evolutionary trees | Small distance reconstruction; character calculation; clock assumption |
| 8 | 16–18 | Expression data, PCA, multiple testing | Centering/projection check; interpret variation; small multiple-testing calculation |
| 9 | 19–20 | Clustering and embeddings | Compare distance choices; perturb an input; separate a plot from a biological claim |
| 10 | 21–22 | Gene regulation | Check a steady state; compare simulation settings; interpret stability |

PS6 should not require students to implement Viterbi, forward, backward, and model training from scratch in one set. Supply working code for some components and assess a small change plus an independent check. Baum-Welch remains optional enrichment.

For PS8, use library SVD rather than asking for an SVD implementation. Keep the multiple-testing part to a hand-checkable table or a short supplied simulation. Neither set should become several problem sets packed under one number.

## First-set deliverable

The student-facing first set is assignments/ps01-scoring-kmers.md. Its point allocation totals 100. The questions use inline toy strings so setup and data loading do not consume the one-hour budget.

The small public reference in data/phix174-NC_001422.1.fasta is available for later word-counting or assembly exercises. Use the provenance and circular-boundary conventions in data/README.md. Label reads sampled from that reference as synthetic; they are not experimental measurements.

## Project

A feasible project uses one course method, a small public or synthetic dataset, an independent validation case, and a reproducible explanation of what the result does and does not establish.

A useful starter gives students the biological question, input format, baseline task, and an attainable validation check. Alignment, assembly, hidden-region models, evolutionary trees, and expression analysis are suitable tracks. Gene-regulation projects need more scaffolding because that material arrives late.

Release starter choices during the probability block so students can explore before spring break. Discuss scope after the second midterm. Project presentations take place April 21 and April 26. For about 25 students working mainly in pairs, plan for 12 or 13 groups, with eight minutes presenting and two minutes for questions per group. Final submission remains in finals week; Canvas supplies the exact deadline and any checkpoint dates.

A project repository should contain runnable code, appropriate data provenance, a small check with a known answer, actual results, and a README understandable to a peer. Keep restricted data, student records, instructor solutions, and private reference extracts out of the public course repository.

## Scope controls

Provide file parsing or plotting support when those steps are not the learning objective. Prefer tiny examples for correctness and a small real dataset for interpretation. Do not require full raw-read processing, large-genome assembly, deep-learning training, or building a web application.

Assessment weights and collaboration requirements belong in the course assessment policy. This document does not change those weights or impose a new group-work requirement.
