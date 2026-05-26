# Assessment Model

## Recommendation

Use two shorter in-class exams and no final exam. Make the final project the finals-week summative assessment.

This fits the course better than one midterm plus a cumulative final because the course is not a linear march toward one final theorem. It is a sequence of mathematical lenses for biological data:

1. Sequence algorithms and graph assembly.
2. Probabilistic sequence models.
3. Phylogenetic trees.
4. Omics matrices and high-throughput statistics.
5. Mechanistic models and synthesis.

The final project is a better cumulative instrument than a final exam because it asks students to choose a method, apply it to data, validate it, document assumptions, and explain failure modes. That is closer to the course objective than asking them to reproduce late-semester derivations under time pressure.

## Proposed Grade Weights

Recommended:

| Component | Weight |
|---|---:|
| Problem sets | 40% |
| Short Midterm 1 | 15% |
| Short Midterm 2 | 15% |
| Final project | 30% |

Why this works:

- The total exam weight drops from 45% to 30%.
- The project rises from 15% to 30%, matching its role as the cumulative assessment.
- Problem sets stay at 40%, preserving weekly technical accountability.

Slightly more conservative alternative:

| Component | Weight |
|---|---:|
| Problem sets | 40% |
| Short Midterm 1 | 17.5% |
| Short Midterm 2 | 17.5% |
| Final project | 25% |

Use the conservative version if the department wants exams to remain a visibly large part of the course grade.

## Exam Placement

### Short Midterm 1

Best placement:

- Around instructional Week 5, ideally after the assembly module.

Coverage:

- Algorithmic complexity.
- Log-odds scoring.
- Needleman-Wunsch.
- Smith-Waterman.
- Affine gaps at the recurrence level.
- Seed-and-extend intuition.
- de Bruijn graphs and Eulerian paths.

What the exam should test:

- Derive or complete a recurrence.
- Fill a small DP table.
- Explain a traceback or tie.
- Construct a small de Bruijn graph.
- Diagnose a repeat/error failure mode.
- Do a short complexity analysis.

What it should not test:

- Long code.
- Memorized names of tools.
- Large biological vocabulary.

### Short Midterm 2

Best placement:

- Around instructional Week 10 or early Week 11.

Cleanest coverage:

- Markov chains.
- CpG likelihood-ratio classification.
- HMM definitions.
- Viterbi.
- Forward/backward/posterior decoding.
- Conceptual HMM training.
- Phylogenetics: UPGMA and parsimony.

Optional coverage if scheduled after the PCA introduction:

- Centering an omics matrix.
- PCA/SVD interpretation.
- Scores, loadings, variance explained.

What the exam should test:

- Compute a likelihood in log space.
- Trace one or two steps of Viterbi or forward.
- Explain why Viterbi and posterior decoding can disagree.
- Run one UPGMA merge or one Fitch parsimony calculation.
- Interpret a model assumption.

What it should not test:

- Baum-Welch algebra in full detail.
- UMAP/t-SNE internals.
- Full differential-expression workflows.
- ODE numerical methods.

## Why No Final Exam Is Better Here

The late-semester material is better assessed through written and computational work:

- PCA/SVD and FDR need interpretation of data products, not just timed computation.
- Clustering and embeddings need skepticism about plots and parameter choices.
- ODEs are included to contrast mechanistic modeling with data-driven bioinformatics, not to become a full differential equations unit.
- The final project naturally forces cumulative method choice and communication.

A cumulative final would either:

- Overweight early DP/HMM mechanics students have already been tested on, or
- Become a broad survey exam with shallow questions on omics, embeddings, ODEs, and projects.

Neither is as aligned as a serious final project.

## Consequences for Homework

Two exams plus a larger project means 12 full problem sets may be too much unless each one is short. There are three viable models.

### Model A: Keep 12 Short Problem Sets

Use if weekly rhythm is important.

Design:

- Keep all 12.
- Make PS1 and PS12 lighter.
- Drop the lowest individual problem-set score.
- Keep the triad HMM set.

Risk:

- Grading load and student load are high when the project ramps up.

### Model B: Use 10 Full Problem Sets

Use if project quality is the priority.

Design:

- Keep PS1-PS8 as technical foundations.
- Combine PCA and FDR into one larger omics problem set.
- Combine clustering/ODE/synthesis into one final short problem set.
- Keep two project checkpoints.

Risk:

- Slightly less weekly accountability late in the semester.

### Model C: 10 Problem Sets Plus 2 Project Workshops

This is probably the best version.

Design:

- Ten standard problem sets.
- Project Workshop 1 replaces a late weekly problem set: proposal, repo skeleton, validation plan.
- Project Workshop 2 replaces another late weekly problem set: preliminary result, README draft, limitation/failure-mode audit.

This makes the project real without adding hidden workload.

Recommended late-semester sequence under Model C:

| Week | Main assessment |
|---|---|
| 10 | Short Midterm 2 |
| 11 | Omics/PCA/FDR problem set plus Project Workshop 1 |
| 12 | Clustering/embedding interpretation problem set |
| 13 | ODE/synthesis short problem set |
| 14 | Project Workshop 2 |
| Finals Week | Final project repository |

## Suggested Revised Grade Weights Under Model C

| Component | Weight |
|---|---:|
| Standard problem sets, 10 total | 35% |
| Project workshops/checkpoints, 2 total | 5% |
| Short Midterm 1 | 15% |
| Short Midterm 2 | 15% |
| Final project repository | 30% |

This keeps weekly mathematical work at 40% total, keeps exams meaningful but not dominant, and gives the final project enough weight to justify the no-final-exam structure.

## Syllabus Language Draft

Possible replacement for the exam section:

> This course has two shorter in-class midterm exams rather than a cumulative final exam. The exams check individual mastery of the core mathematical and algorithmic tools at two natural points in the semester. The final project serves as the cumulative assessment: students apply one course method to a documented biological-data problem, validate the result, and explain assumptions and limitations in a reproducible repository.

Possible grade table:

| Component | Weight |
|---|---:|
| Problem sets and project workshops | 40% |
| Short Midterm Exam 1 | 15% |
| Short Midterm Exam 2 | 15% |
| Final project | 30% |

## Bottom Line

Yes: two shorter midterms plus a final project is better aligned with this course.

The strongest structure is:

- Midterm 1 after sequence algorithms and assembly.
- Midterm 2 after probabilistic models and phylogenetics, optionally including PCA if scheduled early enough.
- No final exam.
- Finals-week final project as the cumulative assessment.

