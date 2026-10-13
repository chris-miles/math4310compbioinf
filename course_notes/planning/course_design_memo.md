# Course Design Memo

## Working Thesis

MATH 4310 should be a theory-first bioinformatics course, not a general data science course and not an applied genomics tools course. The course should teach students why core bioinformatics methods work by building small versions of the algorithms and models from scratch.

The right center of gravity is:

1. Dynamic programming for biological sequence comparison.
2. Probabilistic models for sequence annotation.
3. Graph algorithms for genome assembly.
4. Phylogenetic algorithms as tree inference from sequence data.
5. Linear algebra and statistics for high-dimensional omics data.
6. A compact systems-biology module showing how mechanistic models differ from data-driven models.

This is appropriate for a mathematics department because the recurring objects are recurrences, likelihoods, graphs, Markov chains, matrix factorizations, multiple testing procedures, and dynamical systems. Biology supplies the motivating data and the failure modes.

## Audience

For some students, this may be their only course in computational bioinformatics methods. When choosing examples and allocating time, favor ideas and technical practice that can carry to other biological problems as technologies change: representing data, formulating a model or objective, deriving a computation, and checking what its result supports. A specific biological problem gives those ideas something concrete to explain.

Where useful, connect that practice to working with AI-assisted analysis: students still need to understand the inputs, examine assumptions, validate outputs, and judge biological evidence when software generates code or proposes a method. These are guiding considerations, not required themes for every lesson or promises that a particular algorithm will remain current.

Expected students:

- Bioinformatics BS students with basic cell biology and programming.
- Mathematics or applied mathematics students interested in biological data.
- CS students who want a mathematically grounded introduction to computational biology.
- Biology students with enough programming maturity to write small scripts.

The course should not assume students already know command-line bioinformatics workflows, Bioconductor, scikit-learn, BLAST internals, RNA-seq pipelines, or machine learning. Those can be referenced, but they should not be the instructional core. The BIOL 3150 syllabus suggests that local students may already see IGV, BLAST, Colab/Python, pandas, variant analysis, RNA-seq, GWAS, and visualization in an applied genomics setting; MATH 4310 should complement that course by emphasizing the mathematical mechanisms rather than repeating the workflow training.

## Prerequisites and scaffolding

The student-facing prerequisites are MATH 2270, CS 1410, and BIOL 2020 or BIOL 2021. Plan for basic computational linear algebra and programming, basic cell biology, and little or no bioinformatics. Formal probability, statistics, algorithms, and genetics are not assumed.

Introduce the biological vocabulary and probability or statistical ideas needed for each method locally. Students already know matrix multiplication and loops; use those skills rather than reteaching them. A short readiness check can identify who needs help with functions, string indexing, dictionaries, and table indexing before the first dynamic-programming assignment. It is a support tool, not an additional graded assessment.

## Topic Balance

The allocation keeps foundational methods central while making room for biological comparison and noisy dynamics.

Keep as central:

- Alignment and scoring.
- Dynamic programming.
- de Bruijn graphs and genome assembly.
- Markov chains and HMMs.
- Phylogenetics.
- PCA/SVD.
- Multiple testing and FDR.

Expand slightly:

- Seed-and-extend / BLAST intuition.
- Read mapping as a modern extension of alignment.
- Count matrices, normalization, and batch effects.
- FDR as a core omics method, not a final-week afterthought.

Dynamics scope:

- Retain three dynamics lessons: production/degradation ODEs, biological switches and cell states, and biological dynamics with noisy measurements. Separate process and measurement noise; introduce a state-space bridge without Kalman filtering.

Mention but do not deeply teach:

- UMAP and t-SNE internals.
- Single-cell RNA-seq pipelines.
- Maximum-likelihood phylogenetics.
- Profile HMMs.
- Deep learning.
- Pangenome graphs.

These are useful synthesis or optional enrichment topics. They should not displace the foundational methods.

## Module Order

The approved lesson order is:

1. Sequence algorithms and graph assembly.
2. Probabilistic sequence models.
3. Phylogenetics.
4. Omics linear algebra and statistics.
5. Compact systems-biology modeling.
6. Modern Frontiers and a Flex meeting.

This order lets the first half build a coherent "sequences as structured strings" story. The second half then moves from trees to matrices to models.

## What Students Should Be Able To Do

By the end, a strong student should be able to:

- Derive a dynamic-programming recurrence from an alignment objective.
- Prove a recurrence correct by explaining its optimal-substructure claim.
- Implement global/local alignment on small examples.
- Explain how scoring matrices are log-odds summaries of evolutionary assumptions.
- Build a de Bruijn graph from reads and identify why repeats break assembly.
- Define an HMM and run Viterbi / forward decoding in log space.
- Explain why posterior decoding and Viterbi decoding can disagree.
- Construct a simple phylogenetic tree by distance or parsimony.
- Center an expression matrix and interpret PCA through SVD.
- Compare biological conditions using effect size, independent replication, and a permutation null.
- Apply Benjamini-Hochberg and explain what FDR controls.
- Read an embedding or clustering plot skeptically.
- Simulate a simple ODE model and interpret stability and biological switches.
- Distinguish biological process variability from noisy measurement in a dynamical model.
- Choose an appropriate method for a small biological-data problem and name likely failure modes.

## Main Design Risks

Risk: the course becomes a survey.

Control: every week should have one mathematical object, one algorithm/model, one implementation, and one failure mode. Avoid detached biological storytelling without a computational method attached.

Risk: programming overwhelms the math.

Control: keep implementations short, with starter I/O and tests. Grade algorithmic logic, not software polish.

Risk: biology students lack proof experience.

Control: use proof templates for dynamic programming and HMM identities. Ask for structured explanations rather than abstract theorem-proof style early in the term.

Risk: math students lack biological intuition.

Control: every module starts with a biological question and a toy dataset. Failure-mode reflections should force students to connect assumptions to real data.

Risk: an optional project, if adopted, sprawls.

Control if adopted: bound the project to one course method, one documented dataset, one validation check, and one limitation. Define grading and workload before assigning it.

## Current lesson allocation and assessment decisions

There are 24 lessons: alignment/scoring 5, assembly 2, Markov chains/HMMs 5, phylogenetics 2, omics 6, biological dynamics 3, and Modern Frontiers 1. Two midterms and an April 26 Flex day complete the 27 meetings.

Lesson 17, Comparing Biological Conditions, uses expression measurements to teach one-feature effect size, biological replication and variability, a null hypothesis, and label permutation. Lesson 18 begins where that ends: simultaneous discoveries, FDP/FDR, and Benjamini–Hochberg. Neither lesson requires a full differential-expression assay workflow.

Lesson 21 develops gene-regulation ODEs. Lesson 22 connects nonlinear feedback, stability, and switches to cell states. Lesson 23 distinguishes process noise and measurement noise, briefly separates intrinsic from extrinsic biological variability, and connects hidden dynamics to noisy observations through a state-space formulation. Kalman filtering is outside scope. Lesson 24 retains Modern Frontiers.

Ten problem sets and two midterms are scheduled. Midterm 1 covers Lessons 1–7 and Midterm 2 covers Lessons 8–14. Projects, assessment weights, and any final assessment remain undecided; no presentation days are reserved. See `lesson_calendar.md` and `assessment_model.md`.
