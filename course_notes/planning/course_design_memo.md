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

Expected students:

- Bioinformatics BS students with genetics and some programming.
- Mathematics or applied mathematics students interested in biological data.
- CS students who want a mathematically grounded introduction to computational biology.
- Biology students with enough programming maturity to write small scripts.

The course should not assume students already know command-line bioinformatics workflows, Bioconductor, scikit-learn, BLAST internals, RNA-seq pipelines, or machine learning. Those can be referenced, but they should not be the instructional core. The BIOL 3150 syllabus suggests that local students may already see IGV, BLAST, Colab/Python, pandas, variant analysis, RNA-seq, GWAS, and visualization in an applied genomics setting; MATH 4310 should complement that course by emphasizing the mathematical mechanisms rather than repeating the workflow training.

## Prerequisite Recommendation

The current proposed prerequisites are:

- MATH 2270 Linear Algebra.
- CS 2420 Algorithms and Data Structures.
- BIOL 2030 Genetics.

That is academically clean but probably too restrictive if the goal is to serve the new bioinformatics major and math students. A better catalog stance would be:

Required:

- MATH 2270, or equivalent linear algebra.
- BIOL 2030, or equivalent genetics / molecular biology.
- CS 1410, or equivalent programming experience.

Recommended:

- CS 2420.
- BIOL 3150.
- One prior probability or statistics course.

Rationale:

- MATH 2270 should stay required. PCA/SVD and high-dimensional geometry need real linear algebra.
- BIOL 2030 should stay required or be replaceable by instructor consent. Students need sequence, gene, transcription, mutation, and inheritance vocabulary.
- CS 2420 is helpful but not essential if the course scaffolds dynamic programming, graph traversal, and dictionary/list-based implementations carefully.
- A probability/statistics prerequisite would be useful but may block too many students. The plan includes a probability ramp in Weeks 5-7 and a statistics ramp in Weeks 10-11.

If CS 2420 is relaxed, the course needs a Week 1 programming and algorithmic-thinking diagnostic. The diagnostic should check whether students can:

- Write a function.
- Loop over a string.
- Use a dictionary or map.
- Build and index a 2D table.
- Explain the difference between O(n), O(nm), and exponential search.

Students who cannot do these should be directed to a bootcamp handout or office-hour support before the first dynamic-programming assignment.

This diagnostic should be short. It is a readiness screen and support trigger, not a replacement for BIOL 3150's introductory genomic data science material.

## Topic Balance

The current syllabus has the right major topics. The main revision is weighting.

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

Compress:

- Gene regulation ODEs, Hill functions, and toggle switches. Keep one week. This is valuable mathematics, but it is closer to mathematical biology / systems biology than to the core bioinformatics theory spine.

Mention but do not deeply teach:

- UMAP and t-SNE internals.
- Single-cell RNA-seq pipelines.
- Maximum-likelihood phylogenetics.
- Profile HMMs.
- Deep learning.
- Pangenome graphs.

These are good synthesis or project topics. They should not displace the foundational methods.

## Recommended Module Order

The current syllabus order places ODEs before phylogenetics and omics. I recommend this order instead:

1. Sequence algorithms and graph assembly.
2. Probabilistic sequence models.
3. Phylogenetics.
4. Omics linear algebra and statistics.
5. Compact systems-biology modeling.
6. Synthesis and project work.

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
- Apply Benjamini-Hochberg and explain what FDR controls.
- Read an embedding or clustering plot skeptically.
- Simulate a simple ODE model and identify where mechanistic assumptions enter.
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

Risk: final projects sprawl.

Control: require one course method, one documented dataset, one validation check, and one limitation. Keep project grading focused on reproducibility and explanation.

## Assessment Model Update

The better assessment structure is probably two shorter midterms plus a final project, with no cumulative final exam.

Recommended weights:

- Problem sets and project workshops: 40%.
- Short Midterm 1: 15%.
- Short Midterm 2: 15%.
- Final project: 30%.

This is better aligned than a midterm plus final exam because the final project is the real cumulative task. Students must choose a mathematical object, apply a course method to biological data, validate the result, and communicate assumptions and failure modes. A timed cumulative final would likely either overemphasize early DP/HMM mechanics or become a shallow survey of late-semester topics.

Recommended exam placement:

- Short Midterm 1 after the sequence-algorithms and assembly block.
- Short Midterm 2 after probabilistic models and phylogenetics, optionally including PCA if it has already been introduced.

See `assessment_model.md` for details.
