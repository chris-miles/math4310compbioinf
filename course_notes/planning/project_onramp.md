# Final Project On-Ramp

Status: exploratory internal planning. This document gives options for making the final project workable, but the student-facing syllabus should stay more flexible until the project format settles.

## The Problem

A 30% final project is only fair if students can begin before the last third of the course. But students do not know HMMs, PCA, FDR, phylogenetics, or ODEs at the start. If the project begins only after those topics are taught, the project becomes rushed. If it begins too early, students are choosing methods they do not yet understand.

The solution is to separate project work into stages:

1. Early: understand the menu of possible problem types.
2. Early-middle: practice matching biological questions to mathematical objects.
3. Middle: choose a project track with a known baseline method.
4. Late: add the specific method details after the relevant course module.

Students do not need to know the whole course to start. They need a bounded project track, a curated dataset, and a baseline deliverable whose mathematical method is taught by the time serious implementation begins.

## Core Design Principle

Projects should be track-based, not open-ended research proposals.

Each project starter should specify:

- Biological question.
- Curated small dataset.
- Course method.
- What students can do before learning the method.
- Baseline deliverable.
- Validation check.
- Allowed extensions.
- Explicit out-of-scope work.

This lets students start with problem/data orientation while the technical method is still being taught.

## Project Timeline

### Week 1: Project Gallery, Not Project Selection

Goal:

- Let students see where the course is going without asking them to commit.

In class:

- Show 6-8 one-slide project examples.
- For each example, name the mathematical object: string pair, graph, hidden-state model, tree, matrix, test list, ODE system.
- Emphasize that final projects will use curated starter datasets and course methods.

Deliverable:

- None, or a one-question exit ticket: "Which project type looked most interesting, and why?"

### Week 3 or 4: Method Card 1

Goal:

- After alignment and assembly, students can already understand some viable projects.

Short assignment or in-class activity:

- Students complete one "method card" for either alignment or assembly.

Method card fields:

- What is the input?
- What is the output?
- What mathematical object does the method optimize or construct?
- What assumption can fail?
- What would be a small validation check?

Deliverable:

- Ungraded or 1-point completion.

### Week 6 or 7: Project Starter Menu Released

Goal:

- Give students enough information to form preferences before teams form.

Release:

- A curated list of project tracks.
- One-page starter sheet for each track.
- Tiny example data for each track.
- Expected baseline and extension options.

Deliverable:

- Individual preference form ranking three tracks.

### Week 8: Individual Project Triage Memo

Goal:

- Make students think before teams form.

Prompt:

- Pick two possible project tracks.
- For each, identify the input data, method, baseline output, validation check, and likely failure mode.
- Say what background you still need from the course.

Deliverable:

- One-page individual memo.
- Counts lightly as project-preparation credit or participation.

### Spring Break / Week 9: Team Formation

Goal:

- Teams form with some project preferences already visible.

Process:

- Instructor forms teams using preference forms, or students self-form within project tracks.
- Teams choose one track plus a backup.
- Teams do not yet need a polished research question.

Deliverable:

- Team and track declaration.

### Week 10: Project Contract

Goal:

- Convert the selected track into a concrete scope.

In class:

- Short project clinic after Short Midterm Exam 2 or during the next meeting.
- Each team fills a project contract.

Contract fields:

- Dataset.
- Course method.
- Baseline deliverable.
- Validation check.
- Division of labor.
- One likely limitation.
- One optional extension, explicitly marked optional.

Deliverable:

- Contract draft, used to prepare Checkpoint 1.

### Week 11: Checkpoint 1

Goal:

- Lock project scope.

Deliverable:

- One-page proposal.
- Repository skeleton.
- Data source/provenance.
- Planned validation check.
- Minimal run instructions.

### Week 12 or 13: Micro-Demo

Goal:

- Prevent silent failure.

Deliverable:

- A tiny run of the method on toy data or the first 5-10 samples/sequences.
- One screenshot, table, or output file.
- One paragraph: what failed or changed.

This can be an in-class check rather than a formal graded checkpoint.

### Week 14: Checkpoint 2

Goal:

- Ensure preliminary results exist before finals week.

Deliverable:

- Preliminary result.
- README draft.
- Validation check result or validation plan.
- Limitation/failure-mode paragraph.

### Finals Week: Final Repository

Goal:

- Cumulative assessment.

Deliverable:

- Reproducible repository.
- README.
- Code.
- Data provenance.
- Validation check.
- Team contributions paragraph.
- Revision memo.

## Recommended Project Tracks

### Track A: Alignment and Homology

When students can start:

- Week 3.

Course method:

- Needleman-Wunsch, Smith-Waterman, affine gaps, or seed-and-extend.

What students can do early:

- Choose a small gene/protein family from curated FASTA files.
- Run simple pairwise comparisons.
- Define scoring schemes.

Baseline deliverable:

- From-scratch aligner applied to selected sequences.
- Compare at least two scoring choices.

Validation check:

- Hand-check a tiny alignment.
- Compare one result to an established tool qualitatively, without using that tool as the implementation.

Good extension:

- Add affine gaps or a seed index.

Out of scope:

- Full BLAST clone.
- Large database search.

### Track B: de Bruijn Assembly

When students can start:

- Week 4.

Course method:

- de Bruijn graph construction and Eulerian paths.

What students can do early:

- Inspect reads.
- Count k-mers.
- Build graph summaries.

Baseline deliverable:

- Toy assembler for short reads.
- Compare assemblies across k values.

Validation check:

- Use a simulated genome with known answer.
- Check whether reconstructed contigs match the source string.

Good extension:

- Add simple error filtering by k-mer count.

Out of scope:

- Real-genome assembly from raw reads.
- Large FASTQ processing.

### Track C: HMM Sequence Annotation

When students can start:

- Week 7, after Viterbi and forward-backward.

Course method:

- HMM and Viterbi/posterior decoding.

What students can do early:

- Choose hidden states and emissions from a starter sheet.
- Define biological meaning of states.
- Prepare toy sequences.

Baseline deliverable:

- Viterbi decoder for a two- or three-state HMM.
- Apply to CpG-like, exon/intron-like, or protein-region toy data.

Validation check:

- Hand-solvable sequence with known best path.
- Deterministic emission test.

Good extension:

- Posterior decoding comparison.
- Simple supervised parameter estimation from labeled examples.

Out of scope:

- Full gene finder.
- Unconstrained Baum-Welch project.

### Track D: Phylogenetics

When students can start:

- Week 8.

Course method:

- UPGMA, neighbor-joining conceptually, or parsimony on fixed trees.

What students can do early:

- Pick a curated aligned sequence set.
- Compute pairwise distances.
- Inspect metadata.

Baseline deliverable:

- Build a distance-based tree.
- Compare with a parsimony score on one or more candidate trees.

Validation check:

- Four-taxon hand-check.
- Compare major clades to known labels in the dataset.

Good extension:

- Bootstrap resampling at a conceptual/small-data level.

Out of scope:

- Full maximum-likelihood phylogenetics.
- Alignment-free tree inference unless carefully scoped.

### Track E: Omics PCA and FDR

When students can start:

- Week 6 for data orientation.
- Week 10 or 11 for serious method work.

Course method:

- PCA/SVD, multiple testing, FDR.

What students can do early:

- Inspect sample metadata.
- Load a preprocessed count or expression matrix.
- Identify sample groups.
- Write the biological comparison question.

Baseline deliverable:

- PCA on a curated expression matrix.
- Benjamini-Hochberg analysis of provided or computed test statistics.
- Interpret batch/group structure and discoveries.

Validation check:

- Recover a known sample grouping.
- Simulate null p-values and verify BH behavior.
- Reproduce a small instructor-provided expected output.

Good extension:

- Compare results under two normalization/filtering choices.

Out of scope:

- Raw FASTQ to count matrix.
- Full DESeq2/edgeR workflow as a black box.
- Large single-cell dataset.

### Track F: Clustering and Embedding Stability

When students can start:

- Week 10 after PCA begins.

Course method:

- PCA, distance choice, kNN graph, clustering/embedding interpretation.

What students can do early:

- Load curated reduced matrix.
- Examine metadata and labels.

Baseline deliverable:

- Compare clustering or nearest-neighbor structure under two distance/preprocessing choices.
- Explain what changes and what does not.

Validation check:

- Synthetic data with known clusters.
- Stability under small perturbations.

Good extension:

- Compare PCA visualization to a provided UMAP/t-SNE output.

Out of scope:

- Implement UMAP or t-SNE from scratch.
- Claim cell types from visual clusters alone.

### Track G: Gene-Regulation ODEs

When students can start:

- Week 13, so this should not be a default project track unless heavily scaffolded.

Course method:

- Production/degradation ODEs, Hill functions, toggle switch.

What students can do early:

- Read a short model description.
- Identify variables and parameters.

Baseline deliverable:

- Simulate a provided model.
- Map parameter changes to qualitative behavior.

Validation check:

- Reproduce an instructor-provided phase portrait or steady-state case.

Good extension:

- Sensitivity analysis over two parameters.

Out of scope:

- Parameter inference from real time-course data.
- Large network modeling.

Recommendation:

- Keep ODE projects available only as instructor-approved projects, not as a common starter track.

## Starter Sheet Template

Each project starter should be one page.

Fields:

- Title.
- Biological question.
- Dataset and provenance.
- Mathematical object.
- Course method.
- What you can start now.
- What you will learn later.
- Baseline deliverable.
- Validation check.
- Optional extension.
- Out of scope.
- Suggested first three tasks.

## Example Starter Sheet

Title:

- How does k affect toy genome assembly?

Biological question:

- When short reads come from a repeated genome region, how does k-mer length change assembly ambiguity?

Dataset:

- Instructor-provided simulated reads from a 300-base toy genome with one repeat.

Mathematical object:

- de Bruijn graph.

Course method:

- k-mer graph construction and Eulerian path reasoning.

What you can start now:

- Count k-mers.
- Inspect read coverage.
- Build graph nodes and edges.

What you will learn later:

- More systematic validation and failure-mode analysis.

Baseline deliverable:

- Build graphs for k = 5, 9, and 15.
- Report number of nodes, edges, branching nodes, and contigs.
- Explain which k gives the most faithful reconstruction and why.

Validation check:

- Compare contigs to the known simulated genome.

Optional extension:

- Introduce 1% read errors and test simple k-mer count filtering.

Out of scope:

- Real bacterial genome assembly.
- FASTQ quality-score modeling.

Suggested first three tasks:

1. Load the reads and count k-mers for k = 5.
2. Build the de Bruijn graph and print branching nodes.
3. Run on the no-repeat control dataset to confirm the implementation works.

## Practical Advice

Make early project work individual, then make final project work team-based.

Reason:

- Early individual triage prevents students from joining a team with no idea what they prefer.
- Teams form better when students have already seen the menu and tried method cards.

Do not ask for original project ideas too early.

Better prompt:

- "Choose one starter track and adapt one parameter, dataset subset, or validation question."

Avoid projects whose first step is data cleaning.

Better:

- Provide curated data in the course repository.
- Let students document provenance and limitations without spending two weeks parsing file formats.

Require a baseline before extension.

Every project should have:

- A must-do baseline that earns a solid grade if done well.
- Optional extensions that distinguish excellent projects.

Allow controlled pivots.

Checkpoint 1 should permit a team to switch from one track to a closely related fallback if the data or implementation is failing. After Checkpoint 1, pivots should require instructor approval.

## Bottom Line

The project should start early as orientation and scoping, not as full implementation.

Students can start before knowing the whole course if:

- Project choices are track-based.
- Datasets are curated.
- Baseline deliverables are explicit.
- Early tasks focus on inputs, outputs, assumptions, and validation.
- Teams form only after students have completed individual project triage.
