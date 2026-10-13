# Optional Project On-Ramp

Status: optional project ideas. Projects remain undecided; this public document establishes no deliverables, grade weight, team requirement, or timetable.

## The Problem

If a substantial final project is adopted, it is only fair if students can begin before the last third of the course. But students do not know HMMs, PCA, FDR, phylogenetics, or ODEs at the start. If the project begins only after those topics are taught, the project becomes rushed. If it begins too early, students are choosing methods they do not yet understand.

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

## Timing if projects are adopted

Orientation can begin before a method is taught, but implementation must follow its lesson. No project checkpoints, presentation meetings, or finals-week submission dates are scheduled. April 26 remains Flex.

## Possible Project Tracks

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

- Week 8, after Viterbi and forward-backward.

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

- After Lessons 13–14 (March 3 and March 15).

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
- Lessons 15–17 for PCA and single-feature comparisons; Lesson 18 on March 29 before an FDR analysis.

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

- After Lesson 19 (April 5) for clustering.

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

- Week 14 for ODEs and switches; Week 15 for noisy measurements. This late material requires substantial scaffolding if used in a project.

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

If a team project is adopted, individual exploration before team formation can help students choose a track. Individual and team formats remain options.

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

If checkpoints are adopted, an early scope review can permit a switch to a closely related fallback when data or implementation fails. Decide any later scope-change policy with the project format.

## Design criteria if adopted

The project should start early as orientation and scoping, not as full implementation.

Students can start before knowing the whole course if:

- Project choices are track-based.
- Datasets are curated.
- Baseline deliverables are explicit.
- Early tasks focus on inputs, outputs, assumptions, and validation.
- If teams are used, students have a chance to explore preferences before forming them.
