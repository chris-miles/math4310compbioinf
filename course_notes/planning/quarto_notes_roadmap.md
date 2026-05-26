# Quarto Notes Roadmap

This is a proposed future structure for student-facing Quarto course notes. The files do not need to be created until the planning documents stabilize.

## Proposed Folder Structure

```text
course_notes/
  index.qmd
  notes/
    00_setup.qmd
    01_sequences_and_scores.qmd
    02_global_alignment.qmd
    03_local_alignment_and_blast.qmd
    04_de_bruijn_assembly.qmd
    05_markov_chains.qmd
    06_hidden_markov_models.qmd
    07_forward_backward.qmd
    08_phylogenetics.qmd
    09_omics_matrices_pca.qmd
    10_multiple_testing_fdr.qmd
    11_clustering_embeddings.qmd
    12_gene_regulation_odes.qmd
    13_synthesis.qmd
  assignments/
    ps01_scoring_kmers.qmd
    ps02_global_alignment.qmd
    ps03_local_affine_seed.qmd
    ps04_assembly_graphs.qmd
    ps05_markov_cpg.qmd
    ps06_hmm_viterbi_triad.qmd
    ps07_forward_backward.qmd
    ps08_phylogenetics.qmd
    ps09_pca_svd.qmd
    ps10_fdr.qmd
    workshop01_project_proposal.qmd
    workshop02_project_preliminary_results.qmd
    optional_clustering_embeddings.qmd
    optional_odes_synthesis.qmd
  data/
    toy/
  code/
    starters/
    solutions_private/
```

The `solutions_private/` folder should not be committed to a public course site.

## Standard Chapter Template

Each notes chapter should follow this pattern:

1. Biological question.
2. Mathematical object.
3. Model assumptions.
4. Algorithm or derivation.
5. Small worked example.
6. Short implementation.
7. Validation check.
8. Failure modes.
9. Exercises.

This keeps chapters from becoming either pure biology background or pure math detached from data.

For chapters with code, the notes should explicitly connect the loop or function to the recurrence, graph, model, matrix operation, or statistical procedure it implements. Code examples should model the course policy: tools are allowed, but output must be verified.

## Chapter Sketches

### 00 Setup

Purpose:

- Establish tools, reproducibility expectations, and readiness checks.
- Point students to optional programming support without reteaching the applied-genomics workflow material covered in BIOL 3150.

Include:

- How to run starter notebooks/scripts.
- How to submit code and written work.
- How to disclose AI/tool use.
- How to verify AI-assisted code with hand-checkable tests.
- Basic test-writing examples.
- A short diagnostic: functions, loops over strings, dictionaries/maps, and simple 2D tables.

Avoid:

- A full Python, pandas, Colab, Unix, BLAST, IGV, FASTQ, variant-calling, RNA-seq, or GWAS bootcamp.

### 01 Sequences and Scores

Purpose:

- Introduce biological sequences as finite strings.
- Explain scoring as model choice.

Core math:

- Alphabets.
- k-mers.
- Log-odds ratios.
- Expected score.

Code:

- k-mer counter.
- Toy log-odds matrix.

### 02 Global Alignment

Purpose:

- Teach dynamic programming through Needleman-Wunsch.

Core math:

- Recurrence.
- Boundary conditions.
- Correctness proof.
- Complexity.

Code:

- Score table.
- Traceback.

### 03 Local Alignment and BLAST

Purpose:

- Extend DP to local alignment and realistic search.

Core math:

- Smith-Waterman recurrence.
- Affine gap state expansion.
- Seed-and-extend heuristic.

Code:

- Local alignment.
- Tiny seed index.

### 04 de Bruijn Assembly

Purpose:

- Show genome assembly as graph theory.

Core math:

- k-mer graph construction.
- Eulerian paths.
- Repeats and graph ambiguity.

Code:

- Build graph.
- Find or verify path.

### 05 Markov Chains

Purpose:

- Introduce probabilistic sequence models.

Core math:

- Transition matrices.
- Sequence likelihood.
- Log-likelihood ratio.

Code:

- CpG classifier.

### 06 Hidden Markov Models

Purpose:

- Teach hidden-state sequence models and Viterbi.

Core math:

- Joint probability.
- Conditional independence.
- Max-product recurrence.

Code:

- Viterbi in log space.

### 07 Forward-Backward

Purpose:

- Distinguish best path from total probability and posterior state probability.

Core math:

- Sum-product recurrence.
- Forward and backward variables.
- Posterior decoding.
- Log-sum-exp.

Code:

- Forward algorithm.
- Posterior table.

### 08 Phylogenetics

Purpose:

- Infer trees from sequence data and compare assumptions.

Core math:

- Distance matrices.
- UPGMA.
- Parsimony.
- Fitch algorithm.

Code:

- UPGMA on small matrix.
- Parsimony score on fixed tree.

### 09 Omics Matrices and PCA

Purpose:

- Use linear algebra to understand expression data.

Core math:

- Centering.
- Covariance.
- SVD.
- PCA scores and loadings.
- Low-rank approximation.

Code:

- PCA from SVD on a small matrix.

### 10 Multiple Testing and FDR

Purpose:

- Explain high-throughput statistical testing.

Core math:

- Null p-values.
- False discoveries.
- Benjamini-Hochberg.

Code:

- Simulate p-values.
- Apply BH.

### 11 Clustering and Embeddings

Purpose:

- Connect PCA, distances, neighborhood graphs, and modern visualizations.

Core math:

- Distance choice.
- kNN graph.
- Clustering instability.
- Embedding interpretation.

Code:

- kNN graph.
- Parameter sensitivity experiment.

### 12 Gene Regulation ODEs

Purpose:

- Contrast mechanistic modeling with sequence and matrix methods.

Core math:

- Production-degradation ODEs.
- Hill functions.
- Steady states.
- Nullclines.
- Bistability.

Code:

- Simulate a toggle switch.

### 13 Synthesis

Purpose:

- Help students choose methods and articulate assumptions.

Include:

- Method-selection guide.
- Cross-module comparison table.
- Project validation checklist.
- Course-synthesis and project-review problems.

## Style Rules for Student Notes

- Keep each chapter short enough to be read before class.
- Put long derivations behind collapsible sections or appendices if Quarto supports it cleanly.
- Use the same notation across HMM chapters.
- Use toy data first, then a realistic mini dataset.
- Put "Failure modes" in every chapter.
- Include at least one manually checkable test for every implementation.
- Assignment prompts should ask for AI/tool disclosure, validation, and a short explanation of why the code matches the math.
