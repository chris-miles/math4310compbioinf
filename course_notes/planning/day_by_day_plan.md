# Day-by-Day Plan

This is a 14-instructional-week plan with 3 meetings per week. Spring break is not counted as an instructional week. If the registrar calendar labels spring break as Week 9, then instructional Week 9 below happens after break.

Each class period is approximately 50 minutes.

Assessment note: the current syllabus uses two shorter midterms and no cumulative final exam. This detailed daily plan still contains some older single-midterm phrasing in the middle weeks; when converting to final Quarto notes, use the current syllabus and `assessment_model.md` for exam placement, and preserve the daily topic flow as the content scaffold.

## Assignment Rhythm

Default rhythm:

- Problem sets go out after the third meeting of a week.
- Problem sets are due before the third meeting of the following instructional week.
- There are 10 standard problem sets plus two project workshops/checkpoints.
- The short-midterm weeks should keep new technical homework light.
- Project teams form by the end of spring break.

## Week 1: Sequences, Algorithms, and Scoring

### Day 1: What Counts as Bioinformatics Theory?

Core question: What makes a biological-data problem algorithmic or mathematical?

In class:

- Contrast three tasks: run BLAST, explain BLAST, design a toy alignment method.
- Introduce sequences as strings over finite alphabets.
- Define reads, reference sequences, genes, proteins, mutations, homologs.
- Hand activity: find "similar" DNA strings under several informal scoring rules.
- Programming diagnostic preview: strings, loops, dictionaries, 2D tables.

Notes to write:

- Biology vocabulary primer.
- Why "from scratch" means small, theory-revealing implementations.

### Day 2: Complexity, Recurrences, and Tables

Core question: Why is exhaustive alignment impossible?

In class:

- Big-O review using string comparison examples.
- Difference between brute-force search and dynamic programming.
- Recurrences as computational definitions.
- Optimal substructure as the reason a table can replace a tree of possibilities.
- Mini activity: count paths in a grid and connect to alignment paths.

Notes to write:

- Complexity refresher.
- Table-filling as a mathematical proof strategy.

### Day 3: Log-Odds Scores and Biological Assumptions

Core question: What does an alignment score mean?

In class:

- Match/mismatch/gap scores as model assumptions.
- Toy substitution counts.
- Log-odds score = evidence for relatedness versus background.
- Expected score and why arbitrary positive scores are dangerous.
- Introduce gap penalties.

Assignment:

- PS1 out: scoring, complexity, k-mer counting diagnostic.

Reading suggestions:

- Durbin et al., pairwise-alignment scoring sections.
- Compeau and Pevzner, sequence comparison introduction.

## Week 2: Global Alignment

### Day 4: Needleman-Wunsch Recurrence

Core question: How do we align two whole sequences optimally?

In class:

- Define global alignment objective.
- Derive the three-way recurrence: match/mismatch, deletion, insertion.
- Boundary conditions.
- Fill a small matrix by hand.

Notes to write:

- Global alignment recurrence.
- Boundary conditions as biological assumptions.

### Day 5: Correctness and Traceback

Core question: Why does the recurrence actually find an optimum?

In class:

- Optimal substructure proof template.
- Traceback graph and tie handling.
- Complexity O(nm) time and O(nm) memory.
- Hirschberg mention only: memory can be improved, but not central.

Notes to write:

- Proof template for dynamic programming.
- Traceback and multiple optimal alignments.

### Day 6: Implementation Clinic

Core question: How do table definitions become reliable code?

In class:

- Starter implementation skeleton.
- Unit tests on tiny examples.
- Common bugs: off-by-one indexing, wrong boundary signs, traceback ties.
- Discuss how scoring choices change alignments.

Assignment:

- PS1 due.
- PS2 out: Needleman-Wunsch proof and implementation.

## Week 3: Local Alignment, Gaps, and Database Search

### Day 7: Smith-Waterman Local Alignment

Core question: How do we find the best matching region inside longer sequences?

In class:

- Local versus global biological use cases.
- Add zero to the recurrence.
- Traceback stopping rule.
- Interpret local score as "best evidence region."

Notes to write:

- Local alignment recurrence.
- Why zero changes the problem.

### Day 8: Affine Gap Penalties

Core question: Why should one long gap cost less than many short gaps?

In class:

- Biological motivation: insertions and deletions often occur in runs.
- Gap-open and gap-extend penalties.
- Three-matrix recurrence.
- Complexity stays O(nm), but state space changes.

Notes to write:

- Affine gap model.
- State expansion in dynamic programming.

### Day 9: Seed-and-Extend and BLAST Intuition

Core question: How do alignment tools search large databases quickly?

In class:

- Exact DP does not scale to every query against every database sequence.
- Seeds, word hits, extension.
- Sensitivity versus speed.
- E-values and statistical significance at a conceptual level.
- Failure modes: low-complexity regions, repeats, composition bias.

Assignment:

- PS2 due.
- PS3 out: local alignment, affine gaps, and a tiny seed-and-extend experiment.

## Week 4: Genome Assembly as Graph Theory

### Day 10: Reads, k-mers, and Assembly

Core question: How can short reads reconstruct a long genome?

In class:

- Reads, coverage, sequencing errors.
- Overlap-layout-consensus versus de Bruijn graph framing.
- k-mer spectrum.
- Hand activity: reconstruct a short string from k-mers.

Notes to write:

- Reads and coverage.
- k-mers as lossy summaries.

### Day 11: de Bruijn Graphs and Eulerian Paths

Core question: Why does assembly become an Eulerian path problem?

In class:

- Nodes as (k-1)-mers, edges as k-mers.
- Eulerian path/cycle conditions.
- Hierholzer algorithm conceptually.
- Complexity in terms of number of k-mers.

Notes to write:

- de Bruijn graph construction.
- Eulerian path conditions.

### Day 12: Repeats, Errors, and Graph Cleaning

Core question: Where does the elegant graph model fail?

In class:

- Repeats collapse paths.
- Sequencing errors create tips and bubbles.
- Effect of k.
- Why real assemblers need graph simplification.
- Validation: contiguity, correctness, coverage.

Assignment:

- PS3 due.
- PS4 out: de Bruijn graph construction and toy assembly.

Reading suggestions:

- Compeau and Pevzner, genome assembly / graph algorithms chapter.

## Week 5: Probability and Markov Chains

### Day 13: Probability Refresher for Sequences

Core question: What is the probability of a biological sequence under a model?

In class:

- Random variables, conditional probability, likelihood.
- Independence assumptions.
- Background nucleotide model.
- Log-likelihood to avoid underflow.

Notes to write:

- Probability notation used in the course.
- Likelihood versus probability of parameters.

### Day 14: Markov Chains

Core question: How does local dependence change sequence probability?

In class:

- Transition matrix.
- First-order Markov assumption.
- Stationary distribution concept.
- Simulate short sequences.
- Compare iid and Markov sequence likelihoods.

Notes to write:

- Markov chain definition.
- Sequence likelihood under a Markov chain.

### Day 15: CpG Islands as Likelihood Ratios

Core question: How can a probabilistic model classify genomic regions?

In class:

- CpG biological motivation.
- Two Markov models: island and background.
- Log-likelihood ratio classifier.
- Thresholds, sensitivity, specificity.
- Failure modes: genome-specific composition, training bias.

Assignment:

- PS4 due.
- PS5 out: Markov-chain CpG classifier.

## Week 6: Hidden Markov Models and Viterbi

### Day 16: HMM Definitions

Core question: What changes when the state is hidden?

In class:

- Hidden states and observed emissions.
- Transition and emission probabilities.
- Conditional independence assumptions.
- Toy biological examples: CpG status, exon/intron, transmembrane region.
- Joint probability of state path and observations.

Notes to write:

- HMM definition and notation.
- Biological interpretation of hidden states.

### Day 17: Viterbi Algorithm

Core question: What is the most likely hidden path?

In class:

- Derive the max-product recurrence.
- Convert to log-space max-sum recurrence.
- Backpointers.
- Complexity O(nK^2).
- Compare to ordinary dynamic programming for alignment.

Notes to write:

- Viterbi derivation.
- Log-space implementation.

### Day 18: Viterbi Studio

Core question: How do we debug a probabilistic dynamic program?

In class:

- Work through a tiny HMM by hand.
- Implementation planning in triads.
- Tests: one-state HMM, deterministic emissions, short sequence with known best path.
- Discuss contribution logs and individual addenda.

Assignment:

- PS5 due.
- PS6 out: triad Viterbi / toy gene-finding problem.

## Week 7: Forward-Backward and HMM Training

### Day 19: Forward Algorithm

Core question: What is the total probability of the observed sequence?

In class:

- Sum-product versus max-product.
- Derive forward recurrence.
- Likelihood of observations.
- Scaling and log-sum-exp.

Notes to write:

- Forward recurrence.
- Numerical stability.

### Day 20: Backward Algorithm and Posterior Decoding

Core question: How likely is each state at each position?

In class:

- Backward recurrence.
- Posterior state probabilities.
- Viterbi path versus marginal posterior decoding.
- Example where they disagree.

Notes to write:

- Backward recurrence.
- Posterior decoding.

### Day 21: Baum-Welch Conceptually

Core question: Can the model learn from unlabeled sequences?

In class:

- Supervised versus unsupervised HMM estimation.
- EM idea: expected counts then parameter update.
- Baum-Welch conceptually, not full implementation.
- Overfitting and model choice.

Assignment:

- PS6 due.
- PS7 out: forward-backward and posterior decoding.

## Week 8: Integration and Midterm

### Day 22: Review Across Algorithms and Probability

Core question: What is the common structure behind alignment, assembly, and HMMs?

In class:

- Compare DP tables, graph paths, likelihood recurrences.
- Worked mixed problem.
- Midterm review with proof and implementation-adjacent pseudocode.

Assignment:

- PS7 due.

### Day 23: Midterm Exam

Coverage:

- Scoring and log-odds.
- Global/local/affine alignment.
- de Bruijn graphs.
- Markov chains.
- HMMs, Viterbi, forward-backward.

### Day 24: Midterm Debrief and Project Launch

Core question: How do course methods become project methods?

In class:

- Debrief common midterm issues.
- Present curated final-project starter problems.
- Explain repository skeleton and validation-check requirements.
- Start informal team matching.

Spring-break task:

- Teams form by end of break.
- Students skim project starter list and pick two ranked options.

## Week 9: Phylogenetics

### Day 25: Distances and Molecular Clocks

Core question: How do sequence differences become evolutionary distances?

In class:

- Pairwise distances from aligned sequences.
- Correction for multiple hits conceptually.
- Ultrametric assumption.
- Molecular clock as a model, not a fact.

Notes to write:

- Distance matrix construction.
- Tree assumptions.

### Day 26: UPGMA and Neighbor-Joining Intuition

Core question: How do clustering algorithms build trees?

In class:

- UPGMA algorithm step by step.
- Why UPGMA assumes clock-like evolution.
- Neighbor-joining as a less restrictive extension, conceptually.
- Complexity and tie handling.

Notes to write:

- UPGMA algorithm.
- When distance methods fail.

### Day 27: Parsimony and Fitch's Algorithm

Core question: Can we infer a tree by minimizing character changes?

In class:

- Character matrix.
- Parsimony objective.
- Fitch algorithm on a fixed tree.
- Why tree search is hard.
- Contrast distance and character methods.

Assignment:

- PS8 out: UPGMA and parsimony.

Reading suggestions:

- Durbin et al., phylogenetic tree chapter.

## Week 10: Omics Matrices and PCA

### Day 28: Expression Data as a Matrix

Core question: What is an RNA-seq or omics matrix mathematically?

In class:

- Samples by features.
- Counts, library size, log transforms.
- Centering and scaling.
- Batch effects as unwanted structure.
- Why "more features than samples" matters.

Notes to write:

- Expression-matrix notation.
- Normalization as modeling.

### Day 29: PCA and SVD

Core question: What does PCA optimize?

In class:

- Variance maximization view.
- SVD view.
- Scores, loadings, singular values.
- Low-rank approximation.
- Interpretability warnings.

Notes to write:

- PCA derivation from SVD.
- Geometry of projections.

### Day 30: PCA Lab and Failure Modes

Core question: What can a PCA plot tell us, and what can it hide?

In class:

- Work through a small expression matrix.
- Scree plot interpretation.
- Batch effect example.
- Choosing whether to center, scale, filter genes.
- Link to project data exploration.

Assignment:

- PS8 due.
- PS9 out: PCA/SVD on a small omics matrix.

## Week 11: High-Throughput Testing and FDR

### Day 31: Many Tests, Many Mistakes

Core question: What changes when we test thousands of genes?

In class:

- Null p-value distribution.
- Type I errors across many tests.
- Family-wise error rate concept.
- Why naive p < 0.05 fails for omics.

Notes to write:

- Multiple-testing setup.
- Simulated p-values.

### Day 32: Benjamini-Hochberg FDR

Core question: What does FDR control?

In class:

- Define false discovery proportion and FDR.
- Step-up BH procedure.
- Graphical interpretation with sorted p-values.
- Independence / positive-dependence assumptions at a high level.

Notes to write:

- BH algorithm.
- FDR interpretation.

### Day 33: Differential Expression Mini-Case

Core question: How do normalization, variance, and testing interact?

In class:

- Simulated gene-expression experiment.
- Effect size versus p-value.
- Volcano plots conceptually.
- Validation and reproducibility expectations for projects.

Assignment:

- PS9 due.
- PS10 out: FDR simulation and differential-expression reasoning.
- Project Checkpoint 1 due: proposal plus repository skeleton.

Reading suggestions:

- Holmes and Huber, testing and high-throughput count-data chapters.

## Week 12: Clustering and Modern Dimensionality Reduction

### Day 34: Distances, Clustering, and kNN Graphs

Core question: What does it mean for samples or cells to be near each other?

In class:

- Distance choices: Euclidean, correlation, cosine.
- k-means and hierarchical clustering, briefly.
- k-nearest-neighbor graph construction.
- Curse of dimensionality.

Notes to write:

- Distance choice as modeling choice.
- kNN graph basics.

### Day 35: t-SNE and UMAP Concepts

Core question: Why are modern embeddings useful but easy to overinterpret?

In class:

- Embeddings as neighborhood-preserving summaries.
- t-SNE and UMAP at the level of objectives and tradeoffs.
- Hyperparameters and instability.
- Why visual clusters are not automatically biological cell types.

Notes to write:

- Conceptual comparison of PCA, t-SNE, UMAP.
- Embedding failure modes.

### Day 36: Single-Cell Case Study

Core question: How do the course ideas appear in a modern single-cell workflow?

In class:

- Matrix preprocessing.
- PCA before neighborhood graph.
- Clustering and annotation as separate steps.
- Batch effects and confounding.
- Ethical/reproducibility note: cell labels are claims, not raw data.

Assignment:

- PS10 due.
- Optional clustering / embedding stability exercises available for project teams that need them.

## Week 13: Mechanistic Models of Gene Regulation

### Day 37: ODEs for Gene Regulation

Core question: What is a mechanistic model trying to explain?

In class:

- State variables and parameters.
- Production and degradation.
- Steady states.
- Numerical simulation as evidence, not proof.

Notes to write:

- ODE modeling vocabulary.
- Production-degradation model.

### Day 38: Hill Functions and Bistability

Core question: How does cooperativity create switch-like behavior?

In class:

- Hill activation/repression functions.
- Nullclines.
- Stability by slope / phase-line intuition.
- Bistability as memory.

Notes to write:

- Hill functions.
- Bistability examples.

### Day 39: Toggle Switch and Model Failure Modes

Core question: When is a small mechanistic model useful?

In class:

- Two-gene toggle switch equations.
- Phase-plane sketch.
- Parameter sensitivity.
- Contrast mechanistic ODEs with HMMs and matrix methods.

Assignment:

- Optional toggle-switch and synthesis exercises available.

Reading suggestions:

- Alon, gene-circuit and network-motif chapters.

## Week 14: Synthesis, Projects, and Modern Extensions

### Day 40: Method-Selection Synthesis

Core question: Given a biological question, what mathematical object should we build?

In class:

- Alignment: pairs of strings.
- Assembly: graph from k-mers.
- HMM: hidden state sequence.
- Phylogeny: tree.
- Omics: matrix and tests.
- ODE: state variables and mechanisms.
- Decision-tree activity with project-style prompts.

Notes to write:

- Method-selection guide.

### Day 41: Project Clinic and Preliminary Results

Core question: What evidence makes a project result credible?

In class:

- Teams show current repository state.
- Validation check audit.
- README peer feedback.
- Failure-mode check: what assumption is most fragile?

Assignment:

- Project Checkpoint 2 due: preliminary results.

### Day 42: Final Review and Modern Frontiers

Core question: How do modern bioinformatics methods extend the primitives from this course?

In class:

- Long-read assembly as graph plus error model.
- Pangenome graphs as graph generalization.
- Profile HMMs as richer sequence-family models.
- Deep sequence models as learned scoring / representation functions.
- Final project and course-synthesis review.

Assignment:

- Final project repository due in Finals Week.
