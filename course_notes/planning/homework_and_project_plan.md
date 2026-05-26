# Homework and Project Plan

## Overall Homework Design

There are 10 standard problem sets plus two project workshops. Each standard problem set should have:

- One mathematical derivation or proof.
- One small implementation or computational experiment.
- At least one validation check with a known or independently justified answer.
- One interpretation question.
- One failure-mode reflection, 120-180 words.

Implementation expectations:

- Students may use base Python, R, Julia, or MATLAB if the instructor allows multiple languages.
- For grading consistency, Python starter code is the default.
- Use libraries for arrays, plotting, data loading, and standard numerical linear algebra when allowed.
- Do not use libraries that implement the algorithm being assessed.
- Coding tasks should be compact: usually 30-80 lines for the core method.
- The assignment may provide recurrences, model descriptions, high-level pseudocode, function signatures, and starter tests.
- Students are responsible for the assessed inner loop, validation checks, and explanation.

AI/tool expectations:

- AI coding assistance is allowed when disclosed.
- AI output is treated as unverified output.
- Students must verify code with tests, derivations, hand-checkable examples, or comparisons to known answers.
- Students must be able to explain how their implementation matches the mathematical recurrence, graph, model, or procedure.
- Code that runs but cannot be explained is not sufficient evidence of understanding.

Recommended grading split per problem set:

- 40% mathematical derivation / proof.
- 30% implementation.
- 15% verification tests.
- 10% interpretation.
- 5% AI/tool disclosure when relevant.
- 10% failure-mode reflection.

The exact split can change by week. The key is that students should see the same four forms of evidence repeatedly.

For coding-heavy assignments, use the template in `coding_ai_assessment.md`: mathematical object, implementation, verification, explanation, AI/tool disclosure, and failure-mode reflection.

## PS1: Scoring, Complexity, and k-mers

Timing:

- Out Week 1 Day 3.
- Due Week 2 Day 6.

Learning goals:

- Translate biological sequence assumptions into scores.
- Use asymptotic notation for simple sequence algorithms.
- Demonstrate basic programming readiness.

Possible tasks:

- Given toy substitution and background counts, compute log-odds scores.
- Compare runtime of naive substring search and dictionary-based k-mer counting.
- Implement a k-mer counter.
- Explain one biological case where k-mer counts ignore important information.

Difficulty target:

- Mostly ramp-up. The assignment should identify weak preparation without punishing students too harshly.

## PS2: Global Alignment

Timing:

- Out Week 2 Day 6.
- Due Week 3 Day 9.

Learning goals:

- Derive and prove the Needleman-Wunsch recurrence.
- Implement global alignment with traceback.
- Test a dynamic-programming implementation.

Possible tasks:

- Fill a DP table by hand.
- Prove correctness using an optimal-substructure argument.
- Implement score-table construction and traceback.
- Write tests for empty-string and one-character cases.
- Compare alignments under two scoring systems.

Difficulty target:

- First serious proof and implementation.

## PS3: Local Alignment, Affine Gaps, and Search Heuristics

Timing:

- Out Week 3 Day 9.
- Due Week 4 Day 12.

Learning goals:

- Distinguish local and global alignment.
- Understand state expansion for affine gaps.
- Explain why large database search uses heuristics.

Possible tasks:

- Modify a global-aligner recurrence into Smith-Waterman.
- Derive the three-matrix affine-gap recurrence.
- Implement local alignment, but only derive affine gaps if full implementation is too heavy.
- Build a tiny seed-and-extend search over a list of short sequences.
- Identify a case where seed-and-extend misses an alignment.

Difficulty target:

- Mathematically heavier than PS2. Keep code scope controlled.

## PS4: de Bruijn Graph Assembly

Timing:

- Out Week 4 Day 12.
- Due Week 5 Day 15.

Learning goals:

- Construct a de Bruijn graph from k-mers.
- Relate Eulerian paths to genome reconstruction.
- Diagnose repeats and sequencing errors.

Possible tasks:

- Build graph nodes and edges from reads.
- Find an Eulerian path on a small graph.
- Reconstruct a sequence from a path.
- Compare results for two k values.
- Add one erroneous read and describe its graph signature.

Difficulty target:

- Algorithmically concrete. Good for students who struggled with DP.

## PS5: Markov Chains and CpG Islands

Timing:

- Out Week 5 Day 15.
- Due Week 6 Day 18.

Learning goals:

- Compute sequence likelihoods under iid and Markov models.
- Use log-likelihood ratios for classification.
- Interpret threshold tradeoffs.

Possible tasks:

- Estimate transition matrices from labeled toy sequences.
- Compute log-likelihood under two models.
- Classify windows as CpG-like or background.
- Plot or tabulate sensitivity/specificity at several thresholds.
- Reflect on genome-specific bias.

Difficulty target:

- Probability ramp. Keep notation consistent and provide one fully worked example.

## PS6: Triad HMM / Viterbi Problem Set

Timing:

- Out Week 6 Day 18.
- Due Week 7 Day 21.

Learning goals:

- Implement Viterbi in a group setting.
- Document contribution and verification.
- Explain hidden-state assumptions in a gene-finding toy model.

Possible shared tasks:

- Define a small exon/intron HMM.
- Compute the most likely state path for a short DNA sequence.
- Implement Viterbi in log space.
- Validate with a deterministic or hand-solvable case.
- Compare two parameter settings.

Individual addendum:

- Each student explains one proof step, recurrence term, or code design choice in their own words.

Contribution log:

- One paragraph signed by all members.

Difficulty target:

- Collaboration warm-up. Technical challenge is real, but data and state space should be small.

## PS7: Forward-Backward and Posterior Decoding

Timing:

- Out Week 7 Day 21.
- Due Week 8 Day 22.

Learning goals:

- Distinguish best path from total sequence probability.
- Implement or trace forward-backward.
- Explain numerical underflow and scaling/log-space fixes.

Possible tasks:

- Compute forward probabilities on a short sequence.
- Compute posterior state probabilities.
- Give an example where Viterbi and posterior decoding disagree.
- Implement log-sum-exp or a scaling approach.

Difficulty target:

- Conceptually hard. Keep implementation small and provide expected outputs for tests.

## PS8: Phylogenetics

Timing:

- Out Week 9 Day 27.
- Due Week 10 Day 30.

Learning goals:

- Build and critique distance-based trees.
- Use parsimony on fixed trees.
- Compare model assumptions.

Possible tasks:

- Compute pairwise distances from a short multiple alignment.
- Run UPGMA by hand for four taxa.
- Apply Fitch's algorithm to score a fixed tree.
- Show how non-clock-like evolution breaks UPGMA.
- Compare distance and parsimony conclusions.

Difficulty target:

- Moderate. Good post-midterm restart.

## PS9: PCA and SVD for Omics Data

Timing:

- Out Week 10 Day 30.
- Due Week 11 Day 33.

Learning goals:

- Center an expression matrix.
- Compute and interpret PCA/SVD.
- Identify batch effects or dominant unwanted variation.

Possible tasks:

- Work through a small matrix by hand or with starter code.
- Interpret scores and loadings.
- Compare PCA before and after centering/scaling.
- Explain why a high-variance component need not be biologically meaningful.

Difficulty target:

- Mathematically important. Students with MATH 2270 should be expected to reason carefully here.

## PS10: Multiple Testing and FDR

Timing:

- Out Week 11 Day 33.
- Due Week 12 Day 36.

Learning goals:

- Explain why many tests create many false positives.
- Apply Benjamini-Hochberg.
- Interpret FDR in omics terms.

Possible tasks:

- Simulate p-values under all-null and mixed-null settings.
- Apply BH to a table of p-values.
- Compare Bonferroni and BH.
- Interpret effect size versus p-value.
- Explain how dependence among genes complicates interpretation.

Difficulty target:

- Essential and practical. Keep the proof demands light unless a stats prerequisite is added.

## Optional Exercise Bank: Clustering and Embedding Stability

Timing:

- Use as in-class work, optional practice, or as part of a project starter.
- Do not assign as a full eleventh problem set under the recommended no-final-exam assessment model.

Learning goals:

- Understand distance choice and neighborhood graphs.
- Interpret PCA/t-SNE/UMAP-style plots cautiously.
- Evaluate clustering stability.

Possible tasks:

- Build a kNN graph from low-dimensional data.
- Compare clusters under two distance metrics.
- Run PCA and a provided embedding function, if allowed.
- Perturb data or parameters and document what changes.
- Explain why visual separation does not prove biological identity.

Difficulty target:

- Interpretation-heavy. Do not require students to implement UMAP or t-SNE from scratch.

## Optional Exercise Bank: Toggle Switch and Course Synthesis

Timing:

- Use as in-class work, optional practice, or as a heavily scaffolded final-project track.
- Do not assign as a full twelfth problem set under the recommended no-final-exam assessment model.

Learning goals:

- Formulate a simple gene-regulation ODE model.
- Analyze steady states and parameter effects.
- Synthesize method choice across the course.

Possible tasks:

- Simulate a production-degradation model.
- Sketch nullclines for a two-gene toggle switch.
- Identify parameter regimes with one versus multiple steady states.
- Given short project-style biological questions, choose a course method and justify the choice.

Difficulty target:

- Conceptual synthesis. Avoid a long numerical-analysis assignment at the end of term.

## Final Project

Recommended form:

- Teams of three.
- One course method.
- One real or realistic dataset.
- One validation check.
- One reproducible repository.
- One README aimed at a math/CS/biology peer who has not taken the course.
- One revision memo.

Important design point:

- The project should not begin as an open-ended research proposal. Students do not know enough early in the semester to scope one well.
- Instead, the project should begin as a track-based starter project. Each track gives a biological question, curated small dataset, course method, baseline deliverable, validation check, and optional extension.
- Early project work should focus on identifying inputs, outputs, assumptions, and validation checks. Full implementation begins only after the relevant course method has been taught.

Checkpoint 1:

- Due Week 11 Day 33.
- One-page proposal.
- Repository skeleton.
- Data source identified.
- Planned validation check.

Recommended pre-checkpoint on-ramp:

- Week 1: project gallery, no commitment.
- Week 3 or 4: method card for alignment or assembly.
- Week 6 or 7: project starter menu released.
- Week 8: individual project triage memo ranking two or three tracks.
- Spring break: teams form around starter tracks.
- Week 10: project contract converts the chosen track into a concrete plan.

Checkpoint 2:

- Due Week 14 Day 41.
- Preliminary result.
- README draft.
- Evidence that code runs on a small input.

Final submission:

- Due finals week.
- Repository with code, README, data provenance, requirements file, validation check, contribution paragraph, and revision memo.

## Curated Project Starters

Good starter projects:

- Implement Needleman-Wunsch or Smith-Waterman and compare scoring schemes on homologous genes.
- Build a toy de Bruijn assembler and test how k and errors affect contigs.
- Train a two-state Markov or HMM classifier for CpG-like regions.
- Implement Viterbi for a simplified gene-finding or transmembrane-domain model.
- Build UPGMA trees from short viral or mitochondrial sequence alignments.
- Run PCA and FDR analysis on a small public expression dataset.
- Compare clustering stability on a reduced single-cell or synthetic expression matrix.
- Analyze a toggle-switch model and connect parameter changes to qualitative behavior.

Better framing:

- Alignment, assembly, HMM, phylogenetics, and omics projects should be the main starter tracks.
- ODE/toggle-switch projects should be instructor-approved or heavily scaffolded because the ODE module occurs late.

Projects to discourage unless tightly scoped:

- Full RNA-seq pipeline from FASTQ.
- Full single-cell analysis from raw reads.
- Deep-learning model training.
- Large-genome assembly.
- Novel web application.
- Any project requiring large compute or unclear data permissions.

## Workload Controls

Use small data:

- Toy examples for algorithm correctness.
- Curated mini datasets for omics.
- Avoid files larger than 50 MB for required work.

Use starter code for:

- File parsing.
- Plotting.
- Test harnesses.
- Data downloads, if any.

Do not scaffold:

- The recurrence.
- The core loop of the algorithm.
- The mathematical explanation.

Use autograded tests only for mechanical checks. Human grading is needed for proofs, interpretation, and failure-mode reflections.
