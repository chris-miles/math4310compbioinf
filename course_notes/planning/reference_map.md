# Reference Map

This file maps course modules to books and open resources. Links were checked during planning on 2026-05-23.

## Reuse and Attribution Policy

The course notes and assignments should be original public-facing materials, even when they are built from standard textbook and open-course precedents.

- Use books and downloaded references internally to check math, notation, topic order, examples, and difficulty.
- Borrow problem types freely, but rewrite prompts in the course's own voice and use new toy data unless the source license clearly permits direct reuse.
- Cite sources when a problem or exposition is inspired by a specific exercise, section, dataset, or open course resource.
- For open-license sources, follow the stated attribution, noncommercial, and share-alike requirements.
- For copyrighted books, citation does not by itself permit copying substantial prose, figures, tables, datasets, or distinctive problem statements into the public Quarto book.
- Keep private reference PDFs and solution materials out of any public course site or public repository.

## Primary Course Spine

### Compeau and Pevzner, Bioinformatics Algorithms

Use for:

- Student-facing biological motivation.
- Algorithmic exercises.
- Genome assembly.
- Dynamic programming for sequence comparison.
- Rosalind-style coding tasks.

Best fit:

- Weeks 1-4.

Useful chapters / sections:

- Chapter 3: genome assembly and graph algorithms.
- Chapter 5: comparing biological sequences and dynamic programming.
- Earlier algorithmic warmups can supply optional practice for weaker programmers.

Notes:

- The author site says the textbook website contains several free chapters and lecture materials.
- Use this as the closest thing to a student spine for the first month.
- Do not copy proprietary problem text directly unless license/permission is clear. Adapt problem types and toy datasets.

Links:

- https://compeau.cbd.cmu.edu/online-education-projects/bioinformatics-algorithms-an-active-learning-approach/
- https://www.bioinformaticsalgorithms.org/

### Durbin, Eddy, Krogh, and Mitchison, Biological Sequence Analysis

Use for:

- Instructor reference.
- Clean probabilistic framing.
- HMM derivations.
- Pairwise alignment.
- Phylogenetics.

Best fit:

- Weeks 2-7 and Week 9.

Useful chapters:

- Chapter 2: pairwise alignment.
- Chapter 3: Markov chains and hidden Markov models.
- Chapter 7: building phylogenetic trees.
- Chapter 11: probability background.

Notes:

- Excellent but dense for many undergraduates.
- Assign short excerpts or instructor-written notes rather than expecting students to read long sections unaided.

Link:

- https://www.cambridge.org/core/books/biological-sequence-analysis/

### Holmes and Huber, Modern Statistics for Modern Biology

Use for:

- Omics statistics.
- Multivariate analysis.
- Multiple testing.
- Count data.
- Reproducible examples.

Best fit:

- Weeks 10-12.

Useful chapters:

- Chapter 7: testing.
- Chapter 8: multivariate analysis.
- Chapter 9: high-throughput count data and generalized linear models.
- Chapter 10: multivariate methods for heterogeneous data.
- Chapter 11: networks and trees, selectively.

Notes:

- The online book is actively maintained and includes code/data.
- It is R/Bioconductor-oriented, so instructor notes should translate the core ideas into the language used for the course.

Link:

- https://www.huber.embl.de/msmb/

### Alon, An Introduction to Systems Biology

Use for:

- Gene-regulation ODEs.
- Hill functions.
- Network motifs.
- Toggle-switch intuition.

Best fit:

- Week 13 only.

Notes:

- Use selectively. This is a systems-biology book, not the main bioinformatics text.
- The second edition has many exercises and classroom-tested examples.

Links:

- https://www.routledge.com/An-Introduction-to-Systems-Biology-Design-Principles-of-Biological-Circuits/Alon/p/book/9781439837177
- https://www.weizmann.ac.il/mcb/UriAlon/introduction-systems-biology-design-principles-biological-circuits

## Peer Course Precedents

### UC Davis ECS 124: Theory and Practice of Bioinformatics

Why it matters:

- Strong precedent for an upper-division theory/practice course.
- Covers sequence analysis, DP, scoring matrices, BLAST, assembly, systems biology, clustering/classification, networks, and phylogenetics.
- Its course-overlap language explicitly distinguishes mathematical/algorithmic treatment from tool-use courses.

Use for:

- Justifying the topic mix.
- Explaining distinction from BIOL 3150 and MATH 4100.

Link:

- https://cs.ucdavis.edu/schedules-classes/ecs-124-theory-practice-bioinformatics

### Duke COMPSCI 561 / CBB 561: Computational Sequence Biology

Why it matters:

- Good precedent for emphasizing HMMs and probabilistic sequence methods.
- Includes assembly, homology detection, gene/promoter finding, motifs, comparative genomics, phylogenetics, RNA structure, and post-transcriptional regulation.

Use for:

- Justifying HMMs and probabilistic approaches as central.

Link:

- https://cs.duke.edu/courses/computational-sequence-biology

### University of Minnesota Duluth MATH 5233: Mathematical Foundations of Bioinformatics

Why it matters:

- Math-prefixed course precedent.
- Includes Needleman-Wunsch, Smith-Waterman, BLAST, Clustal, PAM/BLOSUM, DNA-sequence statistics, CpG islands, and phylogenetic methods.

Use for:

- Departmental justification that this material belongs in mathematics.

Link:

- https://www.d.umn.edu/~mhampton/m5233s11.html

### UNC COMP 555: Bioalgorithms

Why it matters:

- Classic bioalgorithms course.
- Broad algorithmic topics include motifs, rearrangements, alignments, gene prediction, graph algorithms, sequencing, HMMs, and randomized algorithms.

Use for:

- Optional examples and comparison when considering whether to add motifs or rearrangements.

Link:

- https://csbio.unc.edu/mcmillan/

## Competency and Curriculum References

### ISCB Competency Framework v3

Use for:

- External validation of data-science and bioinformatics competencies.
- Especially relevant to method choice, data complexity, and responsible data management.

Link:

- https://doi.org/10.1093/bioadv/vbae166
- https://competency.ebi.ac.uk/framework/iscb/3.0/

### NIBLSE Bioinformatics Core Competencies

Use for:

- Undergraduate justification.
- Supports computational concepts, algorithms, scripting, tool use, and genomic data literacy.

Link:

- https://doi.org/10.1371/journal.pone.0196878

## Topics Considered But Not Central

### Motif Finding

Argument for inclusion:

- Important classic bioinformatics problem.
- Good source of randomized algorithms and Gibbs sampling.

Argument against inclusion:

- The course already has HMMs and FDR for probabilistic/statistical thinking.
- Motif finding could crowd assembly or omics.

Decision:

- Mention as optional project or enrichment topic.

### Genome Rearrangements

Argument for inclusion:

- Beautiful discrete math.
- Strong Compeau/Pevzner material.

Argument against inclusion:

- Less central for a first modern bioinformatics theory class than alignment, assembly, HMMs, phylogeny, and omics.

Decision:

- Optional project or final-day extension only.

### RNA Secondary Structure

Argument for inclusion:

- Excellent dynamic-programming example.

Argument against inclusion:

- Adds another biological domain and another DP formalism after students already have alignment.

Decision:

- Optional project or replacement topic if local alignment is shortened.

### Machine Learning Classifiers

Argument for inclusion:

- Modern bioinformatics uses ML everywhere.

Argument against inclusion:

- MATH 4100 already covers general data science / ML.
- A shallow ML survey would dilute the bioinformatics-specific math.

Decision:

- Keep ML language around method choice, validation, overfitting, and embeddings. Do not make supervised ML a module.

### Deep Learning / Foundation Models

Argument for inclusion:

- Modern and exciting.

Argument against inclusion:

- Too much background required for a mathematically honest treatment.
- Risk of becoming a tool demo.

Decision:

- Mention on final synthesis day as an extension of scoring functions, embeddings, and sequence models.
