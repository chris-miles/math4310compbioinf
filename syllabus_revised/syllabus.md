---
title: "MATH 4310 Computational Bioinformatics"
subtitle: "Spring Semester 2027"
---

# MATH 4310 Computational Bioinformatics

**Spring Semester 2027**

## Course Description

Provides the theoretical foundations of modern bioinformatics, focusing on the essential algorithms and mathematical models for analyzing biological data. The course covers the principles behind core computational tasks, from classical alignment algorithms to modern large-scale data challenges. Key topics include fundamental approaches in probabilistic inference, mechanistic modeling, and dimensionality reduction. The emphasis is on rigorous, theory-revealing implementations that explain *why* bioinformatics tools work.

[Building these tools from scratch, rather than calling library functions, asks you to do the kind of work bioinformatics careers are built on: deriving correct algorithms, evaluating which method fits a real dataset, communicating reasoning to readers from another discipline, and producing reproducible computational scholarship.]{.underline}

- **Credit Hours:** 3
- **General Education:** This course does not meet a GE requirement.
- **Prerequisites (Required):** A grade of C or better in MATH 2270 (Linear Algebra), CS 2420 (Introduction to Algorithms & Data Structures), and BIOL 2030 (Genetics).
- **Prerequisites (Recommended):** Ideally taken after BIOL 3150.

## Course Outcomes and Objectives

By the end of this course, you will be able to:

- **Analyze** the correctness and computational complexity of dynamic programming algorithms for sequence alignment, and justify the use of scoring matrices like BLOSUM.
- **Implement and interpret** probabilistic models, including Markov Chains and Hidden Markov Models (HMMs), to solve problems like gene finding.
- **Explain** how graph theory (e.g., de Bruijn graphs) provides a scalable solution to the problem of genome assembly from short-read sequencing data.
- **Use** linear algebra to justify dimensionality reduction techniques (PCA/SVD) for transcriptomics, and formulate simple dynamical systems models for gene regulatory networks.

## Durable Skills in This Course

Durable skills are broad abilities that transfer across courses, research, jobs, and civic life: things like solving unfamiliar problems, communicating clearly, working responsibly with others, and learning new tools. The University of Utah's [Durable Skills Framework](https://durableskills.utah.edu/) gives us a shared way to name those skills and help you recognize where you are practicing them in your coursework.

The technical objectives above are the primary content of MATH 4310. This course also explicitly develops five durable skills:

| Durable skill | Where it is practiced and assessed |
|---|---|
| [**Critical Thinking & Problem Solving**]{.underline} | [Correctness proofs, complexity analysis, method trade-offs, and the Failure-Mode Reflection on each problem set.]{.underline} |
| [**Information & Technology Literacy**]{.underline} | [Documenting data sources, package versions, AI-assisted coding, validation checks, and tool choice in problem sets and the final project.]{.underline} |
| [**Communication**]{.underline} | [Short technical explanations for peers from math, CS, or biology backgrounds; extended in the final-project README.]{.underline} |
| [**Collaboration & Teamwork**]{.underline} | [The Week 6 triad problem set and the final group project, each with contribution documentation.]{.underline} |
| [**Career & Self-Development**]{.underline} | [A reproducible final-project repository, two feedback checkpoints, and a revision memo describing how the work changed.]{.underline} |

The schedule at the end of this syllabus notes which one or two of these skills each module foregrounds.

## Course Requirements

Your grade is based on the following components:

| Component | Weight |
|---|---:|
| Problem Sets and Project Workshops (10 problem sets, including one triad set, plus two project checkpoints) | 40% |
| Short Midterm Exam 1 | 15% |
| Short Midterm Exam 2 | 15% |
| Final Project Repository | 30% |

### Problem Sets and Project Workshops (40%)

Ten problem sets, due roughly weekly during the technical core of the course. Each set is a mix of mathematical work (e.g., prove the correctness of a recurrence) and a small implementation task (e.g., write a 40-line Viterbi algorithm). Nine sets are individual; **Problem Set 5 (Week 6) is done in triads** as a low-stakes warm-up for the final project.

Coding is required in this course, but the coding tasks are intentionally small and theory-revealing. The purpose is not to turn MATH 4310 into a software engineering course; the purpose is to make the mathematical models concrete enough that you can inspect, test, and trust an implementation. Problem-set coding tasks will usually include a recurrence, model description, or high-level pseudocode, plus starter support for data loading or tests when appropriate. You remain responsible for the core algorithm, validation checks, and explanation of how the implementation matches the mathematics.

[Each problem set ends with a **Failure-Mode Reflection**: a 120 to 180 word paragraph in which you identify one assumption your solution depends on, imagine that assumption tightening in a realistic biological or computational scenario (noisier reads, missing data, a 100x larger dataset, less compute), describe where your solution fails or degrades, and sketch what you would change. Anchor your reflection in something concrete: a counterexample, a test case, a dataset property, a theorem condition, or a specific tool or data limitation. Write it so a peer who has not seen your code can follow the argument. The reflection is graded on a single 0/1/2 rubric covering both specificity and audience-awareness: *absent, generic, or unclear to a peer reader* / *present but partial: assumption named but failure or fix is generic, or argument is hard to follow* / *specific, concrete, and clear to a peer reader*.]{.underline}

Two project workshops, in Weeks 11 and 14, count in this category. They are checkpoint submissions that keep the final project scoped and reproducible: the first establishes the problem, data source, method, and repository skeleton; the second demonstrates preliminary results and a validation plan.

### Final Project Repository (30%)

[Working in groups of three (formed in Week 9), you will apply one of the algorithms or models from this course to a real bioinformatics problem of your choice. A pre-curated list of starter problems will be provided; you may also propose your own.]{.underline}

[Project preparation will begin before teams form through examples, short planning activities, and discussion of feasible scopes. The exact project format and starter options will be announced in class and on Canvas. The goal is to help teams choose a realistic application of a course method without requiring an open-ended research project from scratch.]{.underline}

[The deliverable is a single Git repository containing three things:]{.underline}

[1. **A from-scratch implementation** (or principled extension) of one course algorithm, applied to documented data, with the reproducibility basics built in: a `requirements.txt` or equivalent, run instructions, data source and provenance, and at least one validation check (e.g., a test on a small input with a known answer, or a comparison against a published reference).]{.underline}
[2. **A project README** (one markdown file, roughly 3 to 4 pages of content) that explains the work to a CS or biology peer who has not taken this course. It should cover the problem, your method, your results, where the method fails, and what you would do with another month, written so a non-specialist can follow. Open with a short **Team and Contributions** paragraph: who did what, who reviewed whose code, and how the team handled coordination or disagreements.]{.underline}
[3. **A revision memo** (one page, separate file): summarize one piece of feedback you received at the checkpoints, the change you made in response, one approach that didn't work and how you adapted, and what you would improve next.]{.underline}

[Two checkpoints structure the project and are graded as project workshops in the Problem Sets and Project Workshops category:]{.underline}

[- **Week 11**: a one-page proposal plus a repository skeleton (README stub, problem statement, planned approach). Submitted for credit and feedback.]{.underline}
[- **Week 14**: preliminary results in the repository. Submitted for credit and feedback.]{.underline}

[The final repository is due in Finals Week. There is no cumulative final exam; the final project is the course's cumulative assessment. Grading uses a six-criterion rubric (implementation applied to a bioinformatics problem; reproducibility; tool/data/AI documentation; README clarity and audience-awareness; team and contributions evidence; revision memo) so that each kind of work, and each durable skill it carries evidence of, is credited distinctly.]{.underline}

### Short Midterm Exams (15% each)

There are two shorter in-class written exams rather than a cumulative final exam.

**Short Midterm Exam 1** covers the sequence-algorithm and graph-assembly core: scoring, dynamic programming for alignment, affine gaps at the recurrence level, seed-and-extend intuition, de Bruijn graphs, and Eulerian paths.

**Short Midterm Exam 2** covers probabilistic sequence models and phylogenetics: Markov chains, likelihood-ratio classification, HMMs, Viterbi, forward-backward/posterior decoding, and distance- or parsimony-based tree methods.

The final-project repository, due in Finals Week, replaces the cumulative final exam.

## Tools, Data, and AI Use

This course treats AI assistants and bioinformatics software the same way we treat any powerful computational tool: useful, but not self-verifying. You may use AI coding assistants (e.g., Copilot, Claude, ChatGPT) on problem-set implementations and on the final project, with three requirements:

1. **Disclose** any AI-assisted step in the relevant problem set or in the project README. A line like "I used Copilot to draft the inner loop of the alignment function and verified it against the small test on page X" is sufficient.
2. **Verify**. AI output is unverified output. Confirm correctness with tests, derivations, hand-checkable examples, or comparisons to known answers before submitting.
3. **Explain**. You may be asked to explain how your implementation matches the recurrence, model, or algorithmic idea used in class. Code that runs but cannot be explained is not sufficient evidence of understanding.

The same standard applies to using existing software libraries on the final project: name what you used, document version and source, and verify the output. AI assistants and libraries are tools to be used responsibly, not shortcuts that exempt you from understanding your own work. The useful modern skill is not pretending these tools do not exist; it is knowing how to check whether their output is trustworthy.

Generative-AI use that is not disclosed, or that is used to bypass the assignment's purpose (e.g., asking an LLM to write a proof you submit verbatim), is academic misconduct under University Policy 6-410.

## Group Work

The triad problem set in Week 6 and the final project use the same conventions:

- **Triads, not pairs or larger.** Triads are large enough to share work and small enough that individual contribution stays visible.
- **Contribution log.** Each shared submission includes a one-paragraph contribution log signed by all members.
- **Individual addenda.** On the triad problem set, each student writes a short addendum (one paragraph) explaining one proof step or code design choice in their own words. This is the unit of grading for individual participation.
- **Documented non-participation.** If a teammate is not contributing, document it in the log; individual grades may be adjusted accordingly.

## Required and Recommended Readings

There is no required textbook. Required readings, including academic papers and technical notes, will be posted on Canvas. A list of recommended supplementary references (Durbin et al., Pevzner & Compeau, Strang) will be available on the course site for students who want a textbook companion.

## Grading Scale

| Letter | Range |
|---|---|
| A | 93 to 100% |
| A- | 90 to 92.9% |
| B+ | 87 to 89.9% |
| B | 83 to 86.9% |
| B- | 80 to 82.9% |
| C+ | 77 to 79.9% |
| C | 73 to 76.9% |
| C- | 70 to 72.9% |
| D+ | 67 to 69.9% |
| D | 63 to 66.9% |
| D- | 60 to 62.9% |
| E | below 60% |

## Workload Expectations

This is a 3-credit upper-division course. Plan for roughly 6 to 9 hours of work per week outside of class meetings: reading, problem sets, implementation, exam preparation, and (in the second half) the final project. Problem sets that include implementation tasks tend to take longer than purely mathematical sets; budget accordingly.

## Attendance

Regular attendance is expected. Lecture material includes worked derivations and live coding that is not always reproduced in posted notes. If you must miss a class, you are responsible for the material; the instructor will post slides and any code demos to Canvas.

## University Policies

Mandatory institutional policies, including those covering the Americans with Disabilities Act, Safety at the U, Addressing Sexual Misconduct, and Academic Misconduct, are maintained at:

<https://cte.utah.edu/instructor-education/syllabus/institutional-policies.php>

Students with disabilities or accommodation needs should contact the Center for Disability and Access (CDA) at 801-581-5020 or [disability.utah.edu](https://disability.utah.edu/).

## Preliminary Course Schedule

| Week | Module / Topic | Due |
|---|---|---|
| 1 | **Module 1.** Algorithmic complexity; scoring matrices and log-odds ratios. | |
| 2 | **Module 1.** Global alignment (Needleman-Wunsch); dynamic programming. | PS 1 |
| 3 | **Module 1.** Local alignment (Smith-Waterman); affine gap penalties; seed-and-extend intuition. | PS 2 |
| 4 | **Module 1.** Graph theory for assembly; de Bruijn graphs and Eulerian paths. | PS 3 |
| 5 | **Short Midterm Exam 1**; **Module 2.** Discrete probability and Markov chains for CpG islands. | PS 4 |
| 6 | **Module 2.** Hidden Markov Models; the Viterbi algorithm for gene finding. | PS 5 (triad) |
| 7 | **Module 2.** Forward-Backward algorithms; HMM training. | PS 6 |
| 8 | **Module 3.** Phylogenetics: distance-based (UPGMA) vs. character-based (parsimony); project preparation. | PS 7 |
| 9 | Spring Break. Project teams form by end of week. | |
| 10 | **Short Midterm Exam 2**; **Module 4.** Linear algebra of 'omics data; PCA and SVD; project preparation. | PS 8 |
| 11 | **Module 4.** Statistics for high-throughput tests; False Discovery Rate. | PS 9; **Project Checkpoint 1** (proposal + repo skeleton) |
| 12 | **Module 4.** Modern dimensionality reduction (UMAP, t-SNE); clustering and neighborhood graphs. | PS 10 |
| 13 | **Module 5.** Modeling gene regulation with ODEs; Hill functions and bistability. | |
| 14 | **Module 5.** Biological networks and toggle switches; project workshop. | **Project Checkpoint 2** (preliminary results) |
| 15 | Course synthesis; modern extensions; final project work session. | |
| Finals Week | **Final Project** repository due. No cumulative final exam. | |

[**Durable skills emphasized by module:**]{.underline}

| Module | Durable skills emphasized |
|---|---|
| [Module 1: Algorithms]{.underline} | [Critical Thinking & Problem Solving; Communication]{.underline} |
| [Module 2: Probabilistic Models]{.underline} | [Critical Thinking & Problem Solving; Information & Technology Literacy; Collaboration & Teamwork]{.underline} |
| [Module 3: Phylogenetics]{.underline} | [Critical Thinking & Problem Solving; Communication]{.underline} |
| [Module 4: 'Omics Math]{.underline} | [Information & Technology Literacy; Communication; Career & Self-Development]{.underline} |
| [Module 5: Mechanistic Modeling and Synthesis]{.underline} | [Critical Thinking & Problem Solving; Information & Technology Literacy; Career & Self-Development]{.underline} |
