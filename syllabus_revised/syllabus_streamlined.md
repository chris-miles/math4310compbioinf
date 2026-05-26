---
title: "MATH 4310 Computational Bioinformatics"
subtitle: "Spring Semester 2027"
---

# MATH 4310 Computational Bioinformatics

**Spring Semester 2027**

## Course Overview

MATH 4310 studies the mathematical and algorithmic foundations of modern bioinformatics. We will build and analyze core methods for biological data, including sequence alignment, probabilistic models, graph-based genome assembly, dynamical systems, and dimensionality reduction for 'omics data. The emphasis is on understanding why the methods work, not just how to run existing software.

- **Credit hours:** 3
- **General Education:** This course does not meet a GE requirement.
- **Prerequisites:** C or better in MATH 2270, CS 2420, and BIOL 2030.
- **Recommended preparation:** BIOL 3150.
- **Section details:** Instructor, meeting time, location, office hours, and contact information will be posted in Canvas and the official course schedule.

## Course Outcomes

By the end of this course, you will be able to:

- Analyze the correctness and computational complexity of sequence-alignment algorithms.
- Implement and interpret Markov chains and Hidden Markov Models for bioinformatics problems.
- Explain how graph theory, especially de Bruijn graphs, supports scalable genome assembly.
- Use linear algebra to justify PCA/SVD methods for transcriptomics and formulate simple dynamical-systems models for gene regulation.

## Durable Skills

Durable skills are broad abilities that transfer across courses, research, jobs, and civic life. The University of Utah's [Durable Skills Framework](https://durableskills.utah.edu/) gives us a shared way to name those skills and notice where they are practiced in coursework.

MATH 4310 explicitly develops five durable skills:

| Durable skill | Where it appears in this course |
|---|---|
| [**Critical Thinking & Problem Solving**]{.underline} | [Proofs, complexity analysis, method comparison, and failure-mode reflections.]{.underline} |
| [**Information & Technology Literacy**]{.underline} | [Data provenance, software/version documentation, validation checks, and responsible AI/tool use.]{.underline} |
| [**Communication**]{.underline} | [Short technical explanations in problem sets and a longer final-project README for a peer audience.]{.underline} |
| [**Collaboration & Teamwork**]{.underline} | [One triad problem set and the group final project, both with contribution documentation.]{.underline} |
| [**Career & Self-Development**]{.underline} | [A reproducible final-project repository, feedback checkpoints, and a revision memo.]{.underline} |

## Graded Work

| Component | Weight |
|---|---:|
| Problem sets and project workshops | 40% |
| Short Midterm Exam 1 | 15% |
| Short Midterm Exam 2 | 15% |
| Final project repository | 30% |

**Problem sets and project workshops.** There will be 10 problem sets plus two project workshops/checkpoints. Problem sets combine mathematical analysis, short implementations, and a brief failure-mode reflection. Coding is required, but tasks are intentionally small and theory-revealing: the goal is to connect algorithms to tests and explanations, not to turn this into a software engineering course. Problem Set 5 is completed in triads; the others are individual unless otherwise announced. The two project workshops structure the final project proposal and preliminary-results stages.

**Final project.** In groups of three, you will apply one course method to a real or realistic bioinformatics problem. Project preparation will begin before teams form through examples, short planning activities, and discussion of feasible scopes. The exact project format and starter options will be announced in class and on Canvas. The final submission is a reproducible Git repository with code, documented data sources, a README for a peer audience, contribution notes, and a short revision memo. Two checkpoints, in Weeks 11 and 14, provide structure and feedback and count as project workshops.

**Exams.** There are two shorter in-class written exams rather than a cumulative final exam. Short Midterm Exam 1 covers sequence algorithms and graph assembly. Short Midterm Exam 2 covers probabilistic sequence models and phylogenetics. The final project repository, due in Finals Week, replaces the cumulative final exam.

Detailed prompts, rubrics, deadlines, and examples will be posted on Canvas.

## Tools, Data, and AI Use

You may use AI coding assistants and existing software tools when an assignment permits them, but you remain responsible for the correctness and interpretation of your work.

- Disclose AI-assisted coding, external libraries, and outside tools when they affect your submitted work.
- Verify outputs with tests, derivations, hand-checkable examples, or comparisons to known answers.
- Be prepared to explain how your implementation matches the mathematical recurrence, model, or algorithmic idea used in class.
- Do not submit AI-generated explanations, proofs, or code you do not understand.

Undisclosed or inappropriate use of AI or other tools may be treated as academic misconduct under University Policy 6-410.

## Readings

There is no required textbook. Required readings, papers, technical notes, and recommended supplementary references will be posted on Canvas.

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

## Workload

This is a 3-credit upper-division course. Expect about 6 to 9 hours per week outside of class for reading, problem sets, implementation, exam preparation, and project work.

## Attendance

Regular attendance is expected. Class meetings include derivations, examples, and live coding that may not be fully reproduced in posted materials. If you miss class, you are responsible for catching up using Canvas materials and classmates' notes.

## University Policies

Mandatory institutional policies, including those covering the Americans with Disabilities Act, Safety at the U, Addressing Sexual Misconduct, and Academic Misconduct, are maintained at:

<https://cte.utah.edu/instructor-education/syllabus/institutional-policies.php>

Students with disabilities or accommodation needs should contact the Center for Disability and Access (CDA) at 801-581-5020 or [disability.utah.edu](https://disability.utah.edu/).

## Preliminary Schedule

| Week | Topic | Due |
|---|---|---|
| 1 | Algorithmic complexity; scoring matrices and log-odds ratios | |
| 2 | Global alignment; Needleman-Wunsch; dynamic programming | PS 1 |
| 3 | Local alignment; Smith-Waterman; affine gap penalties; seed-and-extend intuition | PS 2 |
| 4 | Graph theory for assembly; de Bruijn graphs and Eulerian paths | PS 3 |
| 5 | Short Midterm Exam 1; discrete probability; Markov chains for CpG islands | PS 4 |
| 6 | Hidden Markov Models; Viterbi algorithm for gene finding | PS 5, triad |
| 7 | Forward-Backward algorithms; HMM training | PS 6 |
| 8 | Phylogenetics; distance-based and character-based methods; project preparation | PS 7 |
| 9 | Spring Break; project teams form by end of week | |
| 10 | Short Midterm Exam 2; linear algebra of 'omics data; PCA and SVD; project preparation | PS 8 |
| 11 | High-throughput testing; False Discovery Rate | PS 9; Project Checkpoint 1 |
| 12 | Modern dimensionality reduction; UMAP and t-SNE; clustering and neighborhood graphs | PS 10 |
| 13 | Gene-regulation models; ODEs; Hill functions; bistability | |
| 14 | Biological networks; toggle switches; project workshop | Project Checkpoint 2 |
| 15 | Course synthesis; modern extensions; final project work session | |
| Finals Week | Final project repository due; no cumulative final exam | |
