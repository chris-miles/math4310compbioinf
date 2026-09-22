---
title: "What is this class about?"
subtitle: "Biological questions, computational methods, and knowing what to trust"
toc: true
---

A sequencer produces DNA reads. To find a sequence difference, we need to decide where those reads belong and which differences to believe. An experiment measures gene activity in thousands of cells. To identify groups of cells, we need to choose what similarity means and check whether the groups reflect biology or how the samples were processed.

MATH 4310 connects the mathematics, computing, and biology in the Bioinformatics major. We will use mathematical models and small programs to understand what happens inside an analysis, choose methods that fit the data, and check the conclusions. The same ideas are useful for students coming from mathematics who want to work on biological questions.

## The questions behind the tools

The course starts with strings: how should we compare two DNA sequences, and how can we reconstruct a longer sequence from short fragments? Those questions lead to dynamic programming and graphs. We then use probability models to infer regions of a genome whose labels are hidden, and trees to reason about evolutionary relationships. Later, we work with gene-expression matrices, dimensionality reduction, clustering, and the problem of testing many genes at once. We finish with models of gene regulation and a look at current research.

Biology determines what the mathematics needs to do. An alignment between two complete genes has different requirements from an alignment placing a short read inside a chromosome. A cluster of cells means little until we know what was measured and what the distance between two cells represents. Getting the mathematical answer right includes asking whether we posed a useful biological problem.

Computational work can also change the questions we are able to ask. In *All biology is computational biology*, Florian Markowetz argues that computation helps organize biological knowledge and turn ideas about living systems into testable models [@markowetz2017]. That is a good reason to learn these methods whether the eventual work happens mostly at a computer or in a laboratory.

## Writing code to understand a method

We will write small implementations, work examples by hand, and explain why an algorithm gives the answer it claims to give. The code should make the method visible. A short alignment program lets us inspect the recurrence, change a gap penalty, and predict what should happen before rerunning it.

Most working bioinformaticians use established alignment and assembly software. We should too when doing an analysis that needs those tools. Reimplementing a small version in class gives us a way to understand their objectives and limitations. A proof may establish that a program finds the highest-scoring alignment under a particular rule; deciding whether that rule suits the experiment remains a scientific judgment.

The programming here serves that reasoning. We will spend more time choosing a model, checking an output, or explaining a failure than building a large software system. An independently calculated tiny example can tell us more about correctness than a successful run on a huge dataset.

## Where this fits at Utah

The closest comparisons are BIOL 3150 and MATH 4100. MATH 4310 adds depth in the mathematical methods behind biological analysis. The courses have different emphases, with useful overlap; individual offerings may vary.

| Course | Main emphasis |
|---|---|
| [BIOL 3150: Genomics and Bioinformatics](https://catalog.utah.edu/courses/0197411) | Genome structure and function together with practical analysis of genomic data. |
| [MATH 4100: Introduction to Data Science](https://catalog.utah.edu/courses/0161202/general-aoYks) | Acquiring, cleaning, exploring, and analyzing data; using Python tools and communicating results across applications. |
| MATH 4310: Computational Bioinformatics | Alignment algorithms, hidden Markov models, assembly graphs, the linear algebra of PCA, and gene-regulation models, studied through derivations and short implementations. |

BIOL 3150 connects genome biology with hands-on analysis. Here we can take a familiar task such as sequence comparison and spend time deriving its recurrence, proving what it optimizes, and examining what a fast search heuristic can miss. MATH 4100 develops the broader data science workflow. We concentrate on selected methods in biological settings, including the linear algebra behind PCA and what its directions mean for a gene-expression matrix. The [course information page](course-info.html#prerequisites) lists the preparation expected for MATH 4310.

## Why learn this in the age of AI?

It is becoming easier to ask for an analysis in ordinary language and receive code, a plot, and an explanation. [Claude Science](https://www.anthropic.com/news/claude-science-ai-workbench), for example, combines an AI assistant with scientific tools and computing resources, including ways to inspect the code behind its outputs [@anthropic2026science]. Tools like this can reduce the work of setting up and carrying out an analysis. We still need to specify which analysis would answer the question and what evidence would make its result credible.

Suppose an assistant is asked to find cell types in a gene-expression dataset. It returns a colorful plot with several clusters. To assess it, we need to know how expression values were transformed, what distances the method used, and whether cells separated by sample-processing batch rather than cell type. We can ask the assistant to compare alternatives and run checks, but we need enough understanding to decide which comparisons matter. A polished plot is one output of that process, not evidence that every choice was appropriate.

Goh and colleagues make a related argument in *Rethinking bioinformatics expertise in the era of artificial intelligence*: automation shifts attention toward analysis design, data quality, and interpretation [@goh2026expertise]. We will practice those responsibilities on problems small enough to inspect. We can check an alignment's score directly, construct a graph with two possible assemblies, or compare a probabilistic algorithm with exhaustive enumeration. These habits help us evaluate work produced by a library, a collaborator, or an AI assistant. Assignment-specific permissions still follow the [course AI policy](course-info.html#tools-data-and-ai-use).

## New measurements keep changing the problems

Mature tools solve many routine tasks well. New data types and larger datasets still change the computational problem. Longer reads provide different information about repeats. Measurements combining molecular activity with a cell's location require us to consider spatial relationships. Increasing the number of samples can make an algorithm's memory use the limiting factor.

Even reference sequences are changing: a pangenome graph can represent alternatives from many genomes, creating alignment problems beyond comparing two strings. The Minigraph-Cactus work is one example of methods developed for that setting [@hickey2024pangenome]. Learning about states, recurrences, graphs, and uncertainty gives us a starting point for understanding such methods and for developing new ones when the available tools do not answer the question.

[Browse the lessons](lectures.html) or [return to the schedule](index.html#course-schedule).

::: {#refs}
**References**
:::
