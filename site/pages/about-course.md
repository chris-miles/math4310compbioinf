---
title: "What is this course about?"
toc: false
nocite: |
  @markowetz2017, @goh2026expertise, @anthropic2026science
---

How do we piece together a genome from DNA fragments? Figure out which cells are similar? Infer how species are related?

In this class, you'll see how math helps answer those questions. DNA becomes a string to compare, genome assembly becomes a path through a graph, and measurements from thousands of cells become a problem in linear algebra. We'll study alignment, probability models, evolutionary trees, gene-expression analysis, and gene regulation. You'll get to use mathematics you already know on biological problems, and learn some new mathematics along the way.

## Where this fits at Utah

| Course | Emphasis |
|---|---|
| [BIOL 3150: Genomics and Bioinformatics](https://catalog.utah.edu/courses/0197411) | Genome biology and practical analysis of genomic data. |
| [MATH 4100: Introduction to Data Science](https://catalog.utah.edu/courses/0161202/general-aoYks) | The broader workflow of collecting, cleaning, exploring, and analyzing data. |
| MATH 4310: Computational Bioinformatics | How the methods work, why they work, and where they fail. |

There's overlap. Here we'll spend more time on the mathematics: deriving an alignment algorithm, for example, or understanding what PCA does to a gene-expression matrix.

## Why take this in the age of AI?

Tools like [Claude Science](https://www.anthropic.com/news/claude-science-ai-workbench) make it easier to carry out an analysis. But what should you ask them to do? And how will you know whether the answer makes sense?

Suppose an analysis groups cells into clusters. Are those different cell types, or did the samples just get processed differently? Understanding the method helps you decide what to check before giving those clusters a biological interpretation.

You probably won't need to write your own sequence aligner. Established software handles many routine tasks well. Learning how alignment works helps you choose a tool, understand its settings, and spot results you shouldn't trust. Those skills matter whether you run the analysis yourself or ask an AI to do it.

The problems also keep changing. Longer DNA reads, measurements of cells in their spatial context, and larger datasets create a need for new methods. The algorithms we study give us a starting point for understanding and developing those methods.

By the end, you should be more comfortable opening up an unfamiliar method, figuring out what it assumes, and deciding whether it fits your question. That's useful even when the software does all the computing.

::: {#refs}
**Further reading**
:::
