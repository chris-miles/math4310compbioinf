---
title: "What is this course about?"
toc: false
nocite: |
  @markowetz2017, @goh2026expertise, @anthropic2026science, @hampton
---

How do we piece together a genome from DNA fragments? Figure out which cells are similar? Infer how species are related?

In this class, you'll see how math helps answer those questions. DNA becomes a string to compare, genome assembly becomes a path through a graph, and measurements from thousands of cells become a problem in linear algebra. We'll study alignment, probability models, evolutionary trees, gene-expression analysis, and gene regulation. You'll get to use mathematics you already know on biological problems, and learn some new mathematics along the way.

```{.python .execute}
#| label: setup
import sys

sys.path.insert(0, "../../course_notes/style")
import matplotlib.pyplot as plt

import coursefigs as cf

cf.use_style()
```

## What is bioinformatics?

Bioinformatics uses computation to answer biological questions from data. Most of the data comes from instruments. A DNA sequencer outputs short strings of A, C, G, and T, and RNA sequencing counts how many of those strings come from each gene. One sequencing run produces millions of these strings, far too many to inspect by hand.

The course covers these kinds of data:

| Data | Example | Mathematics | Lessons |
|---|---|---|---|
| Sequencing reads and genomes | a string such as `ACGGTCA` | strings, dynamic programming, graphs | [1–7](lessons/lesson01-first-problem.html) |
| A genome to annotate | a long DNA string whose regions differ | probability, hidden Markov models | [8–12](lessons/lesson08-probability-markov-chains.html) |
| Sequences from several species | aligned rows of letters | distances and trees | [13–14](lessons/lesson13-distances-upgma.html) |
| Gene expression | counts per gene, for each sample or cell | linear algebra, statistics | [15–20](lessons/lesson15-expression-matrices-pca.html) |
| Reporter time courses | one cell's fluorescence over time | differential equations | [21–23](lessons/lesson21-gene-regulation-odes.html) |

Measurements are noisy, so we also use probability to tell real patterns from chance.

## DNA, RNA, and proteins

DNA is a chain of *nucleotides*. Each carries one of four bases: adenine (A), cytosine (C), guanine (G), or thymine (T). Read from its 5′ end to its 3′ end, a DNA strand is a string over {A, C, G, T}. DNA usually has two strands running in opposite directions, with A paired to T and C to G, so either strand determines the other. A *genome* is all the DNA of an organism or virus. The phage phiX174 in the [course data](data.html) has 5,386 bases, *E. coli* about 4.6 million, and the human genome about 3.1 billion.

A *gene* is a stretch of DNA that the cell copies into RNA (*transcription*). RNA has uracil (U) where DNA has T. In most human genes, pieces called *introns* are then cut out and the remaining *exons* joined (*splicing*). How much RNA a cell makes from a gene is the gene's *expression*. RNA sequencing measures it by counting reads.

Ribosomes read the messenger RNA (mRNA) of a protein-coding gene three bases at a time. Each triple, called a *codon*, codes for one of 20 amino acids or for "stop". There are $4^3 = 64$ codons, so most amino acids have several. A protein is a string over a 20-letter alphabet (@fig-about-gene-to-protein). [Most protein sequences in UniProt](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/where-does-the-data-come-from/sequence-data) were translated from DNA, not measured directly.

```{.python .execute}
#| label: fig-about-gene-to-protein
#| code-fold: true
#| code-summary: "Draw the gene, its mRNA, and its protein"
#| fig-cap: "A made-up gene with two exons and one intron (top, both DNA strands), its mRNA (middle), and its protein (bottom). Grey bands carry the exons into the mRNA, which drops the intron and has U where the DNA has T; brackets group the mRNA into codons. The red codon TGG becomes UGG and then tryptophan (W). The right column says how each layer usually becomes data. Real introns are far longer."
#| fig-alt: "Three rows of monospace letters. Top: a DNA coding strand ATGAAG GTAAGCAG TGGGAGTAA, bracketed as exon 1, intron, and exon 2, with its complementary strand below and the 5-prime and 3-prime ends marked. Grey bands lead from the two exons to the middle row, the mRNA AUGAAGUGGGAGUAA. Brackets group the mRNA into five codons above the protein letters M, K, W, E, and the word stop. The codon TGG, the codon UGG, and the letter W have red backgrounds. Labels at the left name DNA, transcription and splicing, mRNA, translation, and protein; labels at the right read: read by DNA sequencing, counted by RNA-seq, usually inferred from the gene."
#| fig-width: 6
#| fig-height: 2.6
exons = ["ATGAAG", "TGGGAGTAA"]
introns = ["GTAAGCAG"]
fig, ax = plt.subplots()
cf.gene_to_protein(ax, exons, introns, codon=2,
                   data=["read by DNA\nsequencing", "counted by\nRNA-seq",
                         "usually inferred\nfrom the gene"])
plt.show()
```

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
**References**
:::
