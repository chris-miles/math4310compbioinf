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

## Bioinformatics and its data

Biologists rarely observe molecules directly. An experiment turns them into records: a DNA sequencer reports strings of letters, and RNA sequencing counts how many sequenced fragments came from each gene's RNA. Bioinformatics draws biological conclusions from records like these.

| Data | One record | Mathematics | Lessons |
|---|---|---|---|
| Sequencing reads and genomes | a string such as `ACGGTCA` | strings, dynamic programming, graphs | [1–7](lessons/lesson01-first-problem.html) |
| A genome to annotate | a long DNA string whose regions differ | probability, hidden Markov models | [8–12](lessons/lesson08-probability-markov-chains.html) |
| Sequences from several species | aligned rows of letters | distances and trees | [13–14](lessons/lesson13-distances-upgma.html) |
| Gene expression | counts per gene, for each sample or cell | linear algebra, statistics | [15–20](lessons/lesson15-expression-matrices-pca.html) |
| Reporter time courses | one cell's fluorescence over time | differential equations | [21–23](lessons/lesson21-gene-regulation-odes.html) |

The records are too many to read by eye: one sequencing run reports millions of reads, and the human genome is about 3.1 billion letters long. Before computing anything, we state the question precisely, for instance as the highest-scoring alignment of a read with a gene. Because measurements are noisy, we also need probability to judge whether a pattern, such as a high score or a group of similar cells, is stronger than chance would produce.

## Molecular biology in one page

DNA is a chain of *nucleotides*, each carrying one of four *bases*: adenine, cytosine, guanine, or thymine. Written from its 5′ end to its 3′ end, a DNA molecule is a string over {A, C, G, T}. In a cell, DNA usually has two strands running in opposite directions, with A paired to T and C to G, so either strand determines the other. A *genome* is the complete DNA of an organism or virus: 5,386 bases for the virus phiX174 in the [course data](data.html), about 4.6 million for the bacterium *E. coli*, and about 3.1 billion for one copy of the human genome.

A *gene* is a stretch of DNA that the cell copies into RNA, a step called *transcription*; RNA uses uracil (U) in place of thymine. In most human genes, *splicing* then cuts out stretches called *introns* and joins the remaining *exons*. How much RNA a cell makes from a gene is the gene's *expression*, which RNA sequencing estimates by counting reads.

Ribosomes translate the messenger RNA (mRNA) of a protein-coding gene three bases at a time. Each triple, a *codon*, names one of 20 amino acids or a stop; with $4^3 = 64$ codons, most amino acids have several. A protein is therefore a string over a 20-letter alphabet (@fig-about-gene-to-protein). [Most protein sequences in UniProt](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/where-does-the-data-come-from/sequence-data), the main protein database, were translated from sequenced genes rather than measured from protein.

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
