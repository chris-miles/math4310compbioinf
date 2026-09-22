---
published: false
title: "Evolutionary Distances and UPGMA"
subtitle: "Lesson 14 · Week 8"
---

## Core question

How do sequence differences become a tree?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Where the data comes from: a multiple alignment of homologous sequences across species.

## The data as a mathematical object

- An $n \times n$ distance matrix; the output is a rooted tree with leaves labeled by sequences.

## The question

- Find a tree whose path lengths reproduce the distances.

## The model

- Distances are additive along branches; UPGMA further assumes a molecular clock (ultrametric distances).

## The technique

- Distance from an alignment; multiple-hit correction conceptually; UPGMA merge step; neighbor-joining mentioned.

## Worked example

- A $4 \times 4$ matrix built to be ultrametric, merged by hand.

## Why it works, and what it costs

- UPGMA recovers the tree exactly when the matrix is ultrametric; $\mathcal{O}(n^3)$ naively.

## Where the model breaks

- Unequal rates break the clock and UPGMA returns the wrong topology.

## Exercises

- Run UPGMA on a given matrix.
- Construct a $3 \times 3$ non-ultrametric matrix where UPGMA fails.

## Schedule

PS6 due. PS7 out.
