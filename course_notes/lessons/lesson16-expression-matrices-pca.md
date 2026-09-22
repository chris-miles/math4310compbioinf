---
published: false
title: "Gene Expression and PCA"
subtitle: "Lesson 16 · 2027-03-24 · Week 11"
---

## Core question

What does PCA optimize, and what does the SVD have to do with it?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Where the data comes from: RNA-seq counts, samples by genes, more genes than samples.

## The data as a mathematical object

- A centered $n \times p$ matrix $X$.

## The question

- Which direction in gene space carries the most variance? Which rank-$r$ matrix is closest to $X$?

## The model

- Low rank plus noise; variance is signal.

## The technique

- Centering; $X = U S V^\top$; scores $US$, loadings $V$, variance explained $s_i^2 / \sum s_j^2$.

## Worked example

- A $4 \times 3$ integer matrix built to be rank one plus a small perturbation.

## Why it works, and what it costs

- Eckart--Young stated with a proof sketch; the first right singular vector maximizes variance.

## Where the model breaks

- Not centering makes the first component point at the mean.

## Exercises

- Compute the SVD of a rank-one matrix by hand.
- Show that scores are uncorrelated.

## Schedule

Self-contained: Midterm 2 precedes it and a weekend follows it.
