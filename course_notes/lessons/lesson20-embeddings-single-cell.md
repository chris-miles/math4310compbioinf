---
published: false
title: "Reading Maps of Single-Cell Data"
subtitle: "Lesson 20 · 2027-04-07 · Week 13"
---

## Core question

Why are modern embeddings useful and easy to overinterpret?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Where the data comes from: droplet single-cell RNA-seq, sparsity, dropout.

## The data as a mathematical object

- A kNN graph; a two-dimensional embedding.

## The question

- What does the picture preserve, and what does it invent?

## The model

- Neighborhood preservation only; distances between clusters carry no meaning.

## The technique

- t-SNE and UMAP at the level of objectives; the standard workflow matrix $\to$ PCA $\to$ kNN $\to$ clusters $\to$ labels.

## Worked example

- The same small dataset embedded under two hyperparameter settings.

## Why it works, and what it costs

- What is optimized, and why the optimum is not unique.

## Where the model breaks

- Cluster size and spacing are artifacts.
- Cell-type labels are claims layered on clusters, not data.

## Exercises

- Explain which steps of the workflow use methods from this course.
- Critique a provided embedding figure.
