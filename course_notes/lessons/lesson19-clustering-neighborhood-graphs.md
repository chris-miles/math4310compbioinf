---
published: false
title: "Distances, Clustering, and Neighborhood Graphs"
subtitle: "Lesson 19 · 2027-04-05 · Week 13"
---

## Core question

What does it mean for two samples to be near each other?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Cells or samples with no labels; we want groups.

## The data as a mathematical object

- Points in $\mathbb{R}^p$; a $k$-nearest-neighbor graph.

## The question

- Which samples belong together?

## The model

- A distance is a modeling choice: Euclidean, correlation, cosine each assume something different.

## The technique

- k-means and hierarchical clustering briefly; kNN graph construction; the curse of dimensionality by a small computation.

## Worked example

- Five points in the plane, kNN graph for $k = 2$, and one distance choice that changes it.

## Why it works, and what it costs

- kNN graph cost; why distances concentrate as $p$ grows.

## Where the model breaks

- Clusters found by every method on data with no clusters.

## Exercises

- Build a kNN graph by hand.
- Give a dataset where Euclidean and correlation distance disagree about the nearest neighbor.

## Schedule

PS8 due. PS9 out.
