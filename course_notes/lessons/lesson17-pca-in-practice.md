---
published: false
title: "PCA in Practice"
subtitle: "Lesson 17 · 2027-03-29 · Week 12"
---

## Core question

What can a PCA plot tell us, and what can it hide?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Two labs, two batches, one beautiful and wrong PCA plot.

## The data as a mathematical object

- Raw counts; the choices that turn them into $X$.

## The question

- Which preprocessing choices change the picture, and how do we know which picture to trust?

## The model

- Same low-rank model, with its assumptions made explicit: library size, log transform, scaling, gene filtering.

## The technique

- Normalization as modeling; scree plots; a batch-effect experiment.

## Worked example

- A small matrix with a planted batch effect, before and after correction.

## Why it works, and what it costs

- Why scaling changes the answer: PCA is not invariant to units.

## Where the model breaks

- Variance is not signal when the largest variance is technical.
- Two clusters in PC1 may be one confounder.

## Exercises

- Recompute PCA with and without scaling and explain the difference.
- Design a check that would reveal a batch effect.

## Schedule

PS8 out.
