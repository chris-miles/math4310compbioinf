---
published: false
title: "Multiple Testing and False Discovery Rates"
subtitle: "Lesson 18 · 2027-03-29 · Week 12"
---

## Core question

What changes when we test ten thousand genes at once?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Twenty thousand null hypotheses at $\alpha = 0.05$ give a thousand discoveries from noise.

## The data as a mathematical object

- A vector of $m$ p-values.

## The question

- Choose a rejection set with a controlled fraction of false discoveries.

## The model

- Null p-values are uniform; independence or positive dependence.

## The technique

- FWER versus FDR; Benjamini--Hochberg step-up; the sorted-p-value picture; a differential-expression mini-case with effect size versus p-value.

## Worked example

- Eight p-values, BH by hand at $q = 0.1$.

## Why it works, and what it costs

- BH controls FDR at $q$ under independence; proof sketch of the load-bearing inequality.

## Where the model breaks

- Dependence between genes.
- Significant is not large; the volcano plot's two axes.

## Exercises

- Apply BH to a given list.
- Simulate the null and count false discoveries with and without BH.

## Schedule

Project checkpoint 1 due.
