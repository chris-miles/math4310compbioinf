---
published: false
title: "Explaining Sequence Changes on a Tree"
subtitle: "Lesson 15 · 2027-03-22 · Week 11"
---

## Core question

Can we infer a tree by minimizing the number of changes?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Distances discard the characters; parsimony keeps them.

## The data as a mathematical object

- A character matrix (species $\times$ sites) and a candidate tree.

## The question

- Minimum number of character changes on a fixed tree; then, over all trees.

## The model

- Each change is equally costly; fewer changes is more plausible.

## The technique

- Fitch's set-union recurrence on a fixed tree, one site at a time.

## Worked example

- Four leaves, one site, scored by hand; then a second site to show additivity.

## Why it works, and what it costs

- Fitch is linear in tree size per site; the number of trees is super-exponential, so tree search is heuristic.

## Where the model breaks

- Long-branch attraction.
- Equal-cost changes ignore transition/transversion bias.

## Exercises

- Score a given tree with Fitch.
- Find a second tree with the same parsimony score.

## Schedule

Project starter list at the end of class; teams form over break.
