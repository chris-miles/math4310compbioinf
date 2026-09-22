---
published: false
title: "Modeling How Gene Activity Changes"
subtitle: "Lesson 21 · 2027-04-12 · Week 14"
---

## Core question

What is a mechanistic model trying to explain?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. The headings below follow the order every lesson uses; the bullets are what the finished lesson will cover.
:::

## Motivation

- Everything so far described data; this model describes a mechanism.
- Where the data comes from: fluorescent reporters and time courses.

## The data as a mathematical object

- A concentration $x(t)$ and parameters.

## The question

- What happens over time, and what happens if we change a parameter?

## The model

- $\dot{x} = \beta - \alpha x$: constant production, first-order degradation; well-mixed and deterministic.

## The technique

- Steady state; time scale $1/\alpha$; numerical simulation with a forward-Euler step.

## Worked example

- $\beta = 2$, $\alpha = 1$, $x(0) = 0$, three Euler steps by hand and the exact solution.

## Why it works, and what it costs

- Stability of the steady state from the sign of the slope.

## Where the model breaks

- Low copy number makes the deterministic model wrong.
- Simulation is evidence, not proof.

## Exercises

- Solve the linear model and check the simulation against it.
- Add a second species and write the new system.

## Schedule

PS9 due. PS10 out.
