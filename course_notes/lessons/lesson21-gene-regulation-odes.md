---
published: false
title: "ODE Models of Gene Regulation"
subtitle: "Lesson 21 · 2027-04-12 · Week 14"
nocite: |
  @rosenfeld2002, @mangelGeneNetworks
---

## Core question

After a signal turns on protein production, how much protein accumulates and how quickly does the cell respond?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. It begins a three-lesson gene-regulation sequence: response dynamics, cell-state switches, and noisy time-course measurements.
:::

## Watching a gene respond

The expression matrices in earlier lessons compared samples at a collection time. Suppose we instead switch on production of a protein and follow its abundance over time. A cell must make new protein molecules while existing molecules are removed or diluted as the cell grows. The eventual abundance and the speed of the response answer different biological questions: a large response can still take a long time to develop. Experiments with engineered gene circuits have used this distinction to investigate how regulation changes response times [@rosenfeld2002].

A *fluorescent reporter* is a protein whose emitted light helps us track gene expression. Repeated images supply records with a cell or culture identifier, a time $t_k$, and a measured intensity $y_k$. The intensity is an observation; our model describes a concentration $x(t)$ that changes continuously between images. For the first calculation, suppose calibration converts intensity to concentration and ignore reporter maturation and measurement error. We want to predict the time course after production changes, then distinguish a change in final abundance from a change in response speed.

## A balance of production and removal

Assume the cell is a well-mixed compartment with an effective constant production rate $\beta\geq0$, measured in concentration per unit time. Suppose removal is proportional to abundance, with rate constant $\alpha>0$ in inverse time. This effective removal can include both degradation and dilution. Over a short interval $\Delta t$, production adds approximately $\beta\Delta t$ and removal subtracts approximately $\alpha x(t)\Delta t$ [@mangelGeneNetworks]. Thus
$$
x(t+\Delta t)-x(t)\approx\bigl(\beta-\alpha x(t)\bigr)\Delta t.
$$
Dividing by $\Delta t$ and taking the limit gives the ordinary differential equation, or ODE,
$$
\dot x(t)=\beta-\alpha x(t),\qquad x(0)=x_0.
$$
The dot denotes a time derivative. The inputs are the rates and initial concentration; the desired output is the trajectory $x(t)$. Constant production and proportional removal are assumptions we can test against a measured response, not consequences of fitting a smooth curve.


## The technique

- Solve for $x^*=\beta/\alpha$ and derive
  $$x(t)=x^*+(x_0-x^*)e^{-\alpha t}.$$
  Interpret $1/\alpha$ as the relaxation time and show that increasing production changes abundance without changing this time scale.
- Approximate the same model with forward Euler, $x_{k+1}=x_k+h(\beta-\alpha x_k)$, where $h$ is the time step.

The exact solutions below isolate response speed. The parameter pairs $(\beta,\alpha)=(2,1)$ and $(1,1/2)$ both give $x^*=2$, but the second response relaxes twice as slowly. A final expression measurement alone would miss that difference.

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

```{.python .execute}
#| label: fig-regulation-response-times
#| code-fold: true
#| fig-cap: "Exact responses from zero for two production/removal pairs with the same steady concentration. Equal final abundance can hide different biological response times."
#| fig-alt: "Two increasing concentration curves approach the same horizontal level of two; the beta two, alpha one curve approaches faster than the beta one, alpha one-half curve."
times = np.linspace(0, 8, 161)
fig, ax = plt.subplots()
for beta, alpha, line_style in [(2.0, 1.0, "-"), (1.0, 0.5, "--")]:
    steady_state = beta / alpha
    response = steady_state * (1 - np.exp(-alpha * times))
    ax.plot(times, response, line_style,
            label=fr"$\beta={beta:g},\ \alpha={alpha:g}$")
ax.axhline(2, color=cf.RULE, linestyle=":", linewidth=1)
ax.set(xlabel="Time (model units)", ylabel="Protein concentration (model units)")
ax.legend()
plt.show()
```

## Worked example

- Use $\beta=2$, $\alpha=1$, $x_0=0$ and $h=1/4$. Calculate three Euler steps, compare them with the exact solution at the same times, and plot a full response.
- Change production and removal separately. Interpret which intervention changes the eventual concentration, the response time, or both.

## Why it works, and what it costs

- For a deviation $e=x-x^*$, derive $\dot e=-\alpha e$. The steady state attracts trajectories because deviations decay exponentially.
- Euler deviations satisfy $e_{k+1}=(1-\alpha h)e_k$. Stability requires $0<\alpha h<2$; a decaying, nonoscillating numerical response requires $0<\alpha h\le1$.
- Compare step sizes against the exact solution. A numerical artifact should not become a claim of biological oscillation.

## Where the model breaks

- Separate omitted regulation, molecular randomness, reporter calibration and numerical error. A deterministic model can be a useful abundance approximation even when individual reaction events are random.
- Constant production cannot explain two persistent expression states under the same conditions. Lesson 22 replaces constant production with a regulatory response; Lesson 23 adds process and measurement noise.

## Exercises

- Solve the linear model and check the simulation against it.
- Find two parameter pairs with the same steady state but different response times.
- Choose a time step that destabilizes Euler even though the continuous model is stable, and explain the discrepancy.

::: {#refs}
**References**
:::
