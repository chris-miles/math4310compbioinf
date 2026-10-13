---
published: false
title: "Biological Dynamics and Noisy Measurements"
subtitle: "Lesson 23 · 2027-04-19 · Week 15"
nocite: |
  @elowitz2002
---

## Core question

Does a noisy gene-expression time course reflect a fluctuating biological state, an imperfect measurement, or both?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. It completes the gene-regulation sequence with a simulation comparison and a state-space formulation. It introduces no filtering algorithm.
:::

## A fluctuation in the cell or in the reading?

Suppose we repeatedly image the fluorescent reporter from Lesson 21 in the same cell. Each record contains a cell identifier, an observation time, and a light intensity. After the response has settled, one reading rises above the expected level. Did the cell briefly make more protein, or did the measurement overestimate an unchanged abundance? The distinction matters for Lesson 22's switches: an apparent excursion in a reporter trace does not by itself establish that the cell entered a different persistent state.

A change in the actual abundance enters the cell's subsequent production-removal dynamics. An error in one camera reading changes the recorded value without changing the cell. We will compare these two possibilities under controlled assumptions, using a known simulated trajectory to separate what the cell does from what the instrument reports. Biological fluctuations are themselves possible even among genetically identical cells [@elowitz2002]; the observation model determines how faithfully a time course reveals them.

## A hidden abundance and recorded measurements

Let $x_k$ be the actual, unobserved abundance at observation time $t_k$, and $y_k$ the recorded reporter value. Use calibrated reporter units for this toy comparison. The input is the ordered series $(t_k,y_k)$ together with a model for abundance changes and measurement error; the question is which features of that series distinguish a disturbance to $x_k$ from a disturbance only to $y_k$. Tracking the same cell preserves its temporal history. Measuring a different cell at every time can describe a population response, but it does not directly show how a fluctuation persists within one cell.

## The deterministic prediction

- Reuse $\dot x=\beta-\alpha x$ with $\beta=2$, $\alpha=1$ and $x_0=0$ for the response from zero. Its steady state is $x^*=2$. For the noise comparison below, start a separate experiment at $x_0=x^*=2$ so the deterministic response has already settled.
- For equally spaced observation times with interval $h=1/4$, use the exact deterministic update
  $$
  F_h(x)=x^*+(x-x^*)e^{-\alpha h}.
  $$
  Using this update keeps Euler error from becoming a third source of variability in the comparison.

## The state and the observation

- Put process noise in the state update and measurement noise in the observation:
  $$
  x_{k+1}=F_h(x_k)+w_k,\qquad y_k=x_k+v_k.
  $$
  The independent disturbances $w_k$ and $v_k$ have mean zero and variances $q$ and $r$, respectively. Explain variance as expected squared deviation. Here $q$ describes accumulated disturbance over the fixed interval $h$, and $r$ describes uncertainty in one reading.
- This pair of equations is a state-space model: a hidden biological state evolves, and the instrument reports a noisy function of it. The earlier HMM stepped along DNA positions; here the index steps through observation times. Abundance is now continuous and has a mechanistic update, while the distinction between hidden state and observation remains.
- For the repeated toy simulations, draw independent disturbances uniformly from a small interval such as $[-0.05,0.05]$, centered at zero. This additive approximation does not provide a molecule-count reaction model; larger disturbances can produce negative modeled abundances. Changing $h$ requires reconsidering the process disturbance distribution.

## A controlled simulation comparison

- Begin near the steady state $x^*=2$. In the process model, suppose one disturbance adds $w_k=0.4$ to the state and subsequent disturbances are zero. The next state is $2.4$; after another interval $h=1/4$, it is $2+0.4e^{-1/4}\approx2.312$, followed by $2+0.4e^{-1/2}\approx2.243$. The perturbation decays through the same biological dynamics as in Lesson 21.
- In the measurement model, keep the actual state at $2$ and add one observation error $v_k=0.4$. The recorded value is $2.4$ once, then returns to $2$ when subsequent observation errors are zero. The instrument's blip does not enter the next state update.
- Plot actual abundance and recorded measurements separately for these two cases. Then simulate repeated time courses with reproducible small disturbances, using the same sampling interval and stated disturbance distributions. Equal disturbance scales need not produce equal variation in the recorded values, because process disturbances accumulate through the dynamics. Compare how long deviations persist rather than judging biology from how jagged one curve looks.
- Add both disturbance types in a third comparison. Ask which features of the hidden trajectory remain visible and which apparent excursions arise entirely in the observation.
- Return to Lesson 22's toggle: perturbing the actual state across an attraction-region boundary can lead to a lasting change; placing only an observed point across that boundary leaves the modeled biological trajectory unchanged. An apparent crossing in a reporter trace therefore needs follow-up measurements.

Both cases below produce the same elevated reading at the disturbance time. Their later readings differ because only the process disturbance enters the next state update. This comparison makes the distinction visible under the stated single-disturbance assumptions; a correlated instrument error would require a different comparison.

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
#| label: fig-noisy-dynamics-disturbance
#| code-fold: true
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "A single 0.4-unit disturbance to abundance or to a recorded value, starting at the steady state of two. Only the abundance disturbance persists through subsequent biological dynamics."
#| fig-alt: "Two time-course panels each show a reading of 2.4 at time one; the process panel then decays toward two, while the measurement panel has a constant true abundance of two and one isolated elevated reading."
times = np.arange(0, 3.01, 0.25)
process_state = np.full(len(times), 2.0)
after_disturbance = times >= 1
process_state[after_disturbance] += 0.4 * np.exp(-(times[after_disturbance] - 1))
measurement_state = np.full(len(times), 2.0)
measurement_readings = measurement_state.copy()
measurement_readings[4] += 0.4
fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
for ax, state, readings, title in [
    (axes[0], process_state, process_state, "A  Process disturbance"),
    (axes[1], measurement_state, measurement_readings, "B  Measurement error"),
]:
    # Repeating the disturbance time draws its immediate change in abundance.
    state_times = np.insert(times, 4, 1.0)
    state_values = np.insert(state, 4, 2.0)
    ax.plot(state_times, state_values, "-", color=cf.SERIES[0],
            label="Actual abundance")
    ax.plot(times, readings, "s", color=cf.SERIES[1], label="Recorded values")
    ax.set(xlabel="Time (model units)", title=title, ylim=(1.96, 2.45))
axes[0].set_ylabel("Calibrated reporter units")
axes[1].legend(loc="upper right")
plt.show()
```

## Biological noise and what the comparison establishes

- Introduce intrinsic noise as randomness in the gene-expression reaction events themselves and extrinsic noise as variation in other cellular components or conditions affecting expression. Two identically regulated reporters can help separate their contributions [@elowitz2002].
- Both can affect the biological process. Intrinsic versus extrinsic is a different distinction from process versus measurement; the additive process term does not identify a molecular cause.
- The simulation distinguishes two specified models. A correlated instrument error could also produce temporal dependence, and one short series cannot uniquely assign every fluctuation to a source. Use calibration measurements and biological replication to strengthen the comparison.

## Where the model breaks

- Extend the observation equation to $y_k=cx_k+b+v_k$ when a reporter has a gain $c$ and background $b$. Explain why unknown calibration can confound an abundance interpretation.
- Discuss sparse sampling, cell division, and omitted changes in production as reasons to refine the state model.

## Exercises

- Calculate the next three states after a single disturbance of size $0.4$, and compare them with three observations following a single measurement error of the same size.
- Simulate the process-only, measurement-only, and combined models. Identify a temporal pattern that helps distinguish the first two under the stated independent-noise assumptions.
- Give an experimental control that estimates measurement variation and explain what it cannot establish about intrinsic versus extrinsic biological noise.
- Draw a state-observation diagram for the two-gene circuit and label where process and measurement disturbances enter.

## Optional extension: matching the amount of variation

- To compare persistence while holding marginal variation fixed, write $a=e^{-\alpha h}$ and $z_k=x_k-x^*$. Then $z_{k+1}=az_k+w_k$. In a stationary regime, independent mean-zero disturbances with variance $q$ give $P=a^2P+q$, hence $P=q/(1-a^2)$ for the state variance $P$.
- Choose a target observation variance $s^2$. A measurement-only model with $r=s^2$ and a process-only model with $q=(1-a^2)s^2$ have the same marginal variance. Their lag-one covariances are respectively $0$ and $as^2$ under the independence assumptions.
- For an approximately stationary simulation comparison, discard the first 30 steps. At $h=1/4$ and $\alpha=1$, the remaining initial variance contribution is multiplied by $a^{60}=e^{-15}$.
- Optional exercise: derive these formulas, then compare repeated simulations. Correlated instrument errors would require a different observation model.

## Bridge to the final lesson

A model's usefulness depends on the inference its data can support. Lesson 24 returns to newer technologies and representations with this distinction in mind: additional information, a more flexible model, and independent validation address different parts of a biological claim.

::: {#refs}
**References**
:::
