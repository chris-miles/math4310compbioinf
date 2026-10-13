---
published: false
title: "Biological Switches and Cell States"
subtitle: "Lesson 22 · 2027-04-14 · Week 14"
nocite: |
  @mangelGeneNetworks, @bhattacharya2011
---

## Core question

How can gene regulation support two persistent cell states under the same external conditions?

::: {.callout-note title="Draft status"}
This lesson is a scaffold. It extends Lesson 21's production-removal model to a regulatory circuit; Lesson 23 asks how noisy measurements and biological fluctuations affect its interpretation.
:::

## Can an expression state persist after a signal is removed?

Lesson 21's constant production-removal model approaches one concentration, regardless of its starting value. Suppose instead that a temporary intervention leaves a cell in a different expression pattern even after the intervention ends. To explain that persistence, we need a mechanism that can maintain different abundances under the same continuing conditions. A *cell state* here means a persistent pattern of gene-product abundances; it is a dynamical claim that a cluster in an expression plot cannot establish by itself.

Consider a synthetic *toggle circuit*: two genes produce proteins that repress one another. Repression means that a protein reduces the other gene's production rate. If protein 1 is abundant, it suppresses protein 2; reduced protein 2 then allows production of protein 1 to continue. The reverse pattern may also sustain itself [@mangelGeneNetworks]. We will ask when this feedback actually supports two persistent states, rather than assume that drawing two inhibitory arrows guarantees a switch.

## Two concentrations and their time courses

The observations would be reporter time courses for both gene products, with initial preparation and intervention times recorded. The modeled state is the pair of nonnegative concentrations $(u(t),v(t))$, rather than the measured light intensities themselves. For the toy circuit, scale each concentration by its repression threshold and time by a common removal time, assuming symmetric thresholds and rates. The resulting $u,v,t$ are dimensionless. We seek the steady states and ask which states nearby trajectories approach, so the output includes both persistent abundance pairs and their stability.

## The model

- Introduce the repression response $H(z)=a/(1+z^n)$, with maximum production $a>0$ and response steepness $n\ge1$. Its decreasing shape captures reduced production as the repressor increases; $n$ is an effective response parameter, not automatically a literal molecular count.
- Retain linear removal and use the symmetric toy system
  $$
  \dot u=H(v)-u,\qquad \dot v=H(u)-v.
  $$
  Explain mutual repression as feedback and distinguish it from a steep response of a single gene to an externally fixed input [@mangelGeneNetworks].

For the worked example $a=3,n=2$, the state portrait below makes persistence visible. Each plotted point represents both protein concentrations in one cell, rather than a different cell in an expression dataset. Two trajectories evolve under the same external conditions but approach different abundance pairs. At either endpoint, production balances removal for both proteins. The blue curve marks $\dot u=0$ and the orange curve marks $\dot v=0$, so their intersections locate these balances. The next section develops this construction and checks which intersections attract nearby states.

::: {.algorithm}
**Algorithm: simulate the two-protein balance by Euler stepping.**

Starting from a concentration pair, calculate both production-removal rates. Add $h$ times those rates to the pair, record the new pair, and repeat. Use the same small step $h=0.01$ for both concentrations; this is Lesson 21's Euler update applied to two coupled variables.
:::

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
#| label: fig-switch-persistent-states
#| code-fold: true
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "State portrait for the symmetric toggle at a = 3 and n = 2; black trajectories start at (2, 0.4) and (0.4, 2), colored curves mark zero rates, and red filled/open circles mark stable/unstable equilibria. The same regulatory conditions support two persistent abundance patterns."
#| fig-alt: "A square concentration plot contains two crossing zero-rate curves, filled points near (2.62, 0.38) and (0.38, 2.62), and an open point near (1.21, 1.21); two black trajectories approach opposite filled points."
def u_rate(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Return protein 1's production minus removal in the toy toggle."""
    return 3 / (1 + v ** 2) - u


def v_rate(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Return protein 2's production minus removal in the toy toggle."""
    return 3 / (1 + u ** 2) - v


trajectories = []
for initial in [(2.0, 0.4), (0.4, 2.0)]:
    trajectory = np.zeros((1001, 2))
    trajectory[0] = initial
    for k in range(1000):
        u, v = trajectory[k]
        rates = np.array([u_rate(u, v), v_rate(u, v)])
        trajectory[k + 1] = trajectory[k] + 0.01 * rates
    trajectories.append(trajectory)
high = (3 + np.sqrt(5)) / 2
low = (3 - np.sqrt(5)) / 2
symmetric = 1.21341166276223  # Positive root of s + s^3 = 3.
fig, ax = plt.subplots()
cf.phase_plane(
    ax, u_rate, v_rate, xlim=(0, 3), ylim=(0, 3),
    fixed_points=[(high, low, True), (low, high, True),
                  (symmetric, symmetric, False)],
    trajectories=trajectories,
    labels=("Protein 1 concentration, u (scaled)",
            "Protein 2 concentration, v (scaled)"),
)
plt.show()
```

## The technique

- Plot the nullclines $u=H(v)$ and $v=H(u)$, where the corresponding derivative vanishes. Their intersections are steady states.
- Use the signs of both derivatives to indicate movement between the nullclines. Introduce local stability through small deviations and the matrix
  $$J=\begin{pmatrix}-1&H'(v)\\H'(u)&-1\end{pmatrix}.$$
  Connect the eigenvalue signs to growth or decay of small perturbations.
- Compare trajectories started on different sides of a state boundary, and ask whether a transient input leaves a persistent change after removal.

## Worked example

- Set $n=2$. For $a=1$, verify the single symmetric steady state; for $a=3$, find three steady states.
- For the asymmetric states at $a=3$, subtract the equilibrium equations to obtain $(u-v)(1-uv)=0$. Hence $uv=1$ and $u+v=3$, giving $(u,v)=((3+\sqrt5)/2,(3-\sqrt5)/2)$ and its reversal.
- Locate the symmetric state from $u=v=s$ and $s+s^3=3$. Compute the local stability matrix at all three states, then compare numerical trajectories from different initial conditions.

## Why it works, and what it costs

- For this symmetric $n=2$ model, the symmetric equilibrium loses stability at $a=2$, where $s=1$. For $a>2$, the two asymmetric equilibria are stable and the symmetric equilibrium is a saddle.
- Explain bistability as two attracting steady states, each with its own set of initial conditions that approach it. Cooperativity alone does not guarantee bistability; the circuit and parameter regime matter.
- Keep claims about this toy circuit distinct from a claim that every differentiated cell type is explained by two mutually repressing genes.

## Where the model breaks

Waddington's landscape pictures persistent cell states as valleys and transitions as movement between them [@bhattacharya2011]. This metaphor can help interpret attracting states and their attraction regions, but drawing valleys does not establish the trajectories of a regulatory circuit. Any quantitative landscape must be justified from the specified dynamics.

- Symmetric rates and the Hill response simplify regulatory chemistry. Report what changes when symmetry is relaxed rather than interpreting the toy parameters as measured constants.
- A deterministic trajectory remains in its attraction region. Biological fluctuations can change the trajectory; measurement noise can change only its apparent position. Lesson 23 separates these possibilities.

## Exercises

- Verify the three equilibria and their stability for $a=3,n=2$.
- Show that the symmetric $n=1$ system has only one nonnegative equilibrium and compare its response with the $n=2$ example.
- Describe an intervention and time-course observation that could distinguish state persistence from a transient response.

::: {#refs}
**References**
:::
