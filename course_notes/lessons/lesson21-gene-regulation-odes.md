---
published: false
title: "ODE Models of Gene Regulation"
subtitle: "Lesson 21 · 2027-04-12 · Week 14"
nocite: |
  @rosenfeld2002
---

## Core question

After a signal turns on production of a protein, how much of it accumulates in the cell, and how quickly?

```{.python .execute}
#| label: setup
#| include: false
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "../style")
import coursefigs as cf

cf.use_style()
```

## Watching a gene respond

The expression matrices of the previous lessons record each sample at one collection time. Many biological questions are about timing instead: when a bacterium meets a new sugar and switches on a gene, does the protein reach a useful level before conditions change again? Answering that requires measuring the same cell repeatedly.

A *fluorescent reporter* makes this possible. Biologists place the coding sequence of a fluorescent protein, such as green fluorescent protein (GFP), under the control of the promoter they want to study, so the cell makes the reporter whenever it would make the gene's own product. A microscope photographs growing cells every few minutes, and image software outlines each cell and records its mean fluorescence intensity over the cell's area. That reading tracks concentration rather than molecule number, so a division does not halve it; this is why the model below describes a concentration. One record is a cell identifier, an observation time $t_k$, and a measured intensity $y_k$ in arbitrary units. The records below are simulated to resemble such an experiment in *E. coli*. A chemical inducer switches on reporter production at time zero; images are taken every ten minutes, following one daughter after each division (@fig-regulation-reporter-setting).

```{.python .execute}
#| label: fig-regulation-reporter-setting
#| echo: false
#| fig-width: 5.2
#| fig-height: 2
#| fig-cap: "The reporter experiment. Left: inside one cell, the inducer switches on the promoter, the gfp gene is copied into mRNA and translated into fluorescent protein (red) at rate β, and protein is removed at rate αx by degradation and by dilution as the cell grows. Right: three images of one lineage, following the black-outlined daughter after each 30-minute division; a division halves each cell's protein but not its concentration. The model tracks one number per cell, the concentration x, set by the balance of production and removal."
#| fig-alt: "Left: a rounded rod-shaped cell containing a DNA strand with a grey promoter box and a gfp gene box; a black diamond inducer from outside the cell sits on the promoter. A blue arrow leads from the gene to a wavy mRNA, a second blue arrow labeled beta leads to a cluster of red dots labeled GFP, x, and an orange arrow labeled alpha x points down to the word removed. The equation x-dot equals beta minus alpha x is below the cell. Right: cartoon cells above a time axis at 0, 30, and 60 minutes; the first has no red dots, the second is two daughters with three dots each, one outlined in black and one in grey, and the third is a cell twice as long with nine dots."
fig, ax = plt.subplots()
cf.reporter_setting(ax)
plt.show()
```

```{.python .execute}
#| label: code-regulation-reporter-records
#| code-fold: true
#| code-summary: "Simulate reporter records for three cells"
np.random.seed(21)
doubling_time = 30.0                      # minutes per cell division
removal_rate = np.log(2) / doubling_time  # dilution by growth, per minute
background = 60.0                         # autofluorescence of an uninduced cell
plateaus = {"c1": 700.0, "c2": 1000.0, "c3": 1300.0}
times_min = np.arange(0, 181, 10)

rows = []
for cell, plateau in plateaus.items():
    level = plateau * (1 - np.exp(-removal_rate * times_min))
    readings = background + level + np.random.normal(0, 30, len(times_min))
    for t, y in zip(times_min, readings):
        rows.append({"cell": cell, "time_min": t, "fluorescence": round(y)})
records = pd.DataFrame(rows)
print(records.head(5).to_string(index=False))
```

```{.python .execute}
#| label: fig-regulation-reporter-data
#| code-fold: true
#| code-summary: "Plot the reporter time courses"
#| fig-cap: "Simulated reporter fluorescence for three E. coli cells after induction (points), with the curves used to simulate them (lines); open circles mark half of each rise. The cells settle at different levels, but each covers half its rise in one 30-minute doubling time (dotted line)."
#| fig-alt: "Three sets of points rise from about 60 units and level off near 760, 1060, and 1360 units by 180 minutes; open circles on a dotted vertical line at 30 minutes mark each trace halfway up its rise."
fig, ax = plt.subplots()
for (cell, plateau), marker in zip(plateaus.items(), ["o", "s", "^"]):
    cell_records = records[records["cell"] == cell]
    line = ax.plot(cell_records["time_min"], cell_records["fluorescence"],
                   marker, label=f"cell {cell}")[0]
    smooth_times = np.linspace(0, 180, 200)
    smooth = background + plateau * (1 - np.exp(-removal_rate * smooth_times))
    ax.plot(smooth_times, smooth, color=line.get_color(), linewidth=1)
    ax.plot(doubling_time, background + plateau / 2, "o", ms=9,
            color=line.get_color(), markerfacecolor="none",
            markeredgewidth=1.2)
ax.axvline(doubling_time, color=cf.RULE, linestyle=":", linewidth=1)
ax.set(xlabel="Time after induction (min)",
       ylabel="Reporter fluorescence (a.u.)")
ax.legend(loc="lower right")
plt.show()
```

Two features of each trace in @fig-regulation-reporter-data carry biological meaning. The *level* at which it settles says how much protein the cell maintains, and the *response time*, for instance the time to cover half the rise, says how quickly the cell adapts. The three cells settle at different levels, yet cover half their rise at about the same time. We want a model that says which cellular process controls each feature, and why measurements of an unregulated reporter in *E. coli* found a response time of about one cell generation [@rosenfeld2002].

The intensity $y_k$ is a measurement; the model describes the protein concentration $x(t)$, which changes continuously between images. For now, assume calibration has removed the background and converted intensity to concentration, and ignore measurement error and the time a new GFP molecule needs to become fluorescent. The forward problem is then the following.

::: {#def-induction-response}
## Response to induction
**Input:** a constant production rate, a removal rate constant, and a starting concentration $x_0$.

**Output:** the concentration $x(t)$ for all $t\geq0$, with its steady level and half-time as functions of the two rates.
:::

The biological inference runs in reverse. Once the forward problem is solved, a measured trace becomes a two-number readout: the removal rate from the half-time and the production rate from the plateau.

## A balance of production and removal

A protein's concentration changes for two reasons: ribosomes make new molecules from the gene's messenger RNA, and existing molecules are removed. Suppose the inducer keeps the promoter fully on, so production proceeds at a constant rate $\beta\geq0$, in concentration per unit time. Removal is different: twice as many molecules lose about twice as many per minute, because each faces the same chance of being degraded. We therefore model removal as proportional to the current concentration, with a rate constant $\alpha>0$ in units of inverse time.

Over a short interval $\Delta t$, production adds about $\beta\,\Delta t$ and removal subtracts about $\alpha x(t)\,\Delta t$:
$$
x(t+\Delta t)-x(t)\approx\bigl(\beta-\alpha x(t)\bigr)\Delta t.
$$
The approximation treats $x$ as nearly constant during the interval, which becomes accurate as $\Delta t$ shrinks. Dividing by $\Delta t$ and letting $\Delta t\to0$ turns the left side into a derivative.

::: {#def-production-removal-model}
## Production-removal model
For a production rate $\beta\geq0$, a removal rate constant $\alpha>0$, and an initial concentration $x_0\geq0$, the *production-removal model* is the ordinary differential equation (ODE) and initial condition
$$
\dot x(t)=\beta-\alpha x(t),\qquad x(0)=x_0 .
$$ {#eq-production-removal}
The dot denotes the derivative with respect to time, $\dot x = dx/dt$. The model is *linear* because its right side is a straight-line function of $x$; this is what makes a solution formula possible.
:::

The ODE does not hand us $x(t)$. It gives the slope of the unknown curve at each height, production minus removal, and solving it means finding the curve with those slopes that starts at $x_0$.

In a growing cell, $\alpha$ has a second contribution that is easy to miss. Concentration is a number of molecules $N(t)$ per cell volume $V(t)$, and a growing cell dilutes its contents even if no molecule is destroyed. Suppose molecules are degraded at rate $\gamma N$ and the volume grows exponentially, $V(t)=V_0e^{\mu t}$, so that $\dot V/V=\mu$ and the cell doubles in time $T_d=\ln2/\mu$. The cell makes $\beta V$ molecules per unit time, proportional to its volume, because its numbers of ribosomes and RNA polymerases grow with the cell. The quotient rule then gives
$$
\dot x=\frac{d}{dt}\frac{N}{V}=\frac{\dot N}{V}-\frac{N}{V}\cdot\frac{\dot V}{V}
=\frac{\beta V-\gamma N}{V}-\mu x=\beta-(\gamma+\mu)\,x .
$$
The result has the form of @eq-production-removal with $\alpha=\gamma+\mu$. Many bacterial proteins are stable, with $\gamma$ much smaller than $\mu$, so for them $\alpha\approx\mu=\ln2/T_d$: cell growth, not degradation, removes the protein.

## The steady state and the flow toward it

Before solving the ODE, we can ask where the concentration stops changing. If $x$ ever stays fixed, its derivative is zero, so production must exactly balance removal.

::: {#def-steady-concentration}
## Steady state
A *steady state* (here a steady concentration; also called an equilibrium) of $\dot x=\beta-\alpha x$ is a value $x^*$ with $\beta-\alpha x^*=0$. For the production-removal model it is unique:
$$
x^*=\frac{\beta}{\alpha}.
$$
A solution started at $x_0=x^*$ remains at $x^*$ for all time.
:::

Doubling production doubles the steady level, and so does halving removal. The three cells in @fig-regulation-reporter-data could differ in $\beta$, for instance because they carry different numbers of copies of the reporter gene, and so settle at different plateaus.

Whether a cell started elsewhere reaches $x^*$ depends on the sign of $\dot x$. In @fig-regulation-rate-balance, production is the horizontal line at $\beta=2$ and removal is the line $\alpha x$ with $\alpha=1$, crossing at $x^*=2$. Left of the crossing, production exceeds removal, so $\dot x>0$ and $x$ increases; to the right, removal wins and $x$ decreases. The arrows below the plot record this: every starting concentration flows toward $x^*$, so the equilibrium is drawn as a filled point.

```{.python .execute}
#| label: fig-regulation-rate-balance
#| code-fold: true
#| code-summary: "Draw production and removal against concentration"
#| fig-width: 5
#| fig-height: 3.4
#| fig-cap: "Production rate β = 2 and removal rate αx with α = 1 against the concentration x, crossing at the steady state x* = 2 (red). The arrows show x increasing where production exceeds removal and decreasing where removal wins."
#| fig-alt: "A horizontal blue line at height two and a rising orange line through the origin cross at x equals two, marked by a red point; a grey line beneath carries a right-pointing arrow left of two and a left-pointing arrow right of two."
beta, alpha = 2.0, 1.0
fig, ax = plt.subplots()
cf.rate_balance(
    ax,
    production=lambda x: beta + 0 * x,
    removal=lambda x: alpha * x,
    xlim=(0, 4),
    fixed_points=[(beta / alpha, True)],
    labels=(r"Production $\beta$", r"Removal $\alpha x$"),
    xlabel="Protein concentration, $x$ (model units)",
    ylabel="Rate (concentration per unit time)",
)
plt.show()
```

## Solving the linear model

The flow picture says where the concentration goes but not how fast; for that we solve the equation. The algebra is simplest in terms of the distance from the steady level. Define the *deviation* $e(t)=x(t)-x^*$. Because $x^*$ is constant, $\dot e=\dot x$, and replacing $\beta$ by $\alpha x^*$ gives
$$
\dot e=\beta-\alpha x=\alpha x^*-\alpha x=-\alpha e .
$$
Measured from the steady level, the deviation obeys pure decay: it shrinks at a rate proportional to its own size.

::: {#thm-linear-response}
## Solution of the production-removal model
The unique solution of @eq-production-removal is
$$
x(t)=x^*+(x_0-x^*)\,e^{-\alpha t},\qquad x^*=\frac{\beta}{\alpha}.
$$ {#eq-linear-response}
:::

::: {.proof}
Let $x$ be any solution and $e=x-x^*$, so $\dot e=-\alpha e$. Multiply the deviation by $e^{\alpha t}$ and differentiate with the product rule:
$$
\frac{d}{dt}\bigl(e^{\alpha t}e(t)\bigr)=e^{\alpha t}\bigl(\dot e(t)+\alpha e(t)\bigr)=0 .
$$
A function with zero derivative is constant, so
$$
e^{\alpha t}e(t)=e(0)=x_0-x^*\quad\Longrightarrow\quad e(t)=(x_0-x^*)\,e^{-\alpha t},
$$
which is @eq-linear-response. The argument applies to every solution, so there is only one.
:::

By @eq-linear-response, the concentration is the steady level plus whatever deviation remains. The initial deviation is multiplied by $e^{-\alpha t}$, which starts at one and decays to zero: after time $1/\alpha$, called the *relaxation time*, the factor is $e^{-1}\approx0.37$. The time at which half the deviation remains is easier to read off a trace. Setting $e^{-\alpha t_{1/2}}=1/2$ gives
$$
t_{1/2}=\frac{\ln 2}{\alpha}.
$$ {#eq-half-time}
For a gene switched on from zero, $x_0=0$ and the solution becomes $x(t)=x^*(1-e^{-\alpha t})$. At $t_{1/2}$ the concentration has covered half the distance to $x^*$, and after $3/\alpha$ it has covered $1-e^{-3}\approx95\%$.

::: {.callout-note title="Key idea"}
Production and removal together set the steady level $x^*=\beta/\alpha$, but removal alone sets the response time $t_{1/2}=\ln2/\alpha$.
:::

::: {#exm-regulation-level-speed}
## Changing production or removal
Start from $\beta=2$, $\alpha=1$, $x_0=0$, in model units, so $x(t)=2(1-e^{-t})$ and $t_{1/2}=\ln2\approx0.69$. Doubling production to $\beta=4$ gives $x(t)=4(1-e^{-t})$: the old response multiplied by two at every time, so the level doubles and the half-time stays at $0.69$.

Instead take $\beta=1$ and $\alpha=1/2$. The steady level is again $x^*=2$, but $t_{1/2}=2\ln2\approx1.39$. The response is twice as slow but ends at the same level, so a final measurement cannot tell the two cells apart; a time course can.
:::

```{.python .execute}
#| label: fig-regulation-response-times
#| code-fold: true
#| code-summary: "Plot the three responses"
#| fig-cap: "Exact responses from x0 = 0 for the worked example's three parameter pairs; points mark half-times. Doubling production doubles the level but keeps the half-time; halving removal keeps the level and doubles the half-time."
#| fig-alt: "Three rising curves: a solid curve approaching two and a dashed curve approaching four both pass their halfway points at time 0.69; a dash-dot curve approaching two passes its halfway point at time 1.39."
times = np.linspace(0, 6, 241)
fig, ax = plt.subplots()
for beta, alpha, style in [(2.0, 1.0, "-"), (4.0, 1.0, "--"), (1.0, 0.5, "-.")]:
    response = beta / alpha * (1 - np.exp(-alpha * times))
    line = ax.plot(times, response, style,
                   label=fr"$\beta={beta:g},\ \alpha={alpha:g}$")[0]
    half_time = np.log(2) / alpha
    ax.plot(half_time, beta / alpha / 2, "o", color=line.get_color())
ax.set(xlabel="Time (model units)",
       ylabel="Protein concentration, $x$ (model units)")
ax.legend(loc="lower right")
plt.show()
```

The dilution calculation turns @fig-regulation-response-times into a biological constraint. For a stable protein $\alpha\approx\ln2/T_d$, so @eq-half-time gives $t_{1/2}\approx T_d$: the cell needs about one generation to reach half its new level, however strongly the promoter is induced. That is the 30-minute half-time shared by all three cells in @fig-regulation-reporter-data.

A cell that must respond faster than it divides therefore cannot do so by making more protein. Active degradation raises $\alpha$, but holding the same level then requires proportionally more production, which costs energy. A transcription factor is a protein that binds DNA to regulate RNA production. One that represses its own promoter can instead start with a strong promoter and shut production down near its target; in *E. coli*, this negative autoregulation cut the time to reach half the final level to about one fifth of a cell cycle [@rosenfeld2002]. Its production depends on $x$, so it lies outside the constant-$\beta$ model; @exr-regulation-autoregulation works out an idealized version.

## Simulating the model: forward Euler

Regulated models, in which production depends on the protein concentration itself, generally have no solution formula. We then simulate the time course in small time steps, and before trusting a simulation on a model we cannot solve, we test it on one we can. The idea reverses the derivation of @eq-production-removal: over a short step of length $h$ the slope barely changes, so we use the slope at the start of the step for the whole step.

::: {#def-forward-euler}
## Forward Euler method
For $\dot x=f(x)$ with $x(0)=x_0$ and a step size $h>0$, the *forward Euler* approximation is the sequence
$$
x_{k+1}=x_k+h\,f(x_k),\qquad k=0,1,2,\ldots,
$$
where $x_k$ approximates $x(kh)$. For the production-removal model, $f(x)=\beta-\alpha x$.
:::

::: {.algorithm}
**Algorithm: forward Euler for production and removal.**
Given $\beta$, $\alpha$, $x_0$, a step $h$, and a number of steps $n$: set $x_0$. For $k=0,\ldots,n-1$, compute the net rate $\beta-\alpha x_k$ and set $x_{k+1}=x_k+h(\beta-\alpha x_k)$. Return $x_0,\ldots,x_n$, approximating $x(0),\ldots,x(nh)$.
:::

::: {#exm-regulation-euler-hand}
## Four Euler steps by hand
Take $\beta=2$, $\alpha=1$, $x_0=0$, and $h=1/2$. The first step uses the rate at zero: $x_1=0+\tfrac12(2-0)=1$. The second uses the rate at $x_1=1$, only $2-1=1$, so $x_2=1+\tfrac12\cdot1=1.5$. Continuing, $x_3=1.75$ and $x_4=1.875$ at $t=2$, against the exact $x(2)=2(1-e^{-2})\approx1.729$. Euler overestimates by about $0.15$ because each step uses the rate at its starting concentration, the largest rate in the step.
:::

```{.python .execute}
#| label: code-regulation-euler
def euler_trajectory(beta: float, alpha: float, x0: float,
                     h: float, n_steps: int) -> np.ndarray:
    """Forward Euler approximations x_0, ..., x_n for dx/dt = beta - alpha*x."""
    x = np.zeros(n_steps + 1)
    x[0] = x0
    for k in range(n_steps):
        net_rate = beta - alpha * x[k]  # production minus removal at x_k
        x[k + 1] = x[k] + h * net_rate
    return x


def exact_response(beta: float, alpha: float, x0: float,
                   t: np.ndarray) -> np.ndarray:
    """Exact solution x(t) of the production-removal model."""
    steady = beta / alpha
    return steady + (x0 - steady) * np.exp(-alpha * t)


exact_at_2 = exact_response(2.0, 1.0, 0.0, np.array([2.0]))[0]
for h in [0.5, 0.25]:
    approx = euler_trajectory(2.0, 1.0, 0.0, h, n_steps=round(2 / h))[-1]
    print(f"h = {h}: Euler x(2) = {approx:.3f}, exact {exact_at_2:.3f}, "
          f"error {approx - exact_at_2:+.3f}")
```

Halving the step to $h=1/4$ and taking eight steps gives $x_8\approx1.800$ at $t=2$, an error of about $0.07$: half the step, half the error. A numerical trajectory is an approximation whose error we can shrink, and rerunning with a smaller step is the basic check that a computed feature belongs to the model rather than to the method.

## When the time step is too large

A large step saves work, so how large can $h$ be? The deviation argument answers this for the Euler sequence too. For the production-removal model, forward Euler gives
$$
x_k=x^*+\lambda^k\,(x_0-x^*),\qquad \lambda=1-\alpha h.
$$ {#eq-euler-closed-form}
If $x_0\neq x^*$, the approximations converge to $x^*$ if and only if $0<\alpha h<2$, and they approach it without changing sides of $x^*$ if and only if $\alpha h\leq1$.

::: {.proof}
Sketch. Write $e_k=x_k-x^*$ and use $\beta=\alpha x^*$:
$$
e_{k+1}=x_k+h(\alpha x^*-\alpha x_k)-x^*=(1-\alpha h)\,e_k=\lambda e_k,
\qquad\text{so}\qquad e_k=\lambda^k e_0 .
$$
The exact solution multiplies the deviation by $e^{-\alpha h}$ over one interval of length $h$; Euler multiplies it by $\lambda$ instead. The powers of $\lambda$ behave as the conditions state:
$$
\lambda^k\to0\iff|1-\alpha h|<1\iff0<\alpha h<2,
\qquad
\lambda^k\ge0\ \text{for all }k\iff\lambda\ge0\iff\alpha h\le1 .
$$
:::

Euler's factor keeps only the first two terms of the Taylor series $e^{-\alpha h}=1-\alpha h+(\alpha h)^2/2-\cdots$, which is accurate when $\alpha h$ is small, but unlike the exact factor it can be negative when $\alpha h$ is not.

With $\beta=2$ and $\alpha=1$, the step $h=3/2$ gives $\lambda=-1/2$, and @eq-euler-closed-form yields $3,\ 1.5,\ 2.25,\ 1.875,\ \ldots$, overshooting $x^*=2$ and oscillating around it with shrinking amplitude. The step $h=3$ gives $\lambda=-2$ and the sequence $6,\ -6,\ 18,\ -30$. @fig-regulation-euler-steps compares these with the exact response.

```{.python .execute}
#| label: fig-regulation-euler-steps
#| code-fold: true
#| code-summary: "Compare Euler step sizes with the exact solution"
#| fig-width: 6
#| fig-height: 3.2
#| fig-cap: "Forward Euler for β = 2, α = 1, x0 = 0 (points) against the exact solution (black). (A) Step 1/2 tracks the curve; step 3/2 oscillates around x* = 2. (B) Step 3 gives negative, growing values. The oscillations come from the step size, not the model."
#| fig-alt: "Two panels. Left: a smooth black curve rising to two, with points close to it for step one-half and points alternating above and below two for step three-halves. Right: points at 0, 6, -6, 18, -30 alternating about a black curve near two, with a grey zero line."
fig, axes = plt.subplots(1, 2)
# Each panel lists (step h, number of steps, marker, series color index).
panels = {"A": (axes[0], [(0.5, 16, "o", 0), (1.5, 5, "s", 1)]),
          "B": (axes[1], [(3.0, 4, "^", 2)])}
for name, (ax, steps) in panels.items():
    end_time = max(h * n for h, n, _, _ in steps)
    curve_times = np.linspace(0, end_time, 200)
    ax.plot(curve_times, exact_response(2.0, 1.0, 0.0, curve_times),
            color=cf.INK, label="Exact")
    for h, n_steps, marker, color in steps:
        values = euler_trajectory(2.0, 1.0, 0.0, h, n_steps)
        ax.plot(h * np.arange(n_steps + 1), values, marker + "-",
                color=cf.SERIES[color], linewidth=0.8, label=f"Euler, $h={h:g}$")
    ax.set_title(name, loc="left", fontweight="bold")
    ax.set_xlabel("Time (model units)")
axes[0].legend(loc="lower right")
axes[1].legend(loc="lower left")
axes[1].axhline(0, color=cf.RULE, linewidth=0.8)
axes[0].set_ylabel("Protein concentration, $x$")
plt.show()
```

Unlike the negative concentrations in panel B, the oscillation in panel A looks plausible, because damped oscillations in protein levels do occur in circuits with delayed negative feedback. Here it is purely an artifact: every exact solution approaches $x^*$ monotonically, and reducing the step to $h\leq1/\alpha$ removes the oscillation, whereas a feature of the model would persist.

## What the linear model leaves out

The model holds production constant, while real promoters respond to transcription factors whose concentrations change, including the gene's own product. When production is a curved function of $x$, the production curve in @fig-regulation-rate-balance can cross the removal line more than once, and which steady state a cell reaches then depends on its history; Lesson 22 develops such feedback. The model is also deterministic and treats $y_k$ as equal to $x(t_k)$, yet a cell may hold only tens of copies of a regulatory protein, made in bursts, and real readings carry background and measurement error, as the scatter in @fig-regulation-reporter-data suggests. Lesson 23 separates fluctuations in the cell's concentration from errors in the recorded values.

## Exercises

::: {#exr-regulation-dilution}
An *E. coli* culture doubles every 40 minutes. Find $\alpha$ and the half-time for a stable reporter. If a degradation tag gives the reporter a 10-minute half-life even without growth, find the new $\alpha$ and half-time, and the factor by which production must rise to keep the same steady level.
:::

::: {#exr-regulation-autoregulation}
In an idealized self-repressing gene with $\alpha=1$, production is $\beta_{\max}=10$ while $x<K=2$; once $x$ reaches $K$, repression lowers production to $\alpha K$ and $x$ stays there. Starting from $x_0=0$, find the times to reach $K/2$ and $K$. An unregulated gene with $\beta=2$ has the same steady level but only approaches $K$; compare the two times to reach $K/2$. What does the self-repressing cell pay for its speed?
:::

::: {#exr-regulation-euler-steps}
With $\beta=2$, $\alpha=1$, and $x_0=0$, use @eq-euler-closed-form to describe the Euler sequences for $h=1$ and $h=2$. Why does $h=1$ reach $x^*$ in one step and stay there, and why does $h=2$ neither converge nor diverge?
:::

::: {#exr-regulation-reporter-fit}
A calibrated reporter trace rises from zero, levels off near 1000 units, and reaches 500 units about 30 minutes after induction. Estimate $\alpha$ and $\beta$, with units. What additional measurement would tell whether the removal comes from degradation or from dilution?
:::

::: {#refs}
**References**
:::
