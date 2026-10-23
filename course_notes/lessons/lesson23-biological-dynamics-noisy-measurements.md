---
published: false
title: "Biological Dynamics and Noisy Measurements"
subtitle: "Lesson 23 · 2027-04-19 · Week 15"
nocite: |
  @elowitz2002
---

## Core question

Does a noisy gene-expression time course reflect a fluctuating biological state, an imperfect measurement, or both?

```{.python .execute}
#| label: setup
#| code-fold: true
#| code-summary: "Setup: imports and figure style"
import sys

sys.path.insert(0, "../style")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import coursefigs as cf

cf.use_style()
```

## Reporter time courses from single cells

Lesson 21 followed a fluorescent reporter after a signal switched on its production. A time-lapse microscope images growing bacteria every few minutes, and software outlines each cell and records its mean fluorescence, which tracks the reporter's concentration. After background subtraction and calibration, each record holds a cell identifier, an image time $t_k$, and a reading $y_k$. The chunk below simulates four cells and prints the first records for cell c1; the model behind it is the subject of this lesson.

```{.python .execute}
#| label: code-noisy-dynamics-records
#| code-fold: true
#| code-summary: "Simulate four reporter time courses and list the first records"
def simulate_reporter(
    beta: float, alpha: float, h: float, n_steps: int,
    x0: float, q: float, r: float, c: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate hidden concentrations x_k and recorded readings y_k."""
    x = np.zeros(n_steps + 1)
    x[0] = x0
    for k in range(n_steps):
        drift = h * (beta - alpha * x[k])  # deterministic Euler step
        disturbance = np.random.normal(0, np.sqrt(q))  # w_k, enters the cell
        x[k + 1] = x[k] + drift + disturbance
    errors = np.random.normal(0, np.sqrt(r), size=n_steps + 1)  # v_k
    y = c * x + errors
    return x, y


np.random.seed(23)
h = 0.25
times = h * np.arange(33)
cells = {}
for name in ["c1", "c2", "c3", "c4"]:
    cells[name] = simulate_reporter(2.0, 1.0, h, 32, x0=0.0, q=0.00875, r=0.02)
records = pd.DataFrame({"cell": "c1", "time": times, "reading": cells["c1"][1]})
print(records.head(5).round(3).to_string(index=False))
```

In @fig-noisy-dynamics-traces, each cell rises toward the steady concentration $x^*=\beta/\alpha=2$, but none follows the smooth curve. Cell c2 stays above the curve for several images near $t=5$, while cell c1 has isolated readings, such as the one at $t=3.5$, that jump away and come straight back. Did c2 briefly make more reporter, or did the camera overestimate it? If cells fluctuate, the spread is a property of gene expression; if the instrument does, it only measures our precision. In a bistable circuit the difference is sharper still: a real fluctuation can push a cell into the other persistent state, while a misread value cannot.

```{.python .execute}
#| label: fig-noisy-dynamics-traces
#| code-fold: true
#| code-summary: "Plot the four reporter time courses"
#| fig-width: 6
#| fig-height: 3.6
#| fig-cap: "Simulated readings for four cells after production switches on, with the noise-free Euler response (dashed). The plot alone does not say whether the cells or the camera produce the scatter."
#| fig-alt: "Four panels, cells c1 to c4; black points rise from zero to about two, scattered around a grey dashed curve."
deterministic = 2.0 * (1 - (1 - h * 1.0) ** np.arange(len(times)))  # Euler steps with no noise
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)
for ax, (name, (x, y)) in zip(axes.flat, cells.items()):
    ax.plot(times, deterministic, "--", color=cf.RULE, linewidth=1)
    ax.plot(times, y, "o-", color=cf.INK, markersize=3, linewidth=0.6)
    ax.set_title(f"Cell {name}", loc="left")
fig.supxlabel("Time (model units)", fontsize=10)
fig.supylabel("Reading $y_k$", fontsize=10)
plt.show()
```

Given an ordered series $(t_k,y_k)$, $k=0,\ldots,n$, with spacing $h$, we want to say which kind of variation it contains, how much of each, and what one trace cannot decide. We keep the rates $\beta=2$, $\alpha=1$ and image every $h=1/4$ time units, so one relaxation time $1/\alpha$ spans four images.

## Two places for noise to enter

The deterministic Euler step $x_{k+1}=x_k+h(\beta-\alpha x_k)$ approximates the concentration at successive image times. In a real cell, production and removal happen through discrete molecular events, such as a promoter firing or a protein being degraded, and over one interval they need not add up to exactly $h(\beta-\alpha x_k)$. We collect the difference into a random *process disturbance* $w_k$, which changes the cell itself. Separately, photon counting and errors in outlining the cell make the reading differ from the concentration. That difference is the *measurement error* $v_k$, which changes only the record.

To say how large each kind of noise is, we use the *variance* of a random quantity $U$: its expected squared deviation from its mean, $\operatorname{Var}(U)=E[(U-E[U])^2]$. For a mean-zero disturbance it is simply $E[w_k^2]$.

::: {#def-state-space-model}
## Linear state-space model
The hidden concentration $x_k$ and the reading $y_k$ at time $t_k$ satisfy
$$
x_{k+1}=x_k+h(\beta-\alpha x_k)+w_k,\qquad y_k=c\,x_k+v_k,
$$
where $c>0$ is the reporter gain. The disturbances $w_0,w_1,\ldots$ and errors $v_0,v_1,\ldots$ all have mean zero; each $w_k$ has variance $q$ and each $v_k$ has variance $r$; and all of them are independent of one another and of the initial concentration $x_0$.
:::

So $q$ measures how far one interval's molecular events typically push the concentration off the deterministic step, and $r$ how far one reading typically lands from $c\,x_k$. We use calibrated units, $c=1$, until the final section.

The model has the shape of the hidden Markov models of Lessons 10 to 12. Its dependencies are
$$
\begin{array}{ccccccc}
x_0 & \xrightarrow{\;+\,w_0\;} & x_1 & \xrightarrow{\;+\,w_1\;} & x_2 & \xrightarrow{\;+\,w_2\;} & \cdots\\[2pt]
\big\downarrow{\scriptstyle +\,v_0} & & \big\downarrow{\scriptstyle +\,v_1} & & \big\downarrow{\scriptstyle +\,v_2} & & \\[2pt]
y_0 & & y_1 & & y_2 & &
\end{array}
$$
The horizontal arrows are the HMM's transitions: given $x_k$, the next concentration does not depend on earlier ones. The vertical arrows are its emissions: given $x_k$, the reading $y_k$ is independent of every other state and reading. The B/I model stepped along DNA through two discrete states; here the index is time, the hidden state is a continuous concentration, and the transition comes from the ODE. Estimating the hidden concentrations from the readings has continuous analogues of the HMM computations for Gaussian noise: the Kalman filter uses readings up to the current time, like the forward pass, and the Kalman smoother uses the whole series, like forward-backward. We do not develop them.

## One disturbance, two fates

The cleanest comparison applies one disturbance of each kind to a cell at steady state. Write $e_k=x_k-x^*$ for the deviation from $x^*=\beta/\alpha$. Because $\beta=\alpha x^*$, the drift is $h(\beta-\alpha x_k)=-\alpha h(x_k-x^*)$, and subtracting $x^*$ from both sides of the state update gives
$$
e_{k+1}=\lambda\,e_k+w_k,\qquad \lambda=1-\alpha h.
$$ {#eq-noisy-deviation}
First the existing deviation shrinks by the factor $\lambda$, then the new disturbance is added. With $\alpha=1$ and $h=1/4$, $\lambda=3/4$.

::: {#exm-single-disturbance}
## A 0.4-unit disturbance to the cell or to the reading
Start at $x_0=2$ with all disturbances and errors zero except one.

*Process.* Let $w_3=0.4$, so a burst of extra production raises the concentration at $t_4=1$. Then $e_4=0.4$, and @eq-noisy-deviation gives $e_5=0.75(0.4)=0.3$, $e_6=0.75(0.3)=0.225$, and $e_7\approx0.169$. The readings are $2.4,\ 2.3,\ 2.225,\ 2.169,\ldots$, decaying back toward 2.

*Measurement.* Let $v_4=0.4$ instead. The concentration stays at 2 throughout, and the readings are $2,\ 2.4,\ 2,\ 2,\ldots$ around $t_4$.
:::

At $t_4$ the two records are identical. They differ afterwards because only $w_k$ enters the state update (@fig-noisy-dynamics-disturbance). The cell forgets a process disturbance at the rate it forgets any perturbation: a fraction $\lambda$ survives each interval, so after one relaxation time (four images) about $0.75^4\approx0.32$ remains, close to $e^{-1}\approx0.37$. The instrument's error is gone at the next image. In a toggle switch, a process disturbance that crosses the boundary between attraction regions changes the state the cell settles into; an equal measurement error leaves one displaced point.

```{.python .execute}
#| label: fig-noisy-dynamics-disturbance
#| code-fold: true
#| code-summary: "Plot the single-disturbance comparison"
#| fig-width: 6
#| fig-height: 2.8
#| fig-cap: "Hidden concentration (line) and readings (points) after a 0.4-unit disturbance at time 1, applied to the cell (A) or to one reading (B). Both records show 2.4 at time 1, but only the process disturbance persists, shrinking by a factor 3/4 per image."
#| fig-alt: "A: line and blue points jump from 2 to 2.4 at time 1 and decay toward 2. B: flat line at 2, one orange point at 2.4."
k = np.arange(13)
process_x = np.where(k >= 4, 2 + 0.4 * 0.75 ** (k - 4), 2.0)  # e_k from the example
measured_y = np.where(k == 4, 2.4, 2.0)
fig, axes = plt.subplots(1, 2, sharey=True)
for ax, x, y, color, title in [
        (axes[0], process_x, process_x, cf.SERIES[0], "A  Process disturbance $w_3=0.4$"),
        (axes[1], np.full(13, 2.0), measured_y, cf.SERIES[1], "B  Measurement error $v_4=0.4$")]:
    ax.plot(h * k, x, color=cf.INK, linewidth=1)
    ax.plot(h * k, y, "o", color=color)
    ax.set_title(title, loc="left")
    ax.set_xlabel("Time (model units)")
axes[0].set_ylabel("Concentration or reading")
plt.show()
```

In a real cell a disturbance arrives at every interval. Unrolling @eq-noisy-deviation from $e_0$ gives
$$
e_k=\lambda^k e_0+\sum_{j=0}^{k-1}\lambda^{\,k-1-j}\,w_j,
$$
so the current deviation is a sum of all past disturbances, each shrunk according to its age. As new disturbances arrive and old ones decay, the spread of $x_k$ settles where the two balance.

## How much the concentration fluctuates

::: {#prp-stationary-variance}
## Stationary variance of the concentration
Under @def-state-space-model with $0<\alpha h<2$, so that $|\lambda|<1$, the variance of $x_k$ converges as $k\to\infty$ to
$$
\sigma_x^2=\frac{q}{1-\lambda^2}.
$$
If $x_0$ already has mean $x^*$ and variance $\sigma_x^2$, every $x_k$ does.
:::

::: {.proof}
Sketch. Two facts about averages do the work. First, averages add: the expectation of a sum is the sum of the expectations. Second, for independent quantities the average of a product is the product of the averages. The deviation $e_k$ is built from $e_0$ and $w_0,\ldots,w_{k-1}$, so knowing $e_k$ tells us nothing about the new disturbance $w_k$, and $E[e_kw_k]=E[e_k]E[w_k]=0$. Squaring @eq-noisy-deviation and averaging term by term,
$$
E[e_{k+1}^2]=\lambda^2E[e_k^2]+2\lambda\,E[e_kw_k]+E[w_k^2]=\lambda^2E[e_k^2]+q.
$$
The fixed point of $V\mapsto\lambda^2V+q$ is $V=q/(1-\lambda^2)$, and the distance from it shrinks by $\lambda^2$ per step. The mean deviation also shrinks, $E[e_{k+1}]=\lambda E[e_k]$, so $E[e_k^2]$ and the variance share this limit.
:::

For our parameters $1-\lambda^2=7/16$, so $\sigma_x^2=16q/7\approx2.29\,q$: several intervals' disturbances are present at any moment. Measurement error does not accumulate; with $q=0$ the readings have variance $r$. With both sources, $y_k-x^*=e_k+v_k$, and the cross term vanishes as in the proof, so at stationarity
$$
\operatorname{Var}(y_k)=\sigma_x^2+r.
$$
The formula also shows the removal rate controlling fluctuations. Disturbances in disjoint sub-intervals of an interval come from separate molecular events and are independent, so their variances add and $q$ is proportional to the interval's length, $q=Dh$. For short intervals $1-\lambda^2\approx2\alpha h$, so $\sigma_x^2\approx D/(2\alpha)$ regardless of how often we image. Faster removal forgets disturbances sooner and keeps the concentration closer to $x^*$.

## Telling the two sources apart

@prp-stationary-variance lets us make the two sources equally large. With $r=0.04$ in a measurement-only cell and $q=(7/16)(0.04)=0.0175$ in a process-only cell, both sets of readings have variance $0.04$ around $x^*=2$ (@fig-noisy-dynamics-three-ways, panels A and B). Their histograms would match, but process-only readings stay above or below 2 for several images in a row, while measurement-only readings jump independently. To measure this, we ask how strongly a reading predicts one taken $\ell$ images later.

The *covariance* of two random quantities, $\operatorname{Cov}(U,V)=E[(U-E[U])(V-E[V])]$, is positive when they tend to deviate from their means in the same direction, and it is zero when they are independent. For a series whose mean and variance do not change with $k$, the *autocorrelation at lag* $\ell$ is
$$
\rho(\ell)=\frac{\operatorname{Cov}(y_k,\,y_{k+\ell})}{\operatorname{Var}(y_k)}.
$$
It equals 1 at lag zero. From one series with average $\bar y$, we estimate it by replacing expectations with averages over the series:
$$
\hat\rho(\ell)=\frac{\sum_{k=0}^{n-\ell}(y_k-\bar y)(y_{k+\ell}-\bar y)}{\sum_{k=0}^{n}(y_k-\bar y)^2}.
$$

::: {#prp-noise-autocorrelation}
## Autocorrelation of a noisy reporter
Under @def-state-space-model with $c=1$ and the concentration at stationarity, for every lag $\ell\ge1$
$$
\rho(\ell)=\frac{\sigma_x^2}{\sigma_x^2+r}\,\lambda^{\ell}.
$$
:::

::: {.proof}
Sketch. Unroll @eq-noisy-deviation for $\ell$ steps: $e_{k+\ell}=\lambda^{\ell}e_k+(\text{terms in }w_k,\ldots,w_{k+\ell-1})$. The covariance of two sums is the sum of the covariances of all pairs of terms, because the product of the two deviations multiplies out term by term. In $\operatorname{Cov}(e_k+v_k,\ e_{k+\ell}+v_{k+\ell})$, two kinds of pairs are independent and contribute zero: $e_k$ with a later disturbance $w_k,\ldots,w_{k+\ell-1}$, and any pair containing a measurement error $v$. The one surviving pair, $e_k$ with $\lambda^\ell e_k$, contributes $\lambda^\ell E[e_k^2]=\lambda^\ell\sigma_x^2$. Dividing by $\operatorname{Var}(y_k)=\sigma_x^2+r$ gives the result.
:::

With measurement error only, $\rho(\ell)=0$ at every positive lag. With process noise only, $\rho(\ell)=\lambda^\ell\approx e^{-\alpha\,\ell h}$: the correlation fades at the removal rate $\alpha$, the clock that sets the deterministic response time. With both, the curve drops sharply from lag zero to lag one, then decays geometrically, and extending the geometric part back to lag zero estimates the share of reading variance that belongs to the cell. With $\sigma_x^2=r=0.02$, for example, $\rho(1)=0.5\times0.75=0.375$ and $\rho(2)\approx0.281$.

::: {.algorithm}
**Algorithm: estimate autocorrelations from one series.**
Subtract the series average from every reading. For each lag $\ell=0,1,\ldots,L$, multiply each deviation by the deviation $\ell$ positions later, sum the products, and divide by the sum of squared deviations. Return the $L+1$ ratios.
:::

The function below applies $\hat\rho(\ell)$ to 2,000 simulated readings in each case, all with reading variance $0.04$; the combined case also produced @fig-noisy-dynamics-traces. Each row prints estimates at lags 1 to 4, then @prp-noise-autocorrelation's values.

```{.python .execute}
#| label: code-noisy-dynamics-autocorrelation
def autocorrelation(y: np.ndarray, max_lag: int) -> np.ndarray:
    """Estimate rho(lag) for lags 0..max_lag from one series."""
    deviations = y - y.mean()
    denominator = np.sum(deviations ** 2)
    values = np.zeros(max_lag + 1)
    for lag in range(max_lag + 1):
        # pair each deviation with the one `lag` images later
        products = deviations[:len(y) - lag] * deviations[lag:]
        values[lag] = np.sum(products) / denominator
    return values


np.random.seed(4)
scenarios = {"measurement only": (0.0, 0.04),
             "process only": (0.0175, 0.0),
             "both": (0.00875, 0.02)}
lam = 1 - 1.0 * h  # lambda = 1 - alpha h
series = {}
for label, (q, r) in scenarios.items():
    series[label] = simulate_reporter(2.0, 1.0, h, 2000, x0=2.0, q=q, r=r)
    sigma2 = q / (1 - lam ** 2)
    theory = sigma2 / (sigma2 + r) * lam ** np.arange(1, 5)
    estimate = autocorrelation(series[label][1], 4)[1:]
    print(f"{label:17s}", np.round(estimate, 3), np.round(theory, 3))
```

With 2,000 readings, each estimate lies within about 0.035 of the formula. Panel C of @fig-noisy-dynamics-three-ways plots the same simulations: the three autocorrelation curves separate cleanly although the traces in panels A and B look equally noisy.

```{.python .execute}
#| label: fig-noisy-dynamics-three-ways
#| code-fold: true
#| code-summary: "Plot the simulated traces and their autocorrelations"
#| fig-width: 6
#| fig-height: 4
#| fig-cap: "Readings (points) and hidden concentrations (lines) with measurement error only (A) or process noise only (B), at equal reading variance. C: autocorrelations of these simulations and of one with both sources (points), with the formula (lines). Only process noise makes readings correlated, shrinking by 3/4 per image."
#| fig-alt: "Left, two trace panels: orange points scattered around a flat line; blue points on a wandering line. Right: autocorrelation against lag; orange falls to zero at lag 1, blue decays geometrically, green drops to about 0.38 then decays."
fig, axes = plt.subplot_mosaic([["meas", "acf"], ["proc", "acf"]])
panels = {"measurement only": (axes["meas"], cf.SERIES[1], "o", "A  Measurement only"),
          "process only": (axes["proc"], cf.SERIES[0], "s", "B  Process only"),
          "both": (None, cf.SERIES[2], "^", "")}
show = h * np.arange(61)
all_lags = np.arange(9)
for label, (q, r) in scenarios.items():
    x, y = series[label]
    ax, color, marker, title = panels[label]
    if ax is not None:
        ax.plot(show, x[:61], color=cf.INK, linewidth=1)
        ax.plot(show, y[:61], marker, color=color, markersize=3)
        ax.set(ylim=(1.3, 2.7), ylabel="Reading")
        ax.set_title(title, loc="left")
    sigma2 = q / (1 - lam ** 2)
    theory = np.where(all_lags == 0, 1.0, sigma2 / (sigma2 + r) * lam ** all_lags)
    axes["acf"].plot(all_lags, theory, color=color, linewidth=1)
    axes["acf"].plot(all_lags, autocorrelation(y, 8), marker, color=color, label=label)
axes["proc"].set_xlabel("Time (model units)")
axes["acf"].set(xlabel=r"Lag $\ell$ (images)", ylabel=r"Autocorrelation $\rho(\ell)$")
axes["acf"].set_title("C  Autocorrelation", loc="left")
axes["acf"].legend()
plt.show()
```

::: {.callout-note title="Key idea"}
Measurement error is forgotten at the next image; a process disturbance becomes part of the cell and fades at the removal rate $\alpha$, making successive readings correlated.
:::

## What one trajectory supports

Under the model's assumptions, @prp-noise-autocorrelation turns one long trajectory into two estimates. The ratio $\hat\rho(\ell+1)/\hat\rho(\ell)$ estimates $\lambda$, hence $\alpha=(1-\lambda)/h$, and the extrapolated lag-zero height estimates the cell's share of the variance. Deviations stay correlated for about $1/\alpha$, so a trace lasting $T$ time units holds only about $\alpha T$ roughly independent stretches, and a trace of a few relaxation times gives uncertain autocorrelations however many images it contains. Imaging fixed cells or fluorescent beads, whose concentration cannot change, estimates $r$ directly.

Several claims remain out of reach of one trajectory. Photobleaching or focus drift produces correlated measurement errors that the autocorrelation cannot distinguish from process noise. One cell cannot reveal whether other cells share its production rate $\beta$. And one excursion is not evidence of a state change: a reading far from $x^*$ may be a measurement error that the next image contradicts, or a process excursion that will decay. Only a deviation that persists well beyond the relaxation time supports the claim that the cell has settled somewhere new.

## Intrinsic and extrinsic biological variability

Process noise has biological causes, and Elowitz and colleagues separated two kinds [@elowitz2002]. Their *E. coli* strains carried two copies of the same promoter, one driving a cyan and one a yellow fluorescent protein, in one cellular environment. *Intrinsic* noise, the randomness of each gene's own transcription and translation, affects the copies independently. *Extrinsic* noise, from variation in other cellular components such as ribosome numbers, moves both copies together. They reported that both kinds contribute substantially to the overall variation.

:::: {.callout-note collapse="true" title="Going deeper: a two-reporter simulation"}
A toy version uses our model twice per cell: cell $j$ has a production rate $\beta_j$ shared by both reporters (extrinsic), and each reporter has its own disturbances (intrinsic). At a snapshot,
$$
A_j=S_j+I^A_j,\qquad B_j=S_j+I^B_j,
$$
where $S_j=\beta_j/\alpha$ is the steady level both copies share and $I^A_j,I^B_j$ are independent deviations with mean zero and the stationary variance $\operatorname{Var}(I)$ of @prp-stationary-variance. The cross terms vanish as in the proof, so
$$
\operatorname{Cov}(A,B)=\operatorname{Var}(S),\qquad \tfrac12E\bigl[(A-B)^2\bigr]=\tfrac12E\bigl[(I^A-I^B)^2\bigr]=\operatorname{Var}(I).
$$
The shared term cancels in the difference, so extrinsic variation spreads cells along the diagonal $A=B$ and intrinsic variation spreads them across it (@fig-noisy-dynamics-two-reporter). The chunk samples cells from this snapshot model and prints both summaries beside the generating variances. Elowitz and colleagues divide the same two quantities by the product of the reporters' mean levels, reporting squared coefficients of variation so that noise is comparable across expression levels.

```{.python .execute}
#| label: fig-noisy-dynamics-two-reporter
#| code-fold: true
#| code-summary: "Sample two identical reporters in 300 cells"
#| fig-width: 6
#| fig-height: 3.2
#| fig-cap: "Two identically regulated reporters in 300 simulated cells, when intrinsic (A) or extrinsic (B) variation dominates. Shared causes stretch the cloud along the diagonal; independent causes widen it across."
#| fig-alt: "Two scatter plots of reporter A against B with a grey diagonal; A is a round cloud at (2, 2), B a thin ellipse along the diagonal."
np.random.seed(6)
fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
settings = [("A  Mostly intrinsic", 0.004, 0.04),
            ("B  Mostly extrinsic", 0.04, 0.004)]
for ax, (title, var_s, var_i) in zip(axes, settings):
    shared = 2.0 + np.random.normal(0, np.sqrt(var_s), size=300)  # S_j
    a = shared + np.random.normal(0, np.sqrt(var_i), size=300)  # S_j + I_j^A
    b = shared + np.random.normal(0, np.sqrt(var_i), size=300)  # S_j + I_j^B
    covariance = np.cov(a, b)[0, 1]
    half_msd = 0.5 * np.mean((a - b) ** 2)
    print(f"{title}: Cov(A,B) = {covariance:.4f} (Var S = {var_s}),",
          f"half mean (A-B)^2 = {half_msd:.4f} (Var I = {var_i})")
    ax.plot([1.2, 2.8], [1.2, 2.8], color=cf.RULE, linewidth=0.8)
    ax.plot(a, b, "o", color=cf.SERIES[0], markersize=3, alpha=0.7)
    ax.set_title(title, loc="left")
    ax.set(xlabel="Reporter A reading", xlim=(1.2, 2.8), ylim=(1.2, 2.8))
    ax.set_aspect("equal")
axes[0].set_ylabel("Reporter B reading")
plt.show()
```
::::

Intrinsic versus extrinsic is a different question from process versus measurement: both are process noise, changes in the cell that the instrument then reports. The two distinctions still interact, because independent measurement error in each color channel enters $A-B$ and inflates the apparent intrinsic noise by $r$ (@exr-noisy-dynamics-two-channel). The two-reporter experiment therefore needs a calibrated $r$ before its spread across the diagonal can be read as biology.

## Where the model breaks

Additive disturbances are an approximation. At low molecule counts the size of fluctuations depends on the count itself, and additive noise can even push a modeled concentration below zero. Uneven splitting of molecules at division adds a jump at division times rather than a steady stream of small disturbances. A real reporter also has a gain and background, $y_k=c\,x_k+b+v_k$. With $c$ unknown, doubling $c$ while halving every concentration leaves the readings unchanged, so absolute concentrations cannot be recovered from readings alone. Autocorrelations, the removal rate, and the cell's share of the variance are unchanged by this rescaling.

## Exercises

::: {#exr-noisy-dynamics-impulse}
Repeat @exm-single-disturbance with $\alpha=1$ and $h=1/2$. Compute the readings for four images after a process disturbance of 0.4 and after a measurement error of 0.4, and say how many images it takes for the process deviation to fall below 0.05.
:::

::: {#exr-noisy-dynamics-variance}
For $\alpha=1$, $h=1/2$, and $q=0.03$, compute $\sigma_x^2$. Then, with $q=Dh$, show that $\sigma_x^2\to D/(2\alpha)$ as $h\to0$, and explain biologically why faster removal reduces the stationary variance.
:::

::: {#exr-noisy-dynamics-diagnose}
A cell imaged every $h=1/4$ time units for a long time gives readings with variance $0.05$ and estimated autocorrelations $\hat\rho(1)=0.45$, $\hat\rho(2)=0.34$, $\hat\rho(3)=0.25$. Assuming @def-state-space-model with $c=1$, estimate $\lambda$, the removal rate $\alpha$, $\sigma_x^2$, and $r$. Which noise source dominates?
:::

::: {#exr-noisy-dynamics-repeats}
A cell is imaged $m$ times so quickly that $x_k$ is fixed, giving readings $x_k+v_k^{(i)}$ with independent errors of variance $r$. Show that their average has measurement variance $r/m$ while its deviation $e_k$ is unchanged, and explain why these technical replicates say nothing about process noise.
:::

::: {#exr-noisy-dynamics-two-channel}
In the two-reporter model, add independent measurement errors with variance $r$ to each channel: $A_j=S_j+I^A_j+v^A_j$ and $B_j=S_j+I^B_j+v^B_j$. Show that $\operatorname{Cov}(A,B)$ is unchanged while $\tfrac12E[(A-B)^2]$ increases by $r$.
:::

## Bridge to the final lesson

The three dynamics lessons built a model, found its persistent states, and asked what noisy readings can and cannot reveal about a hidden state. Lesson 24 asks the same question of sequencing: what longer reads, and references built from many genomes, let us recover that short reads against one reference could not.

::: {#refs}
**References**
:::
