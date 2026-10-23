---
published: false
title: "Biological Switches and Cell States"
subtitle: "Lesson 22 · 2027-04-14 · Week 14"
nocite: |
  @mangelGeneNetworks
---

## Core question

How can genetically identical cells in the same conditions hold two different expression states, and keep whichever state they are in?

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

## A circuit that remembers

In Lesson 21, a protein made at a constant rate $\beta$ and removed at rate $\alpha x$ approached the single steady state $x^*=\beta/\alpha$ from every starting concentration. Any earlier treatment was forgotten within a few relaxation times $1/\alpha$. Yet a differentiated cell keeps its identity long after the signals that set it have disappeared, while neighbors with the same genome keep a different one. Describing that memory needs regulation in which the cell's current state shapes its future production.

Gardner, Cantor, and Collins built such a circuit in the bacterium *Escherichia coli* [@gardner2000]. A plasmid, a small DNA circle separate from the bacterial chromosome, carried genes for two repressors, proteins that bind DNA to reduce transcription: LacI and a temperature-sensitive repressor from bacteriophage λ. Each gene had a promoter blocked by the other repressor. The gene for green fluorescent protein (GFP) is transcribed with the λ repressor, so a green cell is making the λ repressor while LacI production is shut off. Two transient inputs push the circuit: the chemical IPTG disables LacI, and a pulse of heat disables the λ repressor.

The readout came from a *flow cytometer*, which passes cells one at a time through a laser and records each cell's scattered light and fluorescence. One record is one cell; a culture's thousands of records are summarized as a fluorescence histogram, on a logarithmic axis because expression spans orders of magnitude.

Cells grown for 6 hours with IPTG became green. After the IPTG was washed out, they stayed green for the next 22 hours of the experiment, while a control plasmid whose GFP promoter was simply induced by IPTG, with no feedback, returned to the dark state. At an intermediate IPTG concentration near the switching threshold, a single culture split into a dark population and a bright population [@gardner2000]. @fig-switch-flow-histograms shows simulated histograms with the same three shapes.

```{.python .execute}
#| label: fig-switch-flow-histograms
#| code-fold: true
#| code-summary: "Simulate flow-cytometry histograms"
#| fig-width: 5
#| fig-height: 4
#| fig-cap: "Simulated single-cell GFP histograms for a toggle-switch culture never induced (top), induced with IPTG for 6 hours and then grown without it (middle), and held near the switching threshold (bottom). Cells occupy one or two distinct states, not a spread of intermediate levels."
#| fig-alt: "Three stacked histograms on a logarithmic fluorescence axis; the top has one peak near 20, the middle one near 1,000, and the bottom both peaks."
rng = np.random.default_rng(22)
n_cells = 5000
# Log10 fluorescence of the two expression states, with cell-to-cell spread.
low = rng.normal(1.3, 0.18, n_cells)
high = rng.normal(3.0, 0.20, n_cells)
mixed = np.where(rng.random(n_cells) < 0.55, low, high)
cultures = [("Never induced", low), ("IPTG pulse, then washout", high),
            ("Near-threshold IPTG", mixed)]
fig, axes = plt.subplots(3, 1, sharex=True)
bins = np.linspace(0.5, 4, 70)
for ax, (name, values) in zip(axes, cultures):
    ax.hist(values, bins=bins, color=cf.SERIES[0], edgecolor="none", alpha=0.9)
    ax.set_title(name, loc="left")
    ax.set_ylabel("Cells")
axes[-1].set_xticks([1, 2, 3, 4], ["$10^1$", "$10^2$", "$10^3$", "$10^4$"])
axes[-1].set_xlabel("GFP fluorescence (a.u., log scale)")
plt.show()
```

The middle histogram is the memory: the input is gone, yet the culture has not returned to the top histogram. The bottom one shows cells under a single condition occupying one of two states. As a mathematical problem, the input is the regulatory wiring, written as a model $\dot{\mathbf{x}}=f(\mathbf{x})$ for the protein concentrations. The output is its steady states, which of them attract nearby states, and the parameter conditions under which two do. An attracting steady state is a candidate persistent cell state, and the boundary between their attraction regions is the threshold a transient input must push a cell across.

## Hill functions: how a repressor sets production

A repressor lowers production by binding a DNA site next to a promoter and blocking RNA polymerase, so production is proportional to the fraction of time the site is free. Idealize the site as either free or occupied by a complex of $n$ repressor molecules that bind together, with binding fast enough that the site is at equilibrium for the current repressor concentration $r$. Occupation needs all $n$ molecules at the site at once, and each is there with a chance proportional to $r$, so the ratio of occupied to free time scales as $r^n$; write it as $(r/K)^n$ for a constant $K>0$. Requiring all $n$ to bind together is what makes the response sharp. The free fraction is
$$
\frac{\text{free}}{\text{free}+\text{occupied}}=\frac{1}{1+(r/K)^n}.
$$

::: {#def-hill-repression}
## Hill repression function
For a repressor at concentration $r\ge0$, a gene with maximal production rate $\beta>0$, repression threshold $K>0$, and Hill coefficient $n>0$ has production rate
$$
\frac{\beta}{1+(r/K)^n}.
$$
:::

At $r=0$ the gene produces at its full rate $\beta$, and at $r=K$ it produces at half that rate, whatever $n$ is. Real binding is less tidy than the all-or-none picture, so measured curves are fitted with $n$ as a free parameter, which is why fitted values are often not integers. The coefficient then measures effective *cooperativity*, the degree to which repressor molecules act together.

```{.python .execute}
#| label: fig-switch-hill-curves
#| code-fold: true
#| code-summary: "Plot Hill repression curves"
#| fig-width: 5
#| fig-height: 3.2
#| fig-cap: "Fraction of maximal production, 1/(1 + (r/K)^n), for Hill coefficients 1, 2, and 4. All curves pass through one half at r = K; dotted lines mark r = K/2 and r = 2K. Larger n turns a gradual dial into a near on/off switch."
#| fig-alt: "Three decreasing curves from one, all crossing one half at r/K = 1; the n = 4 curve is steepest and the n = 1 curve shallowest."
ratio = np.linspace(0, 3, 301)
fig, ax = plt.subplots()
for n, line_style in [(1, ":"), (2, "--"), (4, "-")]:
    ax.plot(ratio, 1 / (1 + ratio ** n), line_style, label=f"$n={n}$")
for mark in (0.5, 2):
    ax.axvline(mark, color=cf.RULE, linestyle=":", linewidth=0.8)
ax.set(xlabel="Repressor concentration, $r/K$",
       ylabel="Fraction of maximal production")
ax.legend()
plt.show()
```

Compare a fourfold change in repressor, from $r=K/2$ to $r=2K$ (@fig-switch-hill-curves). With $n=1$, production falls from $1/(1+0.5)=0.67$ to $1/(1+2)=0.33$ of maximum, a factor of two. With $n=4$, it falls from $1/(1+0.5^4)=0.94$ to $1/(1+2^4)=0.06$, a factor of sixteen. A cooperative repressor therefore turns a modest change in its own concentration into a nearly all-or-none change in its target's production, and the switches below depend on that steepness. An *activator* that must bind to turn its gene on gives the complementary form $\beta (r/K)^n/(1+(r/K)^n)$, which rises from zero to $\beta$.

## One gene that activates itself

Feedback already creates several steady states with one gene; the simplest case is a protein that activates its own promoter. With Hill coefficient $n=2$ and linear removal as in Lesson 21,
$$
\dot x=\beta\,\frac{(x/K)^2}{1+(x/K)^2}-\alpha x .
$$
Three parameters reduce to one if we measure concentration in units of $K$ and time in units of the relaxation time $1/\alpha$. Substituting $x=K\tilde x$ and $t=\tilde t/\alpha$, so that $dx/dt=K\alpha\,d\tilde x/d\tilde t$, and dividing by $K\alpha$ gives
$$
\frac{d\tilde x}{d\tilde t}=a\,\frac{\tilde x^2}{1+\tilde x^2}-\tilde x,\qquad a=\frac{\beta}{\alpha K}.
$$
The parameter $a$ is Lesson 21's steady state $\beta/\alpha$ at full production, measured in units of the threshold $K$. From here on we drop the tildes, so the model is $\dot x=f(x)$ with
$$
f(x)=P(x)-x,\qquad P(x)=\frac{a\,x^2}{1+x^2}.
$$

A steady state (Lesson 21) is a value $x^*$ with $f(x^*)=0$; here that means $P(x^*)=x^*$, production balancing removal. Graphically, steady states are the crossings of the production curve $y=P(x)$ with the removal line $y=x$. Between crossings, the sign of $f$ gives the direction of motion: where production lies above removal, $x$ increases.

::: {#exm-self-activation}
## Three steady states for $a=3$
Set $a=3$. The equation $3x^2/(1+x^2)=x$ has the root $x=0$. For $x>0$, divide by $x$ and multiply by $1+x^2$ to get $3x=1+x^2$, or $x^2-3x+1=0$, so
$$
x^*\in\Bigl\{0,\ \tfrac{3-\sqrt5}{2},\ \tfrac{3+\sqrt5}{2}\Bigr\}\approx\{0,\ 0.382,\ 2.618\}.
$$
At $x=1$, between the second and third crossings, $f(1)=3/2-1=0.5>0$, so concentration rises toward $2.618$. At $x=0.2$, $f(0.2)=0.12/1.04-0.2\approx-0.085<0$, so concentration falls toward $0$.
:::

A precise test of which steady states attract uses the slope of $f$, extending Lesson 21's deviation argument to a nonlinear rate.

::: {#prp-slope-stability}
## Stability from the slope
Let $f$ be differentiable with $f(x^*)=0$. If $f'(x^*)<0$, small deviations from $x^*$ decay and $x^*$ is *stable*. If $f'(x^*)>0$, small deviations grow and $x^*$ is *unstable*.
:::

::: {.proof}
Sketch. Write $e=x-x^*$ for the deviation. Because $f(x^*)=0$, the tangent-line approximation gives $\dot e=f(x^*+e)\approx f'(x^*)\,e$. This is Lesson 21's linear equation with $-\alpha$ replaced by $f'(x^*)$, so $e(t)\approx e(0)\,e^{f'(x^*)t}$. The deviation decays when the exponent is negative and grows when it is positive; the approximation needs only $e$ small.
:::

Here $f'(x)=6x/(1+x^2)^2-1$. At $0$ the slope is $-1$; at $0.382$ it is $+0.745$; at $2.618$ it is $-0.745$. The off state and the on state are stable, and the middle crossing is unstable. The unstable state is the threshold: below it removal wins, above it production wins (@fig-switch-self-activation, panel A).

That threshold explains memory. Add a transient input that raises production by $s$ for a while, $\dot x=s+P(x)-x$. A strong enough input lifts the production curve until the low crossings vanish, and $x$ climbs. When the input ends, the original equation returns, and what matters is only which side of $0.382$ the concentration has reached. We simulate a pulse with Lesson 21's Euler step, $x_{k+1}=x_k+h\bigl(s_k+P(x_k)-x_k\bigr)$, where $s_k$ is the input level during the pulse and zero otherwise.

```{.python .execute}
#| label: code-switch-pulse
def production(x: np.ndarray, a: float = 3.0) -> np.ndarray:
    """Self-activating production a x^2 / (1 + x^2) in scaled units."""
    return a * x ** 2 / (1 + x ** 2)


def pulse_response(duration: float, level: float = 1.0, start: float = 1.0,
                   t_end: float = 8.0, h: float = 0.01) -> np.ndarray:
    """Euler trajectory from x = 0 with input `level` during one pulse."""
    n_steps = int(round(t_end / h))
    x = np.zeros(n_steps + 1)
    for k in range(n_steps):
        in_pulse = start <= k * h < start + duration
        s = level if in_pulse else 0.0
        x[k + 1] = x[k] + h * (s + production(x[k]) - x[k])
    return x


short_pulse = pulse_response(0.25)
long_pulse = pulse_response(0.5)
# Concentration when each pulse ends (step t / h) and at t = 8.
print(f"short: {short_pulse[125]:.3f} -> {short_pulse[-1]:.3f}")
print(f"long:  {long_pulse[150]:.3f} -> {long_pulse[-1]:.3f}")
```

The short pulse leaves the concentration at $0.234$, below the threshold, and the cell drifts back to the off state. The pulse twice as long reaches $0.480$, just past $0.382$, and the cell continues to the on state at $2.618$ after the input has stopped. Two inputs that differ only slightly produce permanently different cells.

```{.python .execute}
#| label: fig-switch-self-activation
#| code-fold: true
#| code-summary: "Draw the rate balance and the pulse responses"
#| fig-width: 6
#| fig-height: 3
#| fig-cap: "(A) Production 3x^2/(1 + x^2) and removal x; filled and open circles are stable and unstable steady states, and arrows below give the direction of change. (B) Responses to input pulses of level 1 lasting 0.25 or 0.5 (shaded) time units. Only the pulse that carries concentration past the unstable state at 0.382 (dotted) leaves a lasting change."
#| fig-alt: "Left: an S-shaped curve crosses a line at 0, 0.38, and 2.62, with arrows toward 0 and 2.62. Right: a dashed time course rises to 0.23 and decays; a solid one passes a dotted line and climbs to 2.6."
times = np.arange(len(long_pulse)) * 0.01
threshold = (3 - np.sqrt(5)) / 2
fig, axes = plt.subplots(1, 2)
cf.rate_balance(
    axes[0], production, lambda x: x, (0, 3.5),
    fixed_points=[(0.0, True), (threshold, False), ((3 + np.sqrt(5)) / 2, True)],
    xlabel="Concentration, $x$ (scaled)", ylabel="Rate (scaled)",
)
for line in axes[0].lines:
    line.set_clip_on(False)  # keep the markers at x = 0 whole
axes[1].axvspan(1, 1.5, color=cf.LIGHT)
axes[1].plot(times, long_pulse, "-", label="0.5-unit pulse")
axes[1].plot(times, short_pulse, "--", label="0.25-unit pulse")
axes[1].axhline(threshold, color=cf.RULE, linestyle=":", linewidth=1)
axes[1].set(xlabel="Time, $t$ (scaled)", ylabel="Concentration, $x$ (scaled)")
axes[1].legend(loc="upper left", bbox_to_anchor=(0.22, 1.0))
for ax, letter in zip(axes, "AB"):
    ax.set_title(letter, loc="left", fontweight="bold", fontsize=11)
plt.show()
```

## The toggle switch: two repressors in a loop

Gardner's circuit has no self-activating gene, but two negatives make a positive: if protein 1 rises, protein 2's production falls, and less protein 2 lets protein 1's production rise further. Mutual repression is a positive feedback loop, so we expect the same threshold behavior. Let $u$ and $v$ be the two repressor concentrations, scaled as above, and assume the two genes share the same maximal production, threshold, and Hill coefficient. With the scaled Hill repression function
$$
H(z)=\frac{a}{1+z^n},
$$
the model of @gardner2000, in its symmetric form, is
$$
\dot u=H(v)-u,\qquad \dot v=H(u)-v .
$$ {#eq-toggle}
Each line is Lesson 21's balance with production set by the other protein. A point $(u,v)$ describes both concentrations in one cell, and its motion is a trajectory of that cell's state.

::: {#def-bistable}
## Bistability
A system is *bistable* when, for fixed parameters, it has two stable steady states: each attracts every trajectory that starts close enough to it.
:::

Steady states of @eq-toggle are points where $\dot u=\dot v=0$. The set where $\dot u=0$ is the curve $u=H(v)$, and the set where $\dot v=0$ is the curve $v=H(u)$. Steady states are exactly the intersections of these *nullclines*. Off them, the signs of $\dot u$ and $\dot v$ say which way a trajectory moves; to the right of $u=H(v)$, for example, removal of protein 1 exceeds its production and $u$ decreases.

::: {#exm-toggle-steady-states}
## Steady states of the toggle with $a=3$, $n=2$
Steady states satisfy $u(1+v^2)=3$ and $v(1+u^2)=3$. Subtracting the second equation from the first gives $u-v+uv(v-u)=0$, which factors as
$$
(u-v)(1-uv)=0 .
$$
Either the two concentrations are equal or their product is one.

*Equal concentrations.* With $u=v=s$, the condition is $s+s^3=3$. The left side increases with $s$, so there is exactly one root; $1.2+1.2^3=2.928$ and $1.22+1.22^3=3.036$ bracket it, and the root is $s\approx1.213$.

*Product one.* With $v=1/u$, the first equation becomes $u+1/u=3$, so $u+v=3$ and $uv=1$. The two concentrations are the roots of $w^2-3w+1=0$, giving
$$
(u,v)=\Bigl(\tfrac{3+\sqrt5}{2},\tfrac{3-\sqrt5}{2}\Bigr)\approx(2.618,0.382)
$$
and its mirror image $(0.382,2.618)$. In the first, protein 1 is high and holds protein 2 low; in the second, the roles are reversed.
:::

To decide which of the three steady states attract, we extend @prp-slope-stability to two variables. Near a steady state $(u^*,v^*)$, write the deviation as $\mathbf e=(u-u^*,\,v-v^*)$. Replacing each rate by its tangent-plane approximation gives the linear system $\dot{\mathbf e}\approx J\mathbf e$, where $J$ is the *Jacobian*, the matrix of partial derivatives of the two rates:
$$
J=\begin{pmatrix}\partial\dot u/\partial u & \partial\dot u/\partial v\\ \partial\dot v/\partial u & \partial\dot v/\partial v\end{pmatrix}
=\begin{pmatrix}-1 & H'(v^*)\\ H'(u^*) & -1\end{pmatrix}.
$$
Suppose $J\mathbf w=\lambda\mathbf w$, and try a deviation that stays along $\mathbf w$, $\mathbf e(t)=c(t)\,\mathbf w$. Then $\dot c\,\mathbf w=c\,J\mathbf w=\lambda c\,\mathbf w$, so $\dot c=\lambda c$: Lesson 21's decay equation with $-\alpha$ replaced by $\lambda$, giving $c(t)=c(0)e^{\lambda t}$. Any small deviation is a combination of the two eigenvectors, so it decays if both eigenvalues are negative. If one eigenvalue is positive, deviations along its eigenvector grow; a steady state with one positive and one negative eigenvalue is a *saddle*, approached along one direction and left along the other.

The eigenvalues solve
$$
\det(J-\lambda I)=(-1-\lambda)^2-H'(u^*)H'(v^*)=0,
\qquad\text{so}\qquad
\lambda=-1\pm\sqrt{H'(u^*)\,H'(v^*)} .
$$
Both slopes of the repression function are negative, so the product under the root is positive and the eigenvalues are real. With $H'(z)=-2az/(1+z^2)^2$ and $a=3$:

- At the high-$u$ state, $H'(0.382)=-6(0.382)/1.146^2\approx-1.745$ and $H'(2.618)=-6(2.618)/7.854^2\approx-0.255$. Their product is $0.444=0.667^2$, so $\lambda\approx-0.333$ and $-1.667$: stable, and by symmetry so is its mirror image.
- At the symmetric state, $H'(1.213)\approx-1.191$, so $\lambda\approx-1+1.191=0.191$ and $-1-1.191=-2.191$. The state is a saddle.

At the saddle both off-diagonal entries equal $H'(s)$, so
$$
J\begin{pmatrix}1\\-1\end{pmatrix}=\bigl(-1-H'(s)\bigr)\begin{pmatrix}1\\-1\end{pmatrix},
\qquad
J\begin{pmatrix}1\\1\end{pmatrix}=\bigl(-1+H'(s)\bigr)\begin{pmatrix}1\\1\end{pmatrix}.
$$
The growing direction, with eigenvalue $-1-H'(s)=0.191$, is $(1,-1)$: a small difference between the two proteins grows. The decaying direction is the diagonal. In fact the whole diagonal is a trajectory, because by symmetry a point with $u=v$ has $\dot u=\dot v$ and so stays on it. Two trajectories cannot cross, since the model assigns a single direction of motion to each point, so no cell can move from one side of the diagonal to the other. A cell that starts with more protein 1 therefore ends in the high-$u$ state, and one that starts with more protein 2 ends in the high-$v$ state (@fig-switch-toggle-phase-plane, panel A).

```{.python .execute}
#| label: fig-switch-toggle-phase-plane
#| code-fold: true
#| code-summary: "Simulate the toggle and draw its phase planes"
#| fig-width: 6
#| fig-height: 3.3
#| fig-cap: "Phase planes of the symmetric toggle with a = 3, for n = 2 (A) and n = 1 (B), with nullclines, stable (filled) and unstable (open) steady states, and trajectories from four shared starts. With n = 2 each start reaches the stable state on its side of the diagonal; the nullclines nearly touch because a = 3 is just above the threshold a = 2. With n = 1 every start reaches the single steady state."
#| fig-alt: "Two plots of v against u. Left: two labeled curves meet at filled points near (2.62, 0.38) and (0.38, 2.62) and an open point near (1.21, 1.21); trajectories split between the filled points. Right: the curves meet once, near (1.30, 1.30), where all trajectories end."
def simulate_toggle(start: tuple[float, float], a: float, n: float,
                    h: float = 0.01, n_steps: int = 3000) -> np.ndarray:
    """Euler trajectory of the symmetric toggle from one starting pair."""
    path = np.zeros((n_steps + 1, 2))
    path[0] = start
    for k in range(n_steps):
        u, v = path[k]
        rates = np.array([a / (1 + v ** n) - u, a / (1 + u ** n) - v])
        path[k + 1] = path[k] + h * rates
    return path


starts = [(2.0, 0.4), (0.4, 2.0), (0.2, 0.5), (2.9, 2.4)]
high, low = (3 + np.sqrt(5)) / 2, (3 - np.sqrt(5)) / 2
s = 1.2134116627622127  # root of s + s^3 = 3 when n = 2
s_one = (-1 + np.sqrt(13)) / 2  # root of s + s^2 = 3 when n = 1
panels = [(2, [(high, low, True), (low, high, True), (s, s, False)]),
          (1, [(s_one, s_one, True)])]
fig, axes = plt.subplots(1, 2)
for ax, (n, fixed_points), letter in zip(axes, panels, "AB"):
    cf.phase_plane(
        ax,
        lambda u, v, n=n: 3 / (1 + v ** n) - u,
        lambda u, v, n=n: 3 / (1 + u ** n) - v,
        xlim=(0, 3), ylim=(0, 3), fixed_points=fixed_points,
        trajectories=[simulate_toggle(start, 3, n) for start in starts],
        labels=("$u$ (scaled)", "$v$ (scaled)"),
    )
    ax.set_aspect("equal")
    ax.set_title(f"{letter}   $n={n}$", loc="left", fontsize=11)
pointer = dict(arrowstyle="-", color=cf.INK, lw=0.6)
axes[0].annotate(r"$\dot v=0$", xy=(0.75, 1.92), xytext=(1.1, 2.4),
                 fontsize=9, arrowprops=pointer)
axes[0].annotate(r"$\dot u=0$", xy=(0.34, 2.88), xytext=(0.65, 2.78),
                 fontsize=9, arrowprops=pointer)
plt.show()
```

With $u$ as LacI and $v$ as the λ repressor, an IPTG pulse lets $v$ rise until the state crosses the diagonal; after washout the cell keeps flowing to the green high-$v$ state, which is the persistence in the middle histogram of @fig-switch-flow-histograms.

## When is the toggle bistable?

Drawing two inhibitory arrows does not guarantee a switch: with $n=1$ and the same $a=3$, the nullclines cross only once (@fig-switch-toggle-phase-plane, panel B). In the cases we can settle by hand, the switch appears exactly when the symmetric state becomes a saddle, which requires the slope of the repression function there to exceed one in magnitude.

::: {#prp-toggle-bistability}
## Cooperativity is needed for the symmetric toggle
Consider @eq-toggle with $a>0$.

1. If $n=1$, there is exactly one steady state, $(s,s)$, and it is stable.
2. If $n=2$, the system is bistable exactly when $a>2$. Then the stable states are $\bigl((a\pm\sqrt{a^2-4})/2,\,(a\mp\sqrt{a^2-4})/2\bigr)$, with eigenvalues $-1\pm2/a$, and $(s,s)$ is a saddle.
:::

::: {.proof}
A symmetric state satisfies $s(1+s^n)=a$, so $a/(1+s^n)=s$ and $|H'(s)|=a\,n\,s^{n-1}/(1+s^n)^2=n\,s^n/(1+s^n)<n$; its eigenvalues are $-1\pm|H'(s)|$.

For $n=1$, the steady-state equations are $u(1+v)=a$ and $v(1+u)=a$. Subtracting gives $u-v=0$, so every steady state is symmetric, and $s+s^2=a$ has one positive root. There $|H'(s)|<1$, so both eigenvalues are negative.

For $n=2$, the subtraction in @exm-toggle-steady-states works for any $a$: asymmetric states satisfy $uv=1$ and $u+v=a$, so $u$ and $v$ are the roots of $w^2-aw+1=0$. These are real and distinct exactly when $a^2>4$, that is $a>2$. At such a state, $1+v^2=a/u$ and $1+u^2=a/v$, so
$$
H'(u)H'(v)=\frac{4a^2uv}{(1+u^2)^2(1+v^2)^2}=\frac{4a^2uv\cdot u^2v^2}{a^4}=\frac{4}{a^2},
$$
using $uv=1$. The eigenvalues $-1\pm2/a$ are negative when $a>2$. At the symmetric state, $|H'(s)|=2s^2/(1+s^2)$ exceeds one exactly when $s>1$, which, because $s+s^3$ increases with $s$, means $a=s+s^3>2$. For $a\le2$, the symmetric state is the only steady state, and it is stable for $a<2$.
:::

The bound $|H'(s)|<n$ is the heart of the matter: the Hill coefficient caps how steep repression can be at the balance point, and with $n=1$ a larger $a$ moves the single steady state without ever splitting it. With $n=2$, the condition $a=\beta/(\alpha K)>2$ says that each gene at full production must accumulate well above the threshold at which it represses its partner.

## What the switch model explains, and what it leaves out

The model accounts qualitatively for the first two histograms. The near-threshold histogram exposes its main omission: in @eq-toggle, cells with the same starting point follow the same trajectory, so the model alone would put a whole culture on one side of the boundary. Gardner and colleagues attributed the two populations to fluctuations in gene expression that scatter cells near the threshold [@gardner2000]. Mutually repressing transcription factors, such as PU.1 and GATA1 in blood progenitor cells, have been modeled with circuits of this kind as an account of how a cell-fate decision is held in place [@bhattacharya2011].

Fluctuations raise a question the deterministic model cannot answer. A real fluctuation in concentration that carries a cell across the diagonal leaves it in the other persistent state, while an error that moves only the recorded reading changes nothing about the cell. When a reporter trace jumps toward the other state, how can we tell which happened? Lesson 23 separates these two kinds of noise.

## Exercises

::: {#exr-switch-toggle-a10}
For the symmetric toggle with $n=2$ and $a=10$, show that $(2,2)$ is a steady state, compute its eigenvalues, and find the two asymmetric steady states exactly.
:::

::: {#exr-switch-hysteresis}
For the self-activating gene with $a=3$ and constant input $s\ge0$, $\dot x=s+3x^2/(1+x^2)-x$, find numerically the largest $s$ at which a low steady state exists. Explain why raising $s$ above this value and then returning it to zero leaves the cell on.
:::

::: {#exr-switch-unbalanced}
Plot the nullclines of the unbalanced toggle $\dot u=3/(1+v^2)-u$, $\dot v=2/(1+u^2)-v$ and simulate it from the starting points of @fig-switch-toggle-phase-plane. How many steady states are there, and what does this say about repressors of unequal strength?
:::

::: {#exr-switch-memory-test}
A culture shows a bright and a dim population. Describe a perturbation and a sequence of measurements that could distinguish a persistent switch from a slow return to a single steady state.
:::

::: {#refs}
**References**
:::
