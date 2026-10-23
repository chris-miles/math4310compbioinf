---
title: "Notation"
---

Symbols are defined when they first appear in a lesson. Some letters have different meanings in different models; the surrounding definition determines their use.

| Symbol | Meaning |
|---|---|
| $\Sigma$ | Finite alphabet |
| $x=x_1\ldots x_n$, $y=y_1\ldots y_m$ | Observed sequences of lengths $n$ and $m$ |
| $|x|$ | Length of sequence $x$ |
| $x_1\ldots x_i$ | Prefix containing the first $i$ symbols |
| $s(a,b)$ | Score for aligning symbols $a$ and $b$ |
| $g<0$ | Signed score per gap character in a linear-gap model |
| $F(i,j)$ | Best global-alignment score for two prefixes |
| $H(i,j)$ | Best local-alignment score ending at two prefix boundaries, allowing an empty alignment |
| $d_o,d_e>0$ | Gap-open and gap-extension penalties |
| $G(\ell)=d_o+(\ell-1)d_e$ | Cost of a length-$\ell$ gap |
| $Q(i,j)$ | Fitting-alignment score for a read prefix against a reference prefix |
| $h$ | Number of seed hits to verify |
| $M,I_x,I_y$ | Affine-alignment tables for paired letters, a gap in $y$, and a gap in $x$ |
| $\mathcal O(\cdot)$ | Asymptotic upper bound on growth |
| $k$ | Word length in a $k$-mer |
| $V,E$ | Vertex and edge sets of a graph |
| $d^-(v),d^+(v)$ | Incoming and outgoing edge counts at vertex $v$, including multiplicity |
| $X_i$, $x_i$ | Random sequence symbol and its observed value |
| $q_b$ | Base probability, or initial-base probability for a Markov chain |
| $a_{ab}$ | Probability of moving from $a$ to $b$; source indexes the row |
| $\theta$, $L(\theta;x)$ | Model parameters and their likelihood for fixed data $x$ |
| $\ell(\theta;x)$ | Log-likelihood |
| $n_{ab}$ | Observed transition count from $a$ to $b$ |
| $\alpha$ | Pseudocount added to each possible outcome |
| $M_+,M_-$ | Island and background sequence models |
| $S(x)$ | Log-likelihood ratio for a window |
| $\beta(a,b)$ | Log ratio contributed by an adjacent pair in the CpG classifier |
| $Q$, $K=|Q|$ | Hidden state set and number of states |
| $\pi=\pi_1\ldots\pi_n$ | Hidden state path |
| $p_k$ | Initial probability of hidden state $k$ |
| $e_k(b)$ | Probability of emitting symbol $b$ in state $k$ |
| $v_k(i)$ | Largest joint probability of an observation prefix and a path ending in state $k$ |
| $V_k(i)=\log v_k(i)$ | Viterbi score in log space |
| $f_k(i)$ | Joint probability of the observed prefix and state $k$ at $i$ |
| $b_k(i)$ | Probability of the observed suffix after $i$, conditional on state $k$ at $i$ |
| $\gamma_k(i)$ | Posterior state probability $P(\pi_i=k\mid x)$ |
| $\xi_{k\ell}(i)$ | Posterior probability of transition $k\to\ell$ from $i$ to $i+1$ |
| $\operatorname{LSE}(z)$ | Log-sum-exp, $\log\sum_j e^{z_j}$ |
| $f_{ij}$ | Fraction of mismatches between aligned sequences $i,j$ on retained columns |
| $D=(d_{ij})$ | Matrix of estimated pairwise evolutionary distances |
| $d_T(i,j)$ | Sum of branch lengths along the path between tips $i,j$ in tree $T$ |
| $\ell_e$, $h_u$ | Branch length of edge $e$ and height of node $u$ above descendant tips in a clock tree |
| $d(C,D)$ | Mean original distance over all tip pairs from clusters $C,D$ |
| $x_{\ell j}$, $z_v$ | Observed character at leaf $\ell$, site $j$, and an assigned ancestral character at node $v$ |
| $c_j(T)$, $c(T)$ | Minimum equal-cost change count at site $j$ on tree $T$, and its sum over sites |
| $S_v$, $c_v$ | Fitch set of minimum-cost characters at subtree root $v$, and its subtree change count |
| $X\in\mathbb R^{n\times p}$, $Z$ | Transformed sample-by-gene measurements and their gene-column-centered matrix |
| $S=Z^\top Z/(n-1)$ | Sample covariance matrix of gene measurements |
| $v_k$, $t_k=Zv_k$ | PCA loading direction in gene space and its sample-score vector |
| $Z=UDV^\top$ | Singular value decomposition, with $D=\operatorname{diag}(d_k)$ in this context |
| $T=ZV=UD$, $\lambda_k=d_k^2/(n-1)$ | PCA score matrix and component variances |
| $\rho_{ij}$ | Probability that a sampled RNA-seq read in sample i is assigned to gene j |
| $C_{ij}$, $L_i$, $R_{ij}$ | Assigned RNA-seq read count, sample's total assigned count, and counts per million |
| $\sigma_j$, $W=ZG^{-1}$ | Gene's sample standard deviation and standardized matrix, with $G=\operatorname{diag}(\sigma_j)$ |
| $\Delta_j$, $\widehat\Delta_j$ | Population mean contrast between treatment and control for feature $j$, and its observed estimate |
| $\mathcal A$, $S(a)$ | Allowed label assignments and the comparison statistic under assignment $a$ |
| $p_j$ | Calibrated p-value for feature $j$ in the testing lessons |
| $V$, $R$, $\mathrm{FDP}=V/\max(R,1)$ | False rejections, total rejections, and realized false discovery proportion |
| $\mathrm{FDR}=E[\mathrm{FDP}]$ | Expected false discovery proportion |
| $m$, $m_0$, $\mathcal N$ | Number of genes tested, number with a true null hypothesis, and the set of those null genes |
| $q$, $p_{(i)}$ | Target FDR level, and the $i$-th smallest of the $m$ p-values |
| $x(t)$, $\beta$, $\alpha$ | Gene-product concentration, production rate, and removal-rate constant in the scalar ODE |
| $x^*=\beta/\alpha$ | Steady concentration for constant production and removal |
| $K$, $n$, $a=\beta/(\alpha K)$ | Repression threshold, Hill coefficient, and scaled maximal production in the switch lesson |
| $u$, $v$, $H(z)=a/(1+z^n)$, $J$ | Scaled toggle-switch repressor concentrations, scaled repression function, and Jacobian |
| $x_k$, $y_k$ | Hidden concentration and its recorded observation at time $t_k$ |
| $F_h(x)$, $h$ | Deterministic state update across an observation interval, and that interval's duration |
| $w_k$, $v_k$ | Process disturbance and measurement error in the noisy-dynamics lesson |
| $q$, $r$, $c$ | Variances of $w_k$ and $v_k$, and the reporter gain in $y_k=c\,x_k+v_k$ |
| $e_k=x_k-x^*$, $\lambda=1-\alpha h$ | Deviation from steady state and its per-step Euler multiplier |
| $\sigma_x^2=q/(1-\lambda^2)$, $\rho(\ell)$ | Stationary concentration variance and the readings' autocorrelation at lag $\ell$ |

In the toy HMM, state $B$ denotes background and state $I$ denotes an island-like region. Observation $S$ groups C/G and observation $W$ groups A/T. These observation symbols are distinct from the window score $S(x)$.

A logarithm's base is stated when its units matter: base 2 gives bits, and the implementations specify when they use natural logarithms.
