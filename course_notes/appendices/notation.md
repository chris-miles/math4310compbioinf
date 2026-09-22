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
| $d,e>0$ | Gap-open and gap-extension penalties |
| $G(k)=d+(k-1)e$ | Cost of a length-$k$ gap |
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

In the toy HMM, state $B$ denotes background and state $I$ denotes an island-like region. Observation $S$ groups C/G and observation $W$ groups A/T. These observation symbols are distinct from the window score $S(x)$.

A logarithm's base is stated when its units matter: base 2 gives bits, and the implementations specify when they use natural logarithms.
