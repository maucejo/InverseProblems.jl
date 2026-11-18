#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Algorithms])
#metadata((description: "Iteratively reweighted algorithm, resolution strategies and linear solvers of SourceEstimation.jl.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

All the formulations of SourceEstimation.jl lead to the same computational kernel: a sequence of weighted Tikhonov problems. This page describes the iteratively reweighted algorithm, the resolution strategies of space–frequency and space–time problems, and the linear solvers.

= Iteratively reweighted algorithm

The IR algorithm (`IR()`) is derived from the direct application of the first-order optimality condition of the regularized problem @Aucejo2019c @Aucejo2020b. It is the algorithm of the additive and multiplicative formulations, and the $u$-step of the Bayesian MAP.

== Fixed point

At iteration $k + 1$, with $V = L U^((k))$:

+ compute the local weights $W$ and the penalties $P_g = R_g (V_g)$ of the groups;
+ update the regularization parameter(s) — fixed, chosen by a rule, adaptive ($alpha$ or $alpha_g$ of the #link("regularization.html")[multiplicative formulations]) or given by the hyperparameters of the #link("bayesian.html")[Bayesian model];
+ solve the weighted Tikhonov problem
  $ (H^H H + L^H D L) u = H^H y, quad D = lambda W. $ <eq-ir>

The iterations stop when the relative variation $norm(U^((k+1)) - U^((k)))_1 slash norm(U^((k)))_1$ falls below `tol`, or after `maxiter` iterations (100 by default). The initial solution is a Tikhonov solution whose parameter is given by `init`: `Heuristic()` by default for the multiplicative formulations, the parameter (or its choice rule) of the formulation for the additive one.

== Weights from the optimality condition

The weights are such that $partial R slash partial v_(i j)^* = W_(i j) v_(i j)$, which turns the non-quadratic penalty into a sequence of quadratic ones:

$ "Lq": W_i = q/2 max(epsilon, |v_i|)^(q-2), $
$ "Lpq": W_(i j) = q/2 max(epsilon^p, S_i)^(q slash p - 1) max(epsilon, |v_(i j)|)^(p - 2), quad S_i = sum_j |v_(i j)|^p. $

The damping $epsilon$ avoids divisions by zero for $q < 2$; it is set once and for all to the 5 % quantile of $|L U^((0))|$, per group. For quadratic terms, $W = 1$: with the additive formulation, the initial solution is the final one and a single linear solve is performed.

= Resolution strategies

The rows of the unknown are the points and the columns the samples. The strategy is chosen by multiple dispatch on the operator and on the coupling induced by the regularization terms:

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Case*], [*Strategy*], [*Resolution*]),
  [space–frequency, `Lq` only], [per sample], [one independent problem per frequency, in parallel (own $alpha$, $lambda$ and convergence)],
  [space–frequency, `Lpq` without `Lc`], [joint (blocks)], [common weights and parameter, but $H^H H + D$ block-diagonal: frequency-by-frequency solves in parallel],
  [space–time (dense), without `Lc`], [joint (dense)], [dense system of size $n_r n_t$, operator $I times.o L$],
  [`Lpq` with `Lc`], [assembled], [automatically vectorized problem (sparse block-diagonal $H$, operator $L_c times.o L$), sparse Cholesky],
)

Without `Lc`, the block-diagonal structure is always exploited and the frequencies are only coupled through the weights. With `Lc`, the term $L^H D L$ links neighbouring frequencies (block-tridiagonal system for a 3-point operator) and the problem is solved as a whole with a sparse factorization.

= Linear systems

== Regime and Woodbury identity

The regime of the problem is a type parameter of the solver cache, and the method is selected by dispatch:

#table(
  columns: (auto, auto, 1fr, auto),
  align: left,
  table.header([*Regime*], [*`L`*], [*Method*], [*System size*]),
  [$n_m >= n_r$], [any], [normal equations @eq-ir], [$n_r$],
  [$n_m < n_r$], [$I$], [Woodbury identity], [$n_m$],
  [$n_m < n_r$], [square, invertible], [Woodbury on $tilde(H) = H L^(-1)$], [$n_m$],
  [$n_m < n_r$], [rectangular], [normal equations ($L^H D L$ singular)], [$n_r$],
)

In the underdetermined regime, the push-through and Woodbury identities replace the $n_r times n_r$ system by an $n_m times n_m$ one:

$ (H^H H + D)^(-1) H^H y = D^(-1) H^H (I + H D^(-1) H^H)^(-1) y, $
$ (H^H H + L^H D L)^(-1) = L^(-1) (tilde(H)^H tilde(H) + D)^(-1) L^(-H). $

For sparse operators (assembled problems), CHOLMOD is used and the symbolic analysis is reused from one iteration to the next. No inverse is ever formed explicitly: $A^(-1)$ is always applied through a factorization.

== Spectral decompositions

The parameter choice rules require the SVD of $H W^(-1 slash 2)$ ($L = I$) or the GSVD of $(H, W^(1 slash 2) L)$. The components of the null space of $L$ are not regularized; they are accounted for in the criteria. The decompositions are computed on the standard-form problem, so that the cost per evaluation of a criterion is negligible.

== Memory and allocations

A linear system object is built once per resolution. It precomputes everything that does not depend on the weights ($H^H H$, transpose of the global smoothing operator, standard form $H L^(-1)$) and preallocates its buffers; the fixed-point loop then works in place. The solvers are compiled for the concrete types of the problem (regime, smoothing operator, real or complex field) thanks to function barriers.

#bibliography("/refs.bib", title: [References], style: "ieee")
