#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Release notes])
#metadata((description: "Changes between the versions of SourceEstimation.jl.")) <website-metadata>

#title()

This page lists the changes between the versions of SourceEstimation.jl. The package follows #link("https://semver.org")[semantic versioning]: before version 1.0, minor versions may contain breaking changes, which are always mentioned here.

= v0.1.0 — in development

#badges(
  badge[first release],
  badge(color: julia.purple)[unregistered],
)

First public version of the package. It gathers, behind a single interface `solve(problem, formulation, algorithm)`, the methods developed in the papers listed in the #link("../theory/references.html")[bibliography].

== Problems

- `InverseProblem` for single-frequency, space–frequency (vector of FRF matrices or 3D array), space–time (block convolution matrix) and assembled (sparse) problems.
- Real and complex data; known diagonal noise covariance (`σ² * I`, vector or matrix of variances, one variance per measurement point), handled by whitening.

== Regularization

- Regularization terms `Lq` ($0 < q <= 2$) and mixed norms `Lpq`, with an operator `L` on the rows and, for `Lpq`, an operator `Lc` on the columns (frequencies or time steps).
- Region-wise regularization `Groups`, with one smoothing operator per group.
- Finite-difference operators: `fd_weights`, `diff_operator` (1D grids and structured 2D/3D point clouds), `space_sample_operator`, `sample_rows`.
- Formulations: `Additive` (generalized Tikhonov, multi-parameter), `Multiplicative` (OMR, MPMR with `multiparameter = true`) and `Bayesian` (hierarchical model, each hyperparameter `Estimated()` or fixed).
- Parameter choice rules `LCurve`, `GCV`, `BayesEstimator`, computed on the SVD or GSVD of the weighted problem and refined by Brent's method; initialization rules `Heuristic` and `LambdaMax`.

== Algorithms

- `IR`: iteratively reweighted algorithm of the additive and multiplicative formulations; a single linear solve for quadratic additive problems.
- `MAP`: optimal Bayesian regularization; `ElementWise`: relevance vector regression.
- Linear solvers selected by dispatch: normal equations, Woodbury identity in the underdetermined regime (also on the standard form $H L^(-1)$), sparse Cholesky (CHOLMOD) for assembled problems, multithreaded frequency-by-frequency or block-wise resolution.

== Uncertainty quantification

- `solve(...; precision = true)` returns the posterior precision (rather than the covariance, which is never formed) and the posterior model.
- `posterior_std`, `posterior_covariance` (marginal blocks or full covariance), `posterior_samples` (Cholesky factor or perturbation–optimization), `credible_interval`.

== Documentation

- This website: getting started, five tutorials, the theory of the methods with their references, and the API reference generated from the docstrings.
