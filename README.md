# InverseProblems.jl

[![Generic badge](https://img.shields.io/badge/Version-0.1.0-cornflowerblue.svg)]()
[![MIT License](https://img.shields.io/badge/License-MIT-forestgreen)](LICENSE)

`InverseProblems.jl` is a Julia package designed to solve regularized linear inverse problems of the form `y = H u + n` in the time and frequency domains. It is primarily intended for force reconstruction and source identification in structural dynamics and acoustics, but it can be used for any linear inverse problem. This package can be used for educational and research purposes. As a teacher, you can use it to illustrate regularization theory and its practical implementation. As a researcher, you can use it to compare regularization strategies on your own data.

> [!WARNING]
> The package is under active development and is developed as part of my research activities on source identification in structural dynamics and acoustics.

## Features

The package provides the following features:

### Problem definition

- **Problem types**
    - Single frequency problem (matrix `H`, vector `y`)
    - Space–frequency problem (one FRF matrix per frequency)
    - Space–time problem (`nm·nt × nr·nt` operator)

- **Data**
    - Real and complex data
    - Known diagonal noise covariance (scalar, diagonal or variance array)

### Regularization

- **Regularization terms (priors)**
    - `Lq` norm, `0 < q ≤ 2`, with an optional smoothing operator `L`
    - `Lpq` mixed norm (space–frequency or space–time sparsity)
    - Group (region-wise) regularization with one smoothing operator per group

- **Formulations**
    - Additive regularization (generalized Tikhonov)
        - Single parameter
        - Multi-parameter
    - Multiplicative regularization
        - Optimal Multiplicative Regularization (OMR)
        - Multi-Parameter Multiplicative Regularization (MPMR)
    - Bayesian regularization
        - Optimal Bayesian regularization (hierarchical model)
        - Element-wise Bayesian regularization (Relevance Vector Regression)
        - Estimated or fixed hyperparameters (noise precision, prior precision, prior shape)

- **Regularization parameter choice rules**
    - L-curve
    - Generalized Cross Validation (GCV)
    - Bayesian estimator
    - Heuristic (initialization)
    - Local refinement using Brent's method

### Resolution algorithms

- **Algorithms**
    - Iteratively reweighted (IR) algorithm
    - Maximum a posteriori (MAP) estimation
    - Element-wise Bayesian algorithm

- **Linear solvers**
    - Normal equations for overdetermined problems
    - Woodbury identity for underdetermined problems
    - SVD / GSVD of the standard-form problem
    - Multithreaded block-wise resolution for space–frequency problems

### Uncertainty quantification

- Posterior covariance (Gaussian approximation at the fixed point)
- Marginal posterior standard deviations
- Posterior samples
- Credible intervals

## Quick start

```julia
using InverseProblems

prob = InverseProblem(H, y)                     # H: matrix, y: vector

# Additive, multiplicative and Bayesian regularizations
sol = solve(prob, Additive(Lq(1.0); λ = BayesEstimator()))
sol = solve(prob, Multiplicative(Groups(inner => Lq(0.5), edges => Lq(2.0; L = D))))
sol = solve(prob, Multiplicative(Groups(inner => Lq(0.5), edges => Lq(2.0)); multiparameter = true))
sol = solve(prob, Bayesian(Lq(1.0)))

# Space–frequency problem with a mixed norm: Hs[j] is the FRF matrix at frequency j
sol = solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5)))

# Known diagonal noise covariance, posterior covariance and credible intervals
sol = solve(prob, Bayesian(Lq(1.0)); R = σ², covariance = true)
lo, up = credible_interval(sol; level = 0.95)

sol.u        # reconstructed field
sol.λ        # regularization parameter(s)
sol.hyper    # other estimated quantities (τn, τs, q for the Bayesian formulation)
```

See [docs/design/api.md](docs/design/api.md) for the design of the API.

## License

MIT licensed

Copyright (C) 2026 Mathieu AUCEJO (maucejo)
