# InverseProblems.jl

Regularized resolution of linear inverse problems `y = H u + n` in the time and frequency
domains (force reconstruction, source identification).

```julia
using InverseProblems

prob = InverseProblem(H, y)                                  # H: matrix, y: vector

# Additive, multiplicative and Bayesian regularizations
sol = solve(prob, Additive(Lq(1.0); λ = BayesEstimator()))
sol = solve(prob, Multiplicative(Groups(inner => Lq(0.5), edges => Lq(2.0; L = D))))  # one smoothing operator per group
sol = solve(prob, Multiplicative(Groups(inner => Lq(0.5), edges => Lq(2.0)); multiparameter = true))
sol = solve(prob, Bayesian(Lq(1.0)))                         # optimal Bayesian regularization

# Space–frequency problem with a mixed norm: Hs[j] is the FRF matrix at frequency j
sol = solve(InverseProblem(Hs, Y), Multiplicative(Lpq(2, 0.5)))

# Known diagonal noise covariance, posterior covariance and credible intervals
sol = solve(prob, Bayesian(Lq(1.0)); R = σ², covariance = true)  # σ²: noise variances
lo, up = credible_interval(sol; level = 0.95)

sol.u        # reconstructed field
sol.λ        # regularization parameter(s)
sol.hyper    # other estimated quantities (τn, τs, q for the Bayesian formulation)
```

See [docs/design/api.md](docs/design/api.md) for the design of the API.
