"""
    InverseProblems

Regularized resolution of linear inverse problems `y = H u + n` in the time and
frequency domains (force reconstruction, source identification, …).

The whole public API boils down to one verb and four nouns:

```julia
prob = InverseProblem(H, y)                         # what to solve
reg  = Multiplicative(Groups(r1 => Lq(0.5),         # how to regularize
                             r2 => Lq(2.0; L = D)))
sol  = solve(prob, reg)                             # default algorithm
sol.u                                               # reconstructed field
```

* Problem:          [`InverseProblem`](@ref)
* Regularization terms (priors): [`Lq`](@ref), [`Lpq`](@ref), [`Groups`](@ref)
* Formulations:     [`Additive`](@ref), [`Multiplicative`](@ref), [`Bayesian`](@ref)
* Parameter choice: [`LCurve`](@ref), [`GCV`](@ref), [`BayesEstimator`](@ref), [`Heuristic`](@ref)
* Algorithms:       [`IR`](@ref), [`MAP`](@ref)
* Result:           [`InverseSolution`](@ref)
"""
module InverseProblems

using LinearAlgebra
using Random: Random, AbstractRNG
using SparseArrays: SparseMatrixCSC, sparse, findnz, rowvals, nonzeros, nzrange
using Statistics: quantile
using SpecialFunctions: loggamma

# Types ------------------------------------------------------------------------------------
include("types/problem.jl")      # forward operators, InverseProblem
export InverseProblem

include("types/priors.jl")       # regularization terms / priors
export Lq, Lpq, Groups

include("types/choices.jl")      # regularization parameter choice rules
export LCurve, GCV, BayesEstimator, Heuristic

include("types/formulations.jl") # additive, multiplicative, Bayesian formulations
export Additive, Multiplicative, Bayesian, Estimated, Fixed

include("types/algorithms.jl")   # resolution algorithms
export IR, MAP, ElementWise

include("types/solution.jl")     # result type
export InverseSolution

# Utilities (internal) ----------------------------------------------------------------------
include("utils/smoothing.jl") # assembled smoothing operator of the groups
include("utils/noise.jl")     # noise covariance, whitening
include("utils/linsolve.jl")  # caches, weighted Tikhonov solves, posterior covariance
include("utils/minimize.jl")         # 1D minimization: Brent, refined grid search
include("utils/spectrum.jl")  # SVD / GSVD of the standard-form problem
include("utils/regparam.jl")  # L-curve, GCV, Bayesian estimator, heuristic
include("utils/norms.jl")     # norm values, damping, IR weights
include("utils/bayes.jl")     # real / complex Bayesian model

include("utils/uncertainty.jl") # posterior samples and credible intervals
export posterior_std, posterior_samples, credible_interval

# Solvers -------------------------------------------------------------------------------------
include("solvers/ir.jl") # IR algorithm: additive and multiplicative formulations
export solve

include("solvers/map.jl") # MAP: hierarchical Bayesian formulation
include("solvers/elementwise.jl") # element-wise Bayesian (RVR)

end
