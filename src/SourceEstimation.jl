"""
    SourceEstimation

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
* Ill-posedness:    [`ill_posedness`](@ref) (discrete Picard condition, conditioning)
* Visualization:    [`picardplot`](@ref) (Makie extension: load a Makie backend first)
* Performance:      [`relative_error`](@ref), [`correlation`](@ref), [`precision_metric`](@ref), [`accuracy_metric`](@ref), [`performance_ratio`](@ref)
"""
module SourceEstimation

using LinearAlgebra
using Random: Random, AbstractRNG
using SparseArrays: SparseArrays, AbstractSparseMatrix, SparseMatrixCSC, sparse, spdiagm, blockdiag, findnz, rowvals, nonzeros, nzrange
using Statistics: quantile
using SpecialFunctions: loggamma

# Types ------------------------------------------------------------------------------------
include("types/problem.jl")      # forward operators, InverseProblem
export InverseProblem

include("types/priors.jl")       # regularization terms / priors
export Lq, Lpq, Groups

include("types/choices.jl")      # regularization parameter choice rules
export LCurve, GCV, BayesEstimator, Heuristic, LambdaMax

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
include("utils/sparsesolve.jl") # sparse operators: sparse Cholesky (CHOLMOD)
include("utils/minimize.jl")         # 1D minimization: Brent, refined grid search
include("utils/spectrum.jl")  # SVD / GSVD of the standard-form problem
include("utils/regparam.jl")  # L-curve, GCV, Bayesian estimator, heuristic
include("utils/weights.jl")   # IR weights (reweight!), damping ε, penalty values
include("utils/bayes.jl")     # real / complex Bayesian model

include("utils/fdiff.jl")    # finite-difference smoothing operators
export fd_weights, diff_operator, space_sample_operator, sample_rows

include("utils/posterior.jl")  # posterior precision, perturbation–optimization, block solves
include("utils/uncertainty.jl") # posterior samples and credible intervals
export posterior_std, posterior_covariance, posterior_samples, credible_interval

include("indicators/performance.jl") # trueness, precision and accuracy of a reconstruction
export relative_error, correlation, precision_metric, accuracy_metric, performance_ratio

include("indicators/ill_posedness.jl") # discrete Picard condition, conditioning, SVD decay
export ill_posedness, IllPosedness

# Visualization (implemented in the Makie extension, ext/SourceEstimationMakieExt.jl) -----------
include("visualization/plots.jl")
export picardplot, picardplot!

# Solvers -------------------------------------------------------------------------------------
include("solvers/ir.jl") # IR algorithm: additive and multiplicative formulations

include("solvers/map.jl") # MAP: hierarchical Bayesian formulation
include("solvers/elementwise.jl") # element-wise Bayesian (RVR)
include("solvers/strategy.jl")    # resolution strategy: joint / frequency by frequency
export solve

function __init__()
    Base.Experimental.register_error_hint(plot_error_hint, MethodError)
end

end
