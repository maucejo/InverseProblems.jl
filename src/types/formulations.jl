# Formulations of the regularized problem ----------------------------------------

"""
    Formulation

How data fidelity and regularization are combined:
[`Additive`](@ref), [`Multiplicative`](@ref) or [`Bayesian`](@ref).
"""
abstract type Formulation end

# Additive ------------------------------------------------------------------------

"""
    Additive(prior; λ = LCurve())

Additive regularization

    û = argmin ‖y - Hu‖₂² + Σ_g λ_g R_g(L_g u_g)

`λ` can be

* a number: fixed regularization parameter;
* a vector (one entry per group): multi-parameter additive regularization;
* a [`ParameterChoice`](@ref) (`LCurve()`, `GCV()`, `BayesEstimator()`): selected at each
  iteration of the reweighting algorithm.

`prior` is a [`RegularizationTerm`](@ref): `Lq(2.0)` gives the Tikhonov regularization,
`Lq(1.0)` the LASSO, [`Groups`](@ref) the group generalized Tikhonov regularization.
"""
struct Additive{R<:RegularizationTerm, P} <: Formulation
    prior::R
    λ::P
end

function Additive(prior::RegularizationTerm; λ = LCurve())
    check_lambda(λ)
    return Additive(prior, λ)
end

check_lambda(::ParameterChoice) = nothing
check_lambda(::Heuristic) = throw(ArgumentError("Heuristic() is only an initialization rule (`init` keyword of `solve`)"))
check_lambda(λ::Real) = λ > 0 || throw(ArgumentError("λ must be positive"))
check_lambda(λ::AbstractVector{<:Real}) = all(>(0), λ) || throw(ArgumentError("λ must be positive"))

# Multiplicative ------------------------------------------------------------------

"""Combination of local regularization terms in the multiplicative formulation."""
abstract type Coupling end
struct SumCoupling <: Coupling end      # ordinary multiplicative regularization
struct ProductCoupling <: Coupling end  # multi-parameter multiplicative regularization

"""
    Multiplicative(prior; multiparameter = false)

Multiplicative regularization — no regularization parameter to tune beforehand:

* `multiparameter = false` (default): ordinary multiplicative regularization (OMR)

      û = argmin ‖y - Hu‖₂² · Σ_g R_g(L_g u_g)

* `multiparameter = true`: multi-parameter multiplicative regularization (MPMR)

      û = argmin ‖y - Hu‖₂² · Π_g R_g(L_g u_g)

  which yields one adaptive regularization parameter per group.
"""
struct Multiplicative{R<:RegularizationTerm, C<:Coupling} <: Formulation
    prior::R
end

Multiplicative(prior::RegularizationTerm; multiparameter::Bool = false) =
    multiparameter ? Multiplicative{typeof(prior), ProductCoupling}(prior) :
                     Multiplicative{typeof(prior), SumCoupling}(prior)

# Bayesian ------------------------------------------------------------------------

"""
    Estimated()

Marks a hyperparameter of the [`Bayesian`](@ref) formulation as optimized jointly with the
solution (MAP of its full conditional distribution).
"""
struct Estimated end

"""
    Fixed()

Keeps the shape parameters `q_g` of the [`Bayesian`](@ref) formulation at the values given
in the prior (`Lq(q_g)`).
"""
struct Fixed end

"""
    Bayesian(prior; τn = default, τs = Estimated(), q = Estimated(), shape_bounds = (0.05, 2.0))

Hierarchical Bayesian formulation: normal likelihood of covariance `R/τn` (`R`: noise
covariance given to [`solve`](@ref)) and, for each group, a generalized Gaussian prior on
`L_g u_g` with shape `q_g` and precision `τs_g`:

    p(u, τn, τs, q | y) ∝ p(y | u, τn) Π_g p(L_g u_g | τs_g, q_g) p(τs_g) p(q_g) p(τn)

Each hyperparameter is either optimized (`Estimated()`) or fixed:

* `τn`: a number or `Estimated()`. `τn = 1` means that `R` is the exact noise covariance;
  with `Estimated()`, `R` is only known up to a scale factor. By default, `τn = 1` when `R`
  is given to [`solve`](@ref), and `τn` is estimated otherwise;
* `τs`: a number (all groups) or a vector (one value per group);
* `q`: `Fixed()` keeps the shapes given in `prior` (`Lq(q_g)`); with `Estimated()`, they are
  the starting point of the search in `shape_bounds`.

Real and complex data lead to different distributions (real or circular complex normal and
generalized Gaussian laws). Weakly informative Gamma hyperpriors are used for the estimated
hyperparameters, as recommended in the papers.
"""
struct Bayesian{R<:RegularizationTerm, N, S, Q<:Union{Estimated,Fixed}} <: Formulation
    prior::R
    τn::N
    τs::S
    q::Q
    shape_bounds::Tuple{Float64,Float64}
end

function Bayesian(prior::RegularizationTerm; τn = nothing, τs = Estimated(),
                  q::Union{Estimated,Fixed} = Estimated(), shape_bounds = (0.05, 2.0))
    τn isa AbstractVector && throw(ArgumentError("τn must be `Estimated()` or a positive number"))
    check_precision(τn)
    check_precision(τs)
    lb, ub = shape_bounds
    0 < lb < ub || throw(ArgumentError("invalid shape bounds $shape_bounds"))
    return Bayesian(prior, τn, τs, q, (Float64(lb), Float64(ub)))
end

check_precision(::Estimated) = nothing
check_precision(::Nothing) = nothing   # τn: default, resolved by `solve` from R
check_precision(τ::Real) = τ > 0 || throw(ArgumentError("precision parameters must be positive"))
check_precision(τ::AbstractVector{<:Real}) =
    all(>(0), τ) || throw(ArgumentError("precision parameters must be positive"))
check_precision(τ) = throw(ArgumentError("a precision parameter must be `Estimated()` or a positive number, got $τ"))
