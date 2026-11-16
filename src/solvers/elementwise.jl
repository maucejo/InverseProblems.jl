# Element-wise Bayesian algorithm (Relevance Vector Regression) ----------------------------
#
# Hierarchical model with one Gaussian prior precision per component of v = L u:
#     p(v_i | τ_i) ∝ exp(-φ τ_i |v_i|²),   τ_i ~ Gamma(a, b)       (φ = 1/2 real, 1 complex)
# The MAP alternates the closed-form updates
#     τ_i = (k/2 + a - 1) / (φ |v_i|² + b),   τn (if estimated),
#     u   = (Hᴴ H + Lᴴ Λ L)⁻¹ Hᴴ y,   Λ = φ τ / (c τn)
# (normal equations or Woodbury identity, see `linsolve.jl`).

"""
    lambda_max(cache)

Initial Tikhonov parameter of the element-wise algorithm, `‖2 Hᴴ y‖∞` (computed per sample
for a block-diagonal operator).

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem

**Outputs**
- `λ` - `Real` (dense operator) or `Vector{<:Real}` (block-diagonal operator, one per sample)
"""
lambda_max(c::DenseCache) = norm(2 .* c.Hty, Inf)
lambda_max(c::BlockCache) = [lambda_max(b) for b in c.blocks]

# Initial weights D (size of L U): one value, or one value per sample (column)
initial_weights(c::DenseCache, reg, λ::Real) = fill(float(λ), reg.nout, c.ns)
initial_weights(c::BlockCache, reg, λ::AbstractVector) = repeat(transpose(float.(λ)), reg.nout, 1)
initial_weights(c::BlockCache, reg, λ::Real) = fill(float(λ), reg.nout, c.ns)

"""
    solve(prob, f::Bayesian, ElementWise(); R = nothing, covariance = false, tol = 1e-8,
          maxiter = 200, init = nothing)

Element-wise Bayesian regularization (Relevance Vector Regression): MAP of a hierarchical
model where each component of `L u` has its own Gaussian prior precision `τ_i`, alternating
the updates of the precisions, of `τn` (if `Estimated()`) and of `u`. Only the smoothing
operators of `f.prior` are used: the prior is Gaussian (`q = 2`), so all its terms must be
`Lq(2.0)`, and the signal precisions must be `Estimated()`.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `f` - `Bayesian`: Bayesian formulation (`τn` fixed or estimated, as for `MAP()`)
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided
- `covariance` - `Bool`: also compute the posterior covariance `(τn (H̃ᴴH̃ + Lᴴ Λ L))⁻¹`
  (`H̃ = R^{-1/2} H`)
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Nothing`, `Real` or `ParameterChoice`: initial Tikhonov parameter — by default
  `‖2 Hᴴ y‖∞` (per sample for a block-diagonal operator)

**Outputs**
- `sol` - `InverseSolution`: MAP estimate, element-wise regularization parameters
  `Λ = φ τ / (c τn)` (layout of `L u`) in `sol.λ`, hyperparameters `(τn, τs)` in `sol.hyper`
  (`τs`: element-wise precisions), convergence information and, if requested, posterior
  covariance
"""
function solve(prob::InverseProblem, f::Bayesian, ::ElementWise; R = nothing,
               covariance::Bool = false, tol = 1e-8, maxiter = 200, init = nothing)
    f.τs isa Estimated ||
        throw(ArgumentError("the element-wise algorithm estimates the signal precisions: use τs = Estimated()"))
    c, reg = setup(prob, f.prior, R)
    all(t -> t isa Lq && t.q == 2, reg.terms) ||
        throw(ArgumentError("the element-wise algorithm uses Gaussian priors: all terms must be Lq(2.0)"))
    F = field(eltype(c))
    N = ndata(c)
    τnspec = noise_spec(f.τn, R)

    U = weighted_solve(c, reg.L, initial_weights(c, reg, initial_lambda(init, c, reg)))
    τn = update_noise_precision(τnspec, F, N, residual_norm2(c, U))
    τs = elementwise_precision(F, apply(reg.L, U))
    D = prior_factor(F, 2.0) .* τs ./ (misfit_factor(F) * τn)

    history = Float64[]
    converged = false
    for _ in 1:maxiter
        τs = elementwise_precision(F, apply(reg.L, U))
        τn = update_noise_precision(τnspec, F, N, residual_norm2(c, U))
        D = prior_factor(F, 2.0) .* τs ./ (misfit_factor(F) * τn)
        Unew = weighted_solve(c, reg.L, D)
        crit = sum(abs, Unew - U) / max(sum(abs, U), floatmin(Float64))
        push!(history, crit)
        U = Unew
        if crit < tol
            converged = true
            break
        end
    end
    Σ = covariance ? scale_covariance(posterior_covariance(c, reg.L, D), τn) : nothing
    hyper = (τn = τn, τs = τs)
    return InverseSolution(output_shape(prob, U), D, hyper, length(history), converged, history, Σ)
end

# Initial Tikhonov parameter: ‖2 Hᴴ y‖∞ by default, a fixed value, or a parameter choice rule
initial_lambda(::Nothing, c, reg) = lambda_max(c)
initial_lambda(λ::Real, c, reg) = float(λ)
initial_lambda(rule::ParameterChoice, c, reg) =
    select_parameter(rule, c, reg.L, ones(real(eltype(c)), reg.nout, c.ns))
