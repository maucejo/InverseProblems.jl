# Element-wise Bayesian algorithm (Relevance Vector Regression) ----------------------------
#
# Hierarchical model with one Gaussian prior precision per component of v = L u:
#     p(v_i | τ_i) ∝ exp(-φ τ_i |v_i|²),   τ_i ~ Gamma(a, b)       (φ = 1/2 real, 1 complex)
# The MAP alternates the closed-form updates
#     τ_i = (k/2 + a - 1) / (φ |v_i|² + b),   τn (if estimated),
#     u   = (Hᴴ H + Lᴴ Λ L)⁻¹ Hᴴ y,   Λ = φ τ / (c τn) = τ / τn
# (normal equations or Woodbury identity, see `linsolve.jl`).


# Initial weights D (size of L U): one value, or one value per sample (column)
initial_weights(c::Union{DenseCache,SparseCache}, reg, λ::Real) = fill(float(λ), reg.nout, c.ns)
initial_weights(c::BlockCache, reg, λ::AbstractVector) = repeat(transpose(float.(λ)), reg.nout, 1)
initial_weights(c::BlockCache, reg, λ::Real) = fill(float(λ), reg.nout, c.ns)

"""
    solve_joint(prob, f::Bayesian, ElementWise(); R = nothing, precision = false, tol = 1e-8, maxiter = 200, init = nothing)

Element-wise Bayesian regularization (Relevance Vector Regression): MAP of a hierarchical model where each component of `L u` has its own Gaussian prior precision `τ_i`, alternating the updates of the precisions, of `τn` (if `Estimated()`) and of `u`. Only the smoothing operators of `f.prior` are used: the prior is Gaussian (`q = 2`), so all its terms must be `Lq(2.0)`, and the signal precisions must be `Estimated()`.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `f` - `Bayesian`: Bayesian formulation (`τn` fixed or estimated, as for `MAP()`)
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided
- `precision` - `Bool`: also return the posterior precision `τn (H̃ᴴH̃ + Lᴴ Λ L)` and the posterior
  model (see [`solve`](@ref))
  (`H̃ = R^{-1/2} H`)
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Nothing`, `Real` or `ParameterChoice`: initial Tikhonov parameter — by default `‖2 Hᴴ y‖∞` (per sample for a block-diagonal operator)

**Outputs**
- `sol` - `InverseSolution`: MAP estimate, element-wise regularization parameters `Λ = τ / τn` (layout of `L u`; MAP `½‖y - Hu‖²_R + ½ Σᵢ Λᵢ |vᵢ|²` for real data, `‖y - Hu‖²_R + Σᵢ Λᵢ |vᵢ|²` for complex data) in `sol.λ`, hyperparameters `(τn, τs)` in `sol.hyper` (`τs`: element-wise precisions), convergence information and, if requested, posterior precision
"""
function solve_joint(prob::InverseProblem, f::Bayesian, ::ElementWise; R = nothing,
               precision = false, tol = 1e-8, maxiter = 200, init = nothing)
    f.τs isa Estimated ||
        throw(ArgumentError("the element-wise algorithm estimates the signal precisions: use τs = Estimated()"))
    c, reg = setup(prob, f.prior, R)
    all(t -> t isa Lq && t.q == 2, reg.terms) ||
        throw(ArgumentError("the element-wise algorithm uses Gaussian priors: all terms must be Lq(2.0)"))
    return elementwise_solve(prob, c, reg, noise_spec(f.τn, R), posterior_request(precision), tol, maxiter, init)
end

# Function barrier (see `ir_solve`)
function elementwise_solve(prob, c, reg, τnspec, request, tol, maxiter, init)
    F = field(eltype(c))
    N = ndata(c)

    sys = linear_system(c, reg.L)
    U = Matrix{eltype(c)}(undef, c.nr, c.ns)
    Unew = similar(U)
    weighted_solve!(U, sys, initial_weights(c, reg, initial_lambda(init, c, reg)))
    Vbuf = smoothing_buffer(reg.L, U)
    τs = zeros(real(eltype(c)), reg.nout, c.ns)
    D = similar(τs)
    τn = 1.0

    history = sizehint!(Float64[], maxiter)
    converged = false
    for _ in 1:maxiter
        elementwise_precision!(τs, F, apply!(Vbuf, reg.L, U))
        τn = update_noise_precision(τnspec, F, N, residual_norm2!(sys, U))
        D .= prior_factor(F, 2.0) .* τs ./ (misfit_factor(F) * τn)
        weighted_solve!(Unew, sys, D)
        crit = relative_change(Unew, U)
        push!(history, crit)
        U, Unew = Unew, U
        if crit < tol
            converged = true
            break
        end
    end
    Q, post = posterior_output(request, c, reg.L, D, τn)
    hyper = (τn = τn, τs = τs)
    return InverseSolution(output_shape(prob, U), D, hyper, length(history), converged, history, Q, post)
end

# Initial Tikhonov parameter: ‖2 Hᴴ y‖∞ by default, a fixed value, or a parameter choice rule
initial_lambda(::Nothing, c, reg) = lambda_max(c)
initial_lambda(λ::Real, c, reg) = float(λ)
initial_lambda(rule::ParameterChoice, c, reg) =
    select_parameter(rule, c, reg.L, ones(real(eltype(c)), reg.nout, c.ns))
