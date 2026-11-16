# Iteratively Reweighted (IR) fixed point, shared by the additive and multiplicative formulations (and by the u-step of the Bayesian MAP).
# At iteration k + 1, with V = L u⁽ᵏ⁾:
#   1. local weights W_g (first-order optimality condition) and penalties P_g = R_g(V_g)
#   2. regularization parameter(s) λ  (dispatch on formulation)
#   3. u⁽ᵏ⁺¹⁾ = (Hᴴ H + Lᴴ λW L)⁻¹ Hᴴ y (dispatch on H and L)

# Regularization parameters -------------------------------------------------------------

"""
    regularization_parameter(formulation, cache, R, U, W, P)

Regularization parameter(s) at the next IR iteration (dispatch on the formulation).

**Inputs**
- `formulation` - `Additive` or `Multiplicative`: formulation of the problem
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `R` - `Regularization`: assembled regularization
- `U` - `AbstractMatrix`: current iterate (`nr × ns`)
- `W` - `AbstractMatrix{<:Real}`: IR weights (size of `L U`)
- `P` - `AbstractVector{<:Real}`: penalties `R_g(L_g U_g)` of the groups

**Outputs**
- `λ` - `Real` or `Vector{<:Real}`: regularization parameter, or one per group (additive
  multi-parameter, MPMR)
"""
regularization_parameter(f::Additive{<:Any,<:Real}, c, R, U, W, P) = f.λ
regularization_parameter(f::Additive{<:Any,<:AbstractVector}, c, R, U, W, P) = f.λ
regularization_parameter(f::Additive{<:Any,<:ParameterChoice}, c, R, U, W, P) =
    select_parameter(f.λ, c, R.L, W)

# OMR: α = ‖y - Hu‖² / Σ_g R_g
regularization_parameter(::Multiplicative{<:Any,SumCoupling}, c, R, U, W, P) =
    residual_norm2(c, U) / sum(P)

# MPMR: α_g = ‖y - Hu‖² / R_g
regularization_parameter(::Multiplicative{<:Any,ProductCoupling}, c, R, U, W, P) =
    residual_norm2(c, U) ./ P

scale_weights!(D, W, λ::Real, R) = (D .= λ .* W)

function scale_weights!(D, W, λ::AbstractVector, R)
    length(λ) == ngroups(R) ||
        throw(DimensionMismatch("$(length(λ)) regularization parameters for $(ngroups(R)) groups"))
    for g in eachindex(R.rows)
        rows = R.rows[g]
        @views D[rows, :] .= λ[g] .* W[rows, :]
    end
    return D
end

# Fixed-point loop ---------------------------------------------------------------------------

"""
    reweighting(cache, R, formulation, U, ε; tol, maxiter)

Iteratively Reweighted (IR) fixed point started from `U`. The stopping criterion is the
relative `ℓ1` variation of the solution between two iterations.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `R` - `Regularization`: assembled regularization
- `formulation` - `Additive` or `Multiplicative`: rule updating the regularization parameter
- `U` - `AbstractMatrix`: starting point (`nr × ns`)
- `ε` - `AbstractVector{<:Real}`: damping parameter of each group
- `tol` - `Real`: tolerance of the stopping criterion
- `maxiter` - `Integer`: maximum number of iterations

**Outputs**
- `st` - `NamedTuple`: `U` (solution), `λ` (regularization parameter(s)), `W` (weights), `D`
  (scaled weights `λ W` of the last linear solve), `converged`, `history`
"""
function reweighting(c, R::Regularization, f::Formulation, U, ε; tol, maxiter)
    maxiter >= 1 || throw(ArgumentError("maxiter must be ≥ 1"))
    T = real(eltype(c))
    W = zeros(T, R.nout, c.ns)
    D = similar(W)
    P = zeros(T, ngroups(R))
    history = Float64[]
    local λ
    converged = false
    for _ in 1:maxiter
        reweight!(W, P, R, apply(R.L, U), ε)
        λ = regularization_parameter(f, c, R, U, W, P)
        scale_weights!(D, W, λ, R)
        Unew = weighted_solve(c, R.L, D)
        crit = sum(abs, Unew - U) / max(sum(abs, U), floatmin(T))
        push!(history, crit)
        U = Unew
        if crit < tol
            converged = true
            break
        end
    end
    return (; U, λ, W, D, converged, history)
end

"""
    initial_solution(cache, R, init)

Tikhonov initial solution `(Hᴴ H + α₀ Lᴴ L)⁻¹ Hᴴ y`.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `R` - `Regularization`: assembled regularization
- `init` - `Real` or `ParameterChoice`: initial regularization parameter, or rule selecting it

**Outputs**
- `U₀` - `Matrix`: initial solution (`nr × ns`)
- `α₀` - `Real`: initial regularization parameter
"""
function initial_solution(c, R::Regularization, init)
    W = ones(real(eltype(c)), R.nout, c.ns)
    α₀ = select_parameter(init, c, R.L, W)
    return weighted_solve(c, R.L, α₀ .* W), α₀
end

default_init(f::Additive) = f.λ
default_init(::Additive{<:RegularizationTerm,<:AbstractVector}) = LCurve()
default_init(::Multiplicative) = Heuristic()

"""
    setup(prob, prior, R)

Setup shared by all solvers: linear cache of the problem whitened by the noise covariance,
and assembled regularization.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `prior` - `RegularizationTerm`: regularization term of the formulation
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided

**Outputs**
- `c` - `DenseCache` or `BlockCache`: linear cache
- `reg` - `Regularization`: assembled regularization
"""
function setup(prob::InverseProblem, prior::RegularizationTerm, R)
    c = linear_cache(prob, R)
    reg = Regularization(groups(prior, c.nr), c.nr)
    return c, reg
end

"""
    solve(prob, formulation[, algorithm]; R = nothing, covariance = false, tol = 1e-4, maxiter = 100, init)

Solve the inverse problem `prob` with the given regularization `formulation`. When the
algorithm is omitted, [`default_algorithm`](@ref) is used.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `formulation` - `Formulation`: `Additive`, `Multiplicative` or `Bayesian`
- `algorithm` - `Algorithm`: resolution algorithm (`IR()`, `MAP()`)
- `R` - `Nothing`, `UniformScaling`, `Diagonal`, `AbstractVector` or `AbstractMatrix`: known
  diagonal noise covariance — not provided (default, equivalent to `I`), `σ²I`, or the noise
  variances (vector, or matrix with the layout of `y`). The data fidelity is `‖y - Hu‖²_{R⁻¹}`
- `covariance` - `Bool`: also compute the posterior covariance of the solution. For the
  additive and multiplicative formulations, `R` is then taken as the exact noise covariance
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Real` or `ParameterChoice`: initial Tikhonov solution — initial regularization
  parameter, or rule selecting it (`Heuristic()`, `LCurve()`, `GCV()`, `BayesEstimator()`)

**Outputs**
- `sol` - `InverseSolution`: solution, regularization parameter(s), hyperparameters,
  convergence information and, if requested, posterior covariance
"""
solve(prob::InverseProblem, f::Formulation; kwargs...) =
    solve(prob, f, default_algorithm(f); kwargs...)

function solve(prob::InverseProblem, f::Union{Additive,Multiplicative}, ::IR;
               R = nothing, covariance::Bool = false, tol = 1e-4, maxiter = 100, init = default_init(f))
    c, reg = setup(prob, f.prior, R)
    U₀, _ = initial_solution(c, reg, init)
    ε = dampings(reg, apply(reg.L, U₀))
    st = reweighting(c, reg, f, U₀, ε; tol, maxiter)
    Σ = covariance ? posterior_covariance(c, reg.L, st.D) : nothing
    return InverseSolution(output_shape(prob, st.U), st.λ, (;),
                           length(st.history), st.converged, st.history, Σ)
end
