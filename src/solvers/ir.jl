# Iteratively Reweighted (IR) fixed point, shared by the additive and multiplicative formulations (and by the u-step of the Bayesian MAP).
# At iteration k + 1, with V = L u⁽ᵏ⁾:
#   1. local weights W_g (first-order optimality condition) and penalties P_g = R_g(V_g)
#   2. regularization parameter(s) λ  (dispatch on formulation)
#   3. u⁽ᵏ⁺¹⁾ = (Hᴴ H + Lᴴ λW L)⁻¹ Hᴴ y (dispatch on H and L)

# Regularization parameters -------------------------------------------------------------

"""
    regularization_parameter(formulation, cache, R, U, W, P)

Regularization parameter(s) at the next iteration, given the current iterate `U`, the
weights `W` and the group penalties `P`.
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

IR fixed point started from `U`. The stopping criterion is the relative `ℓ1` variation of
the solution between two iterations.
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
Tikhonov initial solution `(Hᴴ H + α₀ Lᴴ L)⁻¹ Hᴴ y`, `α₀` given by `init`.
Returns the solution, `α₀` and the unit weights.
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
Setup shared by all solvers: (whitened) linear cache for the noise covariance `R`, and
assembled regularization.
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

# Keywords
* `R`: noise covariance, known and diagonal — not provided (default, equivalent to `I`), `σ²I`, a `Diagonal`, or the
  vector (or matrix with the layout of `y`) of the noise variances. The data fidelity is
  `‖y - Hu‖²_{R⁻¹}`.
* `covariance`: also compute the posterior covariance of the solution (see
  [`posterior_std`](@ref), [`credible_interval`](@ref)). For the additive and multiplicative
  formulations, `R` is then taken as the exact noise covariance.
* `tol`, `maxiter`: stopping criterion (relative `ℓ1` variation of `u`) and iteration limit.
* `init`: initial Tikhonov solution — a number (initial regularization parameter) or a
  [`ParameterChoice`](@ref) (`Heuristic()`, `LCurve()`, `GCV()`, `BayesEstimator()`).
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
