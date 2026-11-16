# Optimal Bayesian regularization — MAP solver (Aucejo & De Smet, MSSP 2019) ------------

# Estimated or fixed hyperparameters ---------------------------------------------------------
#
# Each update dispatches on the specification given in `Bayesian`: `Estimated()` → mode of # the full conditional distribution, number → kept at its value.
# Initialization (MAP of the standard formulation with a single parameter λ0, then most probable precisions given λ0 — marginal likelihood, Pereira et al.): τn = N / yᴴ(I + H (Lᴴλ0WL)⁻¹ Hᴴ)⁻¹ y,   τs_g = λ0 c τn / φ_g
# Default noise precision: estimated when R is not provided, 1 (R is the exact noise covariance) when it is.
noise_spec(::Nothing, ::Nothing) = Estimated()
noise_spec(::Nothing, R) = 1.0
noise_spec(τn, R) = τn

initial_noise_precision(::Estimated, c, reg, λ0, W, N) =
    N / marginal_quadratic(c, weighted_solve(c, reg.L, λ0 .* W))
initial_noise_precision(τ::Real, c, reg, λ0, W, N) = float(τ)

initial_signal_precision(::Estimated, τ0, G) = collect(float.(τ0))
initial_signal_precision(τ::Real, τ0, G) = fill(float(τ), G)
function initial_signal_precision(τ::AbstractVector{<:Real}, τ0, G)
    length(τ) == G || throw(DimensionMismatch("$(length(τ)) signal precisions for $G groups"))
    return float.(collect(τ))
end

update_shape(::Estimated, F, v, τ, q, bounds) = optimal_shape(F, v, τ, bounds)
update_shape(::Fixed, F, v, τ, q, bounds) = q

update_signal_precision(::Estimated, F, v, q, τ) = signal_precision(F, v, q)
# fixed τs (one value for all groups, or one per group): kept at its initial value
update_signal_precision(::Union{Real,AbstractVector{<:Real}}, F, v, q, τ) = τ

update_noise_precision(::Estimated, F, N, r2) = noise_precision(F, N, r2)
update_noise_precision(τ::Real, F, N, r2) = float(τ)

"""
    solve(prob, f::Bayesian, MAP(); R = nothing, covariance = false, tol = 1e-4, maxiter = 100,
          inner_tol = 1e-4, inner_maxiter = 100, init = LCurve())

Optimal Bayesian regularization: MAP of the hierarchical model, alternating the updates of `q_g`, `τs_g`, `τn` (those marked `Estimated()`) and of `u` (IR algorithm).

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `f` - `Bayesian`: Bayesian formulation
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided
- `covariance` - `Bool`: also compute the posterior covariance `(τn (H̃ᴴH̃ + LᴴDL))⁻¹` of the
  Gaussian approximation at the MAP (`H̃ = R^{-1/2} H`)
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u` (outer loop)
- `maxiter` - `Integer`: maximum number of outer iterations
- `inner_tol` - `Real`: tolerance of the IR algorithm (u-step)
- `inner_maxiter` - `Integer`: maximum number of IR iterations
- `init` - `Real` or `ParameterChoice`: rule for the initial single-parameter solution

**Outputs**
- `sol` - `InverseSolution`: MAP estimate, regularization parameters `λ_g = φ_g τs_g/(c τn)` of the u-step (`c` and `φ_g` are the factors of the likelihood `exp(-c τn ‖r‖²)` and of the prior `exp(-φ_g τs_g ‖L_g u_g‖_q^q)`: `c = 1/2`, `φ_g = 1/q_g` for real data, `c = φ_g = 1` for complex data), hyperparameters `(τn, τs, q)` in `sol.hyper`, convergence information and, if requested, posterior covariance
"""
function solve(prob::InverseProblem, f::Bayesian, ::MAP; R = nothing, covariance::Bool = false,
               tol = 1e-4, maxiter = 100, inner_tol = 1e-8, inner_maxiter = 100,
               init = LCurve())
    c, reg = setup(prob, f.prior, R)
    F = field(eltype(c))
    cn = misfit_factor(F)
    q = Float64[shape(t) for t in reg.terms]
    N = ndata(c)

    U0 = initial_solution(c, reg, init)[1]
    ε = dampings(reg, apply(reg.L, U0))
    st = reweighting(c, reg, Additive(f.prior, init), U0, ε; tol = inner_tol, maxiter = inner_maxiter)
    U, λ0, D = st.U, st.λ, st.D
    τnspec = noise_spec(f.τn, R)
    τn = initial_noise_precision(τnspec, c, reg, λ0, st.W, N)
    τs = initial_signal_precision(f.τs, [λ0 * cn * τn / prior_factor(F, qg) for qg in q], ngroups(reg))
    λ = regularization_parameters(F, τs, τn, q)

    history = Float64[]
    converged = false
    for _ in 1:maxiter
        V = apply(reg.L, U)
        for g in eachindex(reg.rows)
            v = view(V, reg.rows[g], :)
            q[g] = update_shape(f.q, F, v, τs[g], q[g], f.shape_bounds)
            τs[g] = update_signal_precision(f.τs, F, v, q[g], τs[g])
        end
        τn = update_noise_precision(τnspec, F, N, residual_norm2(c, U))
        λ = regularization_parameters(F, τs, τn, q)

        regq = with_terms(reg, Tuple(with_shape(t, qg) for (t, qg) in zip(reg.terms, q)))
        inner = reweighting(c, regq, Additive(f.prior, λ), U, ε;
                            tol = inner_tol, maxiter = inner_maxiter)
        crit = sum(abs, inner.U - U) / max(sum(abs, U), floatmin(Float64))
        push!(history, crit)
        U, D = inner.U, inner.D
        if crit < tol
            converged = true
            break
        end
    end
    Σ = covariance ? scale_covariance(posterior_covariance(c, reg.L, D), τn) : nothing
    hyper = (τn = τn, τs = τs, q = q)
    return InverseSolution(output_shape(prob, U), λ, hyper, length(history), converged, history, Σ)
end

scale_covariance(Σ::AbstractMatrix, τn) = Σ / τn
scale_covariance(Σ::AbstractVector{<:AbstractMatrix}, τn) = [S / τn for S in Σ]

# λ_g = φ_g τs_g / (c τn): u-step of the MAP, ‖r‖² + Σ_g λ_g ‖L_g u_g‖_q^q
regularization_parameters(F::Field, τs, τn, q) =
    [prior_factor(F, qg) * τsg / (misfit_factor(F) * τn) for (τsg, qg) in zip(τs, q)]
