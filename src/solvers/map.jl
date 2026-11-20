# Optimal Bayesian regularization — MAP solver (Aucejo & De Smet, MSSP 2019) ------------

# Estimated or fixed hyperparameters ---------------------------------------------------------
#
# Each update dispatches on the specification given in `Bayesian`: `Estimated()` → mode of # the full conditional distribution, number → kept at its value.
# Initialization (MAP of the standard formulation with a single parameter λ0, then most probable precisions given λ0 — marginal likelihood, Pereira et al.): τn = N / yᴴ(I + H (Lᴴλ0WL)⁻¹ Hᴴ)⁻¹ y,   τs_g = λ0 c τn / φ_g
# Default noise precision: estimated when R is not provided, 1 (R is the exact noise covariance) when it is.
noise_spec(::Nothing, ::Nothing) = Estimated()
noise_spec(::Nothing, R) = 1.0
noise_spec(τn, R) = τn

initial_noise_precision(::Estimated, sys, λ0, W, N) =
    N / marginal_quadratic(sys.c, weighted_solve!(Matrix{eltype(sys.c)}(undef, sys.c.nu, sys.c.ns), sys, λ0 .* W))
initial_noise_precision(τ::Real, sys, λ0, W, N) = float(τ)

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
    solve_joint(prob, f::Bayesian, MAP(); R = nothing, precision = false, tol = 1e-8, maxiter = 100, inner_tol = 1e-8, inner_maxiter = 100, init = nothing)

Optimal Bayesian regularization: MAP of the hierarchical model, alternating the updates of `q_g`, `τs_g`, `τn` (those marked `Estimated()`) and of `u` (IR algorithm).

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `f` - `Bayesian`: Bayesian formulation
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided
- `precision` - `Bool`: also return the posterior precision `τn (H̃ᴴH̃ + LᴴDL)` of the Gaussian approximation at the MAP (`H̃ = R^{-1/2} H`) and the posterior model (see [`solve`](@ref))
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u` (outer loop)
- `maxiter` - `Integer`: maximum number of outer iterations
- `inner_tol` - `Real`: tolerance of the IR algorithm (u-step)
- `inner_maxiter` - `Integer`: maximum number of IR iterations
- `init` - `Nothing`, `Real` or `ParameterChoice`: rule for the initial single-parameter solution — by default `LCurve()` (heuristic value for a sparse operator)

**Outputs**
- `sol` - `InverseSolution`: MAP estimate, regularization parameters `λ_g = τs_g/τn` of the MAP in `u` — `½‖y - Hu‖²_R + Σ_g (λ_g/q_g)‖L_g u_g‖_q^q` for real data,`‖y - Hu‖²_R + Σ_g λ_g ‖L_g u_g‖_q^q` for complex data —, hyperparameters `(τn, τs, q)` in `sol.hyper`, convergence information and, if requested, posterior precision
"""
function solve_joint(prob::InverseProblem, f::Bayesian, ::MAP; R = nothing, precision = false, tol = 1e-8, maxiter = 100, inner_tol = 1e-8, inner_maxiter = 100, init = nothing)
    c, reg = setup(prob, f.prior, R)
    return map_solve(prob, f, c, reg, noise_spec(f.τn, R), posterior_request(precision), tol, maxiter, inner_tol, inner_maxiter, map_init(init, c, reg))
end

# Initial single-parameter solution: L-curve by default, heuristic value for sparse operators (no SVD available)
map_init(::Nothing, c, reg) = LCurve()
map_init(::Nothing, c::SparseCache, reg) = select_parameter(Heuristic(), c, reg.L, ones(real(eltype(c)), reg.nout, c.ns))
map_init(init, c, reg) = init

# Function barrier (see `ir_solve`)
function map_solve(prob, f, c, reg, τnspec, request, tol, maxiter, inner_tol, inner_maxiter, init)
    F = field(eltype(c))
    cn = misfit_factor(F)
    q = Float64[shape(t) for t in reg.terms]
    N = ndata(c)

    sys = linear_system(c, reg.L)
    U0 = initial_solution(sys, reg, init)[1]
    ε = dampings(reg, apply(reg.L, U0))
    st = reweighting(sys, reg, Additive(f.prior, init), U0, ε; tol = inner_tol, maxiter = inner_maxiter)
    U, λ0, D = st.U, st.λ, st.D
    τn0 = initial_noise_precision(τnspec, sys, λ0, st.W, N) # not reassigned: no boxing in the closure below
    τs = initial_signal_precision(f.τs, [λ0 * cn * τn0 / prior_factor(F, qg) for qg in q], ngroups(reg))
    τn = τn0
    λu = u_step_parameters(F, τs, τn, q)

    Vbuf = smoothing_buffer(reg.L, U)
    history = sizehint!(Float64[], maxiter)
    converged = false
    for _ in 1:maxiter
        V = apply!(Vbuf, reg.L, U)
        for g in eachindex(reg.rows)
            v = view(V, reg.rows[g], :)
            q[g] = update_shape(f.q, F, v, τs[g], q[g], f.shape_bounds)
            τs[g] = update_signal_precision(f.τs, F, v, q[g], τs[g])
        end
        τn = update_noise_precision(τnspec, F, N, residual_norm2!(sys, U))
        λu = u_step_parameters(F, τs, τn, q)

        regq = with_terms(reg, with_shapes(reg.terms, q))
        inner = reweighting(sys, regq, Additive(f.prior, λu), U, ε;
                            tol = inner_tol, maxiter = inner_maxiter)
        crit = relative_change(inner.U, U)
        push!(history, crit)
        U, D = inner.U, inner.D
        if crit < tol
            converged = true
            break
        end
    end
    Q, post = posterior_output(request, c, reg.L, D, τn)
    hyper = (τn = τn, τs = τs, q = q)
    return InverseSolution(output_shape(prob, U), τs ./ τn, hyper, length(history), converged, history, Q, post)
end


# MAP in u, with λ_g = τs_g / τn:
#     real:    ½‖y - Hu‖²_R + Σ_g (λ_g/q_g) ‖L_g u_g‖_q^q
#     complex:  ‖y - Hu‖²_R + Σ_g  λ_g      ‖L_g u_g‖_q^q
# solved by the IR algorithm as ‖y - Hu‖²_R + Σ_g λu_g ‖L_g u_g‖_q^q, with
#     λu_g = φ_g τs_g / (c τn)  (= 2λ_g/q_g real, λ_g complex)
u_step_parameters(F::Field, τs, τn, q) =
    [prior_factor(F, qg) * τsg / (misfit_factor(F) * τn) for (τsg, qg) in zip(τs, q)]
