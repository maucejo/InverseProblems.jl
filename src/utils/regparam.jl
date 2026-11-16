# Regularization parameter choice rules -------------------------------------------------
#
# All rules act on the spectral description of the weighted problem
#     min ‖y - H u‖² + λ ‖W^{1/2} L u‖²
# given by `spectrum` (SVD or generalized SVD): with the Hansen parameter ℓ (λ = ℓ²), the
# filter factors are f = γ² / (γ² + ℓ²).

"""
    select_parameter(rule, cache, L, W)

Regularization parameter `λ` of `min ‖y - Hu‖² + λ ‖W^{1/2} L u‖²` chosen by `rule`.
"""
function select_parameter(rule::ParameterChoice, c, L::SmoothingOperator, W::AbstractMatrix{<:Real})
    s = spectrum(c, L, W)
    isempty(s.γ) && throw(ArgumentError("no regularized component: the parameter cannot be selected"))
    return select_parameter(rule, s, ndata(c))
end

# A fixed initial value given by the user
select_parameter(λ::Real, c, ::SmoothingOperator, ::AbstractMatrix{<:Real}) = float(λ)

# Search grid on ℓ, between the extreme (generalized) singular values
function lambda_grid(γ; npoints = 200)
    smax = maximum(γ)
    smin = max(minimum(γ), smax * 16eps())
    return exp.(range(log(smax), log(smin); length = npoints))
end

filter_factors(γ, ℓ) = @. γ^2 / (γ^2 + ℓ^2)

"""
    select_parameter(rule, s::NamedTuple, m)

Minimizes the criterion of `rule` w.r.t. `t = log ℓ`: grid search on `lambda_grid`, then
refinement by Brent's method ([`grid_minimize`](@ref)). Returns `λ = ℓ²`.
"""
function select_parameter(rule::ParameterChoice, s::NamedTuple, m)
    J = criterion(rule, s, m)
    t = sort(log.(lambda_grid(s.γ)))
    return exp(2 * grid_minimize(J, t))
end

"""
    criterion(rule, s, m)

Function of `t = log ℓ` minimized by the parameter choice `rule`, for the spectral data `s`
and `m` data.
"""
criterion(::LCurve, s, m) = t -> -lcurve_curvature(s, exp(t))
criterion(::GCV, s, m) = t -> gcv_functional(s, exp(t), m)
criterion(::BayesEstimator, s, m) = t -> bayes_functional(s, exp(t))

"""Curvature of the L-curve in log-log scale at ℓ (Hansen)."""
function lcurve_curvature(s, ℓ)
    (; γ, β2, δ0) = s
    f = filter_factors(γ, ℓ)
    g = @. γ^2 / (γ^2 + ℓ^2)^2                # f²/γ², safe for γ = 0
    η = sum(@. g * β2)                        # ‖W^{1/2} L u‖²
    ρ = sum(@. (1 - f)^2 * β2) + δ0           # ‖y - Hu‖²
    dη = -4 / ℓ * sum(@. (1 - f) * g * β2)    # dη/dℓ
    return 2η * ρ / dη * (ℓ^2 * dη * ρ + 2ℓ * η * ρ + ℓ^4 * η * dη) / (ℓ^2 * η^2 + ρ^2)^(3 / 2)
end

"""Generalized cross validation functional at ℓ."""
function gcv_functional(s, ℓ, m)
    (; γ, β2, δ0, nfree) = s
    f = filter_factors(γ, ℓ)
    return (sum(@. (1 - f)^2 * β2) + δ0) / (m - nfree - sum(f))^2
end

"""Functional `J_MAP` of the Bayesian estimator at ℓ (`reg_param_bayes`)."""
function bayes_functional(s, ℓ)
    (; γ, β2) = s
    M = length(γ)
    d = @. γ^2 + ℓ^2
    return sum(log, d) + (M - 2) * log(sum(β2 ./ d) / M)
end

"""
Heuristic initial parameter (Aucejo & De Smet, MSSP 2017): median of a geometric grid spanning the admissible range of the squared (generalized) singular values, i.e. their geometric mean. When `Hᴴ H + Lᴴ L` restricted to the regularized components is rank deficient, the lower bound is `16 eps` times the upper one.
"""
function select_parameter(::Heuristic, c::Union{DenseCache,BlockCache}, L::SmoothingOperator,
                          W::AbstractMatrix{<:Real})
    s = extreme_eigenvalues(c, L, W)
    smin = max(s.min, s.max * 16eps())
    return sqrt(s.max * smin)
end

function extreme_eigenvalues(c::DenseCache, L, W)
    s = spectrum(c, L, W)
    full = length(s.γ) + s.nfree >= size(c.H, 2)
    return (min = full ? minimum(s.γ)^2 : 0.0, max = maximum(s.γ)^2)
end

function extreme_eigenvalues(c::BlockCache, L, W)
    s = [extreme_eigenvalues(c.blocks[j], L, view(W, :, j:j)) for j in 1:c.ns]
    return (min = minimum(x -> x.min, s), max = maximum(x -> x.max, s))
end
