# Bayesian model (Aucejo & De Smet, MSSP 2019): distributions and full conditionals ------------
#
# Real and complex data lead to different distributions, selected by dispatch on the field
# of the data. With N data and, in group g, M_g components v = L_g u_g:
#
#                     real                                   complex (circular)
#   likelihood        (τn/2π)^{N/2} exp(-τn ‖r‖²/2)          (τn/π)^N exp(-τn ‖r‖²)
#   prior (per comp.) q τs^{1/q} / (2 Γ(1/q)) exp(-τs|v|^q)   q τs^{2/q} / (2π Γ(2/q)) exp(-τs|v|^q)
#
# so that the MAP u-step is ‖r‖² + Σ_g λ_g ‖v_g‖_q^q with λ_g = τs_g / (c τn), c = 1/2 or 1.

"""Field of the data, selecting the real or complex Bayesian model."""
abstract type Field end
struct RealField <: Field end
struct ComplexField <: Field end

field(::Type{<:Real}) = RealField()
field(::Type{<:Complex}) = ComplexField()

# Exponent of τn in the likelihood: N/2 (real), N (complex)
noise_exponent(::RealField, N) = N / 2
noise_exponent(::ComplexField, N) = N

# Factor c of the data misfit in the likelihood exponent: τn c ‖r‖²
misfit_factor(::RealField) = 1 / 2
misfit_factor(::ComplexField) = 1

# Real dimension k of one component: prior normalization τs^{k/q} / Γ(k/q)
prior_dimension(::RealField) = 1
prior_dimension(::ComplexField) = 2

# Weakly informative Gamma hyperpriors (a = 1, b → 0)
const HYPER_A = 1.0
const HYPER_B = 1e-18      # precision parameters
const HYPER_BQ = 1e-4      # shape parameters

shape(t::Lq) = t.q
shape(t::NormTerm) =
    throw(ArgumentError("the Bayesian formulation requires Lq terms (generalized Gaussian priors), got $(typeof(t))"))

"""
Negative log full conditional of the shape parameter `q` of a group `v` (`M` components,
real dimension `k` per component):

    M log Γ(k/q) - (kM/q) log τ - (a + M - 1) log q + b q + τ Σ|vᵢ|^q
"""
function shape_objective(q, τ, logabs, M, k)
    s = sum(x -> exp(q * x), logabs; init = 0.0)
    return M * loggamma(k / q) - (k * M / q) * log(τ) - (HYPER_A + M - 1) * log(q) +
           HYPER_BQ * q + τ * s
end

"""
Optimal shape parameter in `[lb, ub]` (brute-force search followed by a local refinement).
"""
function optimal_shape(F::Field, v, τ, (lb, ub); npoints = 200)
    k = prior_dimension(F)
    logabs = [log(abs(x)) for x in v if !iszero(x)]
    M = length(v)
    obj(q) = shape_objective(q, τ, logabs, M, k)
    grid = range(lb, ub; length = npoints)
    i = argmin(obj.(grid))
    # golden-section refinement on the neighbouring cells
    a, b = grid[max(i - 1, 1)], grid[min(i + 1, npoints)]
    φ = (sqrt(5) - 1) / 2
    for _ in 1:40
        x1, x2 = b - φ * (b - a), a + φ * (b - a)
        obj(x1) < obj(x2) ? (b = x2) : (a = x1)
    end
    return (a + b) / 2
end

"""Mode of the full conditional of the precision of a group (`M` components)."""
signal_precision(F::Field, v, q) =
    (prior_dimension(F) * length(v) / q + HYPER_A - 1) / (sum(x -> abs(x)^q, v) + HYPER_B)

"""Mode of the full conditional of the noise precision (`N` data, squared misfit `r2`)."""
noise_precision(F::Field, N, r2) =
    (noise_exponent(F, N) + HYPER_A - 1) / (misfit_factor(F) * r2 + HYPER_B)
