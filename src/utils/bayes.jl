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
    shape_objective(q, τ, logabs, M, k)

Negative log full conditional of the shape parameter `q` of a group:

    M log Γ(k/q) - (kM/q) log τ - (a + M - 1) log q + b q + τ Σ|vᵢ|^q

**Inputs**
- `q` - `Real`: shape parameter
- `τ` - `Real`: precision of the group
- `logabs` - `AbstractVector{<:Real}`: `log|vᵢ|` of the non-zero components of the group
- `M` - `Integer`: number of components of the group
- `k` - `Integer`: real dimension of a component (1 real, 2 complex)

**Outputs**
- `J` - `Real`: value of the objective
"""
function shape_objective(q, τ, logabs, M, k)
    s = sum(x -> exp(q * x), logabs; init = 0.0)
    return M * loggamma(k / q) - (k * M / q) * log(τ) - (HYPER_A + M - 1) * log(q) +
           HYPER_BQ * q + τ * s
end

"""
    optimal_shape(F, v, τ, (lb, ub); npoints = 200)

Optimal shape parameter of a group (grid search refined by Brent's method).

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `v` - `AbstractArray`: components `L_g u_g` of the group
- `τ` - `Real`: precision of the group
- `(lb, ub)` - `Tuple{Real,Real}`: search interval
- `npoints` - `Integer`: number of points of the search grid

**Outputs**
- `q` - `Real`: mode of the full conditional of the shape parameter
"""
function optimal_shape(F::Field, v, τ, (lb, ub); npoints = 200)
    k = prior_dimension(F)
    logabs = [log(abs(x)) for x in v if !iszero(x)]
    M = length(v)
    return grid_minimize(q -> shape_objective(q, τ, logabs, M, k), range(lb, ub; length = npoints))
end

"""
    signal_precision(F, v, q)

Mode of the full conditional of the precision of a group.

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `v` - `AbstractArray`: components `L_g u_g` of the group
- `q` - `Real`: shape parameter of the group

**Outputs**
- `τs` - `Real`: signal precision `(k M / q + a - 1) / (Σ|vᵢ|^q + b)`
"""
signal_precision(F::Field, v, q) =
    (prior_dimension(F) * length(v) / q + HYPER_A - 1) / (sum(x -> abs(x)^q, v) + HYPER_B)

"""
    noise_precision(F, N, r2)

Mode of the full conditional of the noise precision.

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `N` - `Integer`: number of data
- `r2` - `Real`: squared (whitened) data misfit `‖y - Hu‖²`

**Outputs**
- `τn` - `Real`: noise precision `(ν + a - 1) / (c r2 + b)` (`ν = N/2, c = 1/2` real,
  `ν = N, c = 1` complex)
"""
noise_precision(F::Field, N, r2) =
    (noise_exponent(F, N) + HYPER_A - 1) / (misfit_factor(F) * r2 + HYPER_B)
