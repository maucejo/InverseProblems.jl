# Bayesian model (Aucejo & De Smet, MSSP 2019): distributions and full conditionals ------------
#
# Real and complex data lead to different distributions, selected by dispatch on the field
# of the data. With N data and, in group g, M_g components v = L_g u_g:
#
#                     real                                          complex (circular)
#   likelihood        (τn/2π)^{N/2} exp(-τn ‖r‖²/2)                 (τn/π)^N exp(-τn ‖r‖²)
#   prior (per comp.) q^{1-1/q} τs^{1/q} / (2Γ(1/q)) exp(-τs|v|^q/q)  q τs^{2/q} / (2πΓ(2/q)) exp(-τs|v|^q)
#
# i.e. exp(-c τn ‖r‖²) and exp(-φ τs |v|^q) with c = 1/2, φ = 1/q (real) and c = φ = 1
# (complex); for q = 2, both priors are normal laws of precision τs. With λ_g = τs_g / τn,
# the MAP in u is
#     real:    ½‖r‖² + Σ_g (λ_g/q_g) ‖v_g‖_q^q
#     complex:  ‖r‖² + Σ_g  λ_g      ‖v_g‖_q^q
# i.e. ‖r‖² + Σ_g λu_g ‖v_g‖_q^q with λu_g = φ_g τs_g / (c τn) (= 2λ_g/q_g real, λ_g complex).

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

# Real dimension k of one component: prior normalization (φ τs)^{k/q} / Γ(k/q)
prior_dimension(::RealField) = 1
prior_dimension(::ComplexField) = 2

# Factor φ of the norm in the prior exponent: φ τs ‖v‖_q^q
prior_factor(::RealField, q) = 1 / q
prior_factor(::ComplexField, q) = one(q)

# Weakly informative Gamma hyperpriors (a = 1, b → 0)
const HYPER_A = 1.0
const HYPER_B = 1e-18      # precision parameters
const HYPER_BQ = 1e-4      # shape parameters

shape(t::Lq) = t.q
shape(t::NormTerm) =
    throw(ArgumentError("the Bayesian formulation requires Lq terms (generalized Gaussian priors), got $(typeof(t))"))

"""
    shape_objective(F, q, τ, logabs, M)

Negative log full conditional of the shape parameter `q` of a group:

    M log Γ(k/q) - (kM/q) log(φ τ) - (a + M - 1) log q + b q + φ τ Σ|vᵢ|^q

with `k = 1`, `φ = 1/q` (real) or `k = 2`, `φ = 1` (complex).

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `q` - `Real`: shape parameter
- `τ` - `Real`: precision of the group
- `logabs` - `AbstractVector{<:Real}`: `log|vᵢ|` of the non-zero components of the group
- `M` - `Integer`: number of components of the group

**Outputs**
- `J` - `Real`: value of the objective
"""
function shape_objective(F::Field, q, τ, logabs, M)
    k, φ = prior_dimension(F), prior_factor(F, q)
    s = sum(x -> exp(q * x), logabs; init = 0.0)
    return M * loggamma(k / q) - (k * M / q) * log(φ * τ) - (HYPER_A + M - 1) * log(q) +
           HYPER_BQ * q + φ * τ * s
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
    logabs = [log(abs(x)) for x in v if !iszero(x)]
    M = length(v)
    return grid_minimize(q -> shape_objective(F, q, τ, logabs, M), range(lb, ub; length = npoints))
end

"""
    signal_precision(F, v, q)

Mode of the full conditional of the precision of a group.

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `v` - `AbstractArray`: components `L_g u_g` of the group
- `q` - `Real`: shape parameter of the group

**Outputs**
- `τs` - `Real`: signal precision `(k M / q + a - 1) / (φ Σ|vᵢ|^q + b)` (`k = 1, φ = 1/q`
  real, `k = 2, φ = 1` complex)
"""
signal_precision(F::Field, v, q) =
    (prior_dimension(F) * length(v) / q + HYPER_A - 1) /
    (prior_factor(F, q) * sum(x -> abs(x)^q, v) + HYPER_B)

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

"""
    elementwise_precision(F, V)

Mode of the full conditional of the precision of each component of `V`, for a Gaussian prior
(`q = 2`) with one precision per component (element-wise Bayesian model).

**Inputs**
- `F` - `Field`: `RealField()` or `ComplexField()`
- `V` - `AbstractArray`: components `L u`

**Outputs**
- `τ` - `Array{<:Real}`: precisions `(k/2 + a - 1) / (φ |vᵢ|² + b)`, with the layout of `V`
  (`1/(|vᵢ|² + 2b)` real, `a/(|vᵢ|² + b)` complex)
"""
function elementwise_precision(F::Field, V)
    return elementwise_precision!(similar(V, real(eltype(V))), F, V)
end

# In-place version: τ .= (k/2 + a - 1) / (φ |v|² + b)
function elementwise_precision!(τ, F::Field, V)
    k, φ = prior_dimension(F), prior_factor(F, 2.0)
    @. τ = (k / 2 + HYPER_A - 1) / (φ * abs2(V) + HYPER_B)
    return τ
end
