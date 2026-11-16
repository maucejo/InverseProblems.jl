# Posterior uncertainty ---------------------------------------------------------------------

"""
    posterior_std(sol)

Posterior standard deviation of each entry of the solution, from the diagonal of the
posterior covariance.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; covariance = true)`

**Outputs**
- `σ` - `Array{<:Real}`: standard deviations, with the layout of `sol.u`
"""
posterior_std(sol::InverseSolution) = marginal_std(sol.covariance, sol.u)

marginal_std(::Nothing, u) =
    throw(ArgumentError("no posterior covariance: call `solve` with `covariance = true`"))
marginal_std(Σ::AbstractMatrix, u) = reshape(sqrt.(real(diag(Σ))), size(u))
marginal_std(Σ::AbstractVector{<:AbstractMatrix}, u) =
    reduce(hcat, [sqrt.(real(diag(S))) for S in Σ])

"""
    posterior_samples(sol, nsamples = 10_000; rng = Random.default_rng())

Samples of the Gaussian approximation of the posterior `N(μ, Σ)` (circular complex normal
for complex unknowns), with `μ = sol.u`:

    u_sample = μ + S z,   Σ = S Sᴴ (Cholesky),   z ~ randn(eltype(μ))

For block-diagonal operators, each sample (frequency) is drawn from its own `nr × nr`
covariance.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; covariance = true)`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

**Outputs**
- `X` - `Array`: samples, of size `(size(sol.u)..., nsamples)`
"""
function posterior_samples(sol::InverseSolution, nsamples::Integer = 10_000;
                           rng::AbstractRNG = Random.default_rng())
    nsamples >= 1 || throw(ArgumentError("nsamples must be ≥ 1"))
    return draw(rng, sol.covariance, sol.u, nsamples)
end

draw(rng, ::Nothing, μ, n) =
    throw(ArgumentError("no posterior covariance: call `solve` with `covariance = true`"))

# Global covariance of vec(μ)
function draw(rng, Σ::AbstractMatrix, μ, n)
    S = cholesky(Hermitian(Σ)).L
    X = vec(μ) .+ S * randn(rng, eltype(μ), length(μ), n)
    return reshape(X, size(μ)..., n)
end

# One covariance per column of μ (block-diagonal operator)
function draw(rng, Σ::AbstractVector{<:AbstractMatrix}, μ::AbstractMatrix, n)
    X = similar(μ, size(μ)..., n)
    for j in axes(μ, 2)
        S = cholesky(Hermitian(Σ[j])).L
        X[:, j, :] = μ[:, j] .+ S * randn(rng, eltype(μ), size(μ, 1), n)
    end
    return X
end

"""
    credible_interval(sol; level = 0.95, nsamples = 10_000, rng = Random.default_rng())

Bayesian credible interval of each entry of the solution, estimated by sampling the posterior
([`posterior_samples`](@ref)): empirical quantiles `(1 - level)/2` and `(1 + level)/2` of the
samples, entry by entry. For complex unknowns, the quantiles are computed separately on the
real and imaginary parts.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; covariance = true)`
- `level` - `Real`: credibility level, in `]0, 1[`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

**Outputs**
- `lower` - `Array`: lower bounds, with the layout of `sol.u` (complex for complex unknowns:
  `real(lower) ≤ real(u)`, `imag(lower) ≤ imag(u)`)
- `upper` - `Array`: upper bounds, with the layout of `sol.u`
"""
function credible_interval(sol::InverseSolution; level = 0.95, nsamples::Integer = 10_000,
                           rng::AbstractRNG = Random.default_rng())
    0 < level < 1 || throw(ArgumentError("level must be in ]0, 1["))
    X = posterior_samples(sol, nsamples; rng)
    p = ((1 - level) / 2, (1 + level) / 2)
    d = ndims(X)
    bounds = map(eachslice(X; dims = ntuple(identity, d - 1))) do x
        sample_quantiles(x, p)
    end
    return first.(bounds), last.(bounds)
end

sample_quantiles(x::AbstractVector{<:Real}, p) = Tuple(quantile(x, collect(p)))

function sample_quantiles(x::AbstractVector{<:Complex}, p)
    re, im_ = quantile(real(x), collect(p)), quantile(imag(x), collect(p))
    return (complex(re[1], im_[1]), complex(re[2], im_[2]))
end
