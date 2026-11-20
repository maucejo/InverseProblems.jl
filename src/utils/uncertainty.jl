# Posterior uncertainty ---------------------------------------------------------------------

"""
    posterior_std(sol)

Posterior standard deviation of each entry of the solution: square root of the diagonal of the posterior covariance `Q⁻¹`, obtained from its diagonal blocks (block solves, see [`posterior_covariance`](@ref)) — `Q⁻¹` is never formed.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`

**Outputs**
- `σ` - `Array{<:Real}`: standard deviations, with the layout of `sol.u`
"""
posterior_std(sol::InverseSolution) =
    reshape(reduce(vcat, [sqrt.(real(diag(B))) for B in marginal_covariances(posterior_model(sol))]), size(sol.u))

"""
    posterior_covariance(sol; full = false)

Posterior covariance `Q⁻¹` of the solution, computed from the factorization of the precision.

* `full = false`: marginal covariances of the samples (frequencies or time steps), i.e. the `nu × nu` diagonal blocks `Σ_jj` of the covariance of `vec(u)`. Each block is obtained by applying `A⁻¹` to the `nu` corresponding columns of the identity (block solves with the factorization of the solver: normal equations, Woodbury identity, sparse Cholesky), so that the `n × n` covariance is never formed;
* `full = true`: covariance of `vec(u)`, `(nu·ns) × (nu·ns)`, with the correlations between
  samples (sparse block-diagonal matrix when the samples are independent or block-diagonal).

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `full` - `Bool`: full covariance instead of the marginal blocks

**Outputs**
- `Σ` - `Hermitian` (single sample, or `full = true`), `Vector{<:Hermitian}` (one block per sample) or `SparseMatrixCSC` (`full = true`, independent samples)
"""
function posterior_covariance(sol::InverseSolution; full::Bool = false)
    P = posterior_model(sol)
    full && return full_covariance(P)
    blocks = marginal_covariances(P)
    return length(blocks) == 1 ? only(blocks) : blocks
end

"""
    posterior_samples(sol, nsamples = 10_000; rng = Random.default_rng())

Samples of the Gaussian approximation of the posterior `N(u, Q⁻¹)` (circular complex normal for complex unknowns), with `Q = τn A`, `A = Hᴴ H + Lᴴ D L`. The sampler is chosen by dispatch on the factorization of `A` used by the solver:

* Cholesky factor `A = Rᴴ R` available (normal equations, sparse Cholesky): **standard sampling** `u_sample = u + R⁻¹ z / √τn`, `z ~ N(0, I)` — one triangular solve per sample (`R⁻¹ R⁻ᴴ = A⁻¹`; with a lower factor `A = L Lᴴ`, this reads `L⁻ᴴ z`);
* Woodbury identity (under-determined problems, no Cholesky factor of `A`): **perturbation–optimization**, with `ε₁ ~ N(0, I)`, `ε₂ ~ N(0, I)` of the sizes of the data and of `L u`,

      u_sample = u + A⁻¹ (Hᴴ ε₁ + Lᴴ D^{1/2} ε₂) / √τn,

  whose covariance is `A⁻¹ (Hᴴ H + Lᴴ D L) A⁻¹ / τn = Q⁻¹`, at the cost of one application of `A⁻¹` (Woodbury) per sample, without forming nor factorizing the `n × n` matrix `A`.

Both give exact samples of `N(u, Q⁻¹)`; the joint samples keep the correlations between frequencies or time steps (block-diagonal problems are sampled block by block).

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

**Outputs**
- `X` - `Array`: samples, of size `(size(sol.u)..., nsamples)`
"""
function posterior_samples(sol::InverseSolution, nsamples::Integer = 10_000;
                           rng::AbstractRNG = Random.default_rng())
    nsamples >= 1 || throw(ArgumentError("nsamples must be ≥ 1"))
    return draw_posterior(rng, posterior_model(sol), sol.u, nsamples)
end

posterior_model(sol::InverseSolution) = posterior_model(sol.posterior)
posterior_model(::Nothing) =
    throw(ArgumentError("no posterior: call `solve` with `precision = true`"))
posterior_model(P) = P

"""
    credible_interval(sol; level = 0.95, nsamples = 10_000, rng = Random.default_rng())

Bayesian credible interval of each entry of the solution, estimated by sampling the posterior ([`posterior_samples`](@ref)): empirical quantiles `(1 - level)/2` and `(1 + level)/2` of the samples, entry by entry. For complex unknowns, the quantiles are computed separately on the real and imaginary parts.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `level` - `Real`: credibility level, in `]0, 1[`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

**Outputs**
- `lower` - `Array`: lower bounds, with the layout of `sol.u` (complex for complex unknowns: `real(lower) ≤ real(u)`, `imag(lower) ≤ imag(u)`)
- `upper` - `Array`: upper bounds, with the layout of `sol.u`
"""
function credible_interval(sol::InverseSolution; level = 0.95, nsamples::Integer = 10_000,
                           rng::AbstractRNG = Random.default_rng())
    0 < level < 1 || throw(ArgumentError("level must be in ]0, 1["))
    X = posterior_samples(sol, nsamples; rng)
    p = ((1 - level) / 2, (1 + level) / 2)
    # last dimension of X = samples; the other ones have the layout of sol.u
    lower = similar(X, axes(X)[1:end-1])
    upper = similar(lower)
    for I in CartesianIndices(lower)
        lower[I], upper[I] = sample_quantiles(view(X, I, :), p)
    end
    return output_layout(sol.u, lower), output_layout(sol.u, upper)
end

output_layout(::AbstractVector, B) = vec(B)
output_layout(::AbstractMatrix, B) = B

function sample_quantiles(x::AbstractVector{<:Real}, p)
    q = quantile(x, collect(p))
    return (q[1], q[2])   # explicit 2-tuple: Tuple(q) has a length unknown to inference
end

function sample_quantiles(x::AbstractVector{<:Complex}, p)
    re, im_ = quantile(real(x), collect(p)), quantile(imag(x), collect(p))
    return (complex(re[1], im_[1]), complex(re[2], im_[2]))
end
