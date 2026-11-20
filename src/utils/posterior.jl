# Posterior Gaussian approximation --------------------------------------------------------------
#
# At the fixed point of the algorithms, the posterior is approximated by a Gaussian distribution
# of mean u and precision (inverse covariance)
#
#     Q = τn A,    A = Hᴴ H + Lᴴ D L        (whitened H, D = λ W: weights of the last solve)
#
# (exact for q = 2; τn = 1 for the additive and multiplicative formulations, R being then taken as the exact noise covariance). The precision is returned as is — sparse for an assembled space–frequency problem, one block per sample for independent / block-diagonal samples, dense otherwise —, and the covariance Q⁻¹ is never formed by default:
# * samples, by dispatch on the factorization of A used by the solver:
#   - Cholesky factor available (normal equations, sparse Cholesky): standard sampling
#     δ = R⁻¹ z, A = Rᴴ R (one triangular solve per sample; P'·(…) for CHOLMOD);
#   - Woodbury identity (under-determined problems, no Cholesky factor of A): perturbation–
#     optimization, δ = A⁻¹ (Hᴴ ε₁ + Lᴴ D^{1/2} ε₂) with ε₁ ~ N(0, I) (data size) and
#     ε₂ ~ N(0, I) (size of L u), whose covariance is A⁻¹ (Hᴴ H + Lᴴ D L) A⁻¹ = A⁻¹;
#   z, ε₁, ε₂ are circular complex normal for complex data, and u + δ / √τn is a sample of the
#   posterior;
# * marginal standard deviations and covariances: diagonal blocks of A⁻¹ (one per sample), each computed by applying A⁻¹ to the corresponding columns of the identity (block solves).
"""
Factorization of `A = Hᴴ H + D` (`D` diagonal) by the Woodbury identity, on the standard form
`H̃ = H L⁻¹` for a square invertible smoothing operator:

    A⁻¹ B = D⁻¹ B − D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ H D⁻¹ B

**Fields**
- `KD` - `Matrix`: `H D⁻¹` (`H̃ D⁻¹` for an invertible smoothing operator)
- `d` - `Vector{<:Real}`: diagonal of `D`
- `F` - `Cholesky`: factorization of `I + H D⁻¹ Hᴴ` (size `nm`)
- `S` - `SmoothingOperator`: smoothing operator (`IdentitySmoothing` or `InvertibleSmoothing`)
- `nu` - `Int`: number of rows of the unknown
"""
struct WoodburyFactor{T, R<:Real, F, S<:SmoothingOperator}
    KD::Matrix{T}
    d::Vector{R}
    F::F
    S::S
    nu::Int
end

function WoodburyFactor(K::AbstractMatrix, d::AbstractVector{<:Real}, S::SmoothingOperator, nu)
    KD = K ./ transpose(d)
    return WoodburyFactor(KD, collect(d), cholesky!(Hermitian(KD * K' + I)), S, nu)
end

"""
    posterior_factor(cache, L, D)

Factorization of `A = Hᴴ H + Lᴴ D L` used to apply `A⁻¹` (dispatch on the regime of the problem
and on the smoothing operator, as for [`weighted_solve!`](@ref)).

**Inputs**
- `cache` - `DenseCache`, `SparseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: smoothing operator
- `D` - `AbstractMatrix{<:Real}`: weights of the last linear solve (size of `L U`)

**Outputs**
- `F` - dense `Cholesky`, `WoodburyFactor`, sparse Cholesky factor (CHOLMOD) or vector of
  factorizations (one per sample of a block-diagonal operator)
"""
posterior_factor(c::Overdetermined, S::SmoothingOperator, D) =
    cholesky(Hermitian(c.HtH + weighted_gram(c, S, vec(D))))
posterior_factor(c::Underdetermined, ::IdentitySmoothing, D) =
    WoodburyFactor(c.H, vec(D), IdentitySmoothing(), c.nu)
posterior_factor(c::Underdetermined, S::InvertibleSmoothing, D) =
    WoodburyFactor(standard_form(c, S), vec(D), S, c.nu)
posterior_factor(c::Underdetermined, S::SingularSmoothing, D) =
    cholesky(Hermitian(gram(c) + weighted_gram(c, S, vec(D))))
posterior_factor(c::SparseCache, S::SmoothingOperator, D) =
    cholesky(Hermitian(sparse_system_matrix(c, S, vec(D))))
posterior_factor(c::BlockCache, S::SmoothingOperator, D) =
    [posterior_factor(c.blocks[j], S, view(D, :, j:j)) for j in 1:c.ns]

"""
    inverse_apply(F, B)

`A⁻¹ B` for a matrix of right-hand sides `B` (`n × k`), from a factorization of `A`.

**Inputs**
- `F` - factorization given by [`posterior_factor`](@ref)
- `B` - `AbstractMatrix`: right-hand sides

**Outputs**
- `X` - `Matrix`: `A⁻¹ B`
"""
inverse_apply(F::Factorization, B) = F \ B
function inverse_apply(W::WoodburyFactor, B)
    Bs = to_standard_form(W.S, B, W.nu)                       # L⁻ᴴ B
    X = Bs ./ W.d .- W.KD' * (W.F \ (W.KD * Bs))              # (H̃ᴴH̃ + D)⁻¹ L⁻ᴴ B
    return from_standard_form(W.S, X, W.nu)                   # L⁻¹ (…)
end

to_standard_form(::IdentitySmoothing, B, nu) = B
to_standard_form(S::InvertibleSmoothing, B, nu) = smoothing_ldiv_adjoint(S, Matrix(B), nu)
from_standard_form(::IdentitySmoothing, X, nu) = X
from_standard_form(S::InvertibleSmoothing, X, nu) = smoothing_ldiv(S, X, nu)

"""
    precision_matrix(cache, L, D, τn)

Posterior precision `Q = τn (Hᴴ H + Lᴴ D L)`.

**Inputs**
- `cache` - `DenseCache`, `SparseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: smoothing operator
- `D` - `AbstractMatrix{<:Real}`: weights of the last linear solve
- `τn` - `Real`: noise precision (`1` for the additive and multiplicative formulations)

**Outputs**
- `Q` - `Hermitian` (dense operator), `SparseMatrixCSC` (sparse operator) or vector of
  `Hermitian` (one `nu × nu` block per sample of a block-diagonal operator)
"""
precision_matrix(c::DenseCache, S::SmoothingOperator, D, τn) =
    Hermitian(τn .* (gram(c) + weighted_gram(c, S, vec(D))))
precision_matrix(c::SparseCache, S::SmoothingOperator, D, τn) = τn .* sparse_system_matrix(c, S, vec(D))
precision_matrix(c::BlockCache, S::SmoothingOperator, D, τn) =
    [precision_matrix(c.blocks[j], S, view(D, :, j:j), τn) for j in 1:c.ns]

"""
Posterior Gaussian approximation of a problem solved as a whole.

**Fields**
- `c` - `DenseCache`, `SparseCache` or `BlockCache`: linear cache of the whitened problem
- `S` - `SmoothingOperator`: smoothing operator
- `D` - `Matrix{<:Real}`: weights of the last linear solve (size of `L U`)
- `F` - factorization of `A = Hᴴ H + Lᴴ D L` ([`posterior_factor`](@ref))
- `τn` - `Float64`: noise precision (`Q = τn A`)
- `nb` - `Int`: size of the diagonal blocks of the marginal covariances (number of
  reconstruction points of a sample)
"""
struct Posterior{C, S<:SmoothingOperator, R<:Real, F}
    c::C
    S::S
    D::Matrix{R}
    F::F
    τn::Float64
    nb::Int
end

"""
Posterior of a problem solved sample by sample (independent frequencies): one
[`Posterior`](@ref) per sample.
"""
struct SamplewisePosterior{P<:Posterior}
    posteriors::Vector{P}
end

"""
Request of the posterior in `solve` (`precision = true`): `nb` is the size of the diagonal
blocks (number of reconstruction points of a sample), `nothing` for the rows of the unknown.
"""
struct PosteriorRequest{B<:Union{Nothing,Int}}
    nb::B
end

posterior_request(x::Bool, nb = nothing) = x ? PosteriorRequest(nb) : nothing
posterior_request(r::PosteriorRequest, nb = nothing) = r
posterior_request(::Nothing, nb = nothing) = nothing

block_size(r::PosteriorRequest{Nothing}, c) = c.nu
block_size(r::PosteriorRequest{Int}, c) = r.nb

"""
    posterior_output(request, cache, L, D, τn)

Posterior precision and posterior model returned by the solvers.

**Inputs**
- `request` - `Nothing` or `PosteriorRequest`: `precision` keyword of [`solve`](@ref)
- `cache` - linear cache of the whitened problem
- `L` - `SmoothingOperator`: smoothing operator
- `D` - `AbstractMatrix{<:Real}`: weights of the last linear solve
- `τn` - `Real`: noise precision

**Outputs**
- `(Q, post)` - precision ([`precision_matrix`](@ref)) and [`Posterior`](@ref), or
  `(nothing, nothing)`
"""
posterior_output(::Nothing, c, S, D, τn) = (nothing, nothing)
function posterior_output(r::PosteriorRequest, c, S, D, τn)
    Dm = Matrix(D)
    post = Posterior(c, S, Dm, posterior_factor(c, S, Dm), Float64(τn), block_size(r, c))
    return precision_matrix(c, S, Dm, τn), post
end

# Samples -----------------------------------------------------------------------------------
#
# k samples δ / √τn of the centered posterior N(0, (τn A)⁻¹), as columns, by dispatch on the factorization of A.

# Cholesky factor of A available: standard sampling δ = R⁻¹ z, A = Rᴴ R (dense Cholesky)
centered_samples(rng, c, S, D, F::Cholesky, τn, k) = (F.U \ randn(rng, eltype(c), size(F, 1), k)) ./ sqrt(τn)

# Sparse Cholesky A = (P'L)(L'P) (fill-reducing permutation P): δ = (L'P)⁻¹ z
centered_samples(rng, c, S, D, F::SparseArrays.CHOLMOD.Factor, τn, k) = (F.UP \ randn(rng, eltype(c), size(F, 1), k)) ./ sqrt(τn)

# Woodbury identity (no Cholesky factor of A): perturbation–optimization,
# δ = A⁻¹ (Hᴴ ε₁ + Lᴴ D^{1/2} ε₂), Cov(δ) = A⁻¹ (Hᴴ H + Lᴴ D L) A⁻¹ = A⁻¹
function centered_samples(rng, c, S, D, F::WoodburyFactor, τn, k)
    T = eltype(c)
    E1 = randn(rng, T, size(c.H, 1), k)
    E2 = randn(rng, T, length(D), k)
    B = c.H' * E1 + adjoint_apply(c, S, sqrt.(vec(D)) .* E2)
    return inverse_apply(F, B) ./ sqrt(τn)
end

# Lᴴ X for X with the layout of vec(L U)
adjoint_apply(c, ::IdentitySmoothing, X) = X
adjoint_apply(c, S::MatrixSmoothing, X) = global_operator(c, S)' * X

function draw_posterior(rng, P::Posterior{<:Union{DenseCache,SparseCache}}, μ, k)
    X = centered_samples(rng, P.c, P.S, P.D, P.F, P.τn, k)
    return reshape(vec(μ) .+ X, size(μ)..., k)
end

function draw_posterior(rng, P::Posterior{<:BlockCache}, μ::AbstractMatrix, k)
    X = similar(μ, size(μ)..., k)
    for j in axes(μ, 2)
        X[:, j, :] = μ[:, j] .+ centered_samples(rng, P.c.blocks[j], P.S, view(P.D, :, j:j), P.F[j], P.τn, k)
    end
    return X
end

function draw_posterior(rng, P::SamplewisePosterior, μ::AbstractMatrix, k)
    X = similar(μ, size(μ)..., k)
    for j in axes(μ, 2)
        X[:, j, :] = draw_posterior(rng, P.posteriors[j], μ[:, j], k)
    end
    return X
end

# Marginal covariances (block solves) -----------------------------------------------------------

# Diagonal blocks of A⁻¹ / τn, each obtained by applying A⁻¹ to nb columns of the identity
function marginal_blocks(F, n, nb, τn, T)
    n % nb == 0 || throw(DimensionMismatch("$n unknowns cannot be split into blocks of size $nb"))
    E = zeros(T, n, nb)
    return map(1:(n ÷ nb)) do j
        rows = ((j - 1) * nb + 1):(j * nb)
        fill!(E, zero(T))
        for (k, i) in enumerate(rows)
            E[i, k] = one(T)
        end
        Hermitian(Matrix(inverse_apply(F, E)[rows, :]) ./ τn)
    end
end

marginal_covariances(P::Posterior{<:Union{DenseCache,SparseCache}}) =
    marginal_blocks(P.F, size(P.c.H, 2), P.nb, P.τn, eltype(P.c))
marginal_covariances(P::Posterior{<:BlockCache}) =
    [only(marginal_blocks(P.F[j], P.c.nu, P.c.nu, P.τn, eltype(P.c))) for j in 1:P.c.ns]
marginal_covariances(P::SamplewisePosterior) = reduce(vcat, map(marginal_covariances, P.posteriors))

# Full covariance of vec(u)
function full_covariance(P::Posterior{<:Union{DenseCache,SparseCache}})
    n, T = size(P.c.H, 2), eltype(P.c)
    return Hermitian(Matrix(inverse_apply(P.F, Matrix{T}(I, n, n))) ./ P.τn)
end
full_covariance(P::Union{Posterior{<:BlockCache},SamplewisePosterior}) =      # block-diagonal
    blockdiag(sparse.(parent.(marginal_covariances(P)))...)
