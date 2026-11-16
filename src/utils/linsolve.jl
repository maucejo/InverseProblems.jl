# Weighted Tikhonov solves -----------------------------------------------------------
#
# Every algorithm of the package boils down to repeated solves of
#
#     (Hᴴ H + Lᴴ D L) u = Hᴴ y        (D diagonal, positive)
#
# and, for the posterior covariance, to the inverse of the same matrix. The structure of H (dense / block-diagonal), the regime of the problem and the smoothing operator L are exploited through dispatch:
#
#   regime                    L                   system factorized
#   nm ≥ nr (over/square)     any                 normal equations           nr × nr
#   nm < nr (under)           I                   Woodbury identity          nm × nm
#   nm < nr (under)           square invertible   Woodbury on H̃ = H L⁻¹     nm × nm
#   nm < nr (under)           rectangular         normal equations           nr × nr
#                                                 (Lᴴ D L singular: no Woodbury)
#
# Woodbury identity (D invertible):
#     (Hᴴ H + D)⁻¹ = D⁻¹ - D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ H D⁻¹
#     (Hᴴ H + D)⁻¹ Hᴴ y = D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ y
# and, for L square invertible, (Hᴴ H + Lᴴ D L)⁻¹ = L⁻¹ (H̃ᴴ H̃ + D)⁻¹ L⁻ᴴ.
"""
Cache of the quantities reused across iterations for a dense operator. The type of `HtH` (`Hᴴ H`)
encodes the regime: a matrix when `nm ≥ nr` (over-determined or square problem), `nothing`
when `nm < nr` (under-determined problem, `Hᴴ H` is not needed).
"""
struct DenseCache{T, M<:AbstractMatrix{T}, G<:Union{Matrix{T},Nothing}}
    H::M
    y::Vector{T}
    HtH::G
    Hty::Vector{T}
    nr::Int
    ns::Int
end

"""Over-determined or square problem (`nm ≥ nr`)."""
const Overdetermined = DenseCache{<:Any, <:Any, <:Matrix}

"""Under-determined problem (`nm < nr`)."""
const Underdetermined = DenseCache{<:Any, <:Any, Nothing}

"""
Cache for a block-diagonal operator (one dense cache per sample).
"""
struct BlockCache{C<:DenseCache}
    blocks::Vector{C}
    nr::Int
    ns::Int
end

function dense_cache(H::AbstractMatrix, y::AbstractVecOrMat, nr, ns, σ2 = nothing)
    T = promote_type(eltype(H), eltype(y), Float64)
    Hc = Matrix{T}(H)
    yc = Vector{T}(vec(y))
    whiten!(Hc, yc, σ2)
    HtH = size(Hc, 1) >= size(Hc, 2) ? Hc' * Hc : nothing
    return DenseCache(Hc, yc, HtH, Hc' * yc, nr, ns)
end

"""
    linear_cache(prob, R = nothing)

Cache of the quantities reused across iterations (`H`, `y`, `Hᴴ H`, `Hᴴ y`), for the
problem whitened by the noise covariance `R`.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided

**Outputs**
- `c` - `DenseCache` (dense operator) or `BlockCache` (block-diagonal operator)
"""
linear_cache(prob::InverseProblem{<:DenseOperator}, R = nothing) =
    dense_cache(prob.H.H, prob.y, prob.H.nr, prob.H.ns, noise_variances(R, prob.y))

function linear_cache(prob::InverseProblem{<:BlockDiagonalOperator}, R = nothing)
    nr, ns = unknown_size(prob)
    σ2 = noise_variances(R, prob.y)
    column(j) = σ2 === nothing ? nothing : view(σ2, :, j)
    blocks = [dense_cache(prob.H.blocks[j], view(prob.y, :, j), nr, 1, column(j)) for j in 1:ns]
    return BlockCache(blocks, nr, ns)
end

Base.eltype(::DenseCache{T}) where {T} = T
Base.eltype(c::BlockCache) = eltype(first(c.blocks))
ndata(c::DenseCache) = length(c.y)
ndata(c::BlockCache) = sum(ndata, c.blocks)

# Smoothing operator acting on vec(U) for a dense cache
global_operator(c::DenseCache, S::MatrixSmoothing) =
    c.ns == 1 ? S.L : kron(sparse(I, c.ns, c.ns), S.L)

# Gram matrix Hᴴ H (computed on demand for an under-determined problem)
gram(c::Overdetermined) = c.HtH
gram(c::Underdetermined) = c.H' * c.H

# Lᴴ diag(d) L
weighted_gram(c::DenseCache, ::IdentitySmoothing, d) = Diagonal(d)
function weighted_gram(c::DenseCache, S::MatrixSmoothing, d)
    Lg = global_operator(c, S)
    return Matrix(Lg' * Diagonal(d) * Lg)
end

# Kernels ----------------------------------------------------------------------------------

# (G + M)⁻¹ b, normal equations
normal_solve(G, M, b) = cholesky!(Hermitian(G + M)) \ b

# (Hᴴ H + D)⁻¹ Hᴴ y = D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ y, Woodbury identity
function woodbury_solve(H, d, y)
    HD = H ./ transpose(d)                                     # H D⁻¹
    return HD' * (cholesky!(Hermitian(HD * H' + I)) \ y)
end

# (Hᴴ H + D)⁻¹ = D⁻¹ - D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ H D⁻¹, Woodbury identity
function woodbury_inverse(H, d)
    HD = H ./ transpose(d)
    return Diagonal(1 ./ d) - HD' * (cholesky!(Hermitian(HD * H' + I)) \ HD)
end

# H̃ = H L⁻¹ (standard form for a square invertible L)
standard_form(c::DenseCache, S::InvertibleSmoothing) = Matrix(smoothing_ldiv_adjoint(S, Matrix(c.H'), c.nr)')

"""
    weighted_solve(cache, L, D)

Solution of the weighted Tikhonov system `(Hᴴ H + Lᴴ diag(vec(D)) L) vec(U) = Hᴴ y`. The
normal equations or the Woodbury identity are used depending on the regime of the problem
and on the smoothing operator (see the table at the top of `linsolve.jl`).

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: `IdentitySmoothing()` or `MatrixSmoothing`
- `D` - `AbstractMatrix{<:Real}`: positive diagonal weights, with the size of `L U`

**Outputs**
- `U` - `Matrix`: solution (`nr × ns`)
"""
weighted_solve(c::Overdetermined, S::SmoothingOperator, D::AbstractMatrix{<:Real}) =
    reshape(normal_solve(c.HtH, weighted_gram(c, S, vec(D)), c.Hty), c.nr, c.ns)

weighted_solve(c::Underdetermined, ::IdentitySmoothing, D::AbstractMatrix{<:Real}) =
    reshape(woodbury_solve(c.H, vec(D), c.y), c.nr, c.ns)

function weighted_solve(c::Underdetermined, S::InvertibleSmoothing, D::AbstractMatrix{<:Real})
    v = woodbury_solve(standard_form(c, S), vec(D), c.y)      # v = L u
    return reshape(smoothing_ldiv(S, v, c.nr), c.nr, c.ns)
end

weighted_solve(c::Underdetermined, S::SingularSmoothing, D::AbstractMatrix{<:Real}) =
    reshape(normal_solve(gram(c), weighted_gram(c, S, vec(D)), c.Hty), c.nr, c.ns)

function weighted_solve(c::BlockCache, S::SmoothingOperator, D::AbstractMatrix{<:Real})
    U = Matrix{eltype(c)}(undef, c.nr, c.ns)
    Threads.@threads for j in 1:c.ns
        U[:, j] = weighted_solve(c.blocks[j], S, view(D, :, j:j))
    end
    return U
end

"""
    posterior_covariance(cache, L, D)

Covariance `(Hᴴ H + Lᴴ diag(vec(D)) L)⁻¹` (whitened `H`) of the Gaussian approximation of
the posterior at a fixed point of the IR algorithm (exact for `q = 2`). The inverse is
computed from the normal equations or the Woodbury identity, as in [`weighted_solve`](@ref).

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: `IdentitySmoothing()` or `MatrixSmoothing`
- `D` - `AbstractMatrix{<:Real}`: scaled weights `λ W` of the last linear solve

**Outputs**
- `Σ` - `Hermitian` (dense operator) or `Vector{<:Hermitian}` (block-diagonal operator, one
  `nr × nr` matrix per sample)
"""
posterior_covariance(c::Overdetermined, S::SmoothingOperator, D::AbstractMatrix{<:Real}) =
    Hermitian(inv(cholesky(Hermitian(c.HtH + weighted_gram(c, S, vec(D))))))

posterior_covariance(c::Underdetermined, ::IdentitySmoothing, D::AbstractMatrix{<:Real}) =
    Hermitian(woodbury_inverse(c.H, vec(D)))

function posterior_covariance(c::Underdetermined, S::InvertibleSmoothing, D::AbstractMatrix{<:Real})
    Σv = woodbury_inverse(standard_form(c, S), vec(D))         # covariance of v = L u
    X = smoothing_ldiv(S, Matrix(Σv), c.nr)                    # L⁻¹ Σv
    return Hermitian(Matrix(smoothing_ldiv(S, Matrix(X'), c.nr)'))   # L⁻¹ Σv L⁻ᴴ
end

posterior_covariance(c::Underdetermined, S::SingularSmoothing, D::AbstractMatrix{<:Real}) =
    Hermitian(inv(cholesky(Hermitian(gram(c) + weighted_gram(c, S, vec(D))))))

posterior_covariance(c::BlockCache, S::SmoothingOperator, D::AbstractMatrix{<:Real}) =
    [posterior_covariance(c.blocks[j], S, view(D, :, j:j)) for j in 1:c.ns]

"""
    residual_norm2(cache, U)

Squared (whitened) data misfit.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `U` - `AbstractMatrix`: unknown (`nr × ns`)

**Outputs**
- `r2` - `Real`: `‖y - H vec(U)‖₂²`
"""
residual_norm2(c::DenseCache, U) = sum(abs2, c.y - c.H * vec(U))
residual_norm2(c::BlockCache, U) = sum(residual_norm2(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)

"""
    marginal_quadratic(cache, U)

Quadratic form `yᴴ (I + H M⁻¹ Hᴴ)⁻¹ y = yᴴ (y - H û)` (Woodbury identity) of the marginal
likelihood, used to initialize the noise precision of the Bayesian formulation. Valid even
when `M = Lᴴ D L` is singular.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `U` - `AbstractMatrix`: weighted Tikhonov solution `û = (HᴴH + M)⁻¹ Hᴴ y`

**Outputs**
- `Q` - `Real`: value of the quadratic form
"""
marginal_quadratic(c::DenseCache, U) = real(dot(c.y, c.y - c.H * vec(U)))
marginal_quadratic(c::BlockCache, U) = sum(marginal_quadratic(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)
