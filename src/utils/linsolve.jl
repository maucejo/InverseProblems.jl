# Weighted Tikhonov solves -----------------------------------------------------------
#
# Every algorithm of the package boils down to repeated solves of
#
#     (Hᴴ H + Lᴴ D L) u = Hᴴ y        (D diagonal, positive)
#
# The structure of H (dense / block-diagonal) and of L (identity / matrix) is exploited
# through dispatch. With L = I, the smallest of the two normal systems is factorized:
#
#     nm ≥ nr :  (Hᴴ H + D)⁻¹ Hᴴ y
#     nm < nr :  D⁻¹ Hᴴ (H D⁻¹ Hᴴ + I)⁻¹ y            (push-through identity)

"""
Cache of the quantities reused across iterations for a dense operator.
"""
struct DenseCache{T, M<:AbstractMatrix{T}}
    H::M
    y::Vector{T}
    HᴴH::Matrix{T}
    Hᴴy::Vector{T}
    nr::Int
    ns::Int
end

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
    return DenseCache(Hc, yc, Hc' * Hc, Hc' * yc, nr, ns)
end

"""
    linear_cache(prob, R = nothing)

Cache of the (whitened) problem for the noise covariance `R`.
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

"""
    weighted_solve(cache, L, D)

Solution `U` (`nr × ns`) of `(Hᴴ H + Lᴴ diag(vec(D)) L) vec(U) = Hᴴ y`, where `D` has the
size of `L U`.
"""
function weighted_solve(c::DenseCache, ::IdentitySmoothing, D::AbstractMatrix{<:Real})
    d = vec(D)
    x = if size(c.H, 1) >= size(c.H, 2)
        cholesky!(Hermitian(c.HᴴH + Diagonal(d))) \ c.Hᴴy
    else
        HD = c.H ./ transpose(d)                                # H D⁻¹
        HD' * (cholesky!(Hermitian(HD * c.H' + I)) \ c.y)      # D⁻¹ Hᴴ (H D⁻¹ Hᴴ + I)⁻¹ y
    end
    return reshape(x, c.nr, c.ns)
end

function weighted_solve(c::DenseCache, S::MatrixSmoothing, D::AbstractMatrix{<:Real})
    Lg = global_operator(c, S)
    A = c.HᴴH + Matrix(Lg' * Diagonal(vec(D)) * Lg)
    return reshape(cholesky!(Hermitian(A)) \ c.Hᴴy, c.nr, c.ns)
end

function weighted_solve(c::BlockCache, S::SmoothingOperator, D::AbstractMatrix{<:Real})
    U = Matrix{eltype(c)}(undef, c.nr, c.ns)
    Threads.@threads for j in 1:c.ns
        U[:, j] = weighted_solve(c.blocks[j], S, view(D, :, j:j))
    end
    return U
end

"""
    posterior_covariance(cache, L, D)

`(Hᴴ H + Lᴴ diag(vec(D)) L)⁻¹` (whitened `H`): covariance of the Gaussian approximation of
the posterior at a fixed point of the IR algorithm (exact for `q = 2`). For a block-diagonal
operator, a vector with one `nr × nr` matrix per sample is returned.
"""
posterior_covariance(c::DenseCache, ::IdentitySmoothing, D::AbstractMatrix{<:Real}) =
    Hermitian(inv(cholesky(Hermitian(c.HᴴH + Diagonal(vec(D))))))

function posterior_covariance(c::DenseCache, S::MatrixSmoothing, D::AbstractMatrix{<:Real})
    Lg = global_operator(c, S)
    return Hermitian(inv(cholesky(Hermitian(c.HᴴH + Matrix(Lg' * Diagonal(vec(D)) * Lg)))))
end

posterior_covariance(c::BlockCache, S::SmoothingOperator, D::AbstractMatrix{<:Real}) =
    [posterior_covariance(c.blocks[j], S, view(D, :, j:j)) for j in 1:c.ns]

"""
    residual_norm2(cache, U)

Squared data misfit `‖y - H vec(U)‖₂²`.
"""
residual_norm2(c::DenseCache, U) = sum(abs2, c.y - c.H * vec(U))
residual_norm2(c::BlockCache, U) = sum(residual_norm2(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)

"""
    marginal_quadratic(cache, U)

`yᴴ (I + H M⁻¹ Hᴴ)⁻¹ y = yᴴ (y - H û)` (Woodbury identity), where `û` is the weighted
Tikhonov solution `(HᴴH + M)⁻¹ Hᴴ y`. Used to initialize the noise precision of the Bayesian
formulation (maximum of the marginal likelihood). Valid even when `M = Lᴴ D L` is singular.
"""
marginal_quadratic(c::DenseCache, U) = real(dot(c.y, c.y - c.H * vec(U)))
marginal_quadratic(c::BlockCache, U) = sum(marginal_quadratic(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)
