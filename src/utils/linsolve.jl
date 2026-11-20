# Weighted Tikhonov solves -----------------------------------------------------------
#
# Every algorithm of the package boils down to repeated solves of
#
#     (Hᴴ H + Lᴴ D L) u = Hᴴ y        (D diagonal, positive)
#
# and, for the posterior covariance, to the inverse of the same matrix. The structure of H (dense / block-diagonal), the regime of the problem and the smoothing operator L are exploited through dispatch:
#
#   regime                    L                   system factorized
#   nm ≥ nu (over/square)     any                 normal equations           nu × nu
#   nm < nu (under)           I                   Woodbury identity          nm × nm
#   nm < nu (under)           square invertible   Woodbury on H̃ = H L⁻¹     nm × nm
#   nm < nu (under)           rectangular         normal equations           nu × nu
#                                                 (Lᴴ D L singular: no Woodbury)
#
# Woodbury identity (D invertible):
#     (Hᴴ H + D)⁻¹ = D⁻¹ - D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ H D⁻¹
#     (Hᴴ H + D)⁻¹ Hᴴ y = D⁻¹ Hᴴ (I + H D⁻¹ Hᴴ)⁻¹ y
# and, for L square invertible, (Hᴴ H + Lᴴ D L)⁻¹ = L⁻¹ (H̃ᴴ H̃ + D)⁻¹ L⁻ᴴ.
"""
Cache of the quantities reused across iterations for a dense operator. The type of `HtH` (`Hᴴ H`) encodes the regime: a matrix when `nm ≥ nu` (over-determined or square problem), `nothing` when `nm < nu` (under-determined problem, `Hᴴ H` is not needed).
"""
struct DenseCache{T, M<:AbstractMatrix{T}, G<:Union{Matrix{T},Nothing}}
    H::M
    y::Vector{T}
    HtH::G
    Hty::Vector{T}
    nu::Int
    ns::Int
end

"""
Over-determined or square problem (`nm ≥ nu`).
"""
const Overdetermined = DenseCache{<:Any, <:Any, <:Matrix}

"""
Under-determined problem (`nm < nu`).
"""
const Underdetermined = DenseCache{<:Any, <:Any, Nothing}

"""
Cache for a block-diagonal operator (one dense cache per sample).
"""
struct BlockCache{C<:DenseCache}
    blocks::Vector{C}
    nu::Int
    ns::Int
end

function dense_cache(H::AbstractMatrix, y::AbstractVecOrMat, nu, ns, σ2 = nothing)
    T = promote_type(eltype(H), eltype(y), Float64)
    Hc = Matrix{T}(H)
    yc = Vector{T}(vec(y))
    whiten!(Hc, yc, σ2)
    HtH = size(Hc, 1) >= size(Hc, 2) ? Hc' * Hc : nothing
    return DenseCache(Hc, yc, HtH, Hc' * yc, nu, ns)
end

"""
    linear_cache(prob, R = nothing)

Cache of the quantities reused across iterations (`H`, `y`, `Hᴴ H`, `Hᴴ y`), for the problem whitened by the noise covariance `R`.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `R` - noise covariance (see [`solve`](@ref)), `nothing` if not provided

**Outputs**
- `c` - `DenseCache` (dense operator) or `BlockCache` (block-diagonal operator)
"""
linear_cache(prob::InverseProblem{<:DenseOperator}, R = nothing) =
    dense_cache(prob.H.H, prob.y, prob.H.nu, prob.H.ns, noise_variances(R, prob.y))

function linear_cache(prob::InverseProblem{<:BlockDiagonalOperator}, R = nothing)
    nu, ns = unknown_size(prob)
    σ2 = noise_variances(R, prob.y)
    column(j) = σ2 === nothing ? nothing : view(σ2, :, j)
    block(j) = dense_cache(prob.H.blocks[j], view(prob.y, :, j), nu, 1, column(j))
    return BlockCache(block_caches(block, block(1), ns), nu, ns)
end

# All the blocks have the same size, hence the same regime and type: the vector is typed from
# the first block (function barrier), instead of a comprehension whose element type is only
# known at run time.
function block_caches(block, first_block::C, ns) where {C<:DenseCache}
    blocks = Vector{C}(undef, ns)
    blocks[1] = first_block
    for j in 2:ns
        blocks[j] = block(j)
    end
    return blocks
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

# Lᴴ diag(d) L (allocating version, used once for the posterior covariance)
weighted_gram(c::DenseCache, ::IdentitySmoothing, d) = Diagonal(d)
function weighted_gram(c::DenseCache, S::MatrixSmoothing, d)
    Lg = global_operator(c, S)
    return Matrix(Lg' * Diagonal(d) * Lg)
end


# H̃ = H L⁻¹ (standard form for a square invertible L)
standard_form(c::DenseCache, S::InvertibleSmoothing) = Matrix(smoothing_ldiv_adjoint(S, Matrix(c.H'), c.nu)')

# Linear systems with preallocated workspace ---------------------------------------------------
#
# The iterative algorithms solve the same system with new weights D at each iteration. A
# `LinearSystem` is built once per resolution: it precomputes everything that does not depend on D (Hᴴ H, transposed global smoothing operator, standard form H L⁻¹) and preallocates the buffers, so that `weighted_solve!` and `residual_norm2!` do not allocate large arrays.
"""
Linear system `(Hᴴ H + Lᴴ D L) u = Hᴴ y` of a dense operator, with its workspace. The method (normal equations or Woodbury identity) is encoded in the types of `G` and `K`.

**Fields**
- `c` - `DenseCache`: linear cache of the whitened problem
- `S` - `SmoothingOperator`: smoothing operator
- `G` - `Matrix` or `Nothing`: `Hᴴ H` for the normal equations, `nothing` for Woodbury
- `Lt` - `SparseMatrixCSC` or `Nothing`: transpose of the global smoothing operator (`nothing` for `L = I`)
- `K` - `Matrix` or `Nothing`: operator of the Woodbury identity (`H`, or `H L⁻¹` for a square invertible `L`), `nothing` for the normal equations
- `A` - `Matrix`: system matrix (`n × n` normal equations, `m × m` Woodbury)
- `KD` - `Matrix`: `K D⁻¹` (Woodbury), empty otherwise
- `b` - `Vector`: right-hand side buffer
- `x` - `Vector`: solution `vec(U)`
- `r` - `Vector`: residual `y - H vec(U)`
"""
struct DenseSystem{C<:DenseCache, S<:SmoothingOperator, T, G, Lt, K}
    c::C
    S::S
    G::G
    Lt::Lt
    K::K
    A::Matrix{T}
    KD::Matrix{T}
    b::Vector{T}
    x::Vector{T}
    r::Vector{T}
end

"""
System solved by the normal equations.
"""
const NormalSystem = DenseSystem{<:Any, <:Any, <:Any, <:AbstractMatrix, <:Any, Nothing}

"""
System solved by the Woodbury identity.
"""
const WoodburySystem = DenseSystem{<:Any, <:Any, <:Any, Nothing, <:Any, <:AbstractMatrix}

"""
Linear systems of a block-diagonal operator (one `DenseSystem` per sample, each with its own workspace, so that the samples can be solved in parallel).
"""
struct BlockSystem{C<:BlockCache, B<:DenseSystem}
    c::C
    blocks::Vector{B}
end

# Method of resolution: (Hᴴ H for the normal equations, Woodbury operator)
system_method(c::Overdetermined, ::SmoothingOperator) = (c.HtH, nothing)
system_method(c::Underdetermined, ::IdentitySmoothing) = (nothing, c.H)
system_method(c::Underdetermined, S::InvertibleSmoothing) = (nothing, standard_form(c, S))
system_method(c::Underdetermined, ::SingularSmoothing) = (gram(c), nothing)

transposed_operator(c::DenseCache, ::IdentitySmoothing) = nothing
transposed_operator(c::DenseCache, S::MatrixSmoothing) = sparse(transpose(global_operator(c, S)))

"""
    linear_system(cache, L)

Linear system of the weighted Tikhonov problem, with precomputed operators and preallocated workspace.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: `IdentitySmoothing()` or `MatrixSmoothing`

**Outputs**
- `sys` - `DenseSystem` or `BlockSystem`
"""
function linear_system(c::DenseCache{T}, S::SmoothingOperator) where {T}
    G, K = system_method(c, S)
    m, n = size(c.H)
    p = K === nothing ? n : m
    KD = K === nothing ? Matrix{T}(undef, 0, 0) : similar(K)
    return DenseSystem(c, S, G, transposed_operator(c, S), K, Matrix{T}(undef, p, p), KD,
                       Vector{T}(undef, p), Vector{T}(undef, n), Vector{T}(undef, m))
end

linear_system(c::BlockCache, S::SmoothingOperator) = BlockSystem(c, [linear_system(b, S) for b in c.blocks])

# A += Lᴴ diag(d) L, with Lt = Lᵀ (columns of Lt = rows of L)
add_weighted_gram!(A, ::Nothing, d) = (for i in eachindex(d); A[i, i] += d[i]; end; A)
function add_weighted_gram!(A, Lt::SparseMatrixCSC, d)
    rv, nz = rowvals(Lt), nonzeros(Lt)
    for i in axes(Lt, 2)
        rng = nzrange(Lt, i)
        for k in rng, j in rng
            A[rv[j], rv[k]] += conj(nz[j]) * d[i] * nz[k]
        end
    end
    return A
end

# u = L⁻¹ v after a Woodbury solve on the standard form
from_standard_form!(x, ::SmoothingOperator, nu) = x
from_standard_form!(x, S::InvertibleSmoothing, nu) = (ldiv!(S.F, reshape(x, nu, :)); x)

# vec(U) of the solution, in sys.x
function solve_system!(sys::NormalSystem, d)
    copyto!(sys.A, sys.G)
    add_weighted_gram!(sys.A, sys.Lt, d)
    copyto!(sys.x, sys.c.Hty)
    ldiv!(cholesky!(Hermitian(sys.A)), sys.x)
    return sys.x
end

function solve_system!(sys::WoodburySystem, d)
    K, KD, A = sys.K, sys.KD, sys.A
    KD .= K ./ transpose(d)                                    # K D⁻¹
    mul!(A, KD, K')                                            # K D⁻¹ Kᴴ
    for i in axes(A, 1)
        A[i, i] += 1
    end
    copyto!(sys.b, sys.c.y)
    ldiv!(cholesky!(Hermitian(A)), sys.b)                     # (I + K D⁻¹ Kᴴ)⁻¹ y
    mul!(sys.x, KD', sys.b)                                    # D⁻¹ Kᴴ (…)⁻¹ y
    return from_standard_form!(sys.x, sys.S, sys.c.nu)
end

"""
    weighted_solve!(U, sys, D)

In-place solution of the weighted Tikhonov system `(Hᴴ H + Lᴴ diag(vec(D)) L) vec(U) = Hᴴ y`, using the workspace of `sys`. The normal equations or the Woodbury identity are used depending on the regime of the problem and on the smoothing operator (see the table at the top of `linsolve.jl`).

**Inputs**
- `U` - `AbstractMatrix`: output (`nu × ns`), overwritten
- `sys` - `DenseSystem` or `BlockSystem`: linear system (see [`linear_system`](@ref))
- `D` - `AbstractMatrix{<:Real}`: positive diagonal weights, with the size of `L U`

**Outputs**
- `U` - `AbstractMatrix`: solution
"""
function weighted_solve!(U, sys::DenseSystem, D::AbstractMatrix{<:Real})
    copyto!(U, solve_system!(sys, vec(D)))
    return U
end

function weighted_solve!(U, sys::BlockSystem, D::AbstractMatrix{<:Real})
    Threads.@threads for j in eachindex(sys.blocks)
        copyto!(view(U, :, j), solve_system!(sys.blocks[j], view(D, :, j)))
    end
    return U
end

"""
    weighted_solve(cache, L, D)

Allocating version of [`weighted_solve!`](@ref), for one-off solves.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: `IdentitySmoothing()` or `MatrixSmoothing`
- `D` - `AbstractMatrix{<:Real}`: positive diagonal weights, with the size of `L U`

**Outputs**
- `U` - `Matrix`: solution (`nu × ns`)
"""
weighted_solve(c, S::SmoothingOperator, D::AbstractMatrix{<:Real}) =
    weighted_solve!(Matrix{eltype(c)}(undef, c.nu, c.ns), linear_system(c, S), D)

"""
    residual_norm2!(sys, U)

Squared (whitened) data misfit, computed in the workspace of `sys`.

**Inputs**
- `sys` - `DenseSystem` or `BlockSystem`: linear system
- `U` - `AbstractMatrix`: unknown (`nu × ns`)

**Outputs**
- `r2` - `Real`: `‖y - H vec(U)‖₂²`
"""
function residual_norm2!(sys::DenseSystem, U)
    r = copyto!(sys.r, sys.c.y)
    mul!(r, sys.c.H, vec(U), -1, 1)
    return sum(abs2, r)
end

residual_norm2!(sys::BlockSystem, U) =
    sum(residual_norm2!(sys.blocks[j], view(U, :, j)) for j in eachindex(sys.blocks))

"""
    residual_norm2(cache, U)

Allocating version of [`residual_norm2!`](@ref).

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `U` - `AbstractMatrix`: unknown (`nu × ns`)

**Outputs**
- `r2` - `Real`: `‖y - H vec(U)‖₂²`
"""
residual_norm2(c::DenseCache, U) = sum(abs2, c.y - c.H * vec(U))
residual_norm2(c::BlockCache, U) = sum(residual_norm2(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)

"""
    marginal_quadratic(cache, U)

Quadratic form `yᴴ (I + H M⁻¹ Hᴴ)⁻¹ y = yᴴ (y - H û)` (Woodbury identity) of the marginal likelihood, used to initialize the noise precision of the Bayesian formulation. Valid even when `M = Lᴴ D L` is singular.

**Inputs**
- `cache` - `DenseCache` or `BlockCache`: linear cache of the whitened problem
- `U` - `AbstractMatrix`: weighted Tikhonov solution `û = (HᴴH + M)⁻¹ Hᴴ y`

**Outputs**
- `Q` - `Real`: value of the quadratic form
"""
marginal_quadratic(c::DenseCache, U) = real(dot(c.y, c.y - c.H * vec(U)))
marginal_quadratic(c::BlockCache, U) = sum(marginal_quadratic(c.blocks[j], view(U, :, j:j)) for j in 1:c.ns)
