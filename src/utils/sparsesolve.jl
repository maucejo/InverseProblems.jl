# Sparse operators -------------------------------------------------------------------------
#
# Assembled problems with a sparse forward operator, e.g. a vectorized space–frequency
# problem (block-diagonal H of size (ny·nf) × (nu·nf)) whose frequencies are coupled by the
# regularization. The weighted Tikhonov systems (Hᴴ H + Lᴴ D L) u = Hᴴ y stay sparse and are
# solved with sparse Cholesky factorizations (CHOLMOD); the symbolic analysis is computed at
# the first solve and reused at the following iterations (same sparsity pattern).
#
# The parameter choice rules based on the SVD (L-curve, GCV, Bayesian estimator) are not
# available; the heuristic initial parameter uses cheap estimates of the extreme eigenvalues
# of Hᴴ H.

"""
Cache of the quantities reused across iterations for a sparse operator.

**Fields**
- `H` - `SparseMatrixCSC`: whitened forward operator
- `y` - `Vector`: whitened data
- `HtH` - `SparseMatrixCSC`: `Hᴴ H`
- `Hty` - `Vector`: `Hᴴ y`
- `nr`, `ns` - `Int`: size of the unknown (`nr × ns`)
"""
struct SparseCache{T, M<:SparseMatrixCSC{T}}
    H::M
    y::Vector{T}
    HtH::SparseMatrixCSC{T,Int}
    Hty::Vector{T}
    nr::Int
    ns::Int
end

sparse_whiten(H, y, ::Nothing) = (H, y)
function sparse_whiten(H, y, σ2)
    s = 1 ./ sqrt.(vec(σ2))
    return Diagonal(s) * H, y .* s
end

function sparse_cache(H::AbstractSparseMatrix, y::AbstractVecOrMat, nr, ns, σ2 = nothing)
    T = promote_type(eltype(H), eltype(y), Float64)
    Hc, yc = sparse_whiten(SparseMatrixCSC{T,Int}(H), Vector{T}(vec(y)), σ2)
    return SparseCache(Hc, yc, sparse(Hc' * Hc), Hc' * yc, nr, ns)
end

linear_cache(prob::InverseProblem{<:SparseOperator}, R = nothing) =
    sparse_cache(prob.H.H, prob.y, prob.H.nr, prob.H.ns, noise_variances(R, prob.y))

Base.eltype(::SparseCache{T}) where {T} = T
ndata(c::SparseCache) = length(c.y)

global_operator(c::SparseCache, S::MatrixSmoothing) =
    c.ns == 1 ? S.L : kron(sparse(I, c.ns, c.ns), S.L)

# Lᴴ diag(d) L as a sparse matrix
sparse_weighted_gram(c::SparseCache, ::IdentitySmoothing, d) = spdiagm(0 => d)
function sparse_weighted_gram(c::SparseCache, S::MatrixSmoothing, d)
    Lg = global_operator(c, S)
    return Lg' * spdiagm(0 => d) * Lg
end

# Hᴴ H + Lᴴ diag(d) L
sparse_system_matrix(c::SparseCache, S, d) = c.HtH + sparse_weighted_gram(c, S, d)

"""
Linear system of a sparse operator, with the sparse Cholesky factorization reused from one iteration to the next.

**Fields**
- `c` - `SparseCache`: linear cache of the whitened problem
- `S` - `SmoothingOperator`: smoothing operator
- `F` - `Base.RefValue{Any}`: sparse Cholesky factorization (`nothing` before the first solve)
- `x` - `Vector`: solution `vec(U)`
- `r` - `Vector`: residual `y - H vec(U)`
"""
struct SparseSystem{C<:SparseCache, S<:SmoothingOperator, T}
    c::C
    S::S
    F::Base.RefValue{Any}
    x::Vector{T}
    r::Vector{T}
end

linear_system(c::SparseCache{T}, S::SmoothingOperator) where {T} =
    SparseSystem(c, S, Ref{Any}(nothing), Vector{T}(undef, size(c.H, 2)), Vector{T}(undef, size(c.H, 1)))

# Numeric factorization, reusing the symbolic analysis of the first call
function sparse_factor!(sys::SparseSystem, A)
    if sys.F[] === nothing
        sys.F[] = cholesky(Hermitian(A))
    else
        cholesky!(sys.F[], Hermitian(A))
    end
    return sys.F[]
end

function weighted_solve!(U, sys::SparseSystem, D::AbstractMatrix{<:Real})
    F = sparse_factor!(sys, sparse_system_matrix(sys.c, sys.S, vec(D)))
    sys.x .= F \ sys.c.Hty
    copyto!(U, sys.x)
    return U
end

function residual_norm2!(sys::SparseSystem, U)
    r = copyto!(sys.r, sys.c.y)
    mul!(r, sys.c.H, vec(U), -1, 1)
    return sum(abs2, r)
end

residual_norm2(c::SparseCache, U) = sum(abs2, c.y - c.H * vec(U))
marginal_quadratic(c::SparseCache, U) = real(dot(c.y, c.y - c.H * vec(U)))

"""
    posterior_covariance(cache::SparseCache, L, D)

Posterior covariance `(Hᴴ H + Lᴴ diag(vec(D)) L)⁻¹` for a sparse operator (dense result: its
size is the square of the number of unknowns).

**Inputs**
- `cache` - `SparseCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: smoothing operator
- `D` - `AbstractMatrix{<:Real}`: scaled weights `λ W` of the last linear solve

**Outputs**
- `Σ` - `Hermitian`: posterior covariance of `vec(u)`
"""
function posterior_covariance(c::SparseCache{T}, S::SmoothingOperator, D::AbstractMatrix{<:Real}) where {T}
    F = cholesky(Hermitian(sparse_system_matrix(c, S, vec(D))))
    n = size(c.H, 2)
    return Hermitian(Matrix(F \ Matrix{T}(I, n, n)))
end

spectrum(c::SparseCache, S::SmoothingOperator, W::AbstractMatrix{<:Real}) =
    throw(ArgumentError("the parameter choice rules (L-curve, GCV, Bayesian estimator) need the SVD of H, " *
                        "which is not available for sparse operators: use a fixed λ, the multiplicative " *
                        "or Bayesian formulations, and `init = Heuristic()` or a number"))

"""
    normest(apply, apply_adjoint, x0; tol = 1e-6, maxiter = 100)

Estimate of the 2-norm of a linear operator `M` by the power method on `Mᴴ M` (as `normest`).

**Inputs**
- `apply` - `Function`: `x ↦ M x`
- `apply_adjoint` - `Function`: `x ↦ Mᴴ x`
- `x0` - `AbstractVector`: starting vector
- `tol` - `Real`: tolerance on the relative variation of the estimate
- `maxiter` - `Integer`: maximum number of iterations

**Outputs**
- `ν` - `Real`: estimate of `‖M‖₂` (lower bound)
"""
function normest(apply, apply_adjoint, x0; tol = 1e-6, maxiter = 100)
    ν = norm(x0)
    iszero(ν) && return ν
    x = x0 / ν
    ν0 = zero(ν)
    for _ in 1:maxiter
        abs(ν - ν0) <= tol * ν && break
        ν0 = ν
        Mx = apply(x)
        iszero(norm(Mx)) && (Mx = rand(eltype(Mx), length(Mx)))
        x = apply_adjoint(Mx)
        ν = norm(x) / norm(Mx)
        x ./= norm(x)
    end
    return ν
end

"""
    condest(A)

Estimate of the 2-norm condition number `‖A‖₂ ‖A⁻¹‖₂` of a square sparse matrix (as `condest(A, p = 2)`): both norms are estimated by the power method ([`normest`](@ref)), the inverse being applied through a sparse LU factorization.

**Inputs**
- `A` - `AbstractSparseMatrix`: square matrix

**Outputs**
- `κ` - `Real`: estimate of `cond₂(A)` (`Inf` if `A` is singular)
"""
function condest(A::AbstractSparseMatrix)
    n = LinearAlgebra.checksquare(A)
    n == 0 && return 0.0
    F = lu(A; check = false)
    issuccess(F) || return Inf
    x0 = vec(sum(abs, A; dims = 1))
    Anorm = normest(x -> A * x, x -> A' * x, x0)
    Ainvnorm = normest(x -> F \ x, x -> F' \ x, x0)
    return Anorm * Ainvnorm
end

"""
    select_parameter(::Heuristic, cache::SparseCache, L, W)

Heuristic initial parameter for a sparse operator (MSSP 2017, MSSP 2018 Eqs. (26)–(27)):
geometric mean of `s_max = ‖A‖₂` and `s_min = 1 / ‖A⁻¹‖₂`, the extreme eigenvalues of

* `A = (H L⁻¹)ᴴ (H L⁻¹)` for a square invertible smoothing operator `L` — `A` and `A⁻¹ = L (Hᴴ H)⁻¹ Lᴴ`
  are applied as operators (LU factorizations), not assembled;
* `A = Hᴴ H` for `L = I`, or for a rectangular / singular `L` (no standard form).

Both norms are estimated by the power method ([`normest`](@ref), as `condest` in norm 2);
`s_min = 0` if `Hᴴ H` is singular. The weights are not taken into account.

**Inputs**
- `cache` - `SparseCache`: linear cache of the whitened problem
- `L` - `SmoothingOperator`: smoothing operator
- `W` - `AbstractMatrix{<:Real}`: weights (not used)

**Outputs**
- `α0` - `Real`: initial regularization parameter
"""
function select_parameter(::Heuristic, c::SparseCache, L::SmoothingOperator, W::AbstractMatrix{<:Real})
    apply_A, apply_Ainv, x0 = heuristic_operators(c, L)
    smax = normest(apply_A, apply_A, x0)
    smin = apply_Ainv === nothing ? 0.0 : 1 / normest(apply_Ainv, apply_Ainv, x0)
    return sqrt(smax * max(smin, smax * 16eps()))
end

# A = Hᴴ H and A⁻¹ (nothing if Hᴴ H is singular), as Hermitian operators, and starting vector
function heuristic_operators(c::SparseCache, ::SmoothingOperator)
    A = c.HtH
    F = lu(A; check = false)
    Ainv = issuccess(F) ? (x -> F \ x) : nothing
    return (x -> A * x), Ainv, vec(sum(abs, A; dims = 1))
end

# A = (H L⁻¹)ᴴ (H L⁻¹), A⁻¹ = L (Hᴴ H)⁻¹ Lᴴ
function heuristic_operators(c::SparseCache, S::InvertibleSmoothing)
    Lg = global_operator(c, S)
    F = lu(c.HtH; check = false)
    apply_A(x) = smoothing_ldiv_adjoint(S, c.HtH * smoothing_ldiv(S, x, c.nr), c.nr)
    Ainv = issuccess(F) ? (x -> Lg * (F \ (Lg' * x))) : nothing
    return apply_A, Ainv, ones(real(eltype(c)), size(c.HtH, 1))
end
