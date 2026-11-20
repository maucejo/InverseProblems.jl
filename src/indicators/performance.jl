# Performance indicators ----------------------------------------------------------------------
#
# Quality of a reconstruction `u` with respect to a reference `uref`, in the sense of ISO 5725: trueness (closeness of `u` to `uref`), precision (spread of the posterior distribution) and accuracy (both combined).
#
# The uncertainty-based indicators only need the posterior *variances* of the n = length(u)
# unknowns, i.e. the diagonal of the posterior covariance Σ = Q⁻¹. They are averaged over the
# unknowns, so that the layout of `u` (vector, `nu × ns` matrix) does not matter. The variances are taken from:
#   * an `InverseSolution` obtained with `solve(...; precision = true)`;
#   * a precision matrix Q (default) or a covariance matrix Σ (`covariance = true`): one matrix for all the unknowns, or a vector of blocks, one per sample (frequency or time step) — e.g. `sol.precision`, `posterior_covariance(sol)`, or the covariances P_k of a Kalman filter.

"""
    relative_error(uref, u; p = 1)

Relative error (trueness) of the reconstruction `u` with respect to the reference `uref`, in percent:

    RE = 100 ‖uref - u‖_p / ‖uref‖_p

The norm is taken over all the entries (all the samples for a space–frequency or space–time unknown). Apply it to a column or a row to get a local indicator, e.g. `relative_error(uref[:, k], u[:, k])` at the `k`-th frequency.

**Inputs**
- `uref` - `AbstractArray`: reference field
- `u` - `AbstractArray` or `InverseSolution`: reconstructed field, with the layout of `uref`
- `p` - `Real`: order of the norm (`1` by default, less sensitive to a few large errors than `2`)

**Outputs**
- `RE` - `Real`: relative error, in %
"""
function relative_error(uref::AbstractArray, u::AbstractArray; p::Real = 1)
    check_layout(uref, u)
    return 100 * norm(uref - u, p) / norm(uref, p)
end

relative_error(uref::AbstractArray, sol::InverseSolution; kw...) = relative_error(uref, sol.u; kw...)

"""
    correlation(uref, u)

Correlation (cosine similarity of shape) between the reconstruction `u` and the reference `uref`, in percent:

    C = 100 |⟨uref, u⟩| / (‖uref‖₂ ‖u‖₂)

It equals 100 % when `u` is proportional to `uref`, whatever the scale factor, and is insensitive to a global phase for complex fields.

**Inputs**
- `uref` - `AbstractArray`: reference field
- `u` - `AbstractArray` or `InverseSolution`: reconstructed field, with the layout of `uref`

**Outputs**
- `C` - `Real`: correlation, in % (`0 ≤ C ≤ 100`)
"""
function correlation(uref::AbstractArray, u::AbstractArray)
    check_layout(uref, u)
    return 100 * abs(dot(uref, u)) / (norm(uref) * norm(u))
end

correlation(uref::AbstractArray, sol::InverseSolution) = correlation(uref, sol.u)


"""
    precision_metric(sol)
    precision_metric(Q)
    precision_metric(Σ; covariance = true)

Precision of the estimate, in the sense of ISO 5725 (spread of the estimate, not to be confused with the precision matrix `Q`): root mean posterior variance of the `n` unknowns,

    P = √( tr(Σ) / n ),    Σ = Q⁻¹,

i.e. the quadratic mean of the posterior standard deviations. It has the unit of the unknown; the smaller, the more precise. It does not depend on the layout of the unknown: for a space–frequency or space–time problem, the mean runs over all the entries of `u` (all points, all samples), and the correlations between unknowns are not used.

The matrices are real (real data) or complex Hermitian (complex data); in both cases the variances are real. The variances `diag(Q⁻¹)` are computed from a Cholesky factorization of `Q` (dense, or sparse for a `SparseMatrixCSC`), without forming `Q⁻¹`.

**Inputs**
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)` (variances by block solves, see [`posterior_std`](@ref))
- `Q` - posterior precision matrix of `vec(u)` (`AbstractMatrix`, real symmetric or complex Hermitian positive definite: dense, `Hermitian`/`Symmetric`, `SparseMatrixCSC` or `Diagonal`), or vector of precision blocks, one per sample when the samples are independent (e.g. `sol.precision`)
- `Σ` - same, for a posterior covariance, real or complex (e.g. [`posterior_covariance`](@ref), Kalman filter covariances `P_k`)
- `covariance` - `Bool`: `true` if the matrices are covariances, `false` (default) if they are precisions

**Outputs**
- `P` - `Real`: precision indicator, in the unit of the unknown
"""
precision_metric(sol::InverseSolution) = rms(posterior_variances(sol))
precision_metric(Q; covariance::Bool = false) = rms(posterior_variances(Q, covariance))

"""
    accuracy_metric(uref, sol)
    accuracy_metric(uref, u, Q)
    accuracy_metric(uref, u, Σ; covariance = true)

Accuracy of the estimate, in the sense of ISO 5725 (trueness and precision combined): root mean squared error, per unknown, of the posterior distribution `N(u, Σ)` with respect to the reference,

    A = √( (‖uref - u‖₂² + tr(Σ)) / n ),    Σ = Q⁻¹,

where `n = length(u)`: `‖uref - u‖₂² + tr(Σ)` is the expectation of `‖uref - x‖₂²` for `x ~ N(u, Σ)` (squared bias plus variance). Without uncertainty (`Σ = 0`) it reduces to the RMS error of `u`, and `A² = RMS(uref - u)² + P²` with `P` the [`precision_metric`](@ref). It has the unit of the unknown; the smaller, the more accurate. As for `precision_metric`, the mean runs over all the entries of `u`, whatever its layout.

**Inputs**
- `uref` - `AbstractArray`: reference field (vector, or `nu × ns` matrix)
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `u` - `AbstractArray`: reconstructed field, with the layout of `uref`
- `Q` / `Σ` - posterior precision (default) or covariance (`covariance = true`) of `vec(u)`, real or complex: one `n × n` matrix, or a vector of blocks, one per sample; see [`precision_metric`](@ref)
- `covariance` - `Bool`: `true` if the matrices are covariances, `false` (default) if they are precisions

**Outputs**
- `A` - `Real`: accuracy indicator, in the unit of the unknown
"""
accuracy_metric(uref::AbstractArray, sol::InverseSolution) =
    accuracy_metric(uref, sol.u, posterior_variances(sol))

accuracy_metric(uref::AbstractArray, u::AbstractArray, Q; covariance::Bool = false) =
    accuracy_metric(uref, u, posterior_variances(Q, covariance))

function accuracy_metric(uref::AbstractArray, u::AbstractArray, v::AbstractVector{<:Real})
    check_layout(uref, u)
    length(v) == length(u) ||
        throw(DimensionMismatch("the precision/covariance describes $(length(v)) unknowns, the field has $(length(u))"))
    return sqrt((sum(abs2, uref - u) + sum(v)) / length(u))
end

"""
    performance_ratio(indicators)

Overall performance of several methods compared on several indicators: geometric mean, over the indicators, of the ratio between the indicator of each method and the best one,

    r_m = ( Π_j I_mj / min_m' I_m'j )^(1 / n_ind).

The best possible value is `1`, reached by a method that is the best for every indicator. The indicators must all be such that *smaller is better* (e.g. [`relative_error`](@ref), [`precision_metric`](@ref), [`accuracy_metric`](@ref), or `100 - correlation`) and strictly positive.

**Inputs**
- `indicators` - `AbstractMatrix{<:Real}`: one row per method, one column per indicator

**Outputs**
- `r` - `Vector{<:Real}`: performance ratio of each method (`≥ 1`)

**Example**
```julia
I = [relative_error(uref, s) for s in sols] # one column per indicator
J = [accuracy_metric(uref, s) for s in sols]
performance_ratio([I J])
```
"""
function performance_ratio(indicators::AbstractMatrix{<:Real})
    all(>(0), indicators) || throw(ArgumentError("the indicators must be strictly positive"))
    logratio = log.(indicators) .- log.(minimum(indicators, dims = 1))
    return vec(exp.(sum(logratio, dims = 2) ./ size(indicators, 2)))
end

# Posterior variances (internal) ---------------------------------------------------------------

rms(v::AbstractVector{<:Real}) = sqrt(sum(v) / length(v))

# Variances of all the unknowns, as a vector (order of vec(u))
posterior_variances(sol::InverseSolution) = vec(abs2.(posterior_std(sol)))

posterior_variances(M::AbstractMatrix, covariance::Bool) =
    covariance ? real(diag(M)) : inverse_diagonal(M)

# One block per sample: variances stacked sample by sample, as vec(u)
function posterior_variances(Ms::AbstractVector{<:AbstractMatrix}, covariance::Bool)
    isempty(Ms) && throw(ArgumentError("empty vector of matrices"))
    return reduce(vcat, [posterior_variances(M, covariance) for M in Ms])
end

# diag(Q⁻¹) for a Hermitian positive definite Q, real (symmetric) or complex (Hermitian), without forming Q⁻¹. Wrappers are expanded with `Matrix` / `sparse`, which honour the stored triangle.
inverse_diagonal(Q::Diagonal) = 1 ./ real(Q.diag)
inverse_diagonal(Q::Union{Hermitian{<:Any,<:AbstractSparseMatrix},Symmetric{<:Real,<:AbstractSparseMatrix}}) =
    inverse_diagonal(sparse(Q))

function inverse_diagonal(Q::AbstractMatrix)
    LinearAlgebra.checksquare(Q)
    # Q = Rᴴ R  ⟹  (Q⁻¹)_ii = ‖R⁻ᴴ e_i‖²: forward substitution Rᴴ X = I (one triangular solve, no explicit inverse), then squared norms of the columns of X
    F = cholesky(Hermitian(Matrix(Q)))
    X = ldiv!(F.L, Matrix{eltype(F)}(I, size(F)...))
    return vec(sum(abs2, X; dims = 1))
end

# Sparse precision (assembled space–frequency or space–time problems): sparse Cholesky, then solves with chunks of columns of the identity, keeping only the diagonal entries
function inverse_diagonal(Q::AbstractSparseMatrix; chunk::Int = 256)
    n = LinearAlgebra.checksquare(Q)
    F = cholesky(Hermitian(SparseMatrixCSC(Q)))
    T = eltype(F)
    d = Vector{real(T)}(undef, n)
    for j0 in 1:chunk:n
        cols = j0:min(j0 + chunk - 1, n)
        E = zeros(T, n, length(cols))
        for (k, i) in enumerate(cols)
            E[i, k] = one(T)
        end
        X = F \ E
        for (k, i) in enumerate(cols)
            d[i] = real(X[i, k])
        end
    end
    return d
end

check_layout(uref, u) = size(uref) == size(u) ||
    throw(DimensionMismatch("reference of size $(size(uref)), reconstruction of size $(size(u))"))
