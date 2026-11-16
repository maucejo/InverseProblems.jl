# Noise covariance ---------------------------------------------------------------------
#
# The noise covariance R is known and diagonal. The data fidelity ‖y - Hu‖²_{R⁻¹} is handled once for all by whitening: H ← R^{-1/2} H, y ← R^{-1/2} y.
"""
    noise_variances(R, y)

Diagonal of the known noise covariance `R`, with the layout of the data.

**Inputs**
- `R` - `Nothing`, `UniformScaling`, `Diagonal`, `AbstractVector{<:Real}` or
  `AbstractMatrix{<:Real}`: noise covariance — not provided, `σ²I`, variances of all the
  data (vector of length `length(y)`, matrix of the size of `y`, or diagonal matrix
  `length(y) × length(y)`), or, for matrix data `y` (`nm × ns`), variances of the `nm`
  measurement points common to all the samples (vector of length `nm`, `Diagonal` or diagonal
  matrix `nm × nm`). A matrix of the size of `y` is always read as the variances of all the
  data (use a vector or a `Diagonal` for per-point variances when `nm = ns`)
- `y` - `AbstractVecOrMat`: measured data

**Outputs**
- `σ2` - `Array{<:Real}` or `Nothing`: noise variances with the layout of `y`, `nothing` when `R` is not provided or is `I`
"""
noise_variances(::Nothing, y) = nothing   # R not provided: R = I
noise_variances(R::UniformScaling, y) =
    isone(R.λ) ? nothing : check_variances(fill(float(R.λ), size(y)))
noise_variances(R::Diagonal, y) = noise_variances(R.diag, y)

function noise_variances(R::AbstractVector{<:Real}, y)
    length(R) == length(y) && return check_variances(reshape(float.(collect(R)), size(y)))
    y isa AbstractMatrix && length(R) == size(y, 1) &&             # per measurement point
        return check_variances(repeat(float.(collect(R)), 1, size(y, 2)))
    throw(DimensionMismatch("$(length(R)) noise variances for $(length(y)) data"))
end

function noise_variances(R::AbstractMatrix{<:Real}, y)
    y isa AbstractMatrix && size(R) == size(y) && return check_variances(float.(R))
    (size(R) == (length(y), length(y)) || size(R) == (size(y, 1), size(y, 1))) ||
        throw(DimensionMismatch("noise covariance of size $(size(R)) for $(length(y)) data"))
    isdiag(R) || throw(ArgumentError("the noise covariance R must be diagonal"))
    return noise_variances(diag(R), y)
end

check_variances(σ2) = all(>(0), σ2) ? σ2 : throw(ArgumentError("noise variances must be positive"))

whiten!(H, y, ::Nothing) = nothing
function whiten!(H, y, σ2)
    s = 1 ./ sqrt.(vec(σ2))
    H .*= s
    y .*= s
    return nothing
end
