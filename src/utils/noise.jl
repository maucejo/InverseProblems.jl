# Noise covariance ---------------------------------------------------------------------
#
# The noise covariance R is known and diagonal. The data fidelity ‖y - Hu‖²_{R⁻¹} is handled
# once for all by whitening: H ← R^{-1/2} H, y ← R^{-1/2} y.
"""
    noise_variances(R, y)

Diagonal of the noise covariance `R`, as an array with the layout of `y`, or `nothing` when
`R = I`. `R` can be `I` (or `σ²I`), a `Diagonal`, a vector of variances (`length(y)`), a
matrix of variances with the size of `y`, or a diagonal `length(y) × length(y)` matrix.
"""
noise_variances(::Nothing, y) = nothing   # R not provided: R = I
noise_variances(R::UniformScaling, y) =
    isone(R.λ) ? nothing : check_variances(fill(float(R.λ), size(y)))
noise_variances(R::Diagonal, y) = noise_variances(R.diag, y)

function noise_variances(R::AbstractVector{<:Real}, y)
    length(R) == length(y) ||
        throw(DimensionMismatch("$(length(R)) noise variances for $(length(y)) data"))
    return check_variances(reshape(float.(collect(R)), size(y)))
end

function noise_variances(R::AbstractMatrix{<:Real}, y)
    y isa AbstractMatrix && size(R) == size(y) && return check_variances(float.(R))
    size(R) == (length(y), length(y)) ||
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
