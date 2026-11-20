# Ill-posedness indicators ----------------------------------------------------------------------
#
# Analysis of the singular value decomposition H = Σᵢ σᵢ vᵢ wᵢᴴ of the (whitened) forward operator and of the Fourier coefficients βᵢ = |vᵢᴴ y| of the data: discrete Picard condition, condition number, numerical rank, decay of the singular values. The analysis is performed for each independent problem: one per sample (frequency) for a block-diagonal space–frequency problem, a single one otherwise.

"""
    IllPosedness

Result of [`ill_posedness`](@ref) for one linear system `y = H u + n` (the whitened one when the noise covariance is given).

**Fields**
- `σ` - `Vector{Float64}`: singular values of `H`, in decreasing order (`min(m, n)` values)
- `β` - `Vector{Float64}`: Fourier coefficients of the data `|vᵢᴴ y|`, `vᵢ` being the left singular vectors
- `βm` - `Vector{Float64}`: coefficients smoothed by a moving geometric mean (window `window`)
- `noise_level` - `Float64`: noise level `η` of the coefficients, RMS value of `|vᵢᴴ n|` (`1` when the noise covariance is given, estimated otherwise)
- `noise_estimated` - `Bool`: whether `η` has been estimated from the data
- `picard_index` - `Int`: number `k` of leading components whose smoothed coefficient `β̄ᵢ` exceeds `τ η`, i.e. that emerge from the noise — the components that a regularization can recover
- `picard_slope` - `Float64`: slope `b` of the least-squares fit `log β̄ᵢ ≈ a + b log σᵢ` over the `k` first components (`NaN` if `k < 3`); the coefficients decay as `σᵢᵇ`
- `picard_satisfied` - `Bool`: discrete Picard condition, `b > 1`: above the noise level, the coefficients decay faster than the singular values
- `condition_number` - `Float64`: `σ₁ / σ_min` (`Inf` for a singular operator)
- `numerical_rank` - `Int`: number of singular values above `max(m, n) eps σ₁`
- `nullity` - `Int`: dimension of the null space, `n - numerical_rank`: components of the unknown that do not affect the data at all
- `decay_rate` - `Float64`: degree `α` of the polynomial decay `σᵢ ≈ σ₁ i^(-α)` (least-squares fit over the numerical rank)
- `severity` - `Symbol`: `:mild` (`α ≤ 1`), `:moderate` (`α > 1`) or `:severe` (the singular values decay faster than polynomially: an exponential `σᵢ ≈ σ₁ e^(-μ i)` fits better)
"""
struct IllPosedness
    σ::Vector{Float64}
    β::Vector{Float64}
    βm::Vector{Float64}
    noise_level::Float64
    noise_estimated::Bool
    picard_index::Int
    picard_slope::Float64
    picard_satisfied::Bool
    condition_number::Float64
    numerical_rank::Int
    nullity::Int
    decay_rate::Float64
    severity::Symbol
end

"""
    ill_posedness(prob; R = nothing, window = 5, τ = 2)
    ill_posedness(H, y; R = nothing, window = 5, τ = 2)

Indicators of the ill-posedness of the linear inverse problem `y = H u + n`, computed from the singular value decomposition `H = Σᵢ σᵢ vᵢ wᵢᴴ` and the Fourier coefficients of the data `βᵢ = |vᵢᴴ y|`, whose least-squares solution reads `u_LS = Σᵢ (vᵢᴴ y / σᵢ) wᵢ`.

**Discrete Picard condition.** Writing `vᵢᴴ y = vᵢᴴ H u + vᵢᴴ n`, the coefficients of the exact data `|vᵢᴴ H u| = σᵢ |wᵢᴴ u|` decay with the singular values, whereas those of the noise stay at the level `η = √E|vᵢᴴ n|²`. A meaningful solution can only be recovered if the exact coefficients decay *faster* than the singular values (discrete Picard condition); the components below the noise level must be filtered out by the regularization. The function returns:

* the noise level `η`: `1` when the noise covariance `R` is given (the problem is whitened, so that `E|vᵢᴴ n|² = 1` for real and circular complex noise); otherwise it is estimated from the part of the data outside the range of `H` when `m > rank(H)`, `η² = (‖y‖² - Σᵢ βᵢ²) / (m - rank(H))`, or else from the RMS value of the last quarter of the coefficients (assumed dominated by the noise);
* the Picard index `k`: number of leading components whose coefficients, smoothed by a moving geometric mean over `window` components (`β̄ᵢ`), exceed `τ η`;
* the Picard slope `b`, from the fit `log β̄ᵢ ≈ a + b log σᵢ` over these `k` components: the coefficients decay as `σᵢᵇ`, and the condition holds for `b > 1`.

**Other indicators.** Condition number `σ₁ / σ_min`, numerical rank, dimension of the null space (unknown components invisible in the data, e.g. underdetermined problems), polynomial decay rate `α` of the singular values (`σᵢ ≈ σ₁ i^(-α)`) and the resulting classification of the problem: mildly (`α ≤ 1`), moderately (`α > 1`) or severely ill-posed (faster, exponential decay).

The problem is analyzed as the solver sees it: whitened by `R` (same values as the keyword `R` of [`solve`](@ref)), one analysis per sample for a block-diagonal space–frequency problem, a single one for a space–time or assembled problem (an assembled sparse operator is converted to a dense matrix: the SVD costs `O(m n min(m, n))`). The indicators do not depend on the regularization.

**Inputs**
- `prob` - `InverseProblem`: problem to analyze
- `H`, `y` - `AbstractMatrix`, `AbstractVector`: forward operator and data of a single system (same as `ill_posedness(InverseProblem(H, y))`)
- `R` - noise covariance, as in [`solve`](@ref) (`nothing`: unknown, the noise level is estimated)
- `window` - `Integer`: width of the moving geometric mean of the coefficients (odd, `1` for no smoothing)
- `τ` - `Real`: threshold, in units of the noise level, above which a component emerges from the noise

**Outputs**
- `ip` - [`IllPosedness`](@ref) (single system) or `Vector{IllPosedness}` (one per sample of a block-diagonal problem)

**Example**
```julia
ip = ill_posedness(InverseProblem(H, y); R = σ² * I)
ip.picard_satisfied, ip.picard_index, ip.condition_number
# Picard plot: ip.σ, ip.β and ip.βm ./ ip.σ against the index, on a log scale
```
"""
ill_posedness(H::AbstractMatrix, y::AbstractVector; kw...) = ill_posedness(InverseProblem(H, y); kw...)

function ill_posedness(prob::InverseProblem; R = nothing, window::Integer = 5, τ::Real = 2)
    (window >= 1 && isodd(window)) || throw(ArgumentError("window must be a positive odd integer"))
    τ > 0 || throw(ArgumentError("τ must be positive"))
    σ2 = noise_variances(R, prob.y)
    return analyze_operator(prob.H, prob.y, σ2, window, float(τ))
end

function analyze_operator(op::Union{DenseOperator,SparseOperator}, y, σ2, window, τ)
    H, b = Matrix(op.H), collect(vec(y))
    whiten!(H, b, σ2)
    return analyze_system(H, b, σ2 !== nothing, window, τ)
end

function analyze_operator(op::BlockDiagonalOperator, y, σ2, window, τ)
    return map(eachindex(op.blocks)) do j
        H, b = Matrix(op.blocks[j]), collect(y[:, j])
        whiten!(H, b, σ2 === nothing ? nothing : σ2[:, j])
        analyze_system(H, b, σ2 !== nothing, window, τ)
    end
end

# Indicators of one (whitened) system H u = y
function analyze_system(H::AbstractMatrix, y::AbstractVector, whitened::Bool, window, τ)
    m, n = size(H)
    F = svd(H)
    σ = Vector{Float64}(F.S)
    β = Vector{Float64}(abs.(F.U' * y))
    p = length(σ)

    tol = max(m, n) * eps(Float64) * (p > 0 ? σ[1] : 0.0)
    r = count(>(tol), σ)
    κ = (p > 0 && σ[end] > 0) ? σ[1] / σ[end] : Inf

    η = whitened ? 1.0 : estimate_noise_level(y, β, m, r)
    βm = moving_geometric_mean(β, window)
    k = something(findlast(>(τ * η), view(βm, 1:r)), 0)
    b = k >= 3 ? fit_slope(log.(σ[1:k]), log.(max.(βm[1:k], floatmin()))) : NaN

    α, severity = singular_value_decay(σ, r)
    return IllPosedness(σ, β, βm, η, !whitened, k, b, b > 1, κ, r, n - r, α, severity)
end

# RMS value of |vᵢᴴ n|: from the residual outside the range of H when available, otherwise from the tail of the coefficients
function estimate_noise_level(y, β, m, r)
    if m > r
        δ0 = sum(abs2, y) - sum(abs2, view(β, 1:r))
        δ0 > 0 && return sqrt(δ0 / (m - r))
    end
    tail = view(β, (3 * length(β) ÷ 4 + 1):length(β))
    return sqrt(sum(abs2, tail) / length(tail))
end

# Moving geometric mean over `w` (odd) components, truncated at the ends
function moving_geometric_mean(β, w)
    h, p = w ÷ 2, length(β)
    logβ = log.(max.(β, floatmin()))
    return [exp(sum(view(logβ, max(1, i - h):min(p, i + h))) / (min(p, i + h) - max(1, i - h) + 1)) for i in 1:p]
end

# Slope of the least-squares line y ≈ a + b x
function fit_slope(x, y)
    xm, ym = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- xm) .* (y .- ym)) / sum(abs2, x .- xm)
end

# Polynomial decay rate of the singular values and classification of the ill-posedness
function singular_value_decay(σ, r)
    r < 3 && return (NaN, :undetermined)
    i = 1:r
    logσ = log.(σ[1:r])
    α = -fit_slope(log.(i), logσ)
    # residuals of the polynomial (log σ linear in log i) and exponential (linear in i) fits
    res(x) = (bx = fit_slope(x, logσ); ax = sum(logσ) / r - bx * sum(x) / r; sum(abs2, logσ .- ax .- bx .* x))
    severe = res(collect(Float64, i)) < res(log.(i))
    return α, severe ? :severe : α <= 1 ? :mild : :moderate
end

function Base.show(io::IO, ::MIME"text/plain", ip::IllPosedness)
    println(io, "IllPosedness")
    println(io, "  components:        ", length(ip.σ), " (numerical rank ", ip.numerical_rank, ", nullity ", ip.nullity, ")")
    println(io, "  condition number:  ", ip.condition_number)
    println(io, "  singular values:   ", ip.severity, " decay (α = ", round(ip.decay_rate; digits = 2), ")")
    println(io, "  noise level:       ", ip.noise_level, ip.noise_estimated ? " (estimated)" : " (known)")
    println(io, "  Picard index:      ", ip.picard_index, " component(s) above the noise")
    print(io,   "  Picard condition:  ", ip.picard_satisfied ? "satisfied" : "not satisfied",
          " (slope ", round(ip.picard_slope; digits = 2), ")")
end
