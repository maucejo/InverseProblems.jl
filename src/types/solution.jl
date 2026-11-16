"""
    InverseSolution

Result of [`solve`](@ref). See [`posterior_std`](@ref), [`posterior_samples`](@ref) and
[`credible_interval`](@ref) for the posterior uncertainty.

**Fields**
- `u` - `AbstractVector` or `AbstractMatrix`: reconstructed unknown (vector, or `nr × ns`
  matrix, with the layout of the problem)
- `λ` - `Real` or `Vector{<:Real}`: final regularization parameter(s) — the adaptive one for
  the multiplicative formulations, `φ_g τs_g / (c τn)` for the Bayesian formulation (`c = 1/2`,
  `φ_g = 1/q_g` for real data, `c = φ_g = 1` for complex data)
- `hyper` - `NamedTuple`: other estimated quantities (`τn`, `τs`, `q` for the Bayesian
  formulation, empty otherwise)
- `iterations` - `Int`: number of iterations of the outer loop
- `converged` - `Bool`: whether the stopping criterion was met
- `history` - `Vector{Float64}`: stopping criterion (relative `ℓ1` variation of `u`) at each
  iteration
- `covariance` - `Nothing`, `Hermitian` or `Vector{<:Hermitian}`: posterior covariance of
  `vec(u)` when `solve(...; covariance = true)` — one matrix for a dense operator, one
  `nr × nr` matrix per sample for a block-diagonal operator, `nothing` otherwise
"""
struct InverseSolution{U, L, H<:NamedTuple, C}
    u::U
    λ::L
    hyper::H
    iterations::Int
    converged::Bool
    history::Vector{Float64}
    covariance::C
end

function Base.show(io::IO, ::MIME"text/plain", sol::InverseSolution)
    println(io, "InverseSolution")
    println(io, "  unknown:    ", summary(sol.u))
    println(io, "  λ:          ", sol.λ)
    isempty(sol.hyper) || println(io, "  hyper:      ", sol.hyper)
    sol.covariance === nothing || println(io, "  covariance: ", summary(sol.covariance))
    print(io, "  iterations: ", sol.iterations, sol.converged ? " (converged)" : " (not converged)")
end
