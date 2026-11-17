"""
    InverseSolution

Result of [`solve`](@ref). See [`posterior_std`](@ref), [`posterior_covariance`](@ref), [`posterior_samples`](@ref) and [`credible_interval`](@ref) for the posterior uncertainty.

**Fields**
- `u` - `AbstractVector` or `AbstractMatrix`: reconstructed unknown (vector, or `nr × ns` matrix, with the layout of the problem)
- `λ` - `Real` or `Vector{<:Real}`: final regularization parameter(s) — the adaptive one for the multiplicative formulations, `τs_g / τn` for the Bayesian formulation (MAP `½‖y - Hu‖²_R + Σ_g (λ_g/q_g)‖L_g u_g‖_q^q` for real data, `‖y - Hu‖²_R + Σ_g λ_g‖L_g u_g‖_q^q` for complex data)
- `hyper` - `NamedTuple`: other estimated quantities (`τn`, `τs`, `q` for the Bayesian formulation, empty otherwise)
- `iterations` - `Int`: number of iterations of the outer loop
- `converged` - `Bool`: whether the stopping criterion was met
- `history` - `Vector{Float64}`: stopping criterion (relative `ℓ1` variation of `u`) at each iteration
- `precision` - `Nothing`, `Hermitian`, `SparseMatrixCSC` or `Vector{<:Hermitian}`: posterior
  precision `Q = τn (Hᴴ H + Lᴴ D L)` of `vec(u)` when `solve(...; precision = true)` — dense for
  a dense operator, sparse for an assembled space–frequency problem, one `nr × nr` block per
  sample for independent or block-diagonal samples, `nothing` otherwise
- `posterior` - `Nothing`, `Posterior` or `SamplewisePosterior`: posterior model (factorization
  of the precision) used by [`posterior_std`](@ref), [`posterior_covariance`](@ref) and
  [`posterior_samples`](@ref)
"""
struct InverseSolution{U, L, H<:NamedTuple, Q, P}
    u::U
    λ::L
    hyper::H
    iterations::Int
    converged::Bool
    history::Vector{Float64}
    precision::Q
    posterior::P
end

function Base.show(io::IO, ::MIME"text/plain", sol::InverseSolution)
    println(io, "InverseSolution")
    println(io, "  unknown:    ", summary(sol.u))
    println(io, "  λ:          ", sol.λ)
    isempty(sol.hyper) || println(io, "  hyper:      ", sol.hyper)
    sol.precision === nothing || println(io, "  precision:  ", summary(sol.precision))
    print(io, "  iterations: ", sol.iterations, sol.converged ? " (converged)" : " (not converged)")
end
