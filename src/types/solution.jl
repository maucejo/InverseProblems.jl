"""
    InverseSolution

Result of [`solve`](@ref).

# Fields
* `u`: reconstructed unknown (vector or `nr × ns` matrix, same layout as the problem);
* `λ`: final regularization parameter(s) — the adaptive one for multiplicative formulations,
  `τs ./ (c τn)` for the Bayesian formulation;
* `hyper`: `NamedTuple` of the other quantities (`τn`, `τs`, `q` for the Bayesian formulation);
* `iterations`: number of iterations of the outer loop;
* `converged`: whether the stopping criterion was met;
* `history`: stopping criterion (relative `ℓ1` variation of `u`) at each iteration;
* `covariance`: posterior covariance of `vec(u)` (`solve(...; covariance = true)`), a vector
  of `nr × nr` matrices (one per sample) for block-diagonal operators, `nothing` otherwise.

See [`posterior_std`](@ref) and [`credible_interval`](@ref).
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
