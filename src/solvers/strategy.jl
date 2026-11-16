# Resolution strategy -------------------------------------------------------------------------
#
# A space–frequency problem regularized with Lq terms only is solved frequency by frequency
# (the papers treat each frequency independently); with at least one mixed norm Lpq, the
# frequencies are coupled through the weights and the problem is solved as a whole, block by
# block. With an operator along the samples (Lpq with Lf, MSSP 2018), the samples are coupled
# in the linear systems: the space–frequency or space–time problem is assembled (vectorized,
# sparse block-diagonal H for space–frequency problems) and solved as a whole.
# Single-frequency problems are always solved as a whole.

"""Do the regularization terms couple the samples (frequencies)?"""
abstract type SampleCoupling end

"""`Lq` terms only: the samples are independent."""
struct IndependentSamples <: SampleCoupling end

"""At least one mixed norm `Lpq`: the samples are coupled through the weights."""
struct CoupledSamples <: SampleCoupling end

"""At least one operator along the samples (`Lpq` with `Lf`): the samples are coupled in the
linear systems."""
struct SmoothedSamples <: SampleCoupling end

"""
    sample_coupling(prior)

Coupling of the samples induced by a regularization term (dispatch on its type).

**Inputs**
- `prior` - `RegularizationTerm`: `Lq`, `Lpq` or `Groups`

**Outputs**
- `c` - `IndependentSamples()` (`Lq` terms only), `CoupledSamples()` (at least one `Lpq`
  term) or `SmoothedSamples()` (at least one `Lpq` term with `Lf`)
"""
sample_coupling(::Lq) = IndependentSamples()
sample_coupling(::Lpq) = CoupledSamples()
sample_coupling(::SampleSmoothedLpq) = SmoothedSamples()
sample_coupling(G::Groups) = mapreduce(sample_coupling, strongest, G.terms)

# Strongest coupling of two terms: SmoothedSamples > CoupledSamples > IndependentSamples
strongest(::IndependentSamples, b::SampleCoupling) = b
strongest(a::SmoothedSamples, ::SampleCoupling) = a
strongest(a::CoupledSamples, ::IndependentSamples) = a
strongest(::CoupledSamples, b::Union{CoupledSamples,SmoothedSamples}) = b

"""How a problem is solved: as a whole, or sample by sample."""
abstract type Strategy end
struct Joint <: Strategy end
struct PerSample <: Strategy end
struct Assembled <: Strategy end

strategy(::InverseProblem, f::Formulation) = Joint()
strategy(::InverseProblem{<:BlockDiagonalOperator}, f::Formulation) = strategy(sample_coupling(f.prior))
strategy(::InverseProblem{<:DenseOperator}, f::Formulation) = dense_strategy(sample_coupling(f.prior))
strategy(::IndependentSamples) = PerSample()
strategy(::CoupledSamples) = Joint()
strategy(::SmoothedSamples) = Assembled()
dense_strategy(::SampleCoupling) = Joint()
dense_strategy(::SmoothedSamples) = Assembled()

"""
    solve(prob, formulation[, algorithm]; R = nothing, covariance = false, tol = 1e-4, maxiter = 100, init)

Solve the inverse problem `prob` with the given regularization `formulation`. When the
algorithm is omitted, [`default_algorithm`](@ref) is used.

For a space–frequency problem (block-diagonal operator), the resolution strategy depends on
the regularization terms (dispatch on [`sample_coupling`](@ref)):

* `Lq` terms only: frequency-by-frequency resolution — each frequency is an independent
  problem with its own regularization parameter(s), hyperparameters and stopping criterion,
  solved in parallel;
* at least one `Lpq` term: the mixed norm couples the frequencies, the space–frequency
  problem is solved as a whole.

Single-frequency and space–time problems are always solved as a whole.

**Inputs**
- `prob` - `InverseProblem`: inverse problem
- `formulation` - `Formulation`: `Additive`, `Multiplicative` or `Bayesian`
- `algorithm` - `Algorithm`: resolution algorithm (`IR()` for `Additive` and
  `Multiplicative`, `MAP()` or `ElementWise()` for `Bayesian`)
- `R` - `Nothing`, `UniformScaling`, `Diagonal`, `AbstractVector` or `AbstractMatrix`: known diagonal noise covariance — not provided (default, equivalent to `I`), `σ²I`, or the noise variances (vector, or matrix with the layout of `y`). The data fidelity is `‖y - Hu‖²_{R⁻¹}`
- `covariance` - `Bool`: also compute the posterior covariance of the solution. For the additive and multiplicative formulations, `R` is then taken as the exact noise covariance
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Real` or `ParameterChoice`: initial Tikhonov solution — initial regularization parameter, or rule selecting it (`Heuristic()`, `LCurve()`, `GCV()`, `BayesEstimator()`)

**Outputs**
- `sol` - `InverseSolution`: solution, regularization parameter(s), hyperparameters,
  convergence information and, if requested, posterior covariance. For a
  frequency-by-frequency resolution, the quantities of each frequency are stacked as columns
  (`λ`: vector of length `nf`, or `G × nf` / `nr × nf` matrix; `hyper.τn`, `hyper.τs`,
  `hyper.q` likewise), `hyper` also contains the `iterations`, `converged` and `history` of
  each frequency, `iterations`/`converged`/`history` refer to the slowest frequency / all
  frequencies, and `covariance` is a vector of `nr × nr` matrices
"""
solve(prob::InverseProblem, f::Formulation; kwargs...) =
    solve(prob, f, default_algorithm(f); kwargs...)

solve(prob::InverseProblem, f::Formulation, alg::Algorithm; kwargs...) =
    solve_strategy(strategy(prob, f), prob, f, alg; kwargs...)

solve_strategy(::Joint, prob, f, alg; kwargs...) = solve_joint(prob, f, alg; kwargs...)

# Frequency-by-frequency resolution ---------------------------------------------------------

# Noise variances of sample j (nothing: R not provided)
sample_noise(::Nothing, j) = nothing
sample_noise(σ2::AbstractMatrix, j) = σ2[:, j]

# Per-sample quantities stacked as columns
stack_samples(x::AbstractVector{<:Real}) = x
stack_samples(x::AbstractVector{<:AbstractVecOrMat}) = reduce(hcat, x)
stack_samples(x) = x

function solve_strategy(::PerSample, prob, f, alg; R = nothing, kwargs...)
    σ2 = noise_variances(R, prob.y)
    blocks = prob.H.blocks
    sols = Vector{InverseSolution}(undef, length(blocks))
    Threads.@threads for j in eachindex(blocks)
        sols[j] = solve_joint(InverseProblem(blocks[j], prob.y[:, j]), f, alg;
                              R = sample_noise(σ2, j), kwargs...)
    end
    u = reduce(hcat, [s.u for s in sols])
    λ = stack_samples([s.λ for s in sols])
    hyper = (; (k => stack_samples([s.hyper[k] for s in sols]) for k in keys(sols[1].hyper))...,
             iterations = [s.iterations for s in sols],
             converged = [s.converged for s in sols],
             history = [s.history for s in sols])
    slowest = argmax([s.iterations for s in sols])
    Σ = sols[1].covariance === nothing ? nothing : [s.covariance for s in sols]
    return InverseSolution(u, λ, hyper, sols[slowest].iterations, all(s.converged for s in sols),
                           sols[slowest].history, Σ)
end

# Assembled resolution (operator along the samples) -------------------------------------------
#
# The problem is vectorized: u = vec(U), y = vec(Y), H block-diagonal (sparse) for a
# space–frequency problem or the space–time convolution matrix; each group g of rows r_g becomes
# the group sample_rows(r_g, nr, ns) with the operator kron(Lf_g, L_g) and a mixed norm on the
# (size(L_g, 1) × size(Lf_g, 1)) matrix L_g U_g Lf_gᵀ.

assembled_operator(op::BlockDiagonalOperator) = blockdiag(sparse.(op.blocks)...)
assembled_operator(op::DenseOperator) = op.H

# kron(samples, space), kept as I when both are the identity
vectorized_operator(space::UniformScaling, ::Nothing, n, ns) =
    is_identity(space) ? I : space_sample_operator(n, ns; space)
vectorized_operator(space, ::Nothing, n, ns) = space_sample_operator(n, ns; space)
vectorized_operator(space, Lf::AbstractMatrix, n, ns) = space_sample_operator(n, ns; space, samples = Lf)

sample_count(::Nothing, ns) = ns
sample_count(Lf::AbstractMatrix, ns) = size(Lf, 1)

vectorized(t::Lq, n, ns) = Lq(t.q; L = vectorized_operator(t.L, nothing, n, ns))
vectorized(t::Lpq, n, ns) =
    Lpq(t.p, t.q; L = vectorized_operator(t.L, t.Lf, n, ns), samples = sample_count(t.Lf, ns))

function vectorized(prior::RegularizationTerm, nr, ns)
    G = groups(prior, nr)
    return Groups([sample_rows(rows, nr, ns) for rows in G.rows],
                  Tuple(vectorized(t, length(rows), ns) for (rows, t) in zip(G.rows, G.terms)))
end

with_prior(f::Additive, prior) = Additive(prior, f.λ)
with_prior(::Multiplicative{<:Any,C}, prior) where {C} = Multiplicative{typeof(prior),C}(prior)
with_prior(f::Bayesian, prior) = Bayesian(prior, f.τn, f.τs, f.q, f.shape_bounds)

function solve_strategy(::Assembled, prob, f, alg; R = nothing, kwargs...)
    nr, ns = unknown_size(prob)
    ns > 1 || throw(ArgumentError("an operator along the samples (Lf) requires several frequencies or time samples"))
    σ2 = noise_variances(R, prob.y)
    pv = InverseProblem(assembled_operator(prob.H), vec(prob.y))
    sol = solve_joint(pv, with_prior(f, vectorized(f.prior, nr, ns)), alg;
                      R = σ2 === nothing ? nothing : vec(σ2), kwargs...)
    return InverseSolution(reshape(sol.u, nr, ns), sol.λ, sol.hyper, sol.iterations,
                           sol.converged, sol.history, sol.covariance)
end
