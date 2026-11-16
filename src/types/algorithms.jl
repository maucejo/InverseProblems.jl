# Resolution algorithms -------------------------------------------------------------

"""
    Algorithm

Resolution algorithm, passed as third argument of [`solve`](@ref).
"""
abstract type Algorithm end

"""
    IR()

Iteratively Reweighted algorithm derived from the direct application of the first-order
optimality condition (Aucejo & De Smet, JSV 2020; MSSP 2020). Algorithm of the additive and
multiplicative formulations.
"""
struct IR <: Algorithm end

"""
    MAP()

Maximum a posteriori estimate of the hierarchical Bayesian model — optimal Bayesian
regularization (Aucejo & De Smet, MSSP 2019). Default algorithm for [`Bayesian`](@ref).
"""
struct MAP <: Algorithm end

"""
    default_algorithm(formulation)

Algorithm used by `solve(prob, formulation)`.
"""
default_algorithm(::Additive) = IR()
default_algorithm(::Multiplicative) = IR()
default_algorithm(::Bayesian) = MAP()
