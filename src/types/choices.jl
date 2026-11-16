# Regularization parameter choice rules ------------------------------------------------

"""
    ParameterChoice

Rule for choosing a regularization parameter: [`LCurve`](@ref), [`GCV`](@ref),
[`BayesEstimator`](@ref), [`Heuristic`](@ref) (initialization only).
"""
abstract type ParameterChoice end

"""
    LCurve()

L-curve principle: the parameter maximizing the curvature of the L-curve.
"""
struct LCurve <: ParameterChoice end

"""
    GCV()

Generalized Cross Validation.
"""
struct GCV <: ParameterChoice end

"""
    BayesEstimator()

Bayesian estimator: the parameter maximizing the marginal posterior of the regularization parameter (minimum of the functional `J_MAP`, Pereira et al.).
"""
struct BayesEstimator <: ParameterChoice end

"""
    Heuristic()

Heuristic initial value of the regularization parameter (Aucejo & De Smet, MSSP 2017): the geometric mean of the extreme admissible values derived from the (generalized) singular values of `(H, L)`. Only used to initialize iterative algorithms (`init` keyword).
"""
struct Heuristic <: ParameterChoice end

"""
    LambdaMax()

Initial regularization parameter `λ = ‖2 Hᴴ y‖∞` (whitened `H` and `y`): the smallest value
for which the ℓ1-regularized solution is zero (Kim et al., IEEE 2007). Only used to
initialize iterative algorithms (`init` keyword); default initialization of `ElementWise()`.
"""
struct LambdaMax <: ParameterChoice end
