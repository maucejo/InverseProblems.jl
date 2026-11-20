# Assembled regularization ----------------------------------------------------------------
#
# The group smoothing operators L_g are assembled once into a global operator L acting on the
# rows of the unknown: R(u) = Σ_g R_g((L U)[rows_g, :]). Dispatch on the type of L keeps the
# fast diagonal path when every L_g is the identity (default).

"""Global smoothing operator."""
abstract type SmoothingOperator end

"""All groups use the identity: `L U = U`."""
struct IdentitySmoothing <: SmoothingOperator end

"""
General smoothing operator (sparse, `nout × nu`), with its LU factorization when it is
square and invertible (`F = nothing` otherwise).
"""
struct MatrixSmoothing{M<:AbstractMatrix, F} <: SmoothingOperator
    L::M
    F::F
end

MatrixSmoothing(L::AbstractMatrix) = MatrixSmoothing(L, invertible_factorization(L))

"""Square invertible smoothing operator: the standard-form transformation `H L⁻¹` exists."""
const InvertibleSmoothing = MatrixSmoothing{<:AbstractMatrix, <:Factorization}

"""Rectangular or singular smoothing operator: `Lᴴ D L` is singular."""
const SingularSmoothing = MatrixSmoothing{<:AbstractMatrix, Nothing}

function invertible_factorization(L::AbstractMatrix)
    size(L, 1) == size(L, 2) || return nothing
    F = lu(L; check = false)
    return issuccess(F) ? F : nothing
end

# L⁻¹ X and L⁻ᴴ X, where X acts on vec(U) (U: nu × ns) and L acts on the rows of U
smoothing_ldiv(S::InvertibleSmoothing, X, nu) = reshape(S.F \ reshape(X, nu, :), size(X))
smoothing_ldiv_adjoint(S::InvertibleSmoothing, X, nu) = reshape(S.F' \ reshape(X, nu, :), size(X))

apply(::IdentitySmoothing, U) = U
apply(S::MatrixSmoothing, U) = S.L * U

# In-place application L U → V (no copy for L = I), with its buffer
smoothing_buffer(::IdentitySmoothing, U) = nothing
smoothing_buffer(S::MatrixSmoothing, U) = similar(U, size(S.L, 1), size(U, 2))
apply!(::Nothing, ::IdentitySmoothing, U) = U
apply!(V, S::MatrixSmoothing, U) = mul!(V, S.L, U)

"""
Regularization assembled for the solvers: global smoothing operator, rows of `L U` belonging
to each group, and norm terms.
"""
struct Regularization{O<:SmoothingOperator, T<:Tuple}
    L::O
    rows::Vector{Vector{Int}}
    terms::T
    nout::Int
end

ngroups(R::Regularization) = length(R.terms)
with_terms(R::Regularization, terms) = Regularization(R.L, R.rows, terms, R.nout)

# The identity is detected from the type of the operators (`L = I`, the default of the terms),
# so that the type of the assembled regularization is inferred from the type of the terms.
# Any other operator (`2I`, `1.0I`, matrix) goes through the general path.
is_identity(::UniformScaling{Bool}) = true
is_identity(::Any) = false
all_identity(::Tuple{}) = true
all_identity(terms::Tuple) = is_identity(first(terms).L) && all_identity(Base.tail(terms))

function Regularization(G::Groups, nu::Integer)
    all_identity(G.terms) &&
        return Regularization(IdentitySmoothing(), G.rows, G.terms, nu)

    T = mapreduce(t -> eltype(t.L isa UniformScaling ? t.L.λ : t.L), promote_type, G.terms;
                  init = Float64)
    Is, Js, Vs = Int[], Int[], T[]
    out = Vector{Vector{Int}}()
    offset = 0
    for (rows, t) in zip(G.rows, G.terms)
        n = length(rows)
        Lg = t.L isa UniformScaling ? sparse(t.L, n, n) : sparse(t.L)
        i, j, v = findnz(Lg)
        append!(Is, i .+ offset)
        append!(Js, rows[j])
        append!(Vs, v)
        push!(out, collect(offset .+ (1:size(Lg, 1))))
        offset += size(Lg, 1)
    end
    return Regularization(MatrixSmoothing(sparse(Is, Js, Vs, offset, nu)), out, G.terms, offset)
end
