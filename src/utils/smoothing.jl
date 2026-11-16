# Assembled regularization ----------------------------------------------------------------
#
# The group smoothing operators L_g are assembled once into a global operator L acting on the
# rows of the unknown: R(u) = Σ_g R_g((L U)[rows_g, :]). Dispatch on the type of L keeps the
# fast diagonal path when every L_g is the identity (default).

"""Global smoothing operator."""
abstract type SmoothingOperator end

"""All groups use the identity: `L U = U`."""
struct IdentitySmoothing <: SmoothingOperator end

"""General smoothing operator (sparse, `nout × nr`)."""
struct MatrixSmoothing{M<:AbstractMatrix} <: SmoothingOperator
    L::M
end

apply(::IdentitySmoothing, U) = U
apply(S::MatrixSmoothing, U) = S.L * U

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

is_identity(L::UniformScaling) = isone(L.λ)
is_identity(::AbstractMatrix) = false

function Regularization(G::Groups, nr::Integer)
    all(t -> is_identity(t.L), G.terms) &&
        return Regularization(IdentitySmoothing(), G.rows, G.terms, nr)

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
    return Regularization(MatrixSmoothing(sparse(Is, Js, Vs, offset, nr)), out, G.terms, offset)
end
