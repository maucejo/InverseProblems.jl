# Regularization terms R(u) — equivalently, the shape of the prior distribution.

"""
    RegularizationTerm

A regularization term `R(u)`. In the Bayesian reading, it is the (negative log)
generalized Gaussian prior on the unknown.
"""
abstract type RegularizationTerm end

"""
Norm-like term `‖L u‖` applying to a set of rows of the unknown (can be put in [`Groups`](@ref)).
"""
abstract type NormTerm <: RegularizationTerm end

"""
    Lq(q; L = I)

`ℓq` (quasi-)norm term `R(u) = ‖L u‖_q^q = Σᵢ |(L u)ᵢ|^q`, with `0 < q ≤ 2`.

* `q = 2`: distributed (smooth) sources — Tikhonov;
* `q ≤ 1`: localized (sparse) sources — `q = 1` is the LASSO.

`L` is the smoothing operator (identity by default). For a matrix unknown `U`, it acts on
the spatial dimension (`L U`) and the norm is taken entry-wise.
"""
struct Lq{T<:Real, O} <: NormTerm
    q::T
    L::O
    function Lq(q::T; L = I) where {T<:Real}
        0 < q <= 2 || throw(ArgumentError("q must be in ]0, 2], got $q"))
        check_operator(L)
        return new{float(T), typeof(L)}(q, L)
    end
end

"""
    Lpq(p, q; L = I)

Mixed `ℓp,q` term for matrix unknowns (rows = space, columns = frequency/time):

    R(U) = ‖L U‖_{p,q}^q = Σᵢ (Σⱼ |(L U)ᵢⱼ|^p)^(q/p)

* `q` encodes the prior on the **spatial** distribution (`q ≤ 1` localized, `q = 2` distributed);
* `p` encodes the prior on the **spectral/temporal** content (`p ≤ 1` narrow-band or impulsive,
  `p = 2` broadband or continuous).
"""
struct Lpq{T<:Real, O} <: NormTerm
    p::T
    q::T
    L::O
    function Lpq(p::Real, q::Real; L = I)
        0 < p <= 2 || throw(ArgumentError("p must be in ]0, 2], got $p"))
        0 < q <= 2 || throw(ArgumentError("q must be in ]0, 2], got $q"))
        check_operator(L)
        T = float(promote_type(typeof(p), typeof(q)))
        return new{T, typeof(L)}(p, q, L)
    end
end

check_operator(::UniformScaling) = nothing
check_operator(::AbstractMatrix) = nothing
check_operator(L) = throw(ArgumentError("the smoothing operator must be a matrix or `I`, got $(typeof(L))"))

"""
    Groups(rows₁ => term₁, rows₂ => term₂, ...)

Group (region-wise) regularization: `R(u) = Σ_g R_g(L_g u_g)`, where `u_g` gathers the
rows `rows_g` of the unknown and `L_g` is the smoothing operator of `term_g`
(identity by default). The groups must form a partition of the rows.

```julia
Groups(inner => Lq(0.5), boundary => Lq(2.0; L = D₂))
```
"""
struct Groups{T<:Tuple{Vararg{NormTerm}}} <: RegularizationTerm
    rows::Vector{Vector{Int}}
    terms::T
end

Groups(pairs::Pair{<:Any,<:NormTerm}...) =
    Groups([collect(Int, first(p)) for p in pairs], map(last, pairs))

ngroups(G::Groups) = length(G.terms)

"""
    groups(term, nr)

Normalize any regularization term into a validated [`Groups`](@ref) over `nr` rows, so that
the solvers only have to deal with one representation.
"""
groups(t::NormTerm, nr::Integer) = groups(Groups([collect(1:nr)], (t,)), nr)

function groups(G::Groups, nr::Integer)
    seen = falses(nr)
    for rows in G.rows, i in rows
        1 <= i <= nr || throw(ArgumentError("row index $i out of range 1:$nr"))
        seen[i] && throw(ArgumentError("row $i belongs to several groups"))
        seen[i] = true
    end
    all(seen) || throw(ArgumentError("rows $(findall(!, seen)) do not belong to any group"))
    for (rows, t) in zip(G.rows, G.terms)
        t.L isa AbstractMatrix && size(t.L, 2) != length(rows) &&
            throw(DimensionMismatch("smoothing operator with $(size(t.L, 2)) columns for a group of $(length(rows)) rows"))
    end
    return G
end

# Same terms, other norm parameters (used by the Bayesian solver)
with_shape(t::Lq, q) = Lq(q; L = t.L)
