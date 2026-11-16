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

`ℓq` (quasi-)norm regularization term `R(u) = ‖L u‖_q^q = Σᵢ |(L u)ᵢ|^q` — equivalently, a generalized Gaussian prior of shape `q`. `q = 2` promotes distributed (smooth) sources (Tikhonov), `q ≤ 1` localized (sparse) sources (`q = 1`: LASSO). For a matrix unknown `U`, `L` acts on the spatial dimension (`L U`) and the norm is taken entry-wise.

**Inputs**
- `q` - `Real`: norm (shape) parameter, `0 < q ≤ 2`
- `L` - `AbstractMatrix` or `UniformScaling`: smoothing operator (identity by default)

**Outputs**
- `term` - `Lq`: regularization term
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

Mixed `ℓp,q` regularization term for matrix unknowns (rows = space, columns = frequency or time):

    R(U) = ‖L U‖_{p,q}^q = Σᵢ (Σⱼ |(L U)ᵢⱼ|^p)^(q/p)

**Inputs**
- `p` - `Real`: norm parameter on the spectral/temporal content, `0 < p ≤ 2` (`p ≤ 1`:
  narrow-band or impulsive, `p = 2`: broadband or continuous)
- `q` - `Real`: norm parameter on the spatial distribution, `0 < q ≤ 2` (`q ≤ 1`:
  localized, `q = 2`: distributed)
- `L` - `AbstractMatrix` or `UniformScaling`: smoothing operator acting on the spatial dimension (identity by default)

**Outputs**
- `term` - `Lpq`: regularization term
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

Group (region-wise) regularization `R(u) = Σ_g R_g(L_g u_g)`, where `u_g` gathers the rows `rows_g` of the unknown and `L_g` is the smoothing operator of `term_g`. The groups must form a partition of the rows.

```julia
Groups(inner => Lq(0.5), boundary => Lq(2.0; L = D₂))
```

**Inputs**
- `rows_g => term_g` - `Pair{<:Any, <:NormTerm}`: row indices of the group and its regularization term (`Lq` or `Lpq`, with its own smoothing operator)

**Outputs**
- `G` - `Groups`: group regularization term
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

Normalize any regularization term into a validated [`Groups`](@ref) over `nr` rows, so that the solvers only have to deal with one representation.

**Inputs**
- `term` - `RegularizationTerm`: `Lq`, `Lpq` or `Groups`
- `nr` - `Integer`: number of rows of the unknown

**Outputs**
- `G` - `Groups`: groups forming a partition of `1:nr`, with consistent smoothing operators
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
