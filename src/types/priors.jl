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

`ℓq` (quasi-)norm regularization term `R(u) = ‖L u‖_q^q = Σᵢ |(L u)ᵢ|^q` — equivalently, a
generalized Gaussian prior of shape `q`. `q = 2` promotes distributed (smooth) sources
(Tikhonov), `q ≤ 1` localized (sparse) sources (`q = 1`: LASSO). For a matrix unknown `U`,
`L` acts on the rows (reconstruction points: `L U`) and the norm is taken entry-wise.

**Inputs**
- `q` - `Real`: norm (shape) parameter, `0 < q ≤ 2`
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows of the unknown
  (spatial smoothing operator, identity by default)

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
    Lpq(p, q; L = I, Lc = nothing, samples = nothing)

Mixed `ℓp,q` regularization term for matrix unknowns (rows = reconstruction points, columns =
frequencies or time samples):

    R(U) = ‖L U Lcᵀ‖_{p,q}^q = Σᵢ (Σⱼ |(L U Lcᵀ)ᵢⱼ|^p)^(q/p)

`L` acts on the rows (space, `Ls` in MSSP 2018) and `Lc` on the columns (frequencies / time
samples, `Lf` in MSSP 2018, Eqs. (10)–(15): global operator `kron(Lc, L)` on `vec(U)`), e.g.
`Lc = diff_operator(nf, 1; h = Δf)` to promote the continuity of the spectra. An operator on
the columns couples the samples in the linear systems: the problem is then assembled and
solved as a whole with sparse Cholesky factorizations.

For a **vectorized** problem (`u = vec(U)`, see [`space_sample_operator`](@ref) and
[`sample_rows`](@ref)), `samples` gives the number of columns of `L u`: the vector `L u` of
the group is reshaped (column-major, i.e. points fastest) into a `(length(L u)/samples) ×
samples` matrix before computing the mixed norm. With `L = space_sample_operator(n, ns;
Lc = Ls)`, `samples = size(Ls, 1)`; with `L = I`, `samples = ns`.

**Inputs**
- `p` - `Real`: norm parameter on the columns (spectral/temporal content), `0 < p ≤ 2`
  (`p ≤ 1`: narrow-band or impulsive, `p = 2`: broadband or continuous)
- `q` - `Real`: norm parameter on the rows (spatial distribution), `0 < q ≤ 2` (`q ≤ 1`:
  localized, `q = 2`: distributed)
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows (identity by
  default), or on `vec(U)` for a vectorized problem
- `Lc` - `Nothing` or `AbstractMatrix`: operator acting on the columns (`ns` columns) for a
  space–frequency or space–time problem (`nothing`: none)
- `samples` - `Nothing` or `Integer`: number of columns of `L u` for a vectorized problem
  (`nothing`: matrix unknown); not used with `Lc`

**Outputs**
- `term` - `Lpq`: regularization term
"""
struct Lpq{T<:Real, O, S<:Union{Nothing,Int}, F<:Union{Nothing,AbstractMatrix}} <: NormTerm
    p::T
    q::T
    L::O
    samples::S
    Lc::F
    function Lpq(p::Real, q::Real; L = I, Lc::Union{Nothing,AbstractMatrix} = nothing,
                 samples::Union{Nothing,Integer} = nothing)
        0 < p <= 2 || throw(ArgumentError("p must be in ]0, 2], got $p"))
        0 < q <= 2 || throw(ArgumentError("q must be in ]0, 2], got $q"))
        check_operator(L)
        samples === nothing || samples >= 1 || throw(ArgumentError("samples must be positive"))
        Lc === nothing || samples === nothing ||
            throw(ArgumentError("`Lc` (matrix unknown) and `samples` (vectorized problem) cannot be combined"))
        T = float(promote_type(typeof(p), typeof(q)))
        S = samples === nothing ? nothing : Int(samples)
        return new{T, typeof(L), typeof(S), typeof(Lc)}(p, q, L, S, Lc)
    end
end

"""Mixed norm with an operator on the columns (`Lc`)."""
const SampleSmoothedLpq = Lpq{<:Real, <:Any, <:Union{Nothing,Int}, <:AbstractMatrix}

# Does a regularization term contain an operator on the columns (samples)?
has_sample_operator(::NormTerm) = false
has_sample_operator(::SampleSmoothedLpq) = true

# Number of rows of L u for a group of n rows
output_rows(L::UniformScaling, n) = n
output_rows(L::AbstractMatrix, n) = size(L, 1)

# Consistency of the layout of a vectorized mixed norm
check_layout(t::NormTerm, n) = nothing
function check_layout(t::Lpq{<:Any,<:Any,Int}, n)
    m = output_rows(t.L, n)
    m % t.samples == 0 ||
        throw(DimensionMismatch("$m components of L u cannot be reshaped into $(t.samples) samples"))
    return nothing
end

check_operator(::UniformScaling) = nothing
check_operator(L::UniformScaling{Bool}) =
    L.λ ? nothing : throw(ArgumentError("the smoothing operator must not be zero"))
check_operator(::AbstractMatrix) = nothing
check_operator(L) = throw(ArgumentError("the smoothing operator must be a matrix or `I`, got $(typeof(L))"))

"""
    Groups(rows₁ => term₁, rows₂ => term₂, ...)

Group (region-wise) regularization `R(u) = Σ_g R_g(L_g u_g)`, where `u_g` gathers the rows `rows_g` of the unknown and `L_g` is the operator `L` of `term_g`. The groups must form a partition of the rows.

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
has_sample_operator(G::Groups) = any(has_sample_operator, G.terms)

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
            throw(DimensionMismatch("operator L with $(size(t.L, 2)) columns for a group of $(length(rows)) rows"))
        check_layout(t, length(rows))
    end
    return G
end

# Same terms, other norm parameters (used by the Bayesian solver). The recursion on the
# tuple keeps the types of the terms inferable.
with_shape(t::Lq, q) = Lq(q; L = t.L)
with_shapes(::Tuple{}, q, i = 1) = ()
with_shapes(terms::Tuple, q, i = 1) =
    (with_shape(first(terms), q[i]), with_shapes(Base.tail(terms), q, i + 1)...)
