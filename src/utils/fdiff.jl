# Finite-difference smoothing operators ------------------------------------------------------
#
# Helpers building the smoothing operators L of the regularization terms, e.g.
#     Lq(2.0; L = diff_operator(n, 1))                 # first derivative of a 1D field
#     Lq(2.0; L = diff_operator(coords[rows, :], 2))   # ∂/∂y on the points of a group
# The operators act on the reconstruction points given as columns, in the order of the
# unknown (or of the rows of a group).

"""
    fd_weights(h, npts, order)

Finite-difference weights on a uniform grid (port of `ufdwt`, G. von Winckel). Each row
contains the weights of one stencil of `npts` consecutive points: row `r` approximates the
derivative at the `r`-th point of the stencil (row 1: forward scheme, last row: backward
scheme, middle row: centered scheme when `npts` is odd).

**Inputs**
- `h` - `Real`: grid spacing
- `npts` - `Integer`: number of points of the scheme
- `order` - `Integer`: order of the derivative, `order < npts`

**Outputs**
- `W` - `Matrix{Float64}`: `npts × npts` weight matrix
"""
function fd_weights(h::Real, npts::Integer, order::Integer)
    0 <= order < npts || throw(ArgumentError("the order of the derivative must be in [0, npts - 1]"))
    p = npts - 1
    W = zeros(npts, npts)
    rhs = zeros(npts)
    rhs[order + 1] = 1
    for r in 1:npts
        offsets = ((1 - r):(p + 1 - r)) .* h                    # stencil seen from its r-th point
        M = [x^a / factorial(a) for a in 0:p, x in offsets]     # Taylor expansion matrix
        W[r, :] = M \ rhs
    end
    return W
end

default_npts(order) = isodd(order + 1) ? order + 1 : order + 2      # smallest odd number of points

"""
    diff_operator(n, order; npts = default, h = 1.0)

Finite-difference operator of a 1D uniform grid of `n` points (port of `Dx_1D`).

* minimal scheme (`npts = order + 1`): `(n - order) × n` matrix, one stencil of `npts`
  consecutive points per row, without boundary rows (forward differences for `order = 1`,
  second differences of the interior points for `order = 2`);
* centered scheme (`npts` odd, `npts > order + 1`): `n × n` matrix, the first and last
  points use off-centered stencils of the same size.

**Inputs**
- `n` - `Integer`: number of points
- `order` - `Integer`: order of the derivative (`0` gives the identity)
- `npts` - `Integer`: number of points of the scheme — by default the smallest odd number
  of points (`3` for `order = 1`: centered `n × n`; `3` for `order = 2`: minimal
  `(n - 2) × n`)
- `h` - `Real`: grid spacing

**Outputs**
- `L` - `SparseMatrixCSC{Float64}`: finite-difference operator
"""
function diff_operator(n::Integer, order::Integer; npts::Integer = default_npts(order), h::Real = 1.0)
    order == 0 && return sparse(1.0I, n, n)
    npts > order || throw(ArgumentError("npts must be greater than the order of the derivative"))
    npts - order == 1 || isodd(npts) ||
        throw(ArgumentError("npts must be odd for a centered scheme (npts > order + 1)"))
    n >= npts || throw(ArgumentError("$n points are not enough for a $npts-point scheme"))
    W = fd_weights(h, npts, order)
    Is, Js, Vs = Int[], Int[], Float64[]
    if npts - order == 1                                       # minimal scheme
        for i in 1:(n - order), k in 1:npts
            push!(Is, i); push!(Js, i + k - 1); push!(Vs, W[1, k])
        end
        return sparse(Is, Js, Vs, n - order, n)
    end
    c = (npts + 1) ÷ 2                                         # centered scheme
    for i in 1:n
        r, first = i < c ? (i, 1) :                             # left boundary
                   i > n - c + 1 ? (npts - (n - i), n - npts + 1) :   # right boundary
                   (c, i - c + 1)                               # interior
        for k in 1:npts
            push!(Is, i); push!(Js, first + k - 1); push!(Vs, W[r, k])
        end
    end
    return sparse(Is, Js, Vs, n, n)
end

"""
    diff_operator(coords, dir; order = 1, npts = default, atol = 1e-6 × extent)

Finite-difference operator along the direction `dir` of a structured point cloud (2D or 3D
mesh, e.g. the reconstruction points of a group). The points are gathered in grid lines
along `dir` (same other coordinates, within `atol`); the 1D operator is applied along each
line, with its own spacing, which must be uniform.

**Inputs**
- `coords` - `AbstractMatrix{<:Real}`: `n × d` coordinates of the points, in the order of the
  unknown (or of the rows of a group)
- `dir` - `Integer`: direction of the derivative (`1`: x, `2`: y, `3`: z)
- `order` - `Integer`: order of the derivative
- `npts` - `Integer`: number of points of the scheme (see [`diff_operator(n, order)`](@ref))
- `atol` - `Real`: tolerance used to identify the grid lines and check the uniform spacing

**Outputs**
- `L` - `SparseMatrixCSC{Float64}`: finite-difference operator — `n × n` for a centered
  scheme (row `i`: derivative at point `i`), `(n - order × nlines) × n` for a minimal scheme
  (`npts = order + 1`)
"""
function diff_operator(coords::AbstractMatrix{<:Real}, dir::Integer; order::Integer = 1,
                       npts::Integer = default_npts(order),
                       atol::Real = 1e-6 * maximum(x -> -(reverse(extrema(x))...), eachcol(coords)))
    n, d = size(coords)
    1 <= dir <= d || throw(ArgumentError("direction $dir out of range 1:$d"))
    others = [k for k in 1:d if k != dir]
    lines = Dict{Vector{Int},Vector{Int}}()
    for i in 1:n
        push!(get!(lines, [round(Int, coords[i, k] / atol) for k in others], Int[]), i)
    end
    centered = npts - order > 1
    Is, Js, Vs = Int[], Int[], Float64[]
    offset = 0
    for key in sort!(collect(keys(lines)))
        idx = sort(lines[key]; by = i -> coords[i, dir])
        x = coords[idx, dir]
        steps = diff(x)
        h = sum(steps) / length(steps)
        all(s -> abs(s - h) <= atol, steps) ||
            throw(ArgumentError("non-uniform spacing along direction $dir on a grid line"))
        Ld = diff_operator(length(idx), order; npts, h)
        for (r, k, v) in zip(findnz(Ld)...)
            push!(Is, centered ? idx[r] : offset + r); push!(Js, idx[k]); push!(Vs, v)
        end
        offset += size(Ld, 1)
    end
    return sparse(Is, Js, Vs, centered ? n : offset, n)
end

# Space–frequency / space–time operators (vectorized problems) --------------------------------
#
# For a vectorized problem, the unknown is u = vec(U), U being the nr × ns matrix of the field
# (rows: reconstruction points, columns: frequencies or time samples), so that
# u[i + (j - 1) nr] = U[i, j]. With this ordering,
#     kron(Ls, Lx) * vec(U) = vec(Lx * U * transpose(Ls)).

to_operator(L::UniformScaling, n) = sparse(L, n, n)
to_operator(L::AbstractMatrix, n) =
    size(L, 2) == n ? sparse(L) : throw(DimensionMismatch("operator with $(size(L, 2)) columns for a dimension $n"))

"""
    space_sample_operator(nr, ns; space = I, samples = I)

Smoothing operator of a vectorized space–frequency or space–time problem, acting on
`u = vec(U)` (`U`: `nr × ns`, columns = frequencies or time samples):

    L = kron(samples, space),    L * vec(U) = vec(space * U * transpose(samples))

* `samples = Ls`, `space = I`: smoothing along the frequencies / time samples at each point
  (e.g. `Ls = diff_operator(ns, 1; h = Δf)`);
* `samples = I`, `space = Lx`: spatial smoothing at each frequency / time sample;
* both: product of the two operators (mixed derivative). For a sum of penalties, stack the
  operators: `[space_sample_operator(nr, ns; samples = Ls); space_sample_operator(nr, ns; space = Lx)]`.

**Inputs**
- `nr` - `Integer`: number of reconstruction points
- `ns` - `Integer`: number of frequencies or time samples
- `space` - `AbstractMatrix` or `UniformScaling`: spatial operator (`nr` columns)
- `samples` - `AbstractMatrix` or `UniformScaling`: operator along the samples (`ns` columns)

**Outputs**
- `L` - `SparseMatrixCSC`: operator of size `(size(samples, 1)·size(space, 1)) × (ns·nr)`
"""
space_sample_operator(nr::Integer, ns::Integer; space = I, samples = I) =
    kron(to_operator(samples, ns), to_operator(space, nr))

"""
    sample_rows(rows, nr, ns)

Indices, in a vectorized problem `u = vec(U)` (`U`: `nr × ns`), of the points `rows` at all
the frequencies or time samples (port of `define_group_mixed`). The points vary fastest, so
that `space_sample_operator(length(rows), ns; …)` is the smoothing operator of the group.

**Inputs**
- `rows` - `AbstractVector{<:Integer}`: reconstruction points of the group
- `nr` - `Integer`: number of reconstruction points
- `ns` - `Integer`: number of frequencies or time samples

**Outputs**
- `idx` - `Vector{Int}`: indices `rows[k] + (j - 1) nr`, `k` fastest, `j = 1:ns`
"""
sample_rows(rows::AbstractVector{<:Integer}, nr::Integer, ns::Integer) =
    vec([r + (j - 1) * nr for r in rows, j in 1:ns])
