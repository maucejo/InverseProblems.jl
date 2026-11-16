# Linear operators ------------------------------------------------------------

"""
    LinearOperator

Internal representation of the forward operator `H`. The layout of the unknown is
always a matrix `U` of size `nr × ns`: rows are reconstruction points (space),
columns are samples (frequencies or time steps). A single-snapshot problem has `ns = 1`.
"""
abstract type LinearOperator end

"""
Dense operator acting on `vec(U)` (single snapshot, or space–time convolution matrix).
"""
struct DenseOperator{M<:AbstractMatrix} <: LinearOperator
    H::M
    nr::Int
    ns::Int
end

"""
Block-diagonal operator: `H[j]` maps `U[:, j]` to `y[:, j]` (space–frequency problems).
The block structure is exploited: every frequency is solved independently.
"""
struct BlockDiagonalOperator{M<:AbstractMatrix} <: LinearOperator
    blocks::Vector{M}
end

unknown_size(op::DenseOperator) = (op.nr, op.ns)
unknown_size(op::BlockDiagonalOperator) = (size(first(op.blocks), 2), length(op.blocks))

# Problem ----------------------------------------------------------------------

"""
    InverseProblem(H, y)

Linear inverse problem `y = H u + n`. The layout of the unknown depends on the types of `H`
and `y`:

| `H`                                   | `y`                 | unknown `u`             | use case        |
|:--------------------------------------|:--------------------|:------------------------|:----------------|
| `AbstractMatrix` (`nm × nr`)          | `AbstractVector`    | vector (`nr`)           | one frequency   |
| `AbstractVector` of `nm × nr` matrices| `AbstractMatrix` (`nm × nf`) | matrix `nr × nf` | space–frequency |
| `AbstractArray{T,3}` (`nm × nr × nf`) | `AbstractMatrix` (`nm × nf`) | matrix `nr × nf` | space–frequency |
| `AbstractMatrix` (`nm·nt × nr·nt`)    | `AbstractMatrix` (`nm × nt`) | matrix `nr × nt` | space–time      |

Rows of the unknown always correspond to reconstruction points, so that [`Groups`](@ref),
smoothing operators and mixed norms [`Lpq`](@ref) have the same meaning whatever the domain.

**Inputs**
- `H` - `AbstractMatrix`, `AbstractVector{<:AbstractMatrix}` or `AbstractArray{T,3}`: forward
  operator (transfer functions or convolution matrix)
- `y` - `AbstractVector` or `AbstractMatrix`: measured data

**Outputs**
- `prob` - `InverseProblem`: the inverse problem, to be passed to [`solve`](@ref)
"""
struct InverseProblem{Op<:LinearOperator, Y<:AbstractVecOrMat}
    H::Op
    y::Y
end

function InverseProblem(H::AbstractMatrix, y::AbstractVector)
    size(H, 1) == length(y) ||
        throw(DimensionMismatch("H has $(size(H, 1)) rows but y has length $(length(y))"))
    return InverseProblem(DenseOperator(H, size(H, 2), 1), y)
end

function InverseProblem(H::AbstractMatrix, y::AbstractMatrix)
    ns = size(y, 2)
    size(H, 1) == length(y) ||
        throw(DimensionMismatch("H has $(size(H, 1)) rows but y has $(length(y)) entries"))
    size(H, 2) % ns == 0 ||
        throw(DimensionMismatch("size(H, 2) = $(size(H, 2)) is not a multiple of the number of samples $ns"))
    return InverseProblem(DenseOperator(H, size(H, 2) ÷ ns, ns), y)
end

function InverseProblem(H::AbstractVector{<:AbstractMatrix}, y::AbstractMatrix)
    length(H) == size(y, 2) ||
        throw(DimensionMismatch("$(length(H)) transfer matrices for $(size(y, 2)) samples"))
    nm, nr = size(first(H))
    all(h -> size(h) == (nm, nr), H) ||
        throw(DimensionMismatch("all transfer matrices must have the same size"))
    nm == size(y, 1) ||
        throw(DimensionMismatch("transfer matrices have $nm rows but y has $(size(y, 1))"))
    return InverseProblem(BlockDiagonalOperator(collect(H)), y)
end

InverseProblem(H::AbstractArray{<:Any,3}, y::AbstractMatrix) =
    InverseProblem([H[:, :, j] for j in axes(H, 3)], y)

unknown_size(prob::InverseProblem) = unknown_size(prob.H)

# Shape of the returned unknown: vector for vector data, matrix otherwise.
output_shape(prob::InverseProblem{<:Any,<:AbstractVector}, U::AbstractMatrix) = vec(U)
output_shape(::InverseProblem{<:Any,<:AbstractMatrix}, U::AbstractMatrix) = U

function Base.show(io::IO, prob::InverseProblem)
    nr, ns = unknown_size(prob)
    print(io, "InverseProblem(", length(prob.y), " data, unknown ", nr, "×", ns,
          ", ", nameof(typeof(prob.H)), ")")
end
