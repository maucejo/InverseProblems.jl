# Norm values, damping parameter and IR weights ---------------------------------------

# `penalty(term, V)` is evaluated on V = L_g U_g.

"""
    penalty(term, V)

Value of a regularization term.

**Inputs**
- `term` - `Lq` or `Lpq`: regularization term
- `V` - `AbstractArray`: argument `L U` of the term

**Outputs**
- `p` - `Real`: `‖V‖_q^q` or `‖V‖_{p,q}^q`
"""
penalty(t::Lq, V) = sum(x -> abs(x)^t.q, V)

# Matrix layout (space × samples) of L u for a mixed norm: as is for a matrix unknown,
# reshaped column-major for a vectorized problem (`samples` columns)
mixed_layout(::Lpq{<:Any,<:Any,Nothing}, V) =
    size(V, 2) > 1 ? V :
    throw(ArgumentError("mixed norm on a single column: for a vectorized problem, give the number of samples (`Lpq(p, q; samples = ns)`)"))
function mixed_layout(t::Lpq{<:Any,<:Any,Int}, V)
    size(V, 2) == t.samples && return V
    size(V, 2) == 1 ||
        throw(DimensionMismatch("L u has $(size(V, 2)) columns, expected 1 (vectorized problem) or $(t.samples)"))
    return reshape(V, :, t.samples)
end

function penalty(t::Lpq, V::AbstractMatrix)
    p, q = t.p, t.q
    return sum(eachrow(mixed_layout(t, V))) do row
        sum(x -> abs(x)^p, row)^(q / p)
    end
end

"""
    penalty(R, V)

Total regularization `Σ_g R_g(V[rows_g, :])`.

**Inputs**
- `R` - `Regularization`: assembled regularization
- `V` - `AbstractMatrix`: argument `L U` of the regularization

**Outputs**
- `p` - `Real`: sum of the group penalties
"""
penalty(R::Regularization, V) =
    sum(penalty(R.terms[g], view(V, R.rows[g], :)) for g in eachindex(R.terms))

"""
    damping(V, level = 0.05)

Damping parameter `ε` avoiding infinite IR weights: a fraction `level` of the values of `|V|` are less than or equal to `ε` (cumulative histogram rule used throughout the papers).

**Inputs**
- `V` - `AbstractArray`: argument `L U` of the regularization term (initial solution)
- `level` - `Real`: fraction of the values below `ε`

**Outputs**
- `ε` - `Real`: damping parameter (`eps()` if the quantile is zero)
"""
function damping(V, level = 0.05)
    a = vec(abs.(V))
    ε = quantile(a, level)
    return ε > 0 ? ε : eps(eltype(a))
end

dampings(R::Regularization, V) = [damping(view(V, rows, :)) for rows in R.rows]

# Weights -------------------------------------------------------------------------------
#
# `reweight!(W, term, V, ε)` fills W (same size as V = L_g U_g) with the weights given by the first-order optimality condition, ∂R/∂v̄ = W ⊙ v, and returns R(V).
function reweight!(W, t::Lq, V, ε)
    q = t.q
    @. W = (q / 2) * max(ε, abs(V))^(q - 2)
    return penalty(t, V)
end

function reweight!(W, t::Lpq, V, ε)
    p, q = t.p, t.q
    Vm, Wm = mixed_layout(t, V), mixed_layout(t, W)
    for i in axes(Vm, 1)
        row = view(Vm, i, :)
        ws = (q / 2) * max(ε^p, sum(x -> abs(x)^p, row))^(q / p - 1)
        @. Wm[i, :] = ws * max(ε, abs(row))^(p - 2)
    end
    return penalty(t, V)
end

function reweight!(W, P, R::Regularization, V, ε)
    for g in eachindex(R.terms)
        rows = R.rows[g]
        P[g] = reweight!(view(W, rows, :), R.terms[g], view(V, rows, :), ε[g])
    end
    return W, P
end
