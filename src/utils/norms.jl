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

function penalty(t::Lpq, V::AbstractMatrix)
    p, q = t.p, t.q
    return sum(eachrow(V)) do row
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
    for i in axes(V, 1)
        row = view(V, i, :)
        ws = (q / 2) * max(ε^p, sum(x -> abs(x)^p, row))^(q / p - 1)
        @. W[i, :] = ws * max(ε, abs(row))^(p - 2)
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
