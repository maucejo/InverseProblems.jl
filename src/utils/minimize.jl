# 1D minimization -------------------------------------------------------------------------
#
# Used to refine the grid searches (regularization parameter, shape parameters): a coarse
# grid locates the global minimum, Brent's method refines it in the neighbouring cells.

"""
    brent(f, a, b; rtol = sqrt(eps()), atol = 1e-12, maxiter = 100)

Local minimum of `f` on `[a, b]` by Brent's method (golden-section search with parabolic
interpolation, no derivative).

**Inputs**
- `f` - `Function`: function of one real variable
- `a`, `b` - `Real`: bounds of the search interval
- `rtol`, `atol` - `Real`: relative and absolute tolerances on the minimizer
- `maxiter` - `Integer`: maximum number of iterations

**Outputs**
- `x` - `Real`: minimizer
- `fx` - `Real`: minimum `f(x)`
"""
function brent(f, a, b; rtol = sqrt(eps()), atol = 1e-12, maxiter = 100)
    a, b = minmax(float(a), float(b))
    c = (3 - sqrt(5)) / 2
    x = w = v = a + c * (b - a)
    fx = fw = fv = f(x)
    d = e = zero(x)
    for _ in 1:maxiter
        m = (a + b) / 2
        tol1 = rtol * abs(x) + atol
        tol2 = 2tol1
        abs(x - m) <= tol2 - (b - a) / 2 && break
        golden = true
        if abs(e) > tol1  # try a parabolic step
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2(q - r)
            q > 0 ? (p = -p) : (q = -q)
            if abs(p) < abs(q * e / 2) && p > q * (a - x) && p < q * (b - x)
                e, d = d, p / q
                u = x + d
                (u - a < tol2 || b - u < tol2) && (d = x < m ? tol1 : -tol1)
                golden = false
            end
        end
        if golden  # golden-section step
            e = (x < m ? b : a) - x
            d = c * e
        end
        u = abs(d) >= tol1 ? x + d : x + (d > 0 ? tol1 : -tol1)
        fu = f(u)
        if fu <= fx
            u < x ? (b = x) : (a = x)
            v, fv, w, fw, x, fx = w, fw, x, fx, u, fu
        else
            u < x ? (a = u) : (b = u)
            if fu <= fw || w == x
                v, fv, w, fw = w, fw, u, fu
            elseif fu <= fv || v == x || v == w
                v, fv = u, fu
            end
        end
    end
    return x, fx
end

"""
    grid_minimize(f, grid)

Global minimum of `f` over a grid, refined by [`brent`](@ref) between the neighbours of the
best grid point. The refined point is kept only if it improves on the grid. `NaN` values are
treated as `+Inf`.

**Inputs**
- `f` - `Function`: function of one real variable
- `grid` - `AbstractVector{<:Real}`: sorted search grid

**Outputs**
- `x` - `Real`: minimizer
"""
function grid_minimize(f, grid)
    g(x) = (v = f(x); isnan(v) ? Inf : v)
    vals = g.(grid)
    i = argmin(vals)
    a, b = grid[max(i - 1, firstindex(grid))], grid[min(i + 1, lastindex(grid))]
    a == b && return grid[i]
    x, fx = brent(g, a, b)
    return fx < vals[i] ? x : grid[i]
end
