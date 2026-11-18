# Hidden setup of every documentation page (first chunk of each page, `echo: false`).
#
# Each page runs in its own IJulia kernel, started in the directory of the page: activate the
# documentation environment (SourceEstimation.jl in development mode, StructuralVibration.jl
# for the synthetic data, CairoMakie), load the packages and a few small helpers, and set a
# plotting theme in the Julia colors that reads well in light and dark modes (transparent
# background, mid-gray text).

const DOCS_DIR = dirname(@__DIR__)   # docs/ (this file is docs/_lib/setup.jl)

import Pkg
Pkg.activate(DOCS_DIR; io = devnull)

using SourceEstimation, LinearAlgebra, Random, Printf, Statistics
using CairoMakie
CairoMakie.activate!(type = "svg")

# Synthetic data are generated with StructuralVibration.jl. Both packages export `solve`:
# StructuralVibration is imported under the name SV and only the needed names are brought in.
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn, Sdof, impulse_response

# --- Helpers of the tutorials ------------------------------------------------------------

"""`n` equally spaced interior points of a beam of length `L`."""
grid(n, L = 1.0) = collect(range(0, L, length = n + 2)[2:end-1])

"""Relative reconstruction error (%)."""
relerr(u, uref) = 100 * norm(u - uref) / norm(uref)

"""Variance of the noise added by `agwn(x, snr_dB)`: one value for a vector, one value per
row (measurement point) for a matrix."""
noise_variance(x::AbstractVector, snr_dB) = mean(abs2, x) / 10^(snr_dB / 10)
noise_variance(x::AbstractMatrix, snr_dB) = vec(mean(abs2, x; dims = 2)) ./ 10^(snr_dB / 10)

"""Block lower-triangular (Toeplitz) convolution matrix of size `(nm nt) × (nr nt)` acting on `vec(U)`, `U` being `nr × nt` (rows = points, columns = time steps). The impulse responses are computed mode by mode with StructuralVibration.jl: `h(t) = Σₙ ϕm,ₙ ϕr,ₙᵀ hₙ(t)`."""
function convolution_matrix(ωn, ξ, ϕm, ϕr, t)
    nm, nr, nt = size(ϕm, 1), size(ϕr, 1), length(t)
    hmodes = [impulse_response(Sdof(1.0, ω / 2π, ξ), t) for ω in ωn]
    H = zeros(nm * nt, nr * nt)
    for k in 1:nt
        hk = ϕm * Diagonal([h[k] for h in hmodes]) * transpose(ϕr) .* step(t)
        for l in 1:(nt - k + 1)
            H[(k + l - 2) * nm + 1:(k + l - 1) * nm, (l - 1) * nr + 1:l * nr] = hk
        end
    end
    return H
end

const JULIA_COLORS = (blue = colorant"#4063D8", red = colorant"#CB3C33",
                      green = colorant"#389826", purple = colorant"#9558B2")
const FG = colorant"#7a8291"

set_theme!(Theme(
    fontsize = 14,
    backgroundcolor = :transparent,
    textcolor = FG,
    palette = (color = [JULIA_COLORS.blue, JULIA_COLORS.red, JULIA_COLORS.green, JULIA_COLORS.purple, colorant"#E09F3E"],),
    Axis = (
        backgroundcolor = :transparent,
        xgridcolor = (FG, 0.18), ygridcolor = (FG, 0.18),
        leftspinecolor = (FG, 0.6), bottomspinecolor = (FG, 0.6),
        topspinevisible = false, rightspinevisible = false,
        xtickcolor = (FG, 0.6), ytickcolor = (FG, 0.6),
        titlefont = :bold,
    ),
    Legend = (framevisible = false, backgroundcolor = :transparent),
    Colorbar = (tickcolor = (FG, 0.6), spinewidth = 0),
    Lines = (linewidth = 2.2,),
    Scatter = (markersize = 9,),
))

"""Inline code cell of a Typst table."""
struct Code
    s::String
end

"""Print a Typst table (for chunks with `results: "typst"`)."""
function typst_table(header, rows; align = "left")
    esc(s) = replace(s, r"([\\#*_$@<>\[\]`~])" => s"\\\1")
    cell(x) = x isa Code ? "[`" * x.s * "`]" : "[" * esc(string(x)) * "]"
    io = IOBuffer()
    print(io, "#table(columns: ", length(header), ", align: ", align, ", ")
    print(io, "table.header(", join(("[*" * h * "*]" for h in header), ", "), "), ")
    for r in rows
        print(io, join(cell.(r), ", "), ", ")
    end
    print(io, ")")
    print(String(take!(io)))   # a single output item: evaluated as one Typst snippet
end
