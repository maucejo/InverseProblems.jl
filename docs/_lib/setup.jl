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

"""Variance of the noise added by `agwn(x, snr_dB)`: one value for a vector, one value per
row (measurement point) for a matrix."""
noise_variance(x::AbstractVector, snr_dB) = mean(abs2, x) / 10^(snr_dB / 10)
noise_variance(x::AbstractMatrix, snr_dB) = vec(mean(abs2, x; dims = 2)) ./ 10^(snr_dB / 10)

"""Block lower-triangular (Toeplitz) convolution matrix of size `(nm nt) × (nu nt)` acting on `vec(U)`, `U` being `nu × nt` (rows = points, columns = time steps). The impulse responses are computed mode by mode with StructuralVibration.jl: `h(t) = Σₙ ϕm,ₙ ϕr,ₙᵀ hₙ(t)`."""
function convolution_matrix(ωn, ξ, ϕm, ϕr, t)
    nm, nu, nt = size(ϕm, 1), size(ϕr, 1), length(t)
    hmodes = [impulse_response(Sdof(1.0, ω / 2π, ξ), t) for ω in ωn]
    H = zeros(nm * nt, nu * nt)
    for k in 1:nt
        hk = ϕm * Diagonal([h[k] for h in hmodes]) * transpose(ϕr) .* step(t)
        for l in 1:(nt - k + 1)
            H[(k + l - 2) * nm + 1:(k + l - 1) * nm, (l - 1) * nu + 1:l * nu] = hk
        end
    end
    return H
end

const JULIA_COLORS = (blue = colorant"#4063D8", red = colorant"#CB3C33",
                      green = colorant"#389826", purple = colorant"#9558B2")
const FG = colorant"#7a8291"

# Figures in the font of the site (Inter, static instances in _lib/fonts): the text of the
# SVG figures is embedded as glyph outlines, so the font must be given to CairoMakie
const FONT_DIR = joinpath(DOCS_DIR, "_lib", "fonts")
const SITE_FONTS = (regular = joinpath(FONT_DIR, "Inter-Regular.ttf"),
                    bold = joinpath(FONT_DIR, "Inter-SemiBold.ttf"),
                    italic = joinpath(FONT_DIR, "Inter-Regular.ttf"),
                    bold_italic = joinpath(FONT_DIR, "Inter-SemiBold.ttf"))

set_theme!(Theme(
    fonts = SITE_FONTS,
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
function typst_table(header, rows; align = "center")
    esc(s) = replace(s, r"([\\#*_$@<>\[\]`~])" => s"\\\1")
    cell(x) = x isa Code ? "[`" * x.s * "`]" : "[" * esc(string(x)) * "]"
    io = IOBuffer()
    # Typst's HTML export ignores `align` on table cells: the alignment is applied by the CSS class se-table-<align>
    print(io, "#html.elem(\"div\", attrs: (class: \"se-table se-table-", align, "\"))[#table(columns: ", length(header), ", align: ", align, ", ")
    print(io, "table.header(", join(("[*" * h * "*]" for h in header), ", "), "), ")
    for r in rows
        print(io, join(cell.(r), ", "), ", ")
    end
    print(io, ")]")
    print(String(take!(io)))   # a single output item: evaluated as one Typst snippet
end
