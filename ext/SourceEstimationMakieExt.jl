"""
Makie extension of SourceEstimation.jl: plotting functions, loaded automatically with Makie or
one of its backends (CairoMakie, GLMakie, WGLMakie).
"""
module SourceEstimationMakieExt

using SourceEstimation
using SourceEstimation: IllPosedness, InverseProblem
import SourceEstimation: picardplot, picardplot!
using Makie

# Picard plot ---------------------------------------------------------------------------------

function displayed_components(ip::IllPosedness, components)
    p = length(ip.σ)
    n = components === :rank ? ip.numerical_rank :
        components === :all ? p :
        components isa Integer ? Int(components) :
        throw(ArgumentError("components must be :rank, :all or an integer"))
    return 1:clamp(n, 0, p)
end

# Values that cannot be drawn on a logarithmic scale are skipped (NaN)
positive(x) = [v > 0 && isfinite(v) ? v : NaN for v in x]

function picardplot!(ax::Makie.Axis, ip::IllPosedness; components = :rank, smoothed::Bool = true,
                     noise::Bool = true, index::Bool = true,
                     noisecolor = :gray50, indexcolor = :gray50)
    i = displayed_components(ip, components)
    isempty(i) && return ax
    σ, β = ip.σ[i], ip.β[i]
    ratio = (smoothed ? ip.βm[i] : β) ./ σ
    scatterlines!(ax, i, positive(σ), label = rich("σ", subscript("i")))
    scatterlines!(ax, i, positive(β), label = rich("β", subscript("i"), " = |v", subscript("i"), superscript("H"), " y|"))
    scatterlines!(ax, i, positive(ratio),
                  label = rich("β", subscript("i"), " / σ", subscript("i"), smoothed ? " (smoothed)" : ""))
    noise && hlines!(ax, [ip.noise_level], color = noisecolor, linestyle = :dash, label = "noise level η")
    index && ip.picard_index > 0 &&
        vlines!(ax, [ip.picard_index + 0.5], color = indexcolor, linestyle = :dot, label = "Picard index")
    # Logarithmic scale, set once the limits are those of the (positive) data
    if ax.yscale[] !== log10
        reset_limits!(ax)
        ax.yscale = log10
    end
    return ax
end

function picard_title(ip::IllPosedness)
    slope = isnan(ip.picard_slope) ? "–" : string(round(ip.picard_slope; digits = 2))
    status = ip.picard_satisfied ? "satisfied" : "not satisfied"
    return "Picard condition $status: index $(ip.picard_index), slope $slope, " *
           "condition number $(round(ip.condition_number; sigdigits = 2))"
end

# Integer ticks for the component index, about 10 of them
function index_ticks(n)
    step = max(1, ceil(Int, n / 10))
    step = step <= 2 ? step : step <= 5 ? 5 : 10 * ceil(Int, step / 10)
    return 1:step:max(n, 1)
end

function picardplot(ip::IllPosedness; figure = (;), axis = (;), legend::Bool = true, components = :rank, kwargs...)
    n = length(displayed_components(ip, components))
    fig = Figure(; size = (800, 380), figure...)
    ax = Axis(fig[1, 1]; xlabel = "Index i", xticks = index_ticks(n), title = picard_title(ip), axis...)
    picardplot!(ax, ip; components, kwargs...)
    legend && Legend(fig[2, 1], ax, orientation = :horizontal, framevisible = false)
    return fig
end

function picardplot(prob::InverseProblem; R = nothing, window::Integer = 5, τ::Real = 2,
                    sample = nothing, kwargs...)
    ip = ill_posedness(prob; R, window, τ)
    if ip isa AbstractVector
        sample === nothing &&
            throw(ArgumentError("one analysis per sample (block-diagonal problem): choose one with `sample = j`, 1 ≤ j ≤ $(length(ip))"))
        ip = ip[sample]
    end
    return picardplot(ip; kwargs...)
end

picardplot(H::AbstractMatrix, y::AbstractVector; kwargs...) = picardplot(InverseProblem(H, y); kwargs...)

end
