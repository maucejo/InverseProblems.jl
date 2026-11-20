# Visualization ---------------------------------------------------------------------------------
#
# Plotting functions are declared here and implemented in the package extension
# `SourceEstimationMakieExt` (ext/), loaded automatically with Makie or one of its backends
# (CairoMakie, GLMakie, WGLMakie).

"""
    picardplot(ip::IllPosedness; components = :rank, smoothed = true, figure = (;), axis = (;), legend = true)
    picardplot(prob::InverseProblem; R = nothing, window = 5, τ = 2, sample = nothing, kwargs...)
    picardplot(H, y; R = nothing, kwargs...)

Picard plot of an inverse problem: singular values `σᵢ` of the (whitened) operator, Fourier coefficients of the data `βᵢ = |vᵢᴴ y|` and ratio `β̄ᵢ / σᵢ` of the smoothed coefficients to the singular values (field `βm` of the analysis), against the index `i` (integer ticks), on a logarithmic scale, with the noise level `η` (dashed line) and the Picard index (dotted line). See [`ill_posedness`](@ref) for the computation of these quantities.

The discrete Picard condition holds when, above the noise level, the coefficients `βᵢ` decay faster than the singular values, i.e. when the ratio `β̄ᵢ / σᵢ` decreases; the components beyond the Picard index are dominated by the noise. The title summarizes the analysis: Picard index, Picard slope and condition number.

Requires Makie: load a backend first, e.g. `using CairoMakie`.

**Inputs**
- `ip` - `IllPosedness`: result of [`ill_posedness`](@ref)
- `prob` - `InverseProblem`, or `H`, `y`: problem to analyze, with the keywords `R`, `window` and `τ` of [`ill_posedness`](@ref); `sample` selects the frequency of a block-diagonal space–frequency problem (one analysis per sample)
- `components` - `:rank` (components up to the numerical rank, default), `:all`, or an `Integer`: number of components displayed
- `smoothed` - `Bool`: ratio computed with the smoothed coefficients `β̄ᵢ` (default) or the raw ones `βᵢ`
- `figure`, `axis` - `NamedTuple`: attributes of the `Figure` and of the `Axis`
- `legend` - `Bool`: legend below the axis

**Outputs**
- `fig` - `Makie.Figure`

**Example**
```julia
using CairoMakie
ip = ill_posedness(prob; R = σ² * I)
picardplot(ip)
```
"""
function picardplot end

"""
    picardplot!(ax, ip::IllPosedness; components = :rank, smoothed = true, noise = true, index = true,
                noisecolor = :gray50, indexcolor = :gray50)

Draw the Picard plot of `ip` (see [`picardplot`](@ref)) in an existing Makie `Axis`; its `y` scale is set to `log10` if it is not already logarithmic. The curves carry labels, so that `axislegend(ax)` or `Legend(fig[...], ax)` builds the legend.

**Inputs**
- `ax` - `Makie.Axis`: target axis
- `ip` - `IllPosedness`: result of [`ill_posedness`](@ref)
- `components`, `smoothed` - see [`picardplot`](@ref)
- `noise` - `Bool`: draw the noise level `η`
- `index` - `Bool`: draw the Picard index
- `noisecolor`, `indexcolor`: colors of these two lines

**Outputs**
- `ax` - `Makie.Axis`
"""
function picardplot! end

# Hint when the extension is not loaded
function plot_error_hint(io, exc, argtypes, kwargs)
    if exc.f in (picardplot, picardplot!)
        print(io, "\n`", nameof(exc.f), "` is implemented in the Makie extension of SourceEstimation.jl: load a Makie backend first, e.g. `using CairoMakie`.")
    end
end
