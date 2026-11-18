#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Space–frequency problems])
#metadata((description: "Broadband source identification: frequency-by-frequency resolution versus mixed-norm space-frequency regularization.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

Most experimental identifications are broadband: the sources are reconstructed over a whole frequency band. Solving each frequency independently ignores a valuable piece of prior knowledge — *a source located at a given point is active over many frequencies*. Mixed norms $ell_(p,q)$ exploit it by coupling the frequencies.

#badges(
  badge[space–frequency],
  badge(color: julia.purple)[Lpq],
  badge(color: julia.green)[multithreading],
)

= Building a space–frequency problem

A space–frequency problem is defined by one FRF matrix per frequency and a response matrix `Y` of size `nm × nf` (one column per frequency). The unknown is then a `nr × nf` matrix: rows are points, columns are frequencies. The FRFs of the steel beam of the previous tutorials are computed by StructuralVibration.jl over the whole band, which returns a 3D array `nm × nr × nf`.

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)
ωn, kn = modefreq(beam, 5000.)
xm, xr = grid(40), grid(31)
ϕm, ϕr = modeshape(beam, kn, xm), modeshape(beam, kn, xr)

freqs = range(50, 2000, length = 60)             # Hz
FRF = SV.solve(ModalFRFProblem(ωn, 0.01, freqs, ϕm, ϕr); type = :acc, progress = false).u
size(FRF)
```

Two sources are simulated: a broadband point force and a narrow-band one.

```julia
Uref = zeros(ComplexF64, length(xr), length(freqs))
Uref[9, :] .= 1.0
Uref[22, :] .= 0.5 .* exp.(-((freqs .- 1000) ./ 300) .^ 2)

Y0 = reduce(hcat, [FRF[:, :, j] * Uref[:, j] for j in eachindex(freqs)])
Y = agwn(Y0, 25)
prob = InverseProblem(FRF, Y)
```

`InverseProblem` accepts the 3D array as is (or a vector of `nm × nr` matrices). Internally, the operator is stored as a `BlockDiagonalOperator`, whose blocks are solved in parallel.

= Frequency by frequency

With `Lq` terms only, the frequencies are independent: `solve` solves one problem per frequency, each with its own regularization parameter(s), hyperparameters and stopping criterion, in parallel over the available threads (start Julia with `--threads auto`).

```julia
tik = solve(prob, Additive(Lq(2.0); λ = GCV()))
lq  = solve(prob, Multiplicative(Lq(0.5)))
size(lq.u), size(lq.λ)
```

The regularization parameters (and, for the Bayesian formulation, the hyperparameters) are stacked as columns: `lq.λ` has one entry per frequency.

= Coupling the frequencies with mixed norms

The mixed norm `Lpq(p, q)` is defined on the matrix unknown $U$ by

$ R(U) = norm(U)_(p,q)^q = sum_i (sum_j |U_(i j)|^p)^(q slash p), $

where the inner exponent $p$ acts along the frequencies (spectral content) and the outer one $q$ along the space (spatial distribution). `Lpq(2, 0.5)` thus describes sources that are *localized in space* ($q = 0.5$) and *broadband* ($p = 2$).

```julia
l21  = solve(prob, Multiplicative(Lpq(2, 1)))
l205 = solve(prob, Multiplicative(Lpq(2, 0.5)))
l205.λ, l205.iterations
```

The mixed norm couples the frequencies through the weights of the IR algorithm: a single regularization parameter is shared by all the frequencies, but the linear systems remain block-diagonal and are still solved in parallel.

```julia
#| fig-alt-text: "Space-frequency maps of the force amplitude: reference, Tikhonov, l0.5 frequency by frequency, and l2,0.5 mixed norm."
fig = Figure(size = (840, 620))
maps = ["Reference" => Uref, "Tikhonov (GCV), per frequency" => tik.u,
        "Multiplicative ℓ0.5, per frequency" => lq.u, "Multiplicative ℓ2,0.5" => l205.u]
cmax = 1.1
hm = nothing
for (k, (name, U)) in enumerate(maps)
    ax = Axis(fig[fld1(k, 2), mod1(k, 2)], title = name, xlabel = "Frequency (Hz)",
              ylabel = "x (m)")
    global hm = heatmap!(ax, freqs, xr, abs.(transpose(U)), colormap = :linear_kbc_5_95_c73_n256,
                  colorrange = (0, cmax), highclip = :white)
end
Colorbar(fig[:, 3], hm, label = "|u| (N)")
fig
```

```julia
#| fig-alt-text: "Reconstructed spectra at the two source locations."
fig = Figure(size = (840, 340))
for (k, i) in enumerate((9, 22))
    ax = Axis(fig[1, k], title = @sprintf("Point %d (x = %.2f m)", i, xr[i]),
              xlabel = "Frequency (Hz)", ylabel = "|u| (N)")
    lines!(ax, freqs, abs.(Uref[i, :]), color = FG, linestyle = :dash, label = "Reference")
    lines!(ax, freqs, abs.(tik.u[i, :]), label = "Tikhonov")
    lines!(ax, freqs, abs.(lq.u[i, :]), label = "ℓ0.5")
    lines!(ax, freqs, abs.(l205.u[i, :]), label = "ℓ2,0.5")
end
Legend(fig[2, :], current_axis(), orientation = :horizontal)
fig
```

```julia
#| echo: false
#| results: "typst"
cases = [("Tikhonov (GCV)", Code("Additive(Lq(2.0); λ = GCV())"), "frequency by frequency", tik),
         ("Multiplicative ℓ0.5", Code("Multiplicative(Lq(0.5))"), "frequency by frequency", lq),
         ("Multiplicative ℓ2,1", Code("Multiplicative(Lpq(2, 1))"), "coupled (weights)", l21),
         ("Multiplicative ℓ2,0.5", Code("Multiplicative(Lpq(2, 0.5))"), "coupled (weights)", l205)]
typst_table(["Method", "Call", "Strategy", "Relative error"],
            [(n, c, st, @sprintf("%.1f %%", relerr(s.u, Uref))) for (n, c, st, s) in cases])
```

= Operators along the frequencies

`Lpq` also accepts an operator `Lc` acting on the columns of the unknown (frequencies or time steps), so that the penalty becomes $norm(L U L_c^T)_(p,q)^q$. For instance, a first-order finite-difference operator along the frequencies promotes continuous spectra:

```julia
#| eval: false
Δf = step(freqs)
Multiplicative(Lpq(2, 0.5; Lc = diff_operator(length(freqs), 1; h = Δf)))
```

An operator along the columns couples the frequencies in the linear systems themselves. The problem is then *assembled* automatically (sparse block-diagonal $H$, global operator $L_c times.o L$) and solved with sparse Cholesky factorizations.

#calepin.elements.callout(kind: "warning")[
With `Lc` alone, the penalty only acts on the *variations* of the spectra: fields that are constant along the frequencies are not regularized at all. Combine it with a spatial operator `L`, or use it within `Groups`, so that the regularization term still constrains the spatial distribution.
]

= Resolution strategies

The strategy is selected by multiple dispatch on the regularization terms:

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Terms*], [*Strategy*], [*Resolution*]),
  [`Lq` only], [per sample], [one independent problem per frequency, in parallel],
  [at least one `Lpq` without `Lc`], [joint], [common weights and parameter, block-diagonal systems solved in parallel],
  [at least one `Lpq` with `Lc`], [assembled], [vectorized problem, sparse Cholesky (CHOLMOD)],
)

#next-steps(
  card(icon: [⏱️], title: [Time-domain reconstruction], href: "time-domain.html", color: julia.blue)[The same mixed norms on space–time problems.],
  card(icon: [⚙️], title: [Algorithms], href: "../theory/algorithms.html", color: julia.purple)[IR weights, strategies and linear solvers.],
)
