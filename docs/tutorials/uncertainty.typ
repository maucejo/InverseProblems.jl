#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Uncertainty quantification])
#metadata((description: "Posterior precision, standard deviations, samples and credible intervals of a reconstruction.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

A reconstruction is only half of the answer: how much can it be trusted? When the noise covariance is known, SourceEstimation.jl returns the Gaussian approximation of the posterior distribution of the sources at the solution, from which standard deviations, samples and credible intervals are derived.

#badges(
  badge(color: julia.red)[precision],
  badge[posterior_std],
  badge(color: julia.green)[posterior_samples],
  badge(color: julia.purple)[credible_interval],
)

= Known noise covariance

We reuse the mixed excitation of the #link("groups.html")[group regularization tutorial]. The noise added by `agwn` is Gaussian with a known variance, given by the helper `noise_variance` of the documentation (`σ² = mean(|y₀|²) / SNR`). It is passed to `solve` through the keyword `R`, which accepts `σ² * I`, a vector of variances, a `Diagonal` matrix, or one variance per measurement point for matrix data.

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)
ωn, kn = modefreq(beam, 5000.)
xm, xr = grid(40), grid(31)
ϕm, ϕr = modeshape(beam, kn, xm), modeshape(beam, kn, xr)
f0 = (ωn[5] + ωn[6]) / 4π
H = SV.solve(ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕr); type = :acc, progress = false).u[:, :, 1];

region, others = 17:29, setdiff(1:31, 17:29)
uref = zeros(ComplexF64, 31)
uref[6] = 1.0
uref[region] .= 0.3 .* sin.(π .* (1:13) ./ 14)
y0 = H * uref
y = agwn(y0, 30)
σ² = noise_variance(y0, 30)
prob = InverseProblem(H, y)
```

With `precision = true`, the solution also contains the posterior precision matrix $Q = tau_n (H^H R^(-1) H + L^H D L)$ and a factorization of it, used by all the posterior quantities. The covariance $Q^(-1)$ is never formed unless explicitly requested.

```julia
G = Groups(others => Lq(1.0), region => Lq(2.0))
sol = solve(prob, Bayesian(G); R = σ² * I, precision = true)
```

Since `R` is the exact noise covariance, the noise precision is not estimated ($tau_n = 1$). Use `Bayesian(G; τn = Estimated())` when `R` is only known up to a scale factor.

= Standard deviations and credible intervals

`posterior_std` returns the marginal posterior standard deviation of each entry of the solution, computed from the diagonal blocks of $Q^(-1)$ without forming it. `credible_interval` returns the bounds of the equal-tailed credible intervals, estimated from posterior samples (real and imaginary parts separately for complex data).

```julia
sd = posterior_std(sol)
lo, up = credible_interval(sol; level = 0.95);
```

```julia
#| fig-alt-text: "Real part of the reconstruction with its 95 % credible band and the reference."
fig = Figure(size = (840, 320))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "Re(u) (N)")
band!(ax, xr, real.(lo), real.(up), color = (JULIA_COLORS.purple, 0.25), label = "95 % credible interval")
lines!(ax, xr, real.(uref), color = FG, linestyle = :dash, label = "Reference")
scatterlines!(ax, xr, real.(sol.u), color = JULIA_COLORS.purple, label = "Posterior mean")
axislegend(ax, position = :rt)
fig
```

```julia
inside = count(real.(lo) .<= real.(uref) .<= real.(up))
@printf "%d / %d reference values inside the 95 %% credible intervals\n" inside length(uref)
```

= Posterior samples

`posterior_samples(sol, n)` draws exact samples of the Gaussian approximation $cal(N)(u, Q^(-1))$. Depending on the factorization available, it uses a Cholesky factor of the precision or the *perturbation–optimization* method, which only requires to apply $Q^(-1)$ (see #link("../theory/uncertainty.html")[the theory]).

```julia
S = posterior_samples(sol, 300; rng = Xoshiro(42))
size(S)
```

```julia
#| fig-alt-text: "Posterior samples of the reconstructed field, their empirical standard deviation compared with the exact one."
fig = Figure(size = (840, 320))
ax1 = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)", title = "50 posterior samples")
for k in 1:50
    lines!(ax1, xr, abs.(S[:, k]), color = (JULIA_COLORS.blue, 0.12), linewidth = 1)
end
lines!(ax1, xr, abs.(uref), color = FG, linestyle = :dash)
ax2 = Axis(fig[1, 2], xlabel = "Position x (m)", ylabel = "Standard deviation (N)", title = "Exact vs empirical")
scatterlines!(ax2, xr, sd, label = "posterior_std")
lines!(ax2, xr, vec(sqrt.(sum(abs2, S .- sol.u; dims = 2) ./ size(S, 2))), color = JULIA_COLORS.red,
       linestyle = :dash, label = "from 300 samples")
axislegend(ax2, position = :rt)
fig
```

The standard deviation is much smaller on the sparse part of the beam, where the prior concentrates the solution, than on the distributed load, where the quadratic prior leaves more freedom.

= Space–frequency problems

The same functions apply to space–frequency problems. For independent frequencies, `sol.precision` is a vector of `nr × nr` blocks; with a coupled mixed norm and an operator along the frequencies, it is a sparse matrix of size `nr·nf × nr·nf`, and the joint samples keep the correlations between frequencies.

```julia
freqs = range(50, 2000, length = 40)
FRF = SV.solve(ModalFRFProblem(ωn, 0.01, freqs, ϕm, ϕr); type = :acc, progress = false).u
Uref = zeros(ComplexF64, 31, length(freqs)); Uref[9, :] .= 1
Y0 = reduce(hcat, [FRF[:, :, j] * Uref[:, j] for j in eachindex(freqs)])
Y = agwn(Y0, 25)
Rf = noise_variance(Y0, 25)       # one noise variance per measurement point

solf = solve(InverseProblem(FRF, Y), Multiplicative(Lpq(2, 0.5)); R = Rf, precision = true)
sdf = posterior_std(solf)
size(sdf), length(solf.precision)
```

```julia
#| fig-alt-text: "Posterior standard deviation over space and frequency, larger near the natural frequencies."
fig = Figure(size = (840, 300))
ax = Axis(fig[1, 1], xlabel = "Frequency (Hz)", ylabel = "x (m)", title = "Posterior standard deviation")
hm = heatmap!(ax, freqs, xr, transpose(sdf), colormap = :linear_kbc_5_95_c73_n256)
Colorbar(fig[1, 2], hm, label = "σ (N)")
fig
```

#calepin.elements.callout(kind: "note")[
For the additive and multiplicative formulations, $tau_n = 1$: `R` is taken as the exact noise covariance. The precision is that of the Gaussian approximation at the fixed point of the IR algorithm; it is exact for quadratic priors (`Lq(2.0)`) with fixed hyperparameters.
]

#next-steps(
  card(icon: [📚], title: [Posterior uncertainty], href: "../theory/uncertainty.html", color: julia.purple)[Gaussian approximation, block solves and perturbation–optimization.],
  card(icon: [🔔], title: [Bayesian regularization], href: "../theory/bayesian.html", color: julia.green)[The hierarchical model behind `Bayesian`.],
)
