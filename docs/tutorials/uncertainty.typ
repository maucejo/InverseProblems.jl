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

A reconstruction is only half of the answer: how much can it be trusted? With `solve(...; precision = true)`, `SourceEstimation.jl` also approximates the posterior distribution of the sources by a Gaussian $cal(N)(hat(bold(u)), bold(Q)^(-1))$ centred on the solution. The posterior precision matrix:
$ bold(Q) = tau_n (tilde(bold(H))^H tilde(bold(H)) + bold(L)^H bold(D) bold(L)) $
is built by the solver at convergence, from the forward operator whitened by the noise covariance `R` ($tilde(bold(H)) = bold(R)^(-1 slash 2) bold(H)$) and from $bold(D) = lambda bold(W)$, the weights $bold(W)$ of the last iteratively reweighted solve multiplied by the regularization parameter $lambda$ (one parameter $lambda_g$ per group for a group regularization). It is exact for a quadratic penalty and a local approximation otherwise.

The noise precision $tau_n$ scales the noise covariance `R`. With the Bayesian formulation, $tau_n = 1$ when `R` is given (`R` is then the exact noise covariance) and $tau_n$ is estimated with the other hyperparameters otherwise. The additive and multiplicative formulations always use $tau_n = 1$: without `R`, their uncertainty assumes a noise of unit variance and is only meaningful when `R` is given. The covariance $bold(Q)^(-1)$ is never formed: standard deviations, samples and credible intervals are computed from the factorization of $bold(Q)$ used by the solver.

#badges(
  badge(color: julia.red)[Precision],
  badge[posterior_std],
  badge(color: julia.green)[posterior_samples],
  badge(color: julia.purple)[credible_interval],
  badge(color: julia.red)[accuracy_metric],
)

= Known noise covariance

We reuse the mixed excitation of the #link("groups.html")[group regularization tutorial]. The noise added by `agwn` is Gaussian with a known variance, given by the helper `noise_variance` of the documentation (`σ² = mean(|y₀|²) / SNR`). It is passed to `solve` through the keyword `R`, which accepts `σ² * I`, a vector of variances, a `Diagonal` matrix, or one variance per measurement point for matrix data.

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)
ωn, kn = modefreq(beam, 5000.)
xm, xu = grid(40), grid(31)
ϕm, ϕu = modeshape(beam, kn, xm), modeshape(beam, kn, xu)
f0 = (ωn[5] + ωn[6]) / 4π
H = SV.solve(ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕu); type = :acc, progress = false).u[:, :, 1];

region, others = 17:29, setdiff(1:31, 17:29)
uref = zeros(eltype(H), 31)
uref[6] = 1.0
uref[region] .= 0.3 .* sin.(π .* (1:13) ./ 14)
y0 = H * uref
y = agwn(y0, 30)
σ² = noise_variance(y0, 30)
prob = InverseProblem(H, y)
```

With `precision = true`, the solution also contains the posterior precision matrix $bold(Q) = tau_n (bold(H)^H bold(R)^(-1) bold(H) + bold(L)^H bold(D) bold(L))$ and a factorization of it, used by all the posterior quantities. The covariance $bold(Q)^(-1)$ is never formed unless explicitly requested.

```julia
G = Groups(others => Lq(1.0), region => Lq(2.0))
sol = solve(prob, Bayesian(G); R = σ² * I, precision = true)
```

Since `R` is the exact noise covariance, the noise precision is not estimated ($tau_n = 1$). Use `Bayesian(G; τn = Estimated())` when `R` is only known up to a scale factor.

= Standard deviations and credible intervals

`posterior_std` returns the marginal posterior standard deviation of each entry of the solution, computed from the diagonal blocks of $bold(Q)^(-1)$ without forming it. `credible_interval` returns the bounds of the equal-tailed credible intervals, estimated from posterior samples (real and imaginary parts separately for complex data).

```julia
sd = posterior_std(sol)
lo, up = credible_interval(sol; level = 0.95);
```

```julia
#| fig-alt-text: "Real part of the reconstruction with its 95 % credible band and the reference."
fig = Figure(size = (840, 320))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "Re(u) (N)")
band!(ax, xu, real.(lo), real.(up), color = (JULIA_COLORS.purple, 0.25), label = "95 % credible interval")
lines!(ax, xu, real.(uref), color = FG, linestyle = :dash, label = "Reference")
scatterlines!(ax, xu, real.(sol.u), color = JULIA_COLORS.purple, label = "Posterior mean")
axislegend(ax, position = :rt)
fig
```

```julia
inside = count(real.(lo) .<= real.(uref) .<= real.(up))
@printf "%d / %d reference values inside the 95 %% credible intervals\n" inside length(uref)
```

= Posterior samples

`posterior_samples(sol, n)` draws exact samples of the Gaussian approximation $cal(N)(bold(u)|hat(bold(u)), bold(Q)^(-1))$. Depending on the factorization available, it uses a Cholesky factor of the precision or the *perturbation–optimization* method, which only requires to apply $bold(Q)^(-1)$ (see the #link("../theory/uncertainty.html")[Theory]).

```julia
S = posterior_samples(sol, 300; rng = Xoshiro(42))
size(S)
```

```julia
#| fig-alt-text: "Posterior samples of the reconstructed field, their empirical standard deviation compared with the exact one."
fig = Figure(size = (840, 320))
ax1 = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)", title = "50 posterior samples")
for k in 1:50
    lines!(ax1, xu, abs.(S[:, k]), color = (JULIA_COLORS.blue, 0.12), linewidth = 1)
end
lines!(ax1, xu, abs.(uref), color = FG, linestyle = :dash)
ax2 = Axis(fig[1, 2], xlabel = "Position x (m)", ylabel = "Standard deviation (N)", title = "Exact vs empirical")
scatterlines!(ax2, xu, sd, label = "posterior_std")
lines!(ax2, xu, vec(sqrt.(sum(abs2, S .- sol.u; dims = 2) ./ size(S, 2))), color = JULIA_COLORS.red,
       linestyle = :dash, label = "from 300 samples")
axislegend(ax2, position = :lt)
fig
```

The standard deviation is much smaller on the sparse part of the beam, where the prior concentrates the solution, than on the distributed load, where the quadratic prior leaves more freedom.

= Performance indicators

Since the reference field is known here, the reconstruction can be assessed with the #link("../theory/performance.html")[performance indicators] of the package, in the sense of metrology: *trueness* (`relative_error`, `correlation`), *precision* (`precision_metric`, the spread of the posterior distribution) and *accuracy* (`accuracy_metric`, both combined). Let us compare several regularizations, all solved with the same noise covariance and `precision = true`:

```julia
cases = [
    "Tikhonov (GCV)"          => Additive(Lq(2.0); λ = GCV()),
    "Multiplicative ℓ1"       => Multiplicative(Lq(1.0)),
    "Multiplicative, groups"  => Multiplicative(G),
    "Bayesian, groups"        => Bayesian(G),
]
sols = [name => solve(prob, f; R = σ² * I, precision = true) for (name, f) in cases];
```

`precision_metric` and `accuracy_metric` accept a solution, a posterior precision matrix (`sol.precision`, the default) or a covariance matrix (`covariance = true`). Only the posterior variances are needed; they are computed from a Cholesky factorization of the precision, without forming its inverse:

```julia
s = last(sols[end])
precision_metric(s), precision_metric(s.precision), precision_metric(posterior_covariance(s); covariance = true)
```

```julia
#| echo: false
#| results: "typst"
rms_error(s) = norm(uref - s.u) / sqrt(length(uref))
ind = [[relative_error(uref, s), 100 - correlation(uref, s), accuracy_metric(uref, s)] for (_, s) in sols]
ratio = performance_ratio(permutedims(reduce(hcat, ind)))
typst_table(["Method", "Relative error", "Correlation", "RMS error E", "Precision P", "Accuracy A", "Performance ratio"],
            [(name, @sprintf("%.1f %%", relative_error(uref, s)), @sprintf("%.1f %%", correlation(uref, s)),
              @sprintf("%.3f N", rms_error(s)), @sprintf("%.3f N", precision_metric(s)),
              @sprintf("%.3f N", accuracy_metric(uref, s)), @sprintf("%.2f", r))
             for ((name, s), r) in zip(sols, ratio)])
```

The accuracy combines the RMS error $E$ of the estimate and the precision $P$ of the posterior, $A^2 = E^2 + P^2$. The four posteriors have almost the same spread ($P approx 0.1$ N), set by the noise level and the conditioning of $H$, so that the accuracy is driven by the error: the group priors, which fit the field, are both truer and more accurate. Comparing $E$ and $P$ is also a check of the uncertainty: for the single-prior methods, the actual error is of the order of the predicted spread ($E approx P$), whereas for the group priors it is several times smaller ($E < P$): their posterior is conservative. The opposite situation, $E >> P$, would reveal an overconfident posterior.

The last column summarizes the relative error, the correlation (as $100 - C$, so that smaller is better) and the accuracy with `performance_ratio`: the geometric mean of the ratios of each indicator to its best value among the methods. A value of 1 means that the method is the best for every indicator.

```julia
#| eval: false
I = [[relative_error(uref, s), 100 - correlation(uref, s), accuracy_metric(uref, s)] for (_, s) in sols]
performance_ratio(permutedims(reduce(hcat, I)))   # one row per method, one column per indicator
```

= Space–frequency problems

The same functions apply to space–frequency problems. For independent frequencies, `sol.precision` is a vector of `nu × nu` blocks; with a coupled mixed norm and an operator along the frequencies, it is a sparse matrix of size `nu·nf × nu·nf`, and the joint samples keep the correlations between frequencies.

```julia
freqs = range(50, 2000, length = 40)
FRF = SV.solve(ModalFRFProblem(ωn, 0.01, freqs, ϕm, ϕu); type = :acc, progress = false).u
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
hm = heatmap!(ax, freqs, xu, transpose(sdf), colormap = :linear_kbc_5_95_c73_n256)
Colorbar(fig[1, 2], hm, label = "σ (N)")
fig
```

#calepin.elements.callout(kind: "note")[
For the additive and multiplicative formulations, $tau_n = 1$: `R` is taken as the exact noise covariance. The precision is that of the Gaussian approximation at the fixed point of the IR algorithm; it is exact for quadratic priors (`Lq(2.0)`) with fixed hyperparameters.
]

#next-steps(
  card(icon: [📚], title: [Posterior uncertainty], href: "../theory/uncertainty.html", color: julia.purple)[Gaussian approximation, block solves and perturbation–optimization.],
  card(icon: [🔔], title: [Bayesian regularization], href: "../theory/bayesian.html", color: julia.green)[The hierarchical model behind `Bayesian`.],
  card(icon: [🎯], title: [Performance indicators], href: "../theory/performance.html", color: julia.red)[Trueness, precision and accuracy of a reconstruction.],
)
