#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [SourceEstimation.jl])
#metadata((
  title: "Home",
  layout: "layouts/site-landing.html",
  description: "Regularized resolution of linear inverse problems in Julia: force reconstruction and source identification.",
)) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#hero(
  name: [SourceEstimation.jl],
  text: [Regularized inverse problems, the Julia way],
  tagline: [Reconstruct forces and identify sources from vibration or acoustic measurements. One verb, four nouns: describe the problem, choose a prior and a formulation, and `solve`.],
  actions: (
    (label: [Get started], href: "getting-started.html", primary: true),
    (label: [Tutorials], href: "tutorials/point-force.html"),
    (label: [Theory], href: "theory/inverse-problems.html"),
    (label: [GitHub], href: "https://github.com/maucejo/SourceEstimation.jl"),
  ),
  image: "assets/logo.svg",
  install: [add https://github.com/maucejo/SourceEstimation.jl],
)

#features(
  card(icon: [📡], title: [Time & frequency domains], color: julia.blue)[
    Single or multi-frequency problems with one FRF matrix per frequency, SISO/MIMO convolution problems in the time domain, or any linear operator `y = H u + n`.
  ],
  card(icon: [🎯], title: [Priors for every source], color: julia.red)[
    `Lq` norms from sparse (`q ≤ 1`) to smooth (`q = 2`), space–frequency mixed norms `Lpq`, and region-wise `Groups` with their own smoothing operators.
  ],
  card(icon: [⚖️], title: [Three formulations], color: julia.green)[
    Additive (generalized Tikhonov, multi-parameter), multiplicative (OMR, MPMR — nothing to tune) and hierarchical Bayesian regularization.
  ],
  card(icon: [🧭], title: [Automatic parameter choice], color: julia.purple)[
    L-curve, generalized cross-validation and Bayesian estimator on the (generalized) SVD of the weighted problem, plus heuristic initializations.
  ],
  card(icon: [⚡], title: [Fast & structure-aware], color: julia.blue)[
    Iteratively reweighted solver with normal equations, Woodbury identity, GSVD, sparse Cholesky and multithreaded block-wise resolution.
  ],
  card(icon: [📊], title: [Uncertainty quantification], color: julia.red)[
    Posterior precision, marginal standard deviations, exact posterior samples by perturbation–optimization, and credible intervals.
  ],
)

== The whole API in five lines

#lead[Everything is resolved by multiple dispatch on the types of the problem, the regularization term, the formulation and the algorithm. Switching from Tikhonov to a sparse Bayesian reconstruction is a one-word change.]

```julia
#| eval: false
using SourceEstimation

prob = InverseProblem(H, y)                         # what to solve: y = H u + n
reg  = Multiplicative(Groups(inner => Lq(0.5),      # how to regularize
                             edges => Lq(2.0; L = D)))
sol  = solve(prob, reg)                             # default algorithm (IR)
sol.u, sol.λ, sol.hyper                             # reconstructed field, parameters
```

== What's inside

#split(
  [
    - *Problems* — `InverseProblem` for single-frequency, space–frequency, space–time and assembled (sparse) problems, real or complex, with a known diagonal noise covariance.
    - *Regularization terms* — `Lq`, `Lpq`, `Groups`, finite-difference operators `diff_operator`, `space_sample_operator`, `sample_rows`.
    - *Formulations* — `Additive`, `Multiplicative`, `Bayesian` (with `Estimated()` or `Fixed()` hyperparameters).
  ],
  [
    - *Parameter choice* — `LCurve`, `GCV`, `BayesEstimator`, `Heuristic`, `LambdaMax`.
    - *Algorithms* — `IR` (iteratively reweighted), `MAP` (optimal Bayesian regularization), `ElementWise` (relevance vector regression).
    - *Uncertainty* — `posterior_std`, `posterior_covariance`, `posterior_samples`, `credible_interval`.
  ],
)

== A 30-second example

#lead[A point force acting on a simply supported steel beam is identified from 40 noisy accelerometers (signal-to-noise ratio of 26 dB). The data are simulated with #link("https://maucejo.github.io/StructuralVibration.jl/")[StructuralVibration.jl]. The Tikhonov solution spreads the force all over the beam, a sparse prior pins it down.]

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)  # L, S, I, E, ρ: simply supported steel beam
ωn, kn = modefreq(beam, 5000.)                 # natural frequencies and wavenumbers
xm, xr = grid(40), grid(31)                    # sensors and reconstruction points
ϕm, ϕr = modeshape(beam, kn, xm), modeshape(beam, kn, xr)

f0 = (ωn[5] + ωn[6]) / 4π                      # between the 5th and 6th modes (Hz)
prob_frf = ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕr)
H = SV.solve(prob_frf; type = :acc, progress = false).u[:, :, 1];   # accelerance (40 × 31)

uref = zeros(ComplexF64, 31); uref[9] = 1      # unit point force
y = agwn(H * uref, 26)                         # noisy measurements

prob = InverseProblem(H, y)
tik = solve(prob, Additive(Lq(2.0); λ = GCV()))
sbr = solve(prob, Bayesian(Lq(1.0)));          # q is estimated
```

```julia
#| fig-alt-text: "Reconstructed force distribution along the beam: Tikhonov solution versus sparse Bayesian solution and reference."
fig = Figure(size = (820, 300))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)")
lines!(ax, xr, abs.(uref), color = FG, linestyle = :dash, label = "Reference")
scatterlines!(ax, xr, abs.(tik.u), label = "Tikhonov + GCV")
scatterlines!(ax, xr, abs.(sbr.u), label = "Bayesian, q = $(round(only(sbr.hyper.q), digits = 2))")
axislegend(ax, position = :rt)
fig
```

== Environment

#lead[This documentation was generated with the following Julia installation and package versions. All the code on this site is executed at build time.]

```julia
versioninfo()
```

```julia
#| echo: false
#| results: "typst"
deps = Pkg.dependencies()
info(name) = only(v for v in values(deps) if v.name == name)
se = info("SourceEstimation")
pkgver(name) = (v = info(name).version; v === nothing ? "stdlib" : "v" * string(v))
rows = [(Code("SourceEstimation"), "v" * string(se.version), "the package itself")]
push!(rows, (Code("StructuralVibration"), pkgver("StructuralVibration"), "synthetic data of the tutorials"))
for (name, role) in [("LinearAlgebra", "dense factorizations, SVD"),
                     ("SparseArrays", "sparse operators, CHOLMOD"),
                     ("SpecialFunctions", "log-Gamma of the Bayesian priors"),
                     ("Statistics", "quantiles of the credible intervals"),
                     ("Random", "posterior sampling")]
    push!(rows, (Code(name), pkgver(name), role))
end
push!(rows, (Code("CairoMakie"), pkgver("CairoMakie"), "figures of this documentation"))
typst_table(["Package", "Version", "Used for"], rows)
```

== How to cite

If you use SourceEstimation.jl in your research, please cite the papers describing the methods you use (see the #link("theory/references.html")[references]) and the package itself:

```bibtex
#| eval: false
@software{SourceEstimation.jl,
  author = {Aucejo, Mathieu},
  title  = {SourceEstimation.jl: regularized resolution of linear inverse problems},
  url    = {https://github.com/maucejo/SourceEstimation.jl},
  year   = {2026}
}
```
