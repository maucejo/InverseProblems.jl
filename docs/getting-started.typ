#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Getting started])
#metadata((description: "Install SourceEstimation.jl and solve your first inverse problem.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

`SourceEstimation.jl` solves regularized linear inverse problems of the form

$
  y = H u + n,
$

where $y$ gathers the measurements, $H$ is the forward operator (transfer functions, convolution matrix or any linear model), $u$ is the unknown source field and $n$ is the measurement noise. It is primarily intended for *force reconstruction* and *source identification* in structural dynamics and acoustics, but it can be used for any linear inverse problem.

= Installation

The package requires Julia 1.10 or later. It is not registered yet: install it from its GitHub repository with the package manager. In the Julia REPL, press `]` to enter the `Pkg` mode, then

```julia
#| eval: false
pkg> add https://github.com/maucejo/SourceEstimation.jl
```

or, equivalently, from a script:

```julia
#| eval: false
using Pkg
Pkg.add(url = "https://github.com/maucejo/SourceEstimation.jl")
```

#calepin.elements.callout(kind: "tip")[
The space–frequency problems are solved frequency by frequency (or block by block) in parallel. Start Julia with several threads to benefit from it, e.g. `julia --threads auto`.
]

= The API in a nutshell

The whole public API boils down to one verb, `solve`, and four families of types:

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([*Role*], [*Types*], [*What it describes*]),
  [Problem], [`InverseProblem`], [the forward operator $H$ and the data $y$],
  [Regularization term], [`Lq`, `Lpq`, `Groups`], [the prior knowledge on the sources (sparse, distributed, smooth…)],
  [Formulation], [`Additive`, `Multiplicative`, `Bayesian`], [how data fidelity and regularization are combined],
  [Algorithm], [`IR`, `MAP`, `ElementWise`], [how the problem is solved (optional: a default is chosen)],
)

The regularization parameter of the additive formulation is either fixed or chosen by a rule: `LCurve()`, `GCV()` or `BayesEstimator()`. The multiplicative and Bayesian formulations select it automatically.

= Your first reconstruction

Let us generate a small synthetic problem: 40 noisy measurements of a field of 30 unknowns, in which only three entries are non-zero.

```julia
rng = Xoshiro(2026)
H = randn(rng, 40, 30)                  # forward operator
u_true = zeros(30)
u_true[[5, 14, 23]] = [1.0, -0.6, 0.8]  # sparse source field
y = H * u_true .+ 0.15 .* randn(rng, 40)

prob = InverseProblem(H, y)
```

The problem is then solved with a regularization term and a formulation. A Tikhonov regularization (additive formulation, $ell_2$ norm) with a regularization parameter chosen by generalized cross-validation reads

```julia
sol_tik = solve(prob, Additive(Lq(2.0); λ = GCV()))
```

Since the source field is sparse, an $ell_1$ prior is better suited. With the multiplicative formulation, there is no regularization parameter to tune beforehand:

```julia
sol_l1 = solve(prob, Multiplicative(Lq(1.0)))
```

Finally, the hierarchical Bayesian formulation estimates the noise level, the prior precision _and_ the shape of the prior from the data:

```julia
sol_bayes = solve(prob, Bayesian(Lq(1.0)))
```

```julia
#| fig-alt-text: "Comparison of the Tikhonov, multiplicative l1 and Bayesian reconstructions with the true sparse field."
fig = Figure(size = (820, 320))
ax = Axis(fig[1, 1], xlabel = "Index", ylabel = "u")
scatterlines!(ax, u_true, color = FG, linestyle = :dash, label = "True field")
scatterlines!(ax, sol_tik.u, label = "Tikhonov (GCV)")
scatterlines!(ax, sol_l1.u, label = "Multiplicative ℓ1")
scatterlines!(ax, sol_bayes.u, label = "Bayesian")
axislegend(ax, position = :rb, nbanks = 2)
fig
```

= Anatomy of a solution

`solve` returns an `InverseSolution` whose fields are:

- `u` — the reconstructed field (vector, or `nr × ns` matrix for space–frequency and space–time problems);
- `λ` — the final regularization parameter(s);
- `hyper` — the other estimated quantities (noise precision `τn`, prior precision `τs` and shape `q` for the Bayesian formulation);
- `iterations`, `converged`, `history` — convergence information of the iterative algorithms;
- `precision`, `posterior` — the posterior precision and the posterior model, when `solve(...; precision = true)`.

```julia
sol_bayes.hyper
```

```julia
#| results: "typst"
#| echo: false
err(s) = @sprintf("%.1f %%", relerr(s.u, u_true))
typst_table(["Method", "Call", "Relative error", "Iterations"],
            [("Tikhonov", Code("Additive(Lq(2.0); λ = GCV())"), err(sol_tik), sol_tik.iterations),
             ("Multiplicative ℓ1", Code("Multiplicative(Lq(1.0))"), err(sol_l1), sol_l1.iterations),
             ("Bayesian", Code("Bayesian(Lq(1.0))"), err(sol_bayes), sol_bayes.iterations)])
```

#calepin.elements.callout(kind: "note")[
With `Additive` and quadratic terms only (`Lq(2.0)`), the iteratively reweighted algorithm is not needed: `solve` performs a single linear solve and `sol.iterations == 0`.
]

= Layout of the unknown

Whatever the domain, the *rows* of the unknown are the reconstruction points and the *columns* are the frequencies or time steps. The layout of the problem is deduced from the types of `H` and `y`:

#table(
  columns: (auto, auto, auto, 1fr),
  align: left,
  table.header([*`H`*], [*`y`*], [*unknown `u`*], [*use case*]),
  [matrix `nm × nr`], [vector], [vector `nr`], [single frequency, generic problem],
  [vector of `nm × nr` matrices, or `nm × nr × nf` array], [matrix `nm × nf`], [matrix `nr × nf`], [space–frequency],
  [matrix `nm·nt × nr·nt`], [matrix `nm × nt`], [matrix `nr × nt`], [space–time (convolution)],
  [sparse matrix], [vector or matrix], [as above], [assembled problem],
)

As a consequence, `Groups`, smoothing operators and mixed norms `Lpq` have the same meaning in the time and frequency domains.

The data can be *real* (e.g. time-domain problems) or *complex* (frequency-domain problems). The field is deduced from the element types of `H` and `y`, and the few formulas that depend on it — the probability distributions of the Bayesian formulation — are selected automatically. The #link("theory/inverse-problems.html")[theory pages] use a notation valid in both cases.

= Next steps

#next-steps(
  card(icon: [🎯], title: [Localizing a point force], href: "tutorials/point-force.html", color: julia.red)[Compare Tikhonov, sparse and Bayesian regularizations on a vibrating beam.],
  card(icon: [🧩], title: [Group regularization], href: "tutorials/groups.html", color: julia.green)[Mix localized and distributed sources with region-wise priors.],
  card(icon: [📈], title: [Space–frequency problems], href: "tutorials/space-frequency.html", color: julia.purple)[Couple the frequencies with mixed norms.],
  card(icon: [📚], title: [Theory], href: "theory/inverse-problems.html", color: julia.blue)[Ill-posedness, regularization, Bayesian inference and algorithms.],
)
