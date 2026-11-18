// GENERATED FILE — do not edit by hand. Source: the docstrings of SourceEstimation.jl.
// Regenerate with `julia --project=. _lib/make_api.jl` (from docs/).
#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [API reference])
#metadata((description: "Reference of the public API of SourceEstimation.jl, generated from the docstrings.")) <website-metadata>

#title()

This page documents the public API of SourceEstimation.jl. It is generated from the docstrings of the exported symbols (version 0.1.0). Everything boils down to one verb, #link(<api-solve>)[`solve`], applied to a problem, a regularization formulation and, optionally, an algorithm:

```julia
#| eval: false
prob = InverseProblem(H, y)                         # what to solve
reg  = Multiplicative(Groups(r1 => Lq(0.5),         # how to regularize
                             r2 => Lq(2.0; L = D)))
sol  = solve(prob, reg)                             # default algorithm
sol.u                                               # reconstructed field
```

#calepin.elements.callout(kind: "note")[
Internal functions are not part of the public API and may change without notice. Access them with the `SourceEstimation.` prefix at your own risk.
]

= Index

#table(
  columns: (auto, auto, auto, 1fr),
  align: left,
  table.header([*Name*], [*Kind*], [*Section*], [*Summary*]),
  [#link(<api-InverseProblem>)[`InverseProblem`]], [type], [Problem], [Linear inverse problem `y = H u + n`.],
  [#link(<api-Lq>)[`Lq`]], [type], [Regularization terms], [`ℓq` (quasi-)norm regularization term `R(u) = ‖L u‖_q^q = Σᵢ |(L u)ᵢ|^q` — equivalently, a generalized Gaussian prior of shape `q`.],
  [#link(<api-Lpq>)[`Lpq`]], [type], [Regularization terms], [Mixed `ℓp,q` regularization term for matrix unknowns (rows = reconstruction points, columns = frequencies or time samples):],
  [#link(<api-Groups>)[`Groups`]], [type], [Regularization terms], [Group (region-wise) regularization `R(u) = Σ_g R_g(L_g u_g)`, where `u_g` gathers the rows `rows_g` of the unknown and `L_g` is the operator `L` of `term_g`.],
  [#link(<api-fd_weights>)[`fd_weights`]], [function], [Smoothing operators], [Finite-difference weights on a uniform grid (port of `ufdwt`, G. von Winckel).],
  [#link(<api-diff_operator>)[`diff_operator`]], [function], [Smoothing operators], [Finite-difference operator of a 1D uniform grid of `n` points (port of `Dx_1D`).],
  [#link(<api-space_sample_operator>)[`space_sample_operator`]], [function], [Smoothing operators], [Operator of a vectorized space–frequency or space–time problem, acting on `u = vec(U)` (`U`: `nr × ns`, rows = reconstruction points, columns = frequencies or time samples):],
  [#link(<api-sample_rows>)[`sample_rows`]], [function], [Smoothing operators], [Indices, in a vectorized problem `u = vec(U)` (`U`: `nr × ns`), of the points `rows` at all the frequencies or time samples (port of `define_group_mixed`).],
  [#link(<api-Additive>)[`Additive`]], [type], [Formulations], [Additive regularization],
  [#link(<api-Multiplicative>)[`Multiplicative`]], [type], [Formulations], [Multiplicative regularization — no regularization parameter to tune beforehand:],
  [#link(<api-Bayesian>)[`Bayesian`]], [type], [Formulations], [Hierarchical Bayesian formulation: normal likelihood of covariance `R/τn` (`R`: noise covariance given to #link(<api-solve>)[`solve`]) and, for each group, a generalized Gaussian prior on `L_g u_g` with shape `q_g` and precision `τs_g`:],
  [#link(<api-Estimated>)[`Estimated`]], [type], [Formulations], [Marks a hyperparameter of the #link(<api-Bayesian>)[`Bayesian`] formulation as optimized jointly with the solution (MAP of its full conditional distribution).],
  [#link(<api-Fixed>)[`Fixed`]], [type], [Formulations], [Keeps the shape parameters `q_g` of the #link(<api-Bayesian>)[`Bayesian`] formulation at the values given in the prior (`Lq(q_g)`).],
  [#link(<api-LCurve>)[`LCurve`]], [type], [Parameter choice rules], [L-curve principle: the parameter maximizing the curvature of the L-curve.],
  [#link(<api-GCV>)[`GCV`]], [type], [Parameter choice rules], [Generalized Cross Validation.],
  [#link(<api-BayesEstimator>)[`BayesEstimator`]], [type], [Parameter choice rules], [Bayesian estimator: the parameter maximizing the marginal posterior of the regularization parameter (minimum of the functional `J_MAP`, Pereira et al.).],
  [#link(<api-Heuristic>)[`Heuristic`]], [type], [Parameter choice rules], [Heuristic initial value of the regularization parameter (Aucejo & De Smet, MSSP 2017): the geometric mean of the extreme admissible values derived from the (generalized) singular values of `(H, L)`.],
  [#link(<api-LambdaMax>)[`LambdaMax`]], [type], [Parameter choice rules], [Initial regularization parameter `λ = ‖2 Hᴴ y‖∞` (whitened `H` and `y`): the smallest value for which the ℓ1-regularized solution is zero (Kim et al., IEEE 2007).],
  [#link(<api-IR>)[`IR`]], [type], [Algorithms], [Iteratively Reweighted algorithm derived from the direct application of the first-order optimality condition (Aucejo & De Smet, JSV 2020; MSSP 2020).],
  [#link(<api-MAP>)[`MAP`]], [type], [Algorithms], [Maximum a posteriori estimate of the hierarchical Bayesian model — optimal Bayesian regularization (Aucejo & De Smet, MSSP 2019).],
  [#link(<api-ElementWise>)[`ElementWise`]], [type], [Algorithms], [Element-wise Bayesian algorithm (Relevance Vector Regression) for #link(<api-Bayesian>)[`Bayesian`]: each component of `L u` has its own Gaussian prior precision, estimated jointly with `u` (sparse Bayesian learning).],
  [#link(<api-solve>)[`solve`]], [function], [Resolution], [Solve the inverse problem `prob` with the given regularization `formulation`.],
  [#link(<api-InverseSolution>)[`InverseSolution`]], [type], [Resolution], [Result of #link(<api-solve>)[`solve`].],
  [#link(<api-posterior_std>)[`posterior_std`]], [function], [Uncertainty quantification], [Posterior standard deviation of each entry of the solution: square root of the diagonal of the posterior covariance `Q⁻¹`, obtained from its diagonal blocks (block solves, see #link(<api-posterior_covariance>)[`posterior_covariance`]) — `Q⁻¹` is never formed.],
  [#link(<api-posterior_covariance>)[`posterior_covariance`]], [function], [Uncertainty quantification], [Posterior covariance `Q⁻¹` of the solution, computed from the factorization of the precision.],
  [#link(<api-posterior_samples>)[`posterior_samples`]], [function], [Uncertainty quantification], [Samples of the Gaussian approximation of the posterior `N(u, Q⁻¹)` (circular complex normal for complex unknowns), with `Q = τn A`, `A = Hᴴ H + Lᴴ D L`.],
  [#link(<api-credible_interval>)[`credible_interval`]], [function], [Uncertainty quantification], [Bayesian credible interval of each entry of the solution, estimated by sampling the posterior (#link(<api-posterior_samples>)[`posterior_samples`]): empirical quantiles `(1 - level)/2` and `(1 + level)/2` of the samples, entry by entry.],
)

= Problem

Definition of the inverse problem `y = H u + n`.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `InverseProblem` <api-InverseProblem>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
InverseProblem(H, y)
```

Linear inverse problem `y = H u + n`. The layout of the unknown depends on the types of `H` and `y`:

#table(
  columns: 4,
  align: left,
  table.header([*`H`*], [*`y`*], [*unknown `u`*], [*use case*]),
  [`AbstractMatrix` (`nm × nr`)], [`AbstractVector`], [vector (`nr`)], [one frequency],
  [`AbstractVector` of `nm × nr` matrices], [`AbstractMatrix` (`nm × nf`)], [matrix `nr × nf`], [space–frequency],
  [`AbstractArray{T,3}` (`nm × nr × nf`)], [`AbstractMatrix` (`nm × nf`)], [matrix `nr × nf`], [space–frequency],
  [`AbstractMatrix` (`nm·nt × nr·nt`)], [`AbstractMatrix` (`nm × nt`)], [matrix `nr × nt`], [space–time],
  [`AbstractSparseMatrix`], [`AbstractVector` or `AbstractMatrix`], [as above], [assembled problem (e.g. vectorized space–frequency problem with coupled frequencies)],
)

Rows of the unknown always correspond to reconstruction points, so that #link(<api-Groups>)[`Groups`], smoothing operators and mixed norms #link(<api-Lpq>)[`Lpq`] have the same meaning whatever the domain.

*Inputs*

- `H` - `AbstractMatrix`, `AbstractVector{<:AbstractMatrix}` or `AbstractArray{T,3}`: forward operator (transfer functions or convolution matrix)
- `y` - `AbstractVector` or `AbstractMatrix`: measured data

*Outputs*

- `prob` - `InverseProblem`: the inverse problem, to be passed to #link(<api-solve>)[`solve`]

]

= Regularization terms

Priors on the sources: ℓq norms, mixed norms and region-wise groups.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Lq` <api-Lq>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Lq(q; L = I)
```

`ℓq` (quasi-)norm regularization term `R(u) = ‖L u‖_q^q = Σᵢ |(L u)ᵢ|^q` — equivalently, a generalized Gaussian prior of shape `q`. `q = 2` promotes distributed (smooth) sources (Tikhonov), `q ≤ 1` localized (sparse) sources (`q = 1`: LASSO). For a matrix unknown `U`, `L` acts on the rows (reconstruction points: `L U`) and the norm is taken entry-wise.

*Inputs*

- `q` - `Real`: norm (shape) parameter, `0 < q ≤ 2`
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows of the unknown (spatial smoothing operator, identity by default)

*Outputs*

- `term` - `Lq`: regularization term

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Lpq` <api-Lpq>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Lpq(p, q; L = I, Lc = nothing, samples = nothing)
```

Mixed `ℓp,q` regularization term for matrix unknowns (rows = reconstruction points, columns = frequencies or time samples):

```
R(U) = ‖L U Lcᵀ‖_{p,q}^q = Σᵢ (Σⱼ |(L U Lcᵀ)ᵢⱼ|^p)^(q/p)
```

`L` acts on the rows (space, `Ls` in MSSP 2018) and `Lc` on the columns (frequencies / time samples, `Lf` in MSSP 2018, Eqs. (10)–(15): global operator `kron(Lc, L)` on `vec(U)`), e.g. `Lc = diff_operator(nf, 1; h = Δf)` to promote the continuity of the spectra. An operator on the columns couples the samples in the linear systems: the problem is then assembled and solved as a whole with sparse Cholesky factorizations.

For a *vectorized* problem (`u = vec(U)`, see #link(<api-space_sample_operator>)[`space_sample_operator`] and #link(<api-sample_rows>)[`sample_rows`]), `samples` gives the number of columns of `L u`: the vector `L u` of the group is reshaped (column-major, i.e. points fastest) into a `(length(L u)/samples) × samples` matrix before computing the mixed norm. With `L = space_sample_operator(n, ns; Lc = Ls)`, `samples = size(Ls, 1)`; with `L = I`, `samples = ns`.

*Inputs*

- `p` - `Real`: norm parameter on the columns (spectral/temporal content), `0 < p ≤ 2` (`p ≤ 1`: narrow-band or impulsive, `p = 2`: broadband or continuous)
- `q` - `Real`: norm parameter on the rows (spatial distribution), `0 < q ≤ 2` (`q ≤ 1`: localized, `q = 2`: distributed)
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows (identity by default), or on `vec(U)` for a vectorized problem
- `Lc` - `Nothing` or `AbstractMatrix`: operator acting on the columns (`ns` columns) for a space–frequency or space–time problem (`nothing`: none)
- `samples` - `Nothing` or `Integer`: number of columns of `L u` for a vectorized problem (`nothing`: matrix unknown); not used with `Lc`

*Outputs*

- `term` - `Lpq`: regularization term

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Groups` <api-Groups>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Groups(rows₁ => term₁, rows₂ => term₂, ...)
```

Group (region-wise) regularization `R(u) = Σ_g R_g(L_g u_g)`, where `u_g` gathers the rows `rows_g` of the unknown and `L_g` is the operator `L` of `term_g`. The groups must form a partition of the rows.

```julia
#| eval: false
Groups(inner => Lq(0.5), boundary => Lq(2.0; L = D₂))
```

*Inputs*

- `rows_g => term_g` - `Pair{<:Any, <:NormTerm}`: row indices of the group and its regularization term (`Lq` or `Lpq`, with its own smoothing operator)

*Outputs*

- `G` - `Groups`: group regularization term

]

= Smoothing operators

Finite-difference operators and helpers for vectorized space–frequency or space–time problems.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `fd_weights` <api-fd_weights>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
fd_weights(h, npts, order)
```

Finite-difference weights on a uniform grid (port of `ufdwt`, G. von Winckel). Each row contains the weights of one stencil of `npts` consecutive points: row `r` approximates the derivative at the `r`-th point of the stencil (row 1: forward scheme, last row: backward scheme, middle row: centered scheme when `npts` is odd).

*Inputs*

- `h` - `Real`: grid spacing
- `npts` - `Integer`: number of points of the scheme
- `order` - `Integer`: order of the derivative, `order < npts`

*Outputs*

- `W` - `Matrix{Float64}`: `npts × npts` weight matrix

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `diff_operator` <api-diff_operator>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
diff_operator(n, order; npts = default, h = 1.0)
```

Finite-difference operator of a 1D uniform grid of `n` points (port of `Dx_1D`).

- minimal scheme (`npts = order + 1`): `(n - order) × n` matrix, one stencil of `npts` consecutive points per row, without boundary rows (forward differences for `order = 1`, second differences of the interior points for `order = 2`);
- centered scheme (`npts` odd, `npts > order + 1`): `n × n` matrix, the first and last points use off-centered stencils of the same size.

*Inputs*

- `n` - `Integer`: number of points
- `order` - `Integer`: order of the derivative (`0` gives the identity)
- `npts` - `Integer`: number of points of the scheme — by default the smallest odd number of points (`3` for `order = 1`: centered `n × n`; `3` for `order = 2`: minimal `(n - 2) × n`)
- `h` - `Real`: grid spacing

*Outputs*

- `L` - `SparseMatrixCSC{Float64}`: finite-difference operator

#html.elem("hr", attrs: (class: "se-api-sep"))

```julia
#| eval: false
diff_operator(coords, dir; order = 1, npts = default, atol = 1e-6 × extent)
```

Finite-difference operator along the direction `dir` of a structured point cloud (2D or 3D mesh, e.g. the reconstruction points of a group). The points are gathered in grid lines along `dir` (same other coordinates, within `atol`); the 1D operator is applied along each line, with its own spacing, which must be uniform.

*Inputs*

- `coords` - `AbstractMatrix{<:Real}`: `n × d` coordinates of the points, in the order of the unknown (or of the rows of a group)
- `dir` - `Integer`: direction of the derivative (`1`: x, `2`: y, `3`: z)
- `order` - `Integer`: order of the derivative
- `npts` - `Integer`: number of points of the scheme (see `diff_operator(n, order)`)
- `atol` - `Real`: tolerance used to identify the grid lines and check the uniform spacing

*Outputs*

- `L` - `SparseMatrixCSC{Float64}`: finite-difference operator — `n × n` for a centered scheme (row `i`: derivative at point `i`), `(n - order × nlines) × n` for a minimal scheme (`npts = order + 1`)

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `space_sample_operator` <api-space_sample_operator>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
space_sample_operator(nr, ns; L = I, Lc = I)
```

Operator of a vectorized space–frequency or space–time problem, acting on `u = vec(U)` (`U`: `nr × ns`, rows = reconstruction points, columns = frequencies or time samples):

```julia
#| eval: false
L = kron(Lc, L),    L * vec(U) = vec(L * U * transpose(Lc))
```

- `Lc = Lf`, `L = I`: smoothing along the frequencies / time samples at each point (e.g. `Lf = diff_operator(ns, 1; h = Δf)`);
- `Lc = I`, `L = Lx`: spatial smoothing at each frequency / time sample;
- both: product of the two operators (mixed derivative, MSSP 2018 Eq. (15)). For a sum of penalties, stack the operators: `[space_sample_operator(nr, ns; Lc = Lf); space_sample_operator(nr, ns; L = Lx)]`.

*Inputs*

- `nr` - `Integer`: number of reconstruction points (rows of `U`)
- `ns` - `Integer`: number of frequencies or time samples (columns of `U`)
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows (`nr` columns)
- `Lc` - `AbstractMatrix` or `UniformScaling`: operator acting on the columns (`ns` columns)

*Outputs*

- `L` - `SparseMatrixCSC`: operator of size `(size(Lc, 1)·size(L, 1)) × (ns·nr)`

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `sample_rows` <api-sample_rows>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
sample_rows(rows, nr, ns)
```

Indices, in a vectorized problem `u = vec(U)` (`U`: `nr × ns`), of the points `rows` at all the frequencies or time samples (port of `define_group_mixed`). The points vary fastest, so that `space_sample_operator(length(rows), ns; …)` is the smoothing operator of the group.

*Inputs*

- `rows` - `AbstractVector{<:Integer}`: reconstruction points of the group
- `nr` - `Integer`: number of reconstruction points
- `ns` - `Integer`: number of frequencies or time samples

*Outputs*

- `idx` - `Vector{Int}`: indices `rows[k] + (j - 1) nr`, `k` fastest, `j = 1:ns`

]

= Formulations

How data fidelity and regularization are combined.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Additive` <api-Additive>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Additive(prior; λ = LCurve())
```

Additive regularization

```
û = argmin ‖y - Hu‖²_{R⁻¹} + Σ_g λ_g R_g(L_g u_g)
```

`Lq(2.0)` gives the Tikhonov regularization, `Lq(1.0)` the LASSO, #link(<api-Groups>)[`Groups`] the group generalized Tikhonov regularization.

*Inputs*

- `prior` - `RegularizationTerm`: regularization term (`Lq`, `Lpq` or `Groups`)
- `λ` - `Real`, `AbstractVector{<:Real}` or `ParameterChoice`: fixed regularization   parameter, one parameter per group (multi-parameter regularization), or rule selecting it at each iteration (`LCurve()`, `GCV()`, `BayesEstimator()`)

*Outputs*

- `f` - `Additive`: formulation, to be passed to #link(<api-solve>)[`solve`]

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Multiplicative` <api-Multiplicative>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Multiplicative(prior; multiparameter = false)
```

Multiplicative regularization — no regularization parameter to tune beforehand:

- ordinary multiplicative regularization (OMR): `û = argmin ‖y - Hu‖² · Σ_g R_g(L_g u_g)`
- multi-parameter multiplicative regularization (MPMR): `û = argmin ‖y - Hu‖² · Π_g R_g(L_g u_g)`, which yields one adaptive regularization parameter per group.

*Inputs*

- `prior` - `RegularizationTerm`: regularization term (`Lq`, `Lpq` or `Groups`)
- `multiparameter` - `Bool`: `false` for the OMR (default), `true` for the MPMR

*Outputs*

- `f` - `Multiplicative`: formulation, to be passed to #link(<api-solve>)[`solve`]

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Bayesian` <api-Bayesian>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Bayesian(prior; τn = default, τs = Estimated(), q = Estimated(), shape_bounds = (0.05, 2.0))
```

Hierarchical Bayesian formulation: normal likelihood of covariance `R/τn` (`R`: noise covariance given to #link(<api-solve>)[`solve`]) and, for each group, a generalized Gaussian prior on `L_g u_g` with shape `q_g` and precision `τs_g`:

```
p(u, τn, τs, q | y) ∝ p(y | u, τn) Π_g p(L_g u_g | τs_g, q_g) p(τs_g) p(q_g) p(τn)
```

Each hyperparameter is either optimized (`Estimated()`) or fixed. Real and complex data lead to different distributions (real or circular complex normal and generalized Gaussian laws). Weakly informative Gamma hyperpriors are used for the estimated hyperparameters.

*Inputs*

- `prior` - `RegularizationTerm`: `Lq` terms (possibly in `Groups`) giving the initial shapes `q_g` and the smoothing operators
- `τn` - `Real` or `Estimated`: noise precision. `τn = 1` means that `R` is the exact noise covariance. By default, `τn = 1` when `R` is given to `solve`, estimated otherwise
- `τs` - `Estimated`, `Real` or `AbstractVector{<:Real}`: signal precisions (one value for  all groups or one per group)
- `q` - `Estimated` or `Fixed`: shape parameters estimated, or kept at the values of `prior`
- `shape_bounds` - `Tuple{Real,Real}`: search interval of the estimated shape parameters

*Outputs*

- `f` - `Bayesian`: formulation, to be passed to #link(<api-solve>)[`solve`]

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Estimated` <api-Estimated>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Estimated()
```

Marks a hyperparameter of the #link(<api-Bayesian>)[`Bayesian`] formulation as optimized jointly with the solution (MAP of its full conditional distribution).

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Fixed` <api-Fixed>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Fixed()
```

Keeps the shape parameters `q_g` of the #link(<api-Bayesian>)[`Bayesian`] formulation at the values given in the prior (`Lq(q_g)`).

]

= Parameter choice rules

Selection or initialization of the regularization parameter.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `LCurve` <api-LCurve>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
LCurve()
```

L-curve principle: the parameter maximizing the curvature of the L-curve.

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `GCV` <api-GCV>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
GCV()
```

Generalized Cross Validation.

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `BayesEstimator` <api-BayesEstimator>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
BayesEstimator()
```

Bayesian estimator: the parameter maximizing the marginal posterior of the regularization parameter (minimum of the functional `J_MAP`, Pereira et al.).

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `Heuristic` <api-Heuristic>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
Heuristic()
```

Heuristic initial value of the regularization parameter (Aucejo & De Smet, MSSP 2017): the geometric mean of the extreme admissible values derived from the (generalized) singular values of `(H, L)`. Only used to initialize iterative algorithms (`init` keyword).

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `LambdaMax` <api-LambdaMax>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
LambdaMax()
```

Initial regularization parameter `λ = ‖2 Hᴴ y‖∞` (whitened `H` and `y`): the smallest value for which the ℓ1-regularized solution is zero (Kim et al., IEEE 2007). Only used to initialize iterative algorithms (`init` keyword); default initialization of `ElementWise()`.

]

= Algorithms

Resolution algorithms, passed as third argument of `solve`.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `IR` <api-IR>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
IR()
```

Iteratively Reweighted algorithm derived from the direct application of the first-order optimality condition (Aucejo & De Smet, JSV 2020; MSSP 2020). Algorithm of the additive and multiplicative formulations.

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `MAP` <api-MAP>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
MAP()
```

Maximum a posteriori estimate of the hierarchical Bayesian model — optimal Bayesian regularization (Aucejo & De Smet, MSSP 2019). Default algorithm for #link(<api-Bayesian>)[`Bayesian`].

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `ElementWise` <api-ElementWise>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
ElementWise()
```

Element-wise Bayesian algorithm (Relevance Vector Regression) for #link(<api-Bayesian>)[`Bayesian`]: each component of `L u` has its own Gaussian prior precision, estimated jointly with `u` (sparse Bayesian learning).

]

= Resolution

Solving a problem and inspecting its solution.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `solve` <api-solve>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
solve(prob, formulation[, algorithm]; R = nothing, precision = false, tol = 1e-4, maxiter = 100, init)
```

Solve the inverse problem `prob` with the given regularization `formulation`. When the algorithm is omitted, `default_algorithm` is used.

For a space–frequency problem (block-diagonal operator), the resolution strategy depends on the regularization terms (dispatch on `sample_coupling`):

- `Lq` terms only: frequency-by-frequency resolution — each frequency is an independent problem with its own regularization parameter(s), hyperparameters and stopping criterion, solved in parallel;
- at least one `Lpq` term: the mixed norm couples the frequencies, the space–frequency problem is solved as a whole.

Single-frequency and space–time problems are always solved as a whole.

*Inputs*

- `prob` - `InverseProblem`: inverse problem
- `formulation` - `Formulation`: `Additive`, `Multiplicative` or `Bayesian`
- `algorithm` - `Algorithm`: resolution algorithm (`IR()` for `Additive` and `Multiplicative`, `MAP()` or `ElementWise()` for `Bayesian`)
- `R` - `Nothing`, `UniformScaling`, `Diagonal`, `AbstractVector` or `AbstractMatrix`: known diagonal noise covariance — not provided (default, equivalent to `I`), `σ²I`, or the noise variances (vector, or matrix with the layout of `y`). The data fidelity is `‖y - Hu‖²_{R⁻¹}`
- `precision` - `Bool`: also return the posterior precision `Q = τn (Hᴴ H + Lᴴ D L)` of the Gaussian approximation at the solution (`sol.precision`: dense, sparse for an assembled space–frequency problem, or one `nr × nr` block per frequency for independent or block-diagonal frequencies) and the posterior model (`sol.posterior`) used by #link(<api-posterior_std>)[`posterior_std`], #link(<api-posterior_covariance>)[`posterior_covariance`], #link(<api-posterior_samples>)[`posterior_samples`] and #link(<api-credible_interval>)[`credible_interval`]. For the additive and multiplicative formulations, `τn = 1` (`R` is then taken as the exact noise covariance)
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Real` or `ParameterChoice`: initial Tikhonov solution — initial regularization parameter, or rule selecting it (`Heuristic()`, `LCurve()`, `GCV()`, `BayesEstimator()`)

*Outputs*

- `sol` - `InverseSolution`: solution, regularization parameter(s), hyperparameters, convergence information and, if requested, posterior precision. For a frequency-by-frequency resolution, the quantities of each frequency are stacked as columns (`λ`: vector of length `nf`, or `G × nf` / `nr × nf` matrix; `hyper.τn`, `hyper.τs`, `hyper.q` likewise), `hyper` also contains the `iterations`, `converged` and `history` of each frequency, `iterations`/`converged`/`history` refer to the slowest frequency / all frequencies, and `precision` is a vector of `nr × nr` matrices

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `InverseSolution` <api-InverseSolution>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
InverseSolution
```

Result of #link(<api-solve>)[`solve`]. See #link(<api-posterior_std>)[`posterior_std`], #link(<api-posterior_covariance>)[`posterior_covariance`], #link(<api-posterior_samples>)[`posterior_samples`] and #link(<api-credible_interval>)[`credible_interval`] for the posterior uncertainty.

*Fields*

- `u` - `AbstractVector` or `AbstractMatrix`: reconstructed unknown (vector, or `nr × ns` matrix, with the layout of the problem)
- `λ` - `Real` or `Vector{<:Real}`: final regularization parameter(s) — the adaptive one for the multiplicative formulations, `τs_g / τn` for the Bayesian formulation (MAP `½‖y - Hu‖²_R + Σ_g (λ_g/q_g)‖L_g u_g‖_q^q` for real data, `‖y - Hu‖²_R + Σ_g λ_g‖L_g u_g‖_q^q` for complex data)
- `hyper` - `NamedTuple`: other estimated quantities (`τn`, `τs`, `q` for the Bayesian formulation, empty otherwise)
- `iterations` - `Int`: number of iterations of the outer loop
- `converged` - `Bool`: whether the stopping criterion was met
- `history` - `Vector{Float64}`: stopping criterion (relative `ℓ1` variation of `u`) at each iteration
- `precision` - `Nothing`, `Hermitian`, `SparseMatrixCSC` or `Vector{<:Hermitian}`: posterior precision `Q = τn (Hᴴ H + Lᴴ D L)` of `vec(u)` when `solve(...; precision = true)` — dense for a dense operator, sparse for an assembled space–frequency problem, one `nr × nr` block per sample for independent or block-diagonal samples, `nothing` otherwise
- `posterior` - `Nothing`, `Posterior` or `SamplewisePosterior`: posterior model (factorization of the precision) used by #link(<api-posterior_std>)[`posterior_std`], #link(<api-posterior_covariance>)[`posterior_covariance`] and #link(<api-posterior_samples>)[`posterior_samples`]

]

= Uncertainty quantification

Posterior quantities of a solution obtained with `solve(...; precision = true)`.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `posterior_std` <api-posterior_std>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
posterior_std(sol)
```

Posterior standard deviation of each entry of the solution: square root of the diagonal of the posterior covariance `Q⁻¹`, obtained from its diagonal blocks (block solves, see #link(<api-posterior_covariance>)[`posterior_covariance`]) — `Q⁻¹` is never formed.

*Inputs*

- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`

*Outputs*

- `σ` - `Array{<:Real}`: standard deviations, with the layout of `sol.u`

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `posterior_covariance` <api-posterior_covariance>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
posterior_covariance(sol; full = false)
```

Posterior covariance `Q⁻¹` of the solution, computed from the factorization of the precision.

- `full = false`: marginal covariances of the samples (frequencies or time steps), i.e. the `nr × nr` diagonal blocks `Σ_jj` of the covariance of `vec(u)`. Each block is obtained by applying `A⁻¹` to the `nr` corresponding columns of the identity (block solves with the factorization of the solver: normal equations, Woodbury identity, sparse Cholesky), so that the `n × n` covariance is never formed;
- `full = true`: covariance of `vec(u)`, `(nr·ns) × (nr·ns)`, with the correlations between samples (sparse block-diagonal matrix when the samples are independent or block-diagonal).

*Inputs*

- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `full` - `Bool`: full covariance instead of the marginal blocks

*Outputs*

- `Σ` - `Hermitian` (single sample, or `full = true`), `Vector{<:Hermitian}` (one block per sample) or `SparseMatrixCSC` (`full = true`, independent samples)

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `posterior_samples` <api-posterior_samples>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
posterior_samples(sol, nsamples = 10_000; rng = Random.default_rng())
```

Samples of the Gaussian approximation of the posterior `N(u, Q⁻¹)` (circular complex normal for complex unknowns), with `Q = τn A`, `A = Hᴴ H + Lᴴ D L`. The sampler is chosen by dispatch on the factorization of `A` used by the solver:

- Cholesky factor `A = Rᴴ R` available (normal equations, sparse Cholesky): *standard sampling* `u_sample = u + R⁻¹ z / √τn`, `z ~ N(0, I)` — one triangular solve per sample (`R⁻¹ R⁻ᴴ = A⁻¹`; with a lower factor `A = L Lᴴ`, this reads `L⁻ᴴ z`);
- Woodbury identity (under-determined problems, no Cholesky factor of `A`): *perturbation–optimization*, with `ε₁ ~ N(0, I)`, `ε₂ ~ N(0, I)` of the sizes of the data and of `L u`,
  ```
  u_sample = u + A⁻¹ (Hᴴ ε₁ + Lᴴ D^{1/2} ε₂) / √τn,
  ```
  whose covariance is `A⁻¹ (Hᴴ H + Lᴴ D L) A⁻¹ / τn = Q⁻¹`, at the cost of one application of `A⁻¹` (Woodbury) per sample, without forming nor factorizing the `n × n` matrix `A`.

Both give exact samples of `N(u, Q⁻¹)`; the joint samples keep the correlations between frequencies or time steps (block-diagonal problems are sampled block by block).

*Inputs*

- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

*Outputs*

- `X` - `Array`: samples, of size `(size(sol.u)..., nsamples)`

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `credible_interval` <api-credible_interval>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
credible_interval(sol; level = 0.95, nsamples = 10_000, rng = Random.default_rng())
```

Bayesian credible interval of each entry of the solution, estimated by sampling the posterior (#link(<api-posterior_samples>)[`posterior_samples`]): empirical quantiles `(1 - level)/2` and `(1 + level)/2` of the samples, entry by entry. For complex unknowns, the quantiles are computed separately on the real and imaginary parts.

*Inputs*

- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `level` - `Real`: credibility level, in `]0, 1[`
- `nsamples` - `Integer`: number of samples
- `rng` - `AbstractRNG`: random number generator

*Outputs*

- `lower` - `Array`: lower bounds, with the layout of `sol.u` (complex for complex unknowns: `real(lower) ≤ real(u)`, `imag(lower) ≤ imag(u)`)
- `upper` - `Array`: upper bounds, with the layout of `sol.u`

]

