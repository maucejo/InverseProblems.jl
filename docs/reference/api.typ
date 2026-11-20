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
  [#link(<api-space_sample_operator>)[`space_sample_operator`]], [function], [Smoothing operators], [Operator of a vectorized space–frequency or space–time problem, acting on `u = vec(U)` (`U`: `nu × ns`, rows = reconstruction points, columns = frequencies or time samples):],
  [#link(<api-sample_rows>)[`sample_rows`]], [function], [Smoothing operators], [Indices, in a vectorized problem `u = vec(U)` (`U`: `nu × ns`), of the points `rows` at all the frequencies or time samples (port of `define_group_mixed`).],
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
  [#link(<api-ill_posedness>)[`ill_posedness`]], [function], [Ill-posedness], [Indicators of the ill-posedness of the linear inverse problem `y = H u + n`, computed from the singular value decomposition `H = Σᵢ σᵢ vᵢ wᵢᴴ` and the Fourier coefficients of the data `βᵢ = |vᵢᴴ y|`, whose least-squares solution reads `u_LS = Σᵢ (vᵢᴴ y / σᵢ) wᵢ`.],
  [#link(<api-IllPosedness>)[`IllPosedness`]], [type], [Ill-posedness], [Result of #link(<api-ill_posedness>)[`ill_posedness`] for one linear system `y = H u + n` (the whitened one when the noise covariance is given).],
  [#link(<api-picardplot>)[`picardplot`]], [function], [Visualization], [Picard plot of an inverse problem: singular values `σᵢ` of the (whitened) operator, Fourier coefficients of the data `βᵢ = |vᵢᴴ y|` and ratio `β̄ᵢ / σᵢ` of the smoothed coefficients to the singular values (field `βm` of the analysis), against the index `i` (integer ticks), on a logarithmic scale, with the noise level `η` (dashed line) and the Picard index (dotted line).],
  [#link(<api-picardplot-bang>)[`picardplot!`]], [function], [Visualization], [Draw the Picard plot of `ip` (see #link(<api-picardplot>)[`picardplot`]) in an existing Makie `Axis`; its `y` scale is set to `log10` if it is not already logarithmic.],
  [#link(<api-relative_error>)[`relative_error`]], [function], [Performance indicators], [Relative error (trueness) of the reconstruction `u` with respect to the reference `uref`, in percent:],
  [#link(<api-correlation>)[`correlation`]], [function], [Performance indicators], [Correlation (cosine similarity of shape) between the reconstruction `u` and the reference `uref`, in percent:],
  [#link(<api-precision_metric>)[`precision_metric`]], [function], [Performance indicators], [Precision of the estimate, in the sense of ISO 5725 (spread of the estimate, not to be confused with the precision matrix `Q`): root mean posterior variance of the `n` unknowns,],
  [#link(<api-accuracy_metric>)[`accuracy_metric`]], [function], [Performance indicators], [Accuracy of the estimate, in the sense of ISO 5725 (trueness and precision combined): root mean squared error, per unknown, of the posterior distribution `N(u, Σ)` with respect to the reference,],
  [#link(<api-performance_ratio>)[`performance_ratio`]], [function], [Performance indicators], [Overall performance of several methods compared on several indicators: geometric mean, over the indicators, of the ratio between the indicator of each method and the best one,],
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
  [`AbstractMatrix` (`nm × nu`)], [`AbstractVector`], [vector (`nu`)], [one frequency],
  [`AbstractVector` of `nm × nu` matrices], [`AbstractMatrix` (`nm × nf`)], [matrix `nu × nf`], [space–frequency],
  [`AbstractArray{T,3}` (`nm × nu × nf`)], [`AbstractMatrix` (`nm × nf`)], [matrix `nu × nf`], [space–frequency],
  [`AbstractMatrix` (`nm·nt × nu·nt`)], [`AbstractMatrix` (`nm × nt`)], [matrix `nu × nt`], [space–time],
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
space_sample_operator(nu, ns; L = I, Lc = I)
```

Operator of a vectorized space–frequency or space–time problem, acting on `u = vec(U)` (`U`: `nu × ns`, rows = reconstruction points, columns = frequencies or time samples):

```julia
#| eval: false
L = kron(Lc, L),    L * vec(U) = vec(L * U * transpose(Lc))
```

- `Lc = Lf`, `L = I`: smoothing along the frequencies / time samples at each point (e.g. `Lf = diff_operator(ns, 1; h = Δf)`);
- `Lc = I`, `L = Lx`: spatial smoothing at each frequency / time sample;
- both: product of the two operators (mixed derivative, MSSP 2018 Eq. (15)). For a sum of penalties, stack the operators: `[space_sample_operator(nu, ns; Lc = Lf); space_sample_operator(nu, ns; L = Lx)]`.

*Inputs*

- `nu` - `Integer`: number of reconstruction points (rows of `U`)
- `ns` - `Integer`: number of frequencies or time samples (columns of `U`)
- `L` - `AbstractMatrix` or `UniformScaling`: operator acting on the rows (`nu` columns)
- `Lc` - `AbstractMatrix` or `UniformScaling`: operator acting on the columns (`ns` columns)

*Outputs*

- `L` - `SparseMatrixCSC`: operator of size `(size(Lc, 1)·size(L, 1)) × (ns·nu)`

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `sample_rows` <api-sample_rows>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
sample_rows(rows, nu, ns)
```

Indices, in a vectorized problem `u = vec(U)` (`U`: `nu × ns`), of the points `rows` at all the frequencies or time samples (port of `define_group_mixed`). The points vary fastest, so that `space_sample_operator(length(rows), ns; …)` is the smoothing operator of the group.

*Inputs*

- `rows` - `AbstractVector{<:Integer}`: reconstruction points of the group
- `nu` - `Integer`: number of reconstruction points
- `ns` - `Integer`: number of frequencies or time samples

*Outputs*

- `idx` - `Vector{Int}`: indices `rows[k] + (j - 1) nu`, `k` fastest, `j = 1:ns`

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
- `precision` - `Bool`: also return the posterior precision `Q = τn (Hᴴ H + Lᴴ D L)` of the Gaussian approximation at the solution (`sol.precision`: dense, sparse for an assembled space–frequency problem, or one `nu × nu` block per frequency for independent or block-diagonal frequencies) and the posterior model (`sol.posterior`) used by #link(<api-posterior_std>)[`posterior_std`], #link(<api-posterior_covariance>)[`posterior_covariance`], #link(<api-posterior_samples>)[`posterior_samples`] and #link(<api-credible_interval>)[`credible_interval`]. For the additive and multiplicative formulations, `τn = 1` (`R` is then taken as the exact noise covariance)
- `tol` - `Real`: tolerance on the relative `ℓ1` variation of `u`
- `maxiter` - `Integer`: maximum number of iterations
- `init` - `Real` or `ParameterChoice`: initial Tikhonov solution — initial regularization parameter, or rule selecting it (`Heuristic()`, `LCurve()`, `GCV()`, `BayesEstimator()`)

*Outputs*

- `sol` - `InverseSolution`: solution, regularization parameter(s), hyperparameters, convergence information and, if requested, posterior precision. For a frequency-by-frequency resolution, the quantities of each frequency are stacked as columns (`λ`: vector of length `nf`, or `G × nf` / `nu × nf` matrix; `hyper.τn`, `hyper.τs`, `hyper.q` likewise), `hyper` also contains the `iterations`, `converged` and `history` of each frequency, `iterations`/`converged`/`history` refer to the slowest frequency / all frequencies, and `precision` is a vector of `nu × nu` matrices

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

- `u` - `AbstractVector` or `AbstractMatrix`: reconstructed unknown (vector, or `nu × ns` matrix, with the layout of the problem)
- `λ` - `Real` or `Vector{<:Real}`: final regularization parameter(s) — the adaptive one for the multiplicative formulations, `τs_g / τn` for the Bayesian formulation (MAP `½‖y - Hu‖²_R + Σ_g (λ_g/q_g)‖L_g u_g‖_q^q` for real data, `‖y - Hu‖²_R + Σ_g λ_g‖L_g u_g‖_q^q` for complex data)
- `hyper` - `NamedTuple`: other estimated quantities (`τn`, `τs`, `q` for the Bayesian formulation, empty otherwise)
- `iterations` - `Int`: number of iterations of the outer loop
- `converged` - `Bool`: whether the stopping criterion was met
- `history` - `Vector{Float64}`: stopping criterion (relative `ℓ1` variation of `u`) at each iteration
- `precision` - `Nothing`, `Hermitian`, `SparseMatrixCSC` or `Vector{<:Hermitian}`: posterior precision `Q = τn (Hᴴ H + Lᴴ D L)` of `vec(u)` when `solve(...; precision = true)` — dense for a dense operator, sparse for an assembled space–frequency problem, one `nu × nu` block per sample for independent or block-diagonal samples, `nothing` otherwise
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

- `full = false`: marginal covariances of the samples (frequencies or time steps), i.e. the `nu × nu` diagonal blocks `Σ_jj` of the covariance of `vec(u)`. Each block is obtained by applying `A⁻¹` to the `nu` corresponding columns of the identity (block solves with the factorization of the solver: normal equations, Woodbury identity, sparse Cholesky), so that the `n × n` covariance is never formed;
- `full = true`: covariance of `vec(u)`, `(nu·ns) × (nu·ns)`, with the correlations between samples (sparse block-diagonal matrix when the samples are independent or block-diagonal).

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

= Ill-posedness

Discrete Picard condition and conditioning of a problem, independent of the regularization.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `ill_posedness` <api-ill_posedness>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
ill_posedness(prob; R = nothing, window = 5, τ = 2)
ill_posedness(H, y; R = nothing, window = 5, τ = 2)
```

Indicators of the ill-posedness of the linear inverse problem `y = H u + n`, computed from the singular value decomposition `H = Σᵢ σᵢ vᵢ wᵢᴴ` and the Fourier coefficients of the data `βᵢ = |vᵢᴴ y|`, whose least-squares solution reads `u_LS = Σᵢ (vᵢᴴ y / σᵢ) wᵢ`.

*Discrete Picard condition.* Writing `vᵢᴴ y = vᵢᴴ H u + vᵢᴴ n`, the coefficients of the exact data `|vᵢᴴ H u| = σᵢ |wᵢᴴ u|` decay with the singular values, whereas those of the noise stay at the level `η = √E|vᵢᴴ n|²`. A meaningful solution can only be recovered if the exact coefficients decay _faster_ than the singular values (discrete Picard condition); the components below the noise level must be filtered out by the regularization. The function returns:

- the noise level `η`: `1` when the noise covariance `R` is given (the problem is whitened, so that `E|vᵢᴴ n|² = 1` for real and circular complex noise); otherwise it is estimated from the part of the data outside the range of `H` when `m > rank(H)`, `η² = (‖y‖² - Σᵢ βᵢ²) / (m - rank(H))`, or else from the RMS value of the last quarter of the coefficients (assumed dominated by the noise);
- the Picard index `k`: number of leading components whose coefficients, smoothed by a moving geometric mean over `window` components (`β̄ᵢ`), exceed `τ η`;
- the Picard slope `b`, from the fit `log β̄ᵢ ≈ a + b log σᵢ` over these `k` components: the coefficients decay as `σᵢᵇ`, and the condition holds for `b > 1`.

*Other indicators.* Condition number `σ₁ / σ_min`, numerical rank, dimension of the null space (unknown components invisible in the data, e.g. underdetermined problems), polynomial decay rate `α` of the singular values (`σᵢ ≈ σ₁ i^(-α)`) and the resulting classification of the problem: mildly (`α ≤ 1`), moderately (`α > 1`) or severely ill-posed (faster, exponential decay).

The problem is analyzed as the solver sees it: whitened by `R` (same values as the keyword `R` of #link(<api-solve>)[`solve`]), one analysis per sample for a block-diagonal space–frequency problem, a single one for a space–time or assembled problem (an assembled sparse operator is converted to a dense matrix: the SVD costs `O(m n min(m, n))`). The indicators do not depend on the regularization.

*Inputs*

- `prob` - `InverseProblem`: problem to analyze
- `H`, `y` - `AbstractMatrix`, `AbstractVector`: forward operator and data of a single system (same as `ill_posedness(InverseProblem(H, y))`)
- `R` - noise covariance, as in #link(<api-solve>)[`solve`] (`nothing`: unknown, the noise level is estimated)
- `window` - `Integer`: width of the moving geometric mean of the coefficients (odd, `1` for no smoothing)
- `τ` - `Real`: threshold, in units of the noise level, above which a component emerges from the noise

*Outputs*

- `ip` - #link(<api-IllPosedness>)[`IllPosedness`] (single system) or `Vector{IllPosedness}` (one per sample of a block-diagonal problem)

*Example*

```julia
#| eval: false
ip = ill_posedness(InverseProblem(H, y); R = σ² * I)
ip.picard_satisfied, ip.picard_index, ip.condition_number
# Picard plot: ip.σ, ip.β and ip.βm ./ ip.σ against the index, on a log scale
```

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `IllPosedness` <api-IllPosedness>

#el("span", class: "se-api-kind se-api-kind-type")[Type]
]

```julia
#| eval: false
IllPosedness
```

Result of #link(<api-ill_posedness>)[`ill_posedness`] for one linear system `y = H u + n` (the whitened one when the noise covariance is given).

*Fields*

- `σ` - `Vector{Float64}`: singular values of `H`, in decreasing order (`min(m, n)` values)
- `β` - `Vector{Float64}`: Fourier coefficients of the data `|vᵢᴴ y|`, `vᵢ` being the left singular vectors
- `βm` - `Vector{Float64}`: coefficients smoothed by a moving geometric mean (window `window`)
- `noise_level` - `Float64`: noise level `η` of the coefficients, RMS value of `|vᵢᴴ n|` (`1` when the noise covariance is given, estimated otherwise)
- `noise_estimated` - `Bool`: whether `η` has been estimated from the data
- `picard_index` - `Int`: number `k` of leading components whose smoothed coefficient `β̄ᵢ` exceeds `τ η`, i.e. that emerge from the noise — the components that a regularization can recover
- `picard_slope` - `Float64`: slope `b` of the least-squares fit `log β̄ᵢ ≈ a + b log σᵢ` over the `k` first components (`NaN` if `k < 3`); the coefficients decay as `σᵢᵇ`
- `picard_satisfied` - `Bool`: discrete Picard condition, `b > 1`: above the noise level, the coefficients decay faster than the singular values
- `condition_number` - `Float64`: `σ₁ / σ_min` (`Inf` for a singular operator)
- `numerical_rank` - `Int`: number of singular values above `max(m, n) eps σ₁`
- `nullity` - `Int`: dimension of the null space, `n - numerical_rank`: components of the unknown that do not affect the data at all
- `decay_rate` - `Float64`: degree `α` of the polynomial decay `σᵢ ≈ σ₁ i^(-α)` (least-squares fit over the numerical rank)
- `severity` - `Symbol`: `:mild` (`α ≤ 1`), `:moderate` (`α > 1`) or `:severe` (the singular values decay faster than polynomially: an exponential `σᵢ ≈ σ₁ e^(-μ i)` fits better)

]

= Visualization

Plots, provided by the Makie extension of the package: load a Makie backend first, e.g. `using CairoMakie`.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `picardplot` <api-picardplot>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
picardplot(ip::IllPosedness; components = :rank, smoothed = true, figure = (;), axis = (;), legend = true)
picardplot(prob::InverseProblem; R = nothing, window = 5, τ = 2, sample = nothing, kwargs...)
picardplot(H, y; R = nothing, kwargs...)
```

Picard plot of an inverse problem: singular values `σᵢ` of the (whitened) operator, Fourier coefficients of the data `βᵢ = |vᵢᴴ y|` and ratio `β̄ᵢ / σᵢ` of the smoothed coefficients to the singular values (field `βm` of the analysis), against the index `i` (integer ticks), on a logarithmic scale, with the noise level `η` (dashed line) and the Picard index (dotted line). See #link(<api-ill_posedness>)[`ill_posedness`] for the computation of these quantities.

The discrete Picard condition holds when, above the noise level, the coefficients `βᵢ` decay faster than the singular values, i.e. when the ratio `β̄ᵢ / σᵢ` decreases; the components beyond the Picard index are dominated by the noise. The title summarizes the analysis: Picard index, Picard slope and condition number.

Requires Makie: load a backend first, e.g. `using CairoMakie`.

*Inputs*

- `ip` - `IllPosedness`: result of #link(<api-ill_posedness>)[`ill_posedness`]
- `prob` - `InverseProblem`, or `H`, `y`: problem to analyze, with the keywords `R`, `window` and `τ` of #link(<api-ill_posedness>)[`ill_posedness`]; `sample` selects the frequency of a block-diagonal space–frequency problem (one analysis per sample)
- `components` - `:rank` (components up to the numerical rank, default), `:all`, or an `Integer`: number of components displayed
- `smoothed` - `Bool`: ratio computed with the smoothed coefficients `β̄ᵢ` (default) or the raw ones `βᵢ`
- `figure`, `axis` - `NamedTuple`: attributes of the `Figure` and of the `Axis`
- `legend` - `Bool`: legend below the axis

*Outputs*

- `fig` - `Makie.Figure`

*Example*

```julia
#| eval: false
using CairoMakie
ip = ill_posedness(prob; R = σ² * I)
picardplot(ip)
```

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `picardplot!` <api-picardplot-bang>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
picardplot!(ax, ip::IllPosedness; components = :rank, smoothed = true, noise = true, index = true,
            noisecolor = :gray50, indexcolor = :gray50)
```

Draw the Picard plot of `ip` (see #link(<api-picardplot>)[`picardplot`]) in an existing Makie `Axis`; its `y` scale is set to `log10` if it is not already logarithmic. The curves carry labels, so that `axislegend(ax)` or `Legend(fig[...], ax)` builds the legend.

*Inputs*

- `ax` - `Makie.Axis`: target axis
- `ip` - `IllPosedness`: result of #link(<api-ill_posedness>)[`ill_posedness`]
- `components`, `smoothed` - see #link(<api-picardplot>)[`picardplot`]
- `noise` - `Bool`: draw the noise level `η`
- `index` - `Bool`: draw the Picard index
- `noisecolor`, `indexcolor`: colors of these two lines

*Outputs*

- `ax` - `Makie.Axis`

]

= Performance indicators

Trueness, precision and accuracy of a reconstruction with respect to a reference field.

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `relative_error` <api-relative_error>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
relative_error(uref, u; p = 1)
```

Relative error (trueness) of the reconstruction `u` with respect to the reference `uref`, in percent:

```
RE = 100 ‖uref - u‖_p / ‖uref‖_p
```

The norm is taken over all the entries (all the samples for a space–frequency or space–time unknown). Apply it to a column or a row to get a local indicator, e.g. `relative_error(uref[:, k], u[:, k])` at the `k`-th frequency.

*Inputs*

- `uref` - `AbstractArray`: reference field
- `u` - `AbstractArray` or `InverseSolution`: reconstructed field, with the layout of `uref`
- `p` - `Real`: order of the norm (`1` by default, less sensitive to a few large errors than `2`)

*Outputs*

- `RE` - `Real`: relative error, in %

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `correlation` <api-correlation>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
correlation(uref, u)
```

Correlation (cosine similarity of shape) between the reconstruction `u` and the reference `uref`, in percent:

```
C = 100 |⟨uref, u⟩| / (‖uref‖₂ ‖u‖₂)
```

It equals 100 % when `u` is proportional to `uref`, whatever the scale factor, and is insensitive to a global phase for complex fields.

*Inputs*

- `uref` - `AbstractArray`: reference field
- `u` - `AbstractArray` or `InverseSolution`: reconstructed field, with the layout of `uref`

*Outputs*

- `C` - `Real`: correlation, in % (`0 ≤ C ≤ 100`)

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `precision_metric` <api-precision_metric>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
precision_metric(sol)
precision_metric(Q)
precision_metric(Σ; covariance = true)
```

Precision of the estimate, in the sense of ISO 5725 (spread of the estimate, not to be confused with the precision matrix `Q`): root mean posterior variance of the `n` unknowns,

```
P = √( tr(Σ) / n ),    Σ = Q⁻¹,
```

i.e. the quadratic mean of the posterior standard deviations. It has the unit of the unknown; the smaller, the more precise. It does not depend on the layout of the unknown: for a space–frequency or space–time problem, the mean runs over all the entries of `u` (all points, all samples), and the correlations between unknowns are not used.

The matrices are real (real data) or complex Hermitian (complex data); in both cases the variances are real. The variances `diag(Q⁻¹)` are computed from a Cholesky factorization of `Q` (dense, or sparse for a `SparseMatrixCSC`), without forming `Q⁻¹`.

*Inputs*

- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)` (variances by block solves, see #link(<api-posterior_std>)[`posterior_std`])
- `Q` - posterior precision matrix of `vec(u)` (`AbstractMatrix`, real symmetric or complex Hermitian positive definite: dense, `Hermitian`/`Symmetric`, `SparseMatrixCSC` or `Diagonal`), or vector of precision blocks, one per sample when the samples are independent (e.g. `sol.precision`)
- `Σ` - same, for a posterior covariance, real or complex (e.g. #link(<api-posterior_covariance>)[`posterior_covariance`], Kalman filter covariances `P_k`)
- `covariance` - `Bool`: `true` if the matrices are covariances, `false` (default) if they are precisions

*Outputs*

- `P` - `Real`: precision indicator, in the unit of the unknown

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `accuracy_metric` <api-accuracy_metric>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
accuracy_metric(uref, sol)
accuracy_metric(uref, u, Q)
accuracy_metric(uref, u, Σ; covariance = true)
```

Accuracy of the estimate, in the sense of ISO 5725 (trueness and precision combined): root mean squared error, per unknown, of the posterior distribution `N(u, Σ)` with respect to the reference,

```
A = √( (‖uref - u‖₂² + tr(Σ)) / n ),    Σ = Q⁻¹,
```

where `n = length(u)`: `‖uref - u‖₂² + tr(Σ)` is the expectation of `‖uref - x‖₂²` for `x ~ N(u, Σ)` (squared bias plus variance). Without uncertainty (`Σ = 0`) it reduces to the RMS error of `u`, and `A² = RMS(uref - u)² + P²` with `P` the #link(<api-precision_metric>)[`precision_metric`]. It has the unit of the unknown; the smaller, the more accurate. As for `precision_metric`, the mean runs over all the entries of `u`, whatever its layout.

*Inputs*

- `uref` - `AbstractArray`: reference field (vector, or `nu × ns` matrix)
- `sol` - `InverseSolution`: solution obtained with `solve(...; precision = true)`
- `u` - `AbstractArray`: reconstructed field, with the layout of `uref`
- `Q` / `Σ` - posterior precision (default) or covariance (`covariance = true`) of `vec(u)`, real or complex: one `n × n` matrix, or a vector of blocks, one per sample; see #link(<api-precision_metric>)[`precision_metric`]
- `covariance` - `Bool`: `true` if the matrices are covariances, `false` (default) if they are precisions

*Outputs*

- `A` - `Real`: accuracy indicator, in the unit of the unknown

]

#el("div", class: "se-api-entry")[
#el("div", class: "se-api-head")[
== `performance_ratio` <api-performance_ratio>

#el("span", class: "se-api-kind se-api-kind-function")[Function]
]

```julia
#| eval: false
performance_ratio(indicators)
```

Overall performance of several methods compared on several indicators: geometric mean, over the indicators, of the ratio between the indicator of each method and the best one,

```
r_m = ( Π_j I_mj / min_m' I_m'j )^(1 / n_ind).
```

The best possible value is `1`, reached by a method that is the best for every indicator. The indicators must all be such that _smaller is better_ (e.g. #link(<api-relative_error>)[`relative_error`], #link(<api-precision_metric>)[`precision_metric`], #link(<api-accuracy_metric>)[`accuracy_metric`], or `100 - correlation`) and strictly positive.

*Inputs*

- `indicators` - `AbstractMatrix{<:Real}`: one row per method, one column per indicator

*Outputs*

- `r` - `Vector{<:Real}`: performance ratio of each method (`≥ 1`)

*Example*

```julia
#| eval: false
I = [relative_error(uref, s) for s in sols] # one column per indicator
J = [accuracy_metric(uref, s) for s in sols]
performance_ratio([I J])
```

]

