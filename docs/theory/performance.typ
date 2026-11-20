#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Performance indicators])
#metadata((description: "Trueness, precision and accuracy of a reconstruction, and comparison of several methods.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

When the sources are known — numerical validation, laboratory test with a measured excitation — the quality of a reconstruction $hat(u)$ can be assessed against the reference $u_"ref"$. SourceEstimation.jl provides indicators built on the vocabulary of metrology @ISO5725 @VIM2012, which distinguishes three properties of an estimate:

- *trueness*: closeness of the estimate to the reference value (systematic error, or bias);
- *precision*: spread of the possible estimates around their mean (random error);
- *accuracy*: closeness of the possible estimates to the reference, i.e. trueness and precision combined.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: left,
    table.header([*Property*], [*Function*], [*Needs*]),
    [trueness], [`relative_error`, `correlation`], [$u_"ref"$, $hat(u)$],
    [precision], [`precision_metric`], [posterior precision $Q$ (or covariance $Sigma$)],
    [accuracy], [`accuracy_metric`], [$u_"ref"$, $hat(u)$ and $Q$ (or $Sigma$)],
    [overall ranking], [`performance_ratio`], [a table of indicators],
  ),
  caption: [Performance indicators of SourceEstimation.jl.],
)

In what follows, $hat(u)$ and $u_"ref"$ gather the $n$ unknowns of the problem — all the reconstruction points of all the samples (frequencies or time steps) for space–frequency and space–time problems, stacked as $u = "vec"(U)$ (see the #link("inverse-problems.html#a-common-layout")[common layout]). Real and complex fields are handled alike, with the #link("inverse-problems.html#real-or-complex-data")[notation] of the theory pages.

= Trueness

== Relative error

The relative error compares the reconstruction to the reference in the $ell_p$ norm,

$ "RE" = 100 norm(u_"ref" - hat(u))_p / norm(u_"ref")_p quad (%). $ <eq-re>

`relative_error(uref, u; p = 1)` uses the $ell_1$ norm by default: an error of a few large entries (e.g. a peak at a natural frequency) weighs less than with the $ell_2$ norm, and the indicator better reflects the error on the whole field. Applied to a column $u_(:, k)$ or a row $u_(i, :)$, @eq-re gives a *local* indicator, at a given frequency or at a given point.

== Correlation

The relative error mixes errors on the *shape* and on the *amplitude* of the field. The correlation isolates the former:

$ C = 100 (|u_"ref"^H hat(u)|) / (norm(u_"ref")_2 norm(hat(u))_2) quad (%). $ <eq-corr>

By the Cauchy–Schwarz inequality, $0 <= C <= 100 %$, with $C = 100 %$ if and only if $hat(u) = alpha u_"ref"$ for some scalar $alpha$. A reconstruction underestimated by a constant factor — typical of a too strong regularization — thus has a large relative error but a correlation close to $100 %$. For complex fields, $C$ is insensitive to a global phase ($alpha in bb(C)$): it is the modal assurance criterion (MAC) of the two fields, in percent.

= Precision

The trueness indicators only use the estimate $hat(u)$. When the problem is solved with `precision = true`, the #link("uncertainty.html")[Gaussian approximation of the posterior] $cal(N)(hat(u), Sigma)$, $Sigma = Q^(-1)$, also describes how much the estimate can vary. Its spread is measured by the root mean posterior variance of the unknowns,

$ P = sqrt(1/n tr(Sigma)) = sqrt(1/n sum_(i=1)^n sigma_i^2), quad sigma_i^2 = (Q^(-1))_(i i), $ <eq-precision>

i.e. the quadratic mean of the posterior standard deviations returned by `posterior_std`. $P$ has the unit of the unknown; the smaller, the more precise. It does not use the reference: it can be computed for any reconstruction, including on measured data.

#calepin.elements.callout(kind: "note")[
The word *precision* has two meanings here: the precision of an estimate in the sense of metrology (@eq-precision), and the precision matrix $Q = Sigma^(-1)$ of a Gaussian distribution. A large precision matrix means a small spread, hence a small value of $P$.
]

== Variances from the precision matrix

`precision_metric` and `accuracy_metric` take the posterior precision $Q$ by default (e.g. `sol.precision`), or a covariance $Sigma$ with `covariance = true` (e.g. the covariances $P_k$ of a Kalman filter). They accept one matrix for all the unknowns, or a vector of blocks, one per sample. Only the diagonal of $Sigma$ is needed: the correlations between unknowns are not used, and $Q^(-1)$ is never formed. With the Cholesky factorization $Q = R^H R$ — symmetric for real data, Hermitian for complex data — the variances read

$ (Q^(-1))_(i i) = e_i^H R^(-1) R^(-H) e_i = norm(R^(-H) e_i)_2^2, $

where $e_i$ is the $i$-th column of the identity: they are the squared norms of the columns of $X$, solution of the lower-triangular system $R^H X = I$, obtained by forward substitution. No explicit inverse is computed, and the variances are positive by construction. For a sparse precision (assembled space–frequency or space–time problems), $X$ is computed by blocks of columns with the sparse Cholesky factorization, keeping only the diagonal entries. For an `InverseSolution`, the variances are obtained from the factorization already computed by the solver, by #link("uncertainty.html#marginal-covariances-by-block-solves")[block solves].

= Accuracy

Accuracy combines trueness and precision. A natural measure is the mean squared distance between the reference and the possible estimates $x tilde cal(N)(hat(u), Sigma)$:

$ bb(E)[norm(u_"ref" - x)_2^2] = norm(u_"ref" - hat(u))_2^2 + tr(Sigma), $

which follows from $x = hat(u) + epsilon$, $bb(E)[epsilon] = 0$, $bb(E)[epsilon^H epsilon] = tr(Sigma)$: the squared bias plus the variance. Normalized by the number of unknowns, it defines the accuracy indicator

$ A = sqrt(1/n (norm(u_"ref" - hat(u))_2^2 + tr(Sigma))), quad A^2 = E^2 + P^2, $ <eq-accuracy>

where $E = norm(u_"ref" - hat(u))_2 slash sqrt(n)$ is the root mean square error of $hat(u)$ and $P$ the precision @eq-precision. $A$ has the unit of the unknown; the smaller, the more accurate. It penalizes both a biased estimate with a small uncertainty and an unbiased estimate with a large one. Without uncertainty ($Sigma = 0$), it reduces to the RMS error $E$.

The decomposition @eq-accuracy also gives a simple check of the uncertainty: for a well calibrated posterior, the actual error is of the order of the predicted spread, $E approx P$. An error much larger than $P$ means that the posterior is overconfident — for instance because the noise covariance `R` is underestimated, or because the regularization introduces a bias that the Gaussian approximation does not account for.

= Comparing methods

Several indicators, with different units, are usually computed for several methods. `performance_ratio` summarizes them in a single number per method, the geometric mean of the ratios to the best value of each indicator,

$ r_m = (product_(j=1)^J I_(m j) / (min_(m') I_(m' j)))^(1 slash J), $ <eq-ratio>

where $I_(m j) > 0$ is the $j$-th indicator of method $m$. The best possible value is $r_m = 1$, for a method that is the best for every indicator. Since each indicator is divided by its best value, $r_m$ does not depend on the units nor on the normalization of the indicators (per unknown or per sample, percent or fraction), and the geometric mean gives the same weight to a relative improvement of any indicator. All the indicators must be such that *smaller is better*: use $100 - C$ for the correlation.

```julia
#| eval: false
sols = [solve(prob, f; precision = true) for f in formulations]
I = [relative_error(uref, s) for s in sols]
C = [100 - correlation(uref, s) for s in sols]
A = [accuracy_metric(uref, s) for s in sols]
performance_ratio([I C A])     # one row per method, one column per indicator
```

The #link("../tutorials/uncertainty.html#performance-indicators")[uncertainty quantification tutorial] applies these indicators to compare several regularizations.

#bibliography("/refs.bib", title: [References], style: "ieee")
