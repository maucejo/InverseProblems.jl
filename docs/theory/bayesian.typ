#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Bayesian regularization])
#metadata((description: "Hierarchical Bayesian model, optimal Bayesian regularization and relevance vector regression.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

The Bayesian framework gives a probabilistic reading of regularization: the regularization term is the negative logarithm of a *prior* distribution on the sources, and the regularized solution is the *maximum a posteriori* (MAP) estimate. Going one step further, the hyperparameters — noise level, prior precision, prior shape — can be treated as random variables and estimated from the data @Aucejo2016 @Aucejo2018c @Aucejo2019a.

= Hierarchical model

The formulation `Bayesian(prior; τn, τs, q)` is built on the following hierarchical model.

/ Likelihood: the noise is Gaussian with covariance $R slash tau_n$, where $R$ is the (whitening) covariance given to `solve` and $tau_n$ the *noise precision*.

/ Prior: for each group $g$, the components of $v_g = L_g u_g$ follow a *generalized Gaussian* distribution of shape $q_g$ and precision $tau_(s,g)$. The shape $q_g = 2$ gives a Gaussian prior, $q_g = 1$ a Laplace prior.

/ Hyperpriors: weakly informative Gamma distributions on the precisions and the shapes.

The joint posterior reads

$ p(u, tau_n, tau_s, q | y) prop p(y | u, tau_n) product_g p(L_g u_g | tau_(s,g), q_g) p(tau_(s,g)) p(q_g) p(tau_n). $

Real and complex (circular) data lead to different distributions:

#table(
  columns: (auto, 1fr, 1fr),
  align: left,
  table.header([], [*Real data*], [*Complex data (circular)*]),
  [Likelihood], [$(tau_n slash 2 pi)^(N slash 2) exp(-tau_n norm(r)^2 slash 2)$], [$(tau_n slash pi)^N exp(-tau_n norm(r)^2)$],
  [Prior (per component)], [$prop exp(-(tau_s slash q) |v|^q)$], [$prop exp(-tau_s |v|^q)$],
  [MAP in $u$], [$1/2 norm(y - H u)_R^2 + sum_g lambda_g / q_g norm(L_g u_g)_q^q$], [$norm(y - H u)_R^2 + sum_g lambda_g norm(L_g u_g)_q^q$],
)

with $r = y - H u$ and $lambda_g = tau_(s,g) slash tau_n$. For a fixed set of hyperparameters, the MAP estimate is thus an #link("regularization.html")[additive regularization] whose parameters are the ratios of the prior and noise precisions.

= Optimal Bayesian regularization

The default algorithm, `MAP()`, maximizes the joint posterior by alternating closed-form updates of the hyperparameters and an IR step on $u$ @Aucejo2019a. With $M_g$ components in group $g$ and Gamma hyperpriors of parameters $(a, b)$, the updates read, for real data,

$ tau_n = (N slash 2 + a - 1) / (norm(r)^2 slash 2 + b), quad
  tau_(s,g) = (M_g slash q_g + a - 1) / (sum_i |v_(g,i)|^(q_g) slash q_g + b), $

and for complex data

$ tau_n = (N + a - 1) / (norm(r)^2 + b), quad
  tau_(s,g) = (2 M_g slash q_g + a - 1) / (sum_i |v_(g,i)|^(q_g) + b). $

The shape parameters $q_g$ are obtained by a one-dimensional minimization of the negative log-posterior on the interval `shape_bounds` (by default $[0.05, 2]$). In the real case, the objective is

$ M log Gamma(1/q) - M/q log(tau_s / q) - (a + M - 1) log q + b q + tau_s / q sum_i |v_i|^q. $

Estimating $q$ lets the data decide between a sparse and a distributed reconstruction, which makes the method fully automatic. The weakly informative hyperpriors use $a = 1$ and a small $b$.

== Estimated or fixed hyperparameters

Each hyperparameter is either estimated (`Estimated()`) or fixed:

```julia
#| eval: false
Bayesian(prior)                                     # τs, q estimated; τn = 1 if R is given, estimated otherwise
Bayesian(prior; τn = Estimated())                   # R known up to a scale factor
Bayesian(Groups(...); τs = [τ₁, τ₂], q = Fixed())   # τs fixed per group, q from the prior
Bayesian(Lq(2.0); τn = 1.0, τs = 0.5, q = Fixed())  # all fixed: exact Gaussian posterior
```

With `R` the exact noise covariance, $tau_n = 1$; when only the relative variances are known, $tau_n$ is estimated as a global scale factor.

= Relevance vector regression

The element-wise algorithm, `solve(prob, Bayesian(Lq(2.0)), ElementWise())`, assigns *one Gaussian prior precision per component* of $v = L u$, in the spirit of sparse Bayesian learning and the relevance vector machine @Tipping2001:

$ p(v_i | tau_i) prop exp(-phi tau_i |v_i|^2), quad tau_i tilde "Gamma"(a, b), $

with $phi = 1 slash 2$ for real data and $phi = 1$ for complex data. The MAP alternates

$ tau_i = (k slash 2 + a - 1) / (phi |v_i|^2 + b), quad
  u = (H^H H + L^H Lambda L)^(-1) H^H y, quad Lambda = tau / tau_n, $

(and the update of $tau_n$ if estimated), where $k = 1$ for real and $k = 2$ for complex data. The precisions of the irrelevant components grow without bound, which drives them to zero: the prior is Gaussian component-wise, but the resulting estimate is sparse. The algorithm is initialized with $lambda = norm(2 H^H y)_infinity$ @Kim2007.

#bibliography("/refs.bib", title: [References], style: "ieee")
