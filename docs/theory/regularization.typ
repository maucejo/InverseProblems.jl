#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Regularization])
#metadata((description: "Regularization terms, additive and multiplicative formulations, and regularization parameter choice rules.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

Regularization restores the well-posedness of the #link("inverse-problems.html")[inverse problem] by adding prior knowledge on the sources. In SourceEstimation.jl, a regularization is made of a *regularization term* $R(u)$ — the prior — and a *formulation* combining it with the data fidelity.

= Regularization terms

== $ell_q$ norms

The basic term is the $ell_q$ (quasi-)norm of a linear transform of the unknown,

$ R(u) = norm(L u)_q^q = sum_i |(L u)_i|^q, quad 0 < q <= 2, $

implemented by `Lq(q; L = I)`. The shape parameter $q$ encodes the nature of the sources @Aucejo2014:

- $q = 2$ (Tikhonov) promotes *distributed* fields, with energy spread over many points;
- $q = 1$ (LASSO) promotes *sparse* fields, i.e. localized sources;
- $q < 1$ promotes even sparser fields, at the cost of a non-convex problem.

The operator $L$ is the identity by default. Finite-difference operators (`diff_operator`) penalize the derivatives of the field, i.e. promote smooth distributions. For a matrix unknown, $L$ acts on the rows (space): $L U$.

== Mixed norms <mixed-norms>

For space–frequency or space–time problems, the mixed norm $ell_(p,q)$ (`Lpq(p, q; L = I, Lc = nothing)`) exploits the joint structure of the unknown @Aucejo2018b @Aucejo2020b:

$ R(U) = norm(L U L_c^T)_(p,q)^q = sum_i (sum_j |(L U L_c^T)_(i j)|^p)^(q slash p). $

The inner exponent $p$ acts on each row, i.e. on the spectrum (or time history) of each point, and the outer exponent $q$ on the spatial distribution. `Lpq(2, 0.5)` describes broadband ($p = 2$) and localized ($q = 0.5$) sources, `Lpq(1, 2)` narrow-band distributed ones. The operator $L_c$ acts on the columns, e.g. along the frequencies to promote continuous spectra. On $"vec"(U)$, the global operator is a Kronecker product,

$ (L_c times.o L) "vec"(U) = "vec"(L U L_c^T), $

which is the operator of MSSP 2018 (Eq. 15) with $L = L_s$ and $L_c = L_f$ @Aucejo2018b.

== Groups

When the structure contains regions of different natures, `Groups(rows₁ => term₁, …)` defines a region-wise regularization term @Aucejo2016 @Aucejo2020a,

$ R(u) = sum_g R_g (L_g u_g), $

where $u_g$ gathers the rows of the group $g$ and each group has its own term (`Lq` or `Lpq`) and operator $L_g$. The groups must form a partition of the rows. The operators $L_g$ are assembled once into a global sparse operator; when they are all the identity, no matrix is created.

= Additive formulation

The additive formulation (`Additive(prior; λ)`) is the generalized Tikhonov regularization @Aucejo2014:

$ hat(u) = op("argmin")_u norm(y - H u)_2^2 + sum_g lambda_g R_g (L_g u_g). $ <eq-additive>

With a single group, $lambda$ is a scalar; with several groups and a vector $lambda$, @eq-additive is a *multi-parameter* regularization. The regularization parameter is either fixed or selected by a rule at each iteration of the #link("algorithms.html")[IR algorithm]. For quadratic terms only, @eq-additive has the closed-form solution

$ hat(u) = (H^H H + L^H Lambda L)^(-1) H^H y, $

with $Lambda$ the diagonal matrix of the $lambda_g$: a single linear solve is performed.

= Parameter choice rules

All the rules work on the spectral description of the weighted quadratic problem $min_u norm(y - H u)^2 + lambda norm(W^(1 slash 2) L u)^2$ solved at each iteration: SVD of $H W^(-1 slash 2)$ if $L = I$, generalized SVD (GSVD) of the pair $(H, W^(1 slash 2) L)$ otherwise. Let $gamma_i$ denote the (generalized) singular values and $beta_i$ the coefficients of the data in the corresponding basis. The criteria are minimized in $t = log ell$, with $lambda = ell^2$, by a grid search refined by Brent's method.

/ L-curve (`LCurve()`): the parametric curve $(log norm(y - H u_lambda), log norm(L u_lambda))$ has the shape of an "L"; its corner, the point of maximum curvature, balances the residual and the solution norms @Hansen1993.

/ Generalized cross-validation (`GCV()`): minimizes the predictive error estimate @Golub1979
  $ G(lambda) = norm(y - H u_lambda)^2 / (op("tr")(I - H H_lambda^sharp))^2, $
  where $H_lambda^sharp$ is the regularized inverse.

/ Bayesian estimator (`BayesEstimator()`): maximizes the marginal posterior of the regularization parameter, i.e. minimizes @Pereira2015
  $ J_"MAP" (ell) = sum_(i=1)^M log(gamma_i^2 + ell^2) + (M - 2) log(1/M sum_(i=1)^M (|beta_i|^2) / (gamma_i^2 + ell^2)). $

Two rules are only used to *initialize* the iterative algorithms (keyword `init` of `solve`):

/ Heuristic (`Heuristic()`): the geometric mean of the extreme eigenvalues of $(H L^(-1))^H (H L^(-1))$, i.e. of the extreme admissible values of the parameter @Aucejo2017 @Aucejo2018b. For sparse operators, the eigenvalues are estimated by power iterations without assembling the matrix.

/ $lambda_"max"$ (`LambdaMax()`): $lambda = norm(2 H^H y)_infinity$, the smallest value for which the $ell_1$-regularized solution is zero @Kim2007.

= Multiplicative formulations

Choosing the regularization parameter is the main difficulty of the additive formulation. The *multiplicative* regularization removes it by multiplying the data fidelity and the regularization term @Aucejo2017 @Aucejo2018a:

$ hat(u) = op("argmin")_u J(u) = norm(y - H u)_2^2 dot R(u). $ <eq-omr>

Writing the first-order optimality condition $nabla_(u^*) J = 0$ with $nabla_(u^*) R = L^H W L u$ (see the #link("algorithms.html")[IR weights]) gives

$ (H^H H + alpha L^H W L) u = H^H y, quad alpha = norm(y - H u)_2^2 / R(u). $

The problem is thus equivalent to an additive regularization whose parameter $alpha$ is *adaptive*: it is updated at each iteration from the current residual and penalty, without any user input. This is the *ordinary multiplicative regularization* (OMR), `Multiplicative(prior)`, with $R = sum_g R_g$ for groups. It extends to space–frequency mixed norms @Aucejo2018b @Aucejo2020b and to the space–time domain @Aucejo2019b.

== Multi-parameter multiplicative regularization

With groups, the *multi-parameter multiplicative regularization* (MPMR), `Multiplicative(Groups(…); multiparameter = true)`, multiplies the group penalties @Aucejo2020a:

$ hat(u) = op("argmin")_u norm(y - H u)_2^2 dot product_g R_g (L_g u_g). $

The optimality condition now yields one adaptive parameter per group,

$ alpha_g = norm(y - H u)_2^2 / (R_g (L_g u_g)), $

so that each region is regularized according to its own penalty, without having to tune several parameters by hand.

#calepin.elements.callout(kind: "tip")[
The multiplicative formulations require no parameter choice rule, but they still benefit from a good initial solution. By default, the IR algorithm starts from a Tikhonov solution whose parameter is given by the `init` keyword of `solve`.
]

#bibliography("/refs.bib", title: [References], style: "ieee")
