#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Posterior uncertainty])
#metadata((description: "Gaussian approximation of the posterior, marginal covariances by block solves and sampling by perturbation-optimization.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

This page explains how SourceEstimation.jl quantifies the uncertainty of a reconstruction (`solve(...; precision = true)`), and why it works with the posterior *precision* rather than the covariance.

= Gaussian approximation

At the fixed point of the IR algorithm, the posterior distribution of the sources is approximated by a Gaussian distribution of mean $hat(u)$ and precision

$ Q = tau_n A, quad A = H^H H + L^H D L, $ <eq-precision>

where $H$ is the whitened operator (real or complex, $A^H$ being the transpose or the conjugate transpose, so that $A$ is symmetric or Hermitian positive definite), $D = lambda W$ contains the IR weights of the last linear solve and $tau_n$ is the noise precision ($tau_n = 1$ for the additive and multiplicative formulations, $R$ being then taken as the exact covariance). The approximation is exact for quadratic priors ($q = 2$) with fixed hyperparameters.

= Precision rather than covariance

`sol.precision` returns $Q$, not $Q^(-1)$:

- $Q$ is the matrix *already built* by the solver, so it comes at no extra cost;
- it is sparse for assembled space–frequency problems (blocks plus a band), whereas the covariance $Q^(-1)$ is dense: a few megabytes instead of 1.3 GB for 21 points and 451 frequencies;
- it contains all the information, including the correlations between frequencies.

`sol.posterior` stores a factorization of $A$ chosen by dispatch, as for the resolution: dense Cholesky (normal equations), Woodbury factorization (underdetermined regime, factorization of $I + H D^(-1) H^H$ of size $n_m$), sparse Cholesky (CHOLMOD) or one factorization per block. All the posterior quantities are computed on demand by applying $A^(-1)$ to blocks of right-hand sides.

= Marginal covariances by block solves

The diagonal block of sample $j$ (frequency or time step) is obtained by applying $A^(-1)$ to the $n_r$ corresponding columns $E_j$ of the identity:

$ Sigma_(j j) = 1/tau_n (A^(-1) E_j)_([j], :). $

The exact marginal covariances are thus obtained without forming the $n times n$ covariance matrix (memory $n dot n_r$ during the computation, $n_f n_r^2$ for the result). `posterior_std` takes the square root of their diagonal, `posterior_covariance(sol; full = true)` returns the full covariance.

= Sampling by perturbation–optimization

`posterior_samples` draws

$ u_s = hat(u) + 1/sqrt(tau_n) A^(-1) (H^H epsilon_1 + L^H D^(1 slash 2) epsilon_2), quad
  epsilon_1 tilde cal(N)(0, I_m), quad epsilon_2 tilde cal(N)(0, I_p), $ <eq-po>

with circular complex normal distributions for complex data @Orieux2012 @Papandreou2010. Since $epsilon_1$ and $epsilon_2$ are independent,

$ "Cov"(H^H epsilon_1 + L^H D^(1 slash 2) epsilon_2) = H^H H + L^H D L = A
  quad => quad "Cov"(u_s) = 1/tau_n A^(-1) A A^(-1) = Q^(-1). $

The name comes from the interpretation of @eq-po: $u_s$ is the solution of the same weighted Tikhonov problem whose data and prior have been perturbed,

$ u_s = op("argmin")_x norm(y + epsilon_1 slash sqrt(tau_n) - H x)^2 + norm(D^(1 slash 2) L x - epsilon_2 slash sqrt(tau_n))^2. $

== Is it a standard Gaussian sampler?

In distribution, yes: @eq-po draws exactly from $cal(N)(hat(u), Q^(-1))$, like the standard sampler $u_s = hat(u) + R^(-1) z slash sqrt(tau_n)$, $z tilde cal(N)(0, I)$, where $A = R^H R$ is a Cholesky factorization. The difference lies in the square root of $A$ that is used. Writing $A = B^H B$ with the rectangular matrix, known without any computation,

$ B = mat(H; D^(1 slash 2) L), $

and its QR factorization $B = Q_B R$ (where $R$ is the Cholesky factor of $A$),

$ A^(-1) B^H epsilon = R^(-1) R^(-H) R^H Q_B^H epsilon = R^(-1) (Q_B^H epsilon), $

where $Q_B^H epsilon tilde cal(N)(0, I_n)$ since $Q_B$ has orthonormal columns. Perturbation–optimization *is* the standard sampler, without ever computing $R$ or $Q_B$. This justifies its use in the package:

+ *no additional factorization*: applying $A^(-1)$ with the factorization of the solver is enough. In the underdetermined regime, the solver has no Cholesky factor of $A$ ($n_r times n_r$) but a factorization of $I + H D^(-1) H^H$ ($n_m times n_m$): each sample costs $O(n_m n_r)$ instead of an $O(n_r^3)$ factorization;
+ *no covariance*: neither $Q^(-1)$ nor its square root is formed; for assembled problems, the joint samples keep the correlations between frequencies with the sole sparse factorization of $A$;
+ *a single code* for all the regimes and operators (dense, blocks, sparse).

When a Cholesky factor of $A$ is already available (normal equations, CHOLMOD), the standard sampler is used: one triangular solve per sample. Both methods give the same distribution. Variants with approximate optimization corrected by a Metropolis–Hastings step @Gilavert2015, and a unifying review @Vono2022, are available in the literature.

= Credible intervals

`credible_interval(sol; level = 0.95)` computes the empirical quantiles (1 - level)/2 and (1 + level)/2 of posterior samples, component by component (real and imaginary parts separately for complex data).

= Performance indicators

When a reference field is available, the posterior variances also enter the #link("performance.html")[precision and accuracy indicators] (`precision_metric`, `accuracy_metric`), which measure the spread of the posterior and its distance to the reference.

#bibliography("/refs.bib", title: [References], style: "ieee")
