#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Linear inverse problems])
#metadata((description: "Forward models of source identification problems and why their inversion is ill-posed.")) <website-metadata>

#set math.equation(numbering: "(1)")

#title()

This page introduces the class of problems addressed by SourceEstimation.jl and explains why a naive inversion fails. The next pages describe the #link("regularization.html")[regularization strategies], their #link("bayesian.html")[Bayesian interpretation], the #link("algorithms.html")[algorithms] and the #link("uncertainty.html")[posterior uncertainty].

= Real or complex data

Depending on the problem, the data and the unknowns are *real* (time domain, static or real modal problems) or *complex* (frequency domain). SourceEstimation.jl handles both with the same code: the field is deduced from the element types of `H` and `y` (`Float64` or `ComplexF64`), and the formulas that depend on it — essentially the probability distributions of the #link("bayesian.html")[Bayesian model] — are selected by dispatch. All the pages of this section use the following conventions, valid in both cases.

#table(
  columns: (auto, 1fr, 1fr),
  align: left,
  table.header([*Notation*], [*Real data*], [*Complex data*]),
  [$bb(K)$], [$bb(R)$], [$bb(C)$],
  [$x^*$], [$x$], [complex conjugate],
  [$A^H$], [transpose $A^T$], [conjugate transpose $(A^*)^T$],
  [$|x|$], [absolute value], [modulus],
  [$norm(x)_2^2 = x^H x$], [$sum_i x_i^2$], [$sum_i |x_i|^2$],
  [$partial f slash partial x^*$], [$1/2 partial f slash partial x$], [Wirtinger derivative $1/2 (partial f slash partial a + i partial f slash partial b)$, $x = a + i b$],
  [$cal(N)(mu, Sigma)$], [real Gaussian], [circular complex Gaussian],
)

With these conventions, the normal equations, the IR weights and the posterior precision have the same expression for real and complex data; only the constants of the Bayesian model differ.

= Forward models

Source identification consists in recovering the excitation $u in bb(K)^(n)$ of a system from measurements $y in bb(K)^(m)$ of its response. For a linear system, the measurements are related to the sources by a linear forward operator $H in bb(K)^(m times n)$,

$ y = H u + n, $ <eq-model>

where $n in bb(K)^m$ is the measurement noise. The structure of $H$ depends on the domain.

== Frequency domain

At a circular frequency $omega$, the $n_m$ measured responses and the $n_r$ unknown forces are related by the matrix of frequency response functions (FRFs) $H(omega) in bb(C)^(n_m times n_r)$, obtained experimentally or from a numerical model. For a modal model,

$ H_(i j)(omega) = sum_n (phi_n (x_i) phi_n (x_j)) / (omega_n^2 - omega^2 + 2 i xi_n omega_n omega). $

Over a frequency band $omega_1, dots, omega_(n_f)$, the unknowns form the matrix $U = [u(omega_1), dots, u(omega_(n_f))] in bb(C)^(n_r times n_f)$ and the forward operator is block-diagonal,

$ "vec"(Y) = "diag"(H(omega_1), dots, H(omega_(n_f))) "vec"(U) + "vec"(N). $

Each frequency can be solved independently, unless the prior couples them (see #link("regularization.html#mixed-norms")[mixed norms]) @Aucejo2018b.

== Time domain

In the time domain, the response is the convolution of the forces with the impulse responses $h_(i j)(t)$. With a sampling period $Delta t$, the discretized convolution reads

$ y(t_k) = sum_(l=1)^k h(t_k - t_l) u(t_l) Delta t, $

which defines a block lower-triangular Toeplitz matrix $H in bb(R)^(n_m n_t times n_r n_t)$ acting on $"vec"(U)$, $U in bb(R)^(n_r times n_t)$ @Aucejo2019b.

== A common layout

In both domains, the unknown is a matrix $U in bb(K)^(n_r times n_s)$ whose *rows* are the reconstruction points (space) and whose *columns* are the samples (frequencies or time steps); a single-frequency problem has $n_s = 1$. Global vectors are stacked sample by sample, $u = "vec"(U)$, i.e. $u_(i + (j-1) n_r) = U_(i j)$. This convention gives the same meaning to groups, smoothing operators and mixed norms in all domains.

= Ill-posedness

== Singular value decomposition

Let $H = sum_i sigma_i v_i w_i^H$ be the singular value decomposition of $H$ (real or complex singular vectors $v_i$, $w_i$, real singular values), with $sigma_1 >= sigma_2 >= dots >= sigma_r > 0$. The least-squares solution of @eq-model is

$ hat(u)_"LS" = sum_(i=1)^r (v_i^H y) / sigma_i w_i = u + sum_(i=1)^r (v_i^H n) / sigma_i w_i. $ <eq-ls>

Source identification problems are *discrete ill-posed problems* @Hansen2010: the singular values decay gradually towards zero, so that the condition number $sigma_1 slash sigma_r$ is very large. The noise components $v_i^H n$ do not decay, and the second sum of @eq-ls dominates the solution: the least-squares estimate is useless, even with a very low noise level.

== The discrete Picard condition

A meaningful solution can only be recovered if the Fourier coefficients of the exact data $|v_i^H H u|$ decay faster than the singular values — the discrete Picard condition. Above the noise level, the coefficients $|v_i^H y|$ level off and the corresponding components must be filtered out. This is the role of *regularization*: it replaces @eq-ls by a filtered expansion

$ hat(u) = sum_(i=1)^r f_i (v_i^H y) / sigma_i w_i, $

where the filter factors $f_i in [0, 1]$ damp the components dominated by noise. The Tikhonov regularization, for instance, uses $f_i = sigma_i^2 slash (sigma_i^2 + lambda)$.

= Noise

The noise is assumed Gaussian and centered — real Gaussian for real data, circular complex Gaussian for complex data — with a known *diagonal* covariance $R$ (real, positive) (up to a scale factor, see the #link("bayesian.html")[Bayesian formulation]). The data fidelity is then the weighted least-squares term

$ norm(y - H u)_(R^(-1))^2 = (y - H u)^H R^(-1) (y - H u). $

It is handled once and for all by *whitening*: $H <- R^(-1 slash 2) H$, $y <- R^(-1 slash 2) y$, after which the noise covariance is the identity. In SourceEstimation.jl, $R$ is given by the keyword `R` of `solve`: `σ² * I`, the variances of all the data, or one variance per measurement point shared by all the samples.

#bibliography("/refs.bib", title: [References], style: "ieee")
