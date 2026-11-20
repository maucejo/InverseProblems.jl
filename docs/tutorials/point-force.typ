#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Localizing a point force])
#metadata((description: "Single-frequency force reconstruction on a beam: Tikhonov, sparse and Bayesian regularizations.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

This tutorial reconstructs a harmonic point force acting on a simply supported beam from noisy acceleration measurements in the frequency domain. It compares the classical Tikhonov regularization with sparsity-promoting priors, and shows how the regularization parameter is chosen.

#badges(
  badge[Frequency domain],
  badge(color: julia.red)[Additive],
  badge(color: julia.green)[Multiplicative],
  badge(color: julia.purple)[Bayesian],
)

= The forward model

The synthetic data are generated with #link("https://maucejo.github.io/StructuralVibration.jl/")[StructuralVibration.jl], which provides the beam models, their modal bases, the transfer functions and the measurement noise. Both packages export a function named `solve`: StructuralVibration is therefore imported as `SV`, and its solver is called `SV.solve`.

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn
```

The structure is a simply supported steel beam of length $ell = 1$ m with a rectangular cross-section ($30 times 10$ mm). Its natural frequencies and mass-normalized mode shapes are computed up to 5 kHz, at the 40 sensors and at the 31 reconstruction points:

```julia
L, b, h = 1.0, 0.03, 0.01
beam = Beam(L, b * h, b * h^3 / 12, 2.1e11, 7850.)   # L, S, I, E, ρ
ωn, kn = modefreq(beam, 5000.) # natural frequencies (rad/s), wavenumbers

xm = grid(40) # sensor positions
xu = grid(31) # reconstruction points
ϕm = modeshape(beam, kn, xm)
ϕu = modeshape(beam, kn, xu)
@printf "%d modes below 5 kHz, first natural frequencies (Hz): %s\n" length(ωn) string(round.(ωn[1:6] / 2π, digits = 1))
```

The FRF matrix is the modal sum:

$ H_(i j)(omega) = -omega^2 sum_(n=1)^N (phi_n (x_i) phi_n (x_j)) / (omega_n^2 - omega^2 + 2 i xi omega_n omega), $

computed by `ModalFRFProblem` with a modal damping ratio $xi = 1 %$. We excite the beam between its 5#super[th] and 6#super[th] natural frequencies; `type = :acc` returns accelerances.

```julia
f0 = (ωn[5] + ωn[6]) / 4π  # excitation frequency (Hz)
@printf "Excitation frequency: %.0f Hz\n" f0

prob_frf = ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕu)
H = SV.solve(prob_frf; type = :acc, progress = false).u[:, :, 1]; # 40 × 31
```

The reference force is a unit point force at the 9#super[th] reconstruction point. The measurements are polluted by a complex Gaussian white noise with a signal-to-noise ratio of 26 dB, added by `agwn`.

```julia
uref = zeros(ComplexF64, length(xu))
uref[9] = 1.0
y = agwn(H * uref, 26)
prob = InverseProblem(H, y)
```

= How ill-posed is the problem?

Even with more sensors than unknowns, the problem is ill-posed. `ill_posedness` analyzes the singular value decomposition $bold(H) = sum_i sigma_i bold(v)_i bold(w)_i^H$ of the operator and the coefficients $beta_i = |bold(v)_i^H bold(y)|$ of the data (see the #link("../theory/inverse-problems.html#assessing-the-ill-posedness")[Theory]). With the noise variance, the problem is whitened as in `solve`, so that the noise coefficients have a unit level:

```julia
σ² = noise_variance(H * uref, 26) # variance of the noise added by agwn
ip = ill_posedness(prob; R = σ² * I)
```

The matrix $bold(H)$ is a sum of modal contributions, so its rank cannot exceed the number of modes retained: the numerical rank is 14, and the 17 remaining components of the force (the nullity) have no effect on the measurements at all. Before the numerical zero, the singular values span several orders of magnitude, so that the small ones amplify the noise. The *Picard plot* compares the singular values with the coefficients of the data. It is drawn by `picardplot`, provided by the Makie extension of the package (available as soon as a Makie backend such as CairoMakie is loaded); it accepts an `IllPosedness`, or directly a problem, `picardplot(prob; R = σ² * I)`:

```julia
#| fig-alt-text: "Picard plot: singular values, Fourier coefficients of the data and their ratio, with the noise level and the Picard index."
picardplot(ip)
```

```julia
#| echo: false
@printf "%d components above the noise;\nPicard slope b = %.2f: the coefficients of the data decay as σᵢ^%.2f\n" ip.picard_index ip.picard_slope ip.picard_slope
```

The first nine data coefficients lie above the noise level and therefore contain potentially recoverable information. Beyond the Picard index (dotted vertical line), the coefficients reach a noise-dominated plateau and cannot be reliably distinguished from measurement noise.

Above the noise level, the data coefficients $beta_i$ decay approximately as $sigma_i^1.01$, indicating a Picard slope very close to the limiting value of 1. Consequently, the ratio $beta_i slash sigma_i$, which estimates the magnitude of the force coefficients in the right singular-vector basis, exhibits very weak decay.

This behavior is consistent with the identification of a spatially localized point force, which generally excites a large number of singular vectors. Such a force cannot be accurately represented using only a few dominant singular components. Consequently, regularization methods such as truncated SVD or Tikhonov regularization, which suppress components associated with small singular values, may remove essential spatial information and produce a spatially smeared reconstruction of the point force.

= Tikhonov regularization

The additive formulation minimizes $norm(bold(y) - bold(H space.hair u))_2^2 + lambda space.hair R(bold(u))$. With `Lq(2.0)`, $R(bold(u)) = norm(bold(u))_2^2$, this problem is known as the Tikhonov regularization. Here, the regularization parameter $lambda$ is selected automatically by one of the three rules: the L-curve, generalized cross-validation (GCV) or the Bayesian estimator.

```julia
rules = (LCurve(), GCV(), BayesEstimator())
tik = [solve(prob, Additive(Lq(2.0); λ = r)) for r in rules]
[(nameof(typeof(r)), s.λ) for (r, s) in zip(rules, tik)]
```

```julia
#| fig-alt-text: "Tikhonov reconstructions with the L-curve, GCV and Bayesian estimator: the force is spread around its true position."
fig = Figure(size = (820, 300))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)")
scatterlines!(ax, xu, abs.(uref), color = FG, linestyle = :dash, label = "Reference")
for (r, s) in zip(rules, tik)
    scatterlines!(ax, xu, abs.(s.u), label = string(nameof(typeof(r))))
end
axislegend(ax, position = :rt)
fig
```

The force is correctly localized, but its amplitude is underestimated and it leaks onto the neighboring points: the $ell_2$ norm favors distributed fields.

= Sparsity-promoting priors

A point force is a *sparse* field. The $ell_q$ (quasi-)norm with $q <= 1$ promotes such solutions. Because the regularization term is no longer quadratic, the problem is solved by the iteratively reweighted (IR) algorithm.

```julia
sparse_cases = [
    "Additive ℓ1 (GCV)"  => Additive(Lq(1.0); λ = GCV()),
    "Multiplicative ℓ1"  => Multiplicative(Lq(1.0)),
    "Multiplicative ℓ0.5" => Multiplicative(Lq(0.5)),
    "Bayesian (q estimated)" => Bayesian(Lq(1.0)),
]
sols = [name => solve(prob, f) for (name, f) in sparse_cases];
```

```julia
#| fig-alt-text: "Sparse reconstructions: all of them recover the point force; the Bayesian and l0.5 solutions are almost exact."
fig = Figure(size = (820, 340))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)")
scatterlines!(ax, xu, abs.(uref), color = FG, linestyle = :dash, label = "Reference")
for (name, s) in sols
    scatterlines!(ax, xu, abs.(s.u), label = name)
end
axislegend(ax, position = :rt)
fig
```

The *multiplicative* formulation minimizes $norm(bold(y) - bold(H space.hair u))_2^2 dot R(bold(u))$. Here, the regularization parameter is updated automatically at each iteration, so that nothing has to be tuned. The *Bayesian* formulation goes further and estimates the noise precision $tau_n$, the prior precision $tau_s$ and the shape parameter $q$ of the prior from the data:

```julia
bayes = last(sols[end])
bayes.hyper
```

The estimated shape parameter is close to the lower bound of its search interval (`shape_bounds = (0.01, 2.0)`): the data are best explained by a very sparse field.

= Summary

```julia
#| echo: false
#| results: "typst"
indicators(s) = (@sprintf("%.1f %%", relative_error(uref, s)), @sprintf("%.1f %%", correlation(uref, s)), s.iterations)
rows = [(string(nameof(typeof(r)), " (Tikhonov)"), indicators(s)...) for (r, s) in zip(rules, tik)]
append!(rows, [(name, indicators(s)...) for (name, s) in sols])
typst_table(["Method", "Relative error", "Correlation", "Iterations"], rows)
```

The #link("../theory/performance.html#relative-error")[relative error] (`relative_error`, $ell_1$ norm) measures the overall error, while the #link("../theory/performance.html#correlation")[correlation] (`correlation`) only compares the *shapes* of the fields, whatever their amplitude. The Tikhonov solutions spread the force over the whole beam: the leakage on all the points gives a relative error larger than 100 % and a poor correlation. The multiplicative $ell_1$ solution has the right shape (correlation of 99 %) but an underestimated amplitude, hence its relative error of almost 50 %; the Bayesian solution is right on both counts.

#calepin.elements.callout(kind: "tip")[
The choice between formulations is mostly a matter of prior knowledge:

- use `Additive` when you want to control $lambda$ (or choose it with a well-known rule);
- use `Multiplicative` for a parameter-free reconstruction with a fixed prior shape;
- use `Bayesian` when the shape of the prior is unknown, or when you need a consistent estimate of the noise level and posterior uncertainties.
]

#next-steps(
  card(icon: [🧩], title: [Group regularization], href: "groups.html", color: julia.green)[Localized and distributed sources on the same structure.],
  card(icon: [📚], title: [Regularization theory], href: "../theory/regularization.html", color: julia.blue)[Additive and multiplicative formulations, parameter choice rules.],
)
