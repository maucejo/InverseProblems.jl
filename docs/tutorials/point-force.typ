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

This tutorial reconstructs a harmonic point force acting on a simply supported beam from noisy acceleration measurements. It compares the classical Tikhonov regularization with sparsity-promoting priors, and shows how the regularization parameter is chosen.

#badges(
  badge[single frequency],
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
ωn, kn = modefreq(beam, 5000.)                        # natural frequencies (rad/s), wavenumbers

xm = grid(40)                                         # sensor positions
xr = grid(31)                                         # reconstruction points
ϕm = modeshape(beam, kn, xm)
ϕr = modeshape(beam, kn, xr)
@printf "%d modes below 5 kHz, first natural frequencies (Hz): %s\n" length(ωn) string(round.(ωn[1:6] / 2π, digits = 1))
```

The FRF matrix is the modal sum

$ H_(i j)(omega) = -omega^2 sum_(n=1)^N (phi_n (x_i) phi_n (x_j)) / (omega_n^2 - omega^2 + 2 i xi omega_n omega), $

computed by `ModalFRFProblem` with a modal damping ratio $xi = 1 %$. We excite the beam between its 5#super[th] and 6#super[th] natural frequencies; `type = :acc` returns accelerances.

```julia
f0 = (ωn[5] + ωn[6]) / 4π                             # excitation frequency (Hz)
@printf "Excitation frequency: %.0f Hz\n" f0

prob_frf = ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕr)
H = SV.solve(prob_frf; type = :acc, progress = false).u[:, :, 1];   # 40 × 31
```

The reference force is a unit point force at the 9#super[th] reconstruction point. The measurements are polluted by a complex Gaussian white noise with a signal-to-noise ratio of 26 dB, added by `agwn`.

```julia
uref = zeros(ComplexF64, length(xr))
uref[9] = 1.0
y = agwn(H * uref, 26)
prob = InverseProblem(H, y)
```

Even with more sensors than unknowns, the problem is ill-posed. The matrix $H$ is a sum of modal contributions, so its rank cannot exceed the number of modes retained (14 here): beyond the 14#super[th] singular value, the singular values drop to the numerical zero, and the smaller ones before it amplify the measurement noise.

```julia
#| fig-alt-text: "Singular values of the FRF matrix, decaying over several orders of magnitude."
s = svdvals(H)
fig = Figure(size = (820, 260))
ax = Axis(fig[1, 1], yscale = log10, xlabel = "Index", ylabel = "Singular value",
          title = @sprintf("Rank: %d, numerical condition number: %.1e", count(>(1e-10 * s[1]), s), s[1] / s[end]))
scatterlines!(ax, s)
fig
```

= Tikhonov regularization

The additive formulation minimizes $norm(y - H u)_2^2 + lambda R(u)$. With `Lq(2.0)`, $R(u) = norm(u)_2^2$: this is the Tikhonov regularization. The regularization parameter $lambda$ is selected by a rule.

```julia
rules = (LCurve(), GCV(), BayesEstimator())
tik = [solve(prob, Additive(Lq(2.0); λ = r)) for r in rules]
[(nameof(typeof(r)), s.λ) for (r, s) in zip(rules, tik)]
```

```julia
#| fig-alt-text: "Tikhonov reconstructions with the L-curve, GCV and Bayesian estimator: the force is spread around its true position."
fig = Figure(size = (820, 300))
ax = Axis(fig[1, 1], xlabel = "Position x (m)", ylabel = "|u| (N)")
scatterlines!(ax, xr, abs.(uref), color = FG, linestyle = :dash, label = "Reference")
for (r, s) in zip(rules, tik)
    scatterlines!(ax, xr, abs.(s.u), label = string(nameof(typeof(r))))
end
axislegend(ax, position = :rt)
fig
```

The force is correctly localized, but its amplitude is underestimated and it leaks onto the neighbouring points: the $ell_2$ norm favors distributed fields.

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
scatterlines!(ax, xr, abs.(uref), color = FG, linestyle = :dash, label = "Reference")
for (name, s) in sols
    scatterlines!(ax, xr, abs.(s.u), label = name)
end
axislegend(ax, position = :rt)
fig
```

The *multiplicative* formulation minimizes $norm(y - H u)_2^2 dot R(u)$: the regularization parameter is updated automatically at each iteration, so that nothing has to be tuned. The *Bayesian* formulation goes further and estimates the noise precision $tau_n$, the prior precision $tau_s$ and the shape parameter $q$ of the prior from the data:

```julia
bayes = last(sols[end])
bayes.hyper
```

The estimated shape parameter is close to the lower bound of its search interval (`shape_bounds = (0.05, 2.0)`): the data are best explained by a very sparse field.

= Summary

```julia
#| echo: false
#| results: "typst"
rows = [(string(nameof(typeof(r)), " (Tikhonov)"), @sprintf("%.1f %%", relerr(s.u, uref)), s.iterations)
        for (r, s) in zip(rules, tik)]
append!(rows, [(name, @sprintf("%.1f %%", relerr(s.u, uref)), s.iterations) for (name, s) in sols])
typst_table(["Method", "Relative error", "Iterations"], rows)
```

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
