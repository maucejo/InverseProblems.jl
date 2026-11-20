#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Group regularization])
#metadata((description: "Region-wise priors and multi-parameter regularization for mixed localized and distributed sources.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

Real structures are often excited by sources of different natures: a localized force from a shaker or a mount, together with a distributed pressure from a flow or an acoustic field. A single prior cannot describe both. `Groups` assigns its own regularization term to each region of the structure.

#badges(
  badge(color: julia.blue)[Groups],
  badge(color: julia.green)[multi-parameter regularization],
  badge(color: julia.purple)[Bayesian],
  badge(color: julia.red)[diff_operator],
)

= A mixed excitation

We keep the steel beam of the #link("point-force.html")[previous tutorial] (data generated with `StructuralVibration.jl`) and apply both a point force at the 6#super[th] reconstruction point and a half-sine distributed load on the right part of the beam.

```julia
import StructuralVibration as SV
using StructuralVibration: Beam, ModalFRFProblem, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)
ωn, kn = modefreq(beam, 5000.)
xm, xu = grid(40), grid(31)
ϕm, ϕu = modeshape(beam, kn, xm), modeshape(beam, kn, xu)
f0 = (ωn[5] + ωn[6]) / 4π
H = SV.solve(ModalFRFProblem(ωn, 0.01, f0:f0, ϕm, ϕu); type = :acc, progress = false).u[:, :, 1];

region = 17:29 # rows of the distributed load
others = setdiff(1:31, region) # rest of the beam
uref = zeros(eltype(H), 31)
uref[6] = 1.0
uref[region] .= 0.3 .* sin.(π .* (1:13) ./ 14)

y = agwn(H * uref, 30)
prob = InverseProblem(H, y)
```

= One prior for all: a compromise

Neither the $ell_2$ norm nor the $ell_1$ norm is appropriate everywhere: the former blurs the point force, the latter breaks the distributed load into spikes.

```julia
tik = solve(prob, Additive(Lq(2.0); λ = GCV()))
l1  = solve(prob, Multiplicative(Lq(1.0)));
```

= Region-wise priors

`Groups(rows₁ => term₁, rows₂ => term₂, ...)` defines the regularization term $R(bold(u)) = sum_g R_g (bold(u)_g)$, where the groups must form a partition of the rows of the unknown. Here, a sparse prior on the left part and a quadratic prior on the distributed load:

```julia
G = Groups(others => Lq(1.0), region => Lq(2.0));
```

With the ordinary multiplicative formulation (OMR), a single regularization parameter weights the sum of the group penalties. The *multi-parameter multiplicative regularization* (MPMR) minimizes $norm(bold(y) - bold(H space.hair u))^2 dot.c product_g R_g (bold(u)_g)$ instead, which yields one adaptive regularization parameter per group:

```julia
omr  = solve(prob, Multiplicative(G));
mpmr = solve(prob, Multiplicative(G; multiparameter = true))
mpmr.λ
```

The Bayesian formulation also accepts groups: each group has its own prior precision $tau_(s,g)$ and, if estimated, its own shape parameter $q_g$, initialized by the term of the group.

```julia
bay = solve(prob, Bayesian(G))
bay.hyper
```

```julia
#| fig-alt-text: "Reconstructions of the mixed excitation: the group priors recover both the point force and the distributed load."
fig = Figure(size = (820, 560))
cases = ["Tikhonov (GCV)" => tik, "Multiplicative ℓ1" => l1,
         "Groups, OMR" => omr, "Groups, MPMR" => mpmr, "Groups, Bayesian" => bay]
for (k, (name, s)) in enumerate(cases)
    ax = Axis(fig[fld1(k, 2), mod1(k, 2)], title = name,
              xlabel = "x (m)", ylabel = "|u| (N)", limits = (nothing, (-0.05, 1.15)))
    vspan!(ax, xu[first(region)], xu[last(region)], color = (JULIA_COLORS.green, 0.08))
    lines!(ax, xu, abs.(uref), color = FG, linestyle = :dash)
    scatterlines!(ax, xu, abs.(s.u), color = JULIA_COLORS.blue)
    text!(ax, 0.62, 0.8, text = @sprintf("error: %.0f %%", relative_error(uref, s)),
          space = :relative, fontsize = 12, color = FG)
end
fig
```

```julia
#| echo: false
#| results: "typst"
typst_table(["Method", "Regularization parameter(s)", "Relative error", "Correlation"],
            [(name, s.λ isa AbstractVector ? join((@sprintf("%.3g", v) for v in s.λ), ", ") : @sprintf("%.3g", s.λ),
              @sprintf("%.1f %%", relative_error(uref, s)), @sprintf("%.1f %%", correlation(uref, s))) for (name, s) in cases])
```

#calepin.elements.callout(kind: "note")[
With the additive formulation, a vector of regularization parameters, one per group, gives the *multi-parameter* generalized Tikhonov regularization: `Additive(G; λ = [λ₁, λ₂])`. Choosing several parameters by hand is tedious, which is precisely what the MPMR avoids.
]

= Smoothing operators

Each term accepts an operator `L` acting on the rows of its group, so that the penalty becomes $R_g (bold(L)_g bold(u)_g)$. Finite-difference operators are built with `diff_operator`, on a 1D uniform grid, #sym.dots

```julia
D2 = diff_operator(5, 2; h = 0.1) # second derivative, 3-point scheme, interior points
Matrix(D2)
```

… or along a direction of a structured 2D/3D point cloud, given in any order (here a $7 times 5$ grid of a plate):

```julia
coords = [(x, y) for x in range(0, 0.6, length = 7) for y in range(0, 0.4, length = 5)]
C = [first.(coords) last.(coords)]
Dx = diff_operator(C, 1)            # ∂/∂x along the grid lines
Dy = diff_operator(C, 2; order = 2) # ∂²/∂y²
size(Dx), size(Dy)
```

They are then passed to the terms, for instance to promote a smooth distributed load next to localized sources, or to penalize the gradient of the field:

```julia
#| eval: false
Groups(inner => Lq(0.5), edges => Lq(2.0; L = diff_operator(coords[edges, :], 1)))
Lq(2.0; L = [diff_operator(C, 1); diff_operator(C, 2)]) # gradient (stacked operators)
```

#calepin.elements.callout(kind: "warning")[
The components of the null space of `L` (e.g. constant and linear fields for a second derivative without boundary rows) are not regularized. When this null space is large compared to the information brought by the measurements, prefer an identity or a square, invertible operator (centered scheme `npts = 5`).
]

#next-steps(
  card(icon: [📈], title: [Space–frequency problems], href: "space-frequency.html", color: julia.purple)[Couple the frequencies with mixed norms.],
  card(icon: [📚], title: [Priors and formulations], href: "../theory/regularization.html", color: julia.blue)[The theory behind `Groups` and the MPMR.],
)
