#import "/.calepin/calepin.typ" as calepin
#import "/_lib/ui.typ": *
#show: math-html

#set document(title: [Time-domain reconstruction])
#metadata((description: "Space-time force reconstruction from a block convolution matrix with mixed-norm regularization.")) <website-metadata>

#calepin.setup(fig-width: 100%)

```julia
#| echo: false
#| results: "hide"
#| message: false
#| warning: false
include(isfile("_lib/setup.jl") ? "_lib/setup.jl" : "../_lib/setup.jl")
```

#title()

Transient excitations (impacts, shocks, transient loads) are naturally described in the time domain. The response of a linear structure is then a convolution of the forces with the impulse responses, which, once discretized, is again a linear problem $bold(y) = bold(H space.hair u) + bold(n)$.

#badges(
  badge[Space–time],
  badge(color: julia.red)[MIMO convolution],
  badge(color: julia.purple)[Lpq],
)

= The convolution model

With $n_t$ time steps $t_k = (k-1) Delta t$, the displacement measured by sensor $i$ reads:

$ y_i (t_k) = sum_(l=1)^k sum_j h_(i j)(t_k - t_l) u_j (t_l) Delta t, $

where $h_(i j)$ is the impulse response between point $j$ and sensor $i$. Stacking the time samples gives a block lower-triangular (Toeplitz) matrix $bold(H)$ of size $n_m n_t times n_u n_t$, acting on $"vec"(U)$ where $U$ is the $n_r times n_t$ matrix of the forces. For a modal model, $h_(i j)(t) = sum_n phi_n (x_i) phi_n (x_j) h_n (t)$ where $h_n (t)$ is the impulse response of the $n$-th single-degree-of-freedom mode. The modal basis of the steel beam and the impulse responses $h_n$ come from StructuralVibration.jl (`impulse_response` of a `Sdof`); the helper `convolution_matrix` of the documentation (see `_lib/setup.jl`) assembles them into the matrix $bold(H)$.
```julia
import StructuralVibration as SV
using StructuralVibration: Beam, modefreq, modeshape, agwn

beam = Beam(1.0, 3e-4, 2.5e-9, 2.1e11, 7850.)
t = range(0, 0.06, length = 150)                 # 60 ms, Δt ≈ 0.40 ms
ωn, kn = modefreq(beam, 1 / (2.5 * step(t)))     # modes below fs / 2.5
xm, xu = grid(12), grid(9)
ϕm, ϕu = modeshape(beam, kn, xm), modeshape(beam, kn, xu)

H = convolution_matrix(ωn, 0.01, ϕm, ϕu, t);
@printf "%d modes, H is %d × %d\n" length(ωn) size(H)...
```

The reference excitation combines an impact-like pulse at point 3 and a 200 Hz burst at point 7 (displacement measurements, signal-to-noise ratio of 25 dB):

```julia
Uref = zeros(length(xu), length(t))
Uref[3, :] .= @. exp(-((t - 0.01) / 0.0015)^2)
Uref[7, :] .= @. 0.6 * sin(2π * 200 * t) * (0.02 < t < 0.04)

Y0 = reshape(H * vec(Uref), length(xm), :) # nm × nt measurements
Y = agwn(Y0, 25)
prob = InverseProblem(H, Y)
```

The measurements are simulated with the same discrete model as the one used for the inversion (#emoji.warning _inverse crime_); the impulse responses are however those of the continuous structure, sampled at the time steps. Since `Y` is a matrix and `H` a dense matrix with `size(H, 1) == length(Y)`, the problem is interpreted as a space–time problem: the unknown is a `9 × 150` matrix.

= Reconstructions

The Tikhonov solution is compared with the multiplicative regularization using an $ell_1$ norm and two mixed norms. In the time domain, `Lpq(p, q)` couples the time samples of each point: the sources are localized in space ($q <= 1$) and the inner $ell_2$ norm does not constrain their time history.

```julia
cases = [
    "Tikhonov (GCV)"        => Additive(Lq(2.0); λ = GCV()),
    "Multiplicative ℓ1"     => Multiplicative(Lq(1.0)),
    "Multiplicative ℓ2,1"   => Multiplicative(Lpq(2, 1)),
    "Multiplicative ℓ2,0.5" => Multiplicative(Lpq(2, 0.5)),
]
sols = [name => solve(prob, f) for (name, f) in cases];
```

```julia
#| fig-alt-text: "Reconstructed time histories at the two loaded points, compared with the reference."
fig = Figure(size = (840, 520))
for (k, i) in enumerate((3, 7))
    ax = Axis(fig[k, 1], title = @sprintf("Point %d (x = %.2f m)", i, xu[i]),
              xlabel = k == 2 ? "Time (ms)" : "", ylabel = "Force (N)")
    lines!(ax, 1e3 .* t, Uref[i, :], color = FG, linestyle = :dash, linewidth = 2.5, label = "Reference")
    for (name, s) in sols[[1, 4]]
        lines!(ax, 1e3 .* t, s.u[i, :], label = name)
    end
    k == 1 && axislegend(ax, position = :rt)
end
fig
```

```julia
#| fig-alt-text: "Space-time maps of the force: reference, Tikhonov and l2,0.5 reconstructions."
fig = Figure(size = (840, 280))
maps = ["Reference" => Uref, "Tikhonov (GCV)" => last(sols[1]).u, "Multiplicative ℓ2,0.5" => last(sols[4]).u]
for (k, (name, U)) in enumerate(maps)
    ax = Axis(fig[1, k], title = name, xlabel = "Time (ms)", ylabel = k == 1 ? "x (m)" : "")
    heatmap!(ax, 1e3 .* t, xu, transpose(U), colormap = :balance, colorrange = (-1, 1))
end
fig
```

```julia
#| echo: false
#| results: "typst"
typst_table(["Method", "Relative error", "Correlation", "Iterations"],
            [(name, @sprintf("%.1f %%", relative_error(Uref, s)), @sprintf("%.1f %%", correlation(Uref, s)), s.iterations)
             for (name, s) in sols])
```

The quadratic prior spreads the energy of the forces over all the points, while the mixed norm $ell_(2,0.5)$ concentrates it on the two loaded points.

#calepin.elements.callout(kind: "note")[
A space–time problem is solved as a whole, with a dense system of size $n_r n_t$. Its cost grows quickly with the number of time steps: for long signals, consider a frequency-domain formulation, or a sparse (assembled) operator if the impulse responses are short.
]

#next-steps(
  card(icon: [📊], title: [Uncertainty quantification], href: "uncertainty.html", color: julia.red)[Posterior standard deviations and credible intervals.],
  card(icon: [📚], title: [Inverse problems], href: "../theory/inverse-problems.html", color: julia.blue)[Forward models and ill-posedness.],
)
