# SourceEstimation.jl documentation

Website built with [Calepin](https://vincentarelbundock.github.io/calepin) (Typst + executed
code chunks). Every Julia block of the pages is executed at build time through a Jupyter
(IJulia) kernel.

## Layout

| Path | Content |
|:--|:--|
| `calepin.toml` | site configuration (navigation, sidebar, search, theme) |
| `index.typ` | landing page (features, quick example, `versioninfo()`, dependency versions) |
| `getting-started.typ` | installation and first reconstruction |
| `tutorials/` | point force, groups, space–frequency, time domain, uncertainty |
| `theory/` | inverse problems, regularization, Bayesian model, algorithms, uncertainty, references |
| `refs.bib` | bibliography of the theory pages |
| `theme/` | `calepin` theme + VitePress-like skin in the Julia colors |
| `_lib/ui.typ` | Typst components (hero, cards, badges) |
| `_lib/setup.jl` | hidden setup chunk of every page: environment, Makie theme, small helpers (`grid`, `relerr`, `noise_variance`, `convolution_matrix`) |
| `assets/` | logos and favicon |
| `pyproject.toml` | Python environment (uv): `jupyter-client` + `ipykernel` for the kernel bridge |
| `Project.toml` | Julia environment of the documentation (SourceEstimation in dev mode, StructuralVibration.jl for the synthetic data, CairoMakie) |

## Requirements

- [Calepin](https://vincentarelbundock.github.io/calepin) and [Typst](https://typst.app) ≥ 0.15
- [uv](https://docs.astral.sh/uv/)
- Julia with an IJulia kernel (`using Pkg; Pkg.add("IJulia"); using IJulia; installkernel("Julia")`)

## Build

```sh
cd docs
uv sync                                            # creates .venv (used by calepin.toml)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
calepin compile .                                  # builds _site/
calepin watch . _site --serve --open               # live preview
```

Chunk results are cached: only modified pages are re-executed (`calepin compile . --force`
re-runs everything). Each page starts its own kernel in its directory and includes
`_lib/setup.jl` in a hidden first chunk.
