# Generates reference/api.typ from the docstrings of the exported symbols of SourceEstimation.jl.
#
#   julia --project=docs docs/_lib/make_api.jl      (from the repository root)
#   julia --project=. _lib/make_api.jl              (from docs/)
#
# The page is generated as static Typst source, rather than printed by a chunk with
# `results: "typst"`: code blocks are then ordinary fenced blocks (`#| eval: false`), which
# Calepin highlights with the theme of the site (light and dark modes). Run it again whenever a
# docstring changes; `build.sh` does it before each build.

import Pkg
Pkg.activate(dirname(@__DIR__); io = devnull)

using SourceEstimation
using Markdown

const MOD = SourceEstimation
const OUT = joinpath(dirname(@__DIR__), "reference", "api.typ")

# Sections of the page: title, short introduction (Typst markup), exported symbols
const SECTIONS = [
    ("Problem", "Definition of the inverse problem `y = H u + n`.",
     [:InverseProblem]),
    ("Regularization terms", "Priors on the sources: ℓq norms, mixed norms and region-wise groups.",
     [:Lq, :Lpq, :Groups]),
    ("Smoothing operators", "Finite-difference operators and helpers for vectorized space–frequency or space–time problems.",
     [:fd_weights, :diff_operator, :space_sample_operator, :sample_rows]),
    ("Formulations", "How data fidelity and regularization are combined.",
     [:Additive, :Multiplicative, :Bayesian, :Estimated, :Fixed]),
    ("Parameter choice rules", "Selection or initialization of the regularization parameter.",
     [:LCurve, :GCV, :BayesEstimator, :Heuristic, :LambdaMax]),
    ("Algorithms", "Resolution algorithms, passed as third argument of `solve`.",
     [:IR, :MAP, :ElementWise]),
    ("Resolution", "Solving a problem and inspecting its solution.",
     [:solve, :InverseSolution]),
    ("Uncertainty quantification", "Posterior quantities of a solution obtained with `solve(...; precision = true)`.",
     [:posterior_std, :posterior_covariance, :posterior_samples, :credible_interval]),
]

const DOCUMENTED = Set(s for (_, _, syms) in SECTIONS for s in syms)

# ---------------------------------------------------------------------------- Typst helpers

label(s) = "api-" * string(s)

"""Escape a plain string for Typst markup."""
function esc(s::AbstractString)
    s = replace(s, r"([\\#*_$@<>\[\]`~])" => s"\\\1")
    s = replace(s, "//" => "\\/\\/")          # `//` starts a comment
    return s
end

"""Inline code: backticks, or `#raw` when the code itself contains a backtick."""
inline_code(c) = occursin('`', c) ? "#raw(\"" * escape_string(c) * "\")" : "`" * c * "`"

# ------------------------------------------------------------------- Markdown → Typst

inline(x::AbstractString) = esc(x)
inline(x::Markdown.Code) = inline_code(x.code)
inline(x::Markdown.Bold) = "*" * inline(x.text) * "*"
inline(x::Markdown.Italic) = "_" * inline(x.text) * "_"
inline(x::Markdown.LineBreak) = " \\\n"
inline(x::Markdown.LaTeX) = inline_code(x.formula)
inline(v::AbstractVector) = join(inline.(v))
inline(x) = esc(sprint(Markdown.plain, x))

function inline(x::Markdown.Link)
    text = inline(x.text)
    if x.url == "@ref"
        # [`name`](@ref): link to the entry of the page when it is documented here
        name = Symbol(strip(Markdown.plain(Markdown.Paragraph(x.text))))
        name = Symbol(strip(string(name), ['`', '\n', ' ']))
        return name in DOCUMENTED ? "#link(<" * label(name) * ">)[" * text * "]" : text
    end
    return "#link(\"" * x.url * "\")[" * text * "]"
end

block(x::Markdown.Paragraph; first_code = Ref(false)) = inline(x.content) * "\n"
block(x::Markdown.Header{N}; kw...) where {N} = "*" * inline(x.text) * "*\n"
block(x::Markdown.HorizontalRule; kw...) = "#line(length: 100%)\n"

"""Does the code parse as Julia? (Formulas written in indented blocks do not.)"""
function is_julia(code)
    ex = Meta.parseall(code; filename = "docstring")
    return !any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args)
end

function block(x::Markdown.Code; first_code = Ref(false))
    # Docstrings tag every code block as Julia, including the indented blocks used for
    # formulas. The first block is the signature; the following ones are highlighted as
    # Julia only if they parse as Julia code.
    lang = x.language in ("", "julia", "jldoctest") ? "julia" : x.language
    if lang == "julia" && first_code[] && !is_julia(x.code)
        lang = ""
    end
    first_code[] = true
    if lang == "julia"
        return "```julia\n#| eval: false\n" * rstrip(x.code) * "\n```\n"
    else
        return "```\n" * rstrip(x.code) * "\n```\n"
    end
end

function block(x::Markdown.List; kw...)
    io = IOBuffer()
    for (k, item) in enumerate(x.items)
        marker = x.ordered == -1 ? "- " : "+ "
        parts = [rstrip(block(b; kw...)) for b in item]
        body = join(parts, "\n")
        print(io, marker, replace(body, "\n" => "\n  "), "\n")
    end
    return String(take!(io))
end

function block(x::Markdown.Table; kw...)
    rows = x.rows
    ncol = length(first(rows))
    cell(c) = "[" * inline(c) * "]"
    io = IOBuffer()
    println(io, "#table(")
    println(io, "  columns: ", ncol, ",")
    println(io, "  align: left,")
    println(io, "  table.header(", join(("[*" * inline(c) * "*]" for c in first(rows)), ", "), "),")
    for r in rows[2:end]
        println(io, "  ", join(cell.(r), ", "), ",")
    end
    println(io, ")")
    return String(take!(io))
end

function block(x::Markdown.Admonition; kw...)
    kind = x.category in ("warning", "danger") ? "warning" : x.category == "tip" ? "tip" : "note"
    body = join((block(b; kw...) for b in x.content), "\n")
    return "#calepin.elements.callout(kind: \"" * kind * "\")[\n" * body * "]\n"
end

block(x::Markdown.BlockQuote; kw...) = "#quote(block: true)[\n" * join((block(b; kw...) for b in x.content), "\n") * "]\n"
block(x::Markdown.MD; kw...) = join((block(b; kw...) for b in x.content), "\n")
block(x; kw...) = esc(sprint(Markdown.plain, x)) * "\n"

# A docstring is a `MD` holding one `MD` per method, themselves possibly wrapped in another
# `MD`: return the list of entries, each one with its own blocks.
function entries(md::Markdown.MD)
    if !isempty(md.content) && all(c -> c isa Markdown.MD, md.content)
        return reduce(vcat, entries.(md.content))
    end
    return [md]
end

docentries(sym) = entries(Docs.doc(Docs.Binding(MOD, sym)))

docstring(sym) = join((block(e; first_code = Ref(false)) for e in docentries(sym)), "\n#html.elem(\"hr\", attrs: (class: \"se-api-sep\"))\n\n")

"""First sentence of the description (the paragraph after the signature)."""
function summary(sym)
    for b in first(docentries(sym)).content
        if b isa Markdown.Paragraph
            text = replace(inline(b.content), "\n" => " ")
            m = match(r"^(.*?(?<!\b[A-Z])\.)(\s|$)", text)
            return m === nothing ? text : m.captures[1]
        end
    end
    return ""
end

kind(sym) = (obj = getfield(MOD, sym); obj isa Type ? "type" : obj isa Function ? "function" : "constant")

# ------------------------------------------------------------------------------- page

function write_page(path)
    io = IOBuffer()
    print(io, """
    // GENERATED FILE — do not edit by hand. Source: the docstrings of SourceEstimation.jl.
    // Regenerate with `julia --project=. _lib/make_api.jl` (from docs/).
    #import "/.calepin/calepin.typ" as calepin
    #import "/_lib/ui.typ": *
    #show: math-html

    #set document(title: [API reference])
    #metadata((description: "Reference of the public API of SourceEstimation.jl, generated from the docstrings.")) <website-metadata>

    #title()

    This page documents the public API of SourceEstimation.jl. It is generated from the docstrings of the exported symbols (version $(pkgversion(MOD))). Everything boils down to one verb, #link(<$(label(:solve))>)[`solve`], applied to a problem, a regularization formulation and, optionally, an algorithm:

    ```julia
    #| eval: false
    prob = InverseProblem(H, y)                         # what to solve
    reg  = Multiplicative(Groups(r1 => Lq(0.5),         # how to regularize
                                 r2 => Lq(2.0; L = D)))
    sol  = solve(prob, reg)                             # default algorithm
    sol.u                                               # reconstructed field
    ```

    #calepin.elements.callout(kind: "note")[
    Internal functions are not part of the public API and may change without notice. Access them with the `SourceEstimation.` prefix at your own risk.
    ]

    = Index

    """)
    rows = String[]
    for (title, _, syms) in SECTIONS, s in syms
        push!(rows, "  [#link(<$(label(s))>)[`$(s)`]], [$(kind(s))], [$(esc(title))], [$(summary(s))],")
    end
    print(io, "#table(\n  columns: (auto, auto, auto, 1fr),\n  align: left,\n",
          "  table.header([*Name*], [*Kind*], [*Section*], [*Summary*]),\n", join(rows, "\n"), "\n)\n\n")

    for (title, intro, syms) in SECTIONS
        print(io, "= ", title, "\n\n", intro, "\n\n")   # intros are Typst markup
        for s in syms
            # One rounded box per entry; its header line holds the name (left) and the
            # kind (right). Styles: .se-api-* in theme/css/sourceestimation.css
            print(io, "#el(\"div\", class: \"se-api-entry\")[\n")
            print(io, "#el(\"div\", class: \"se-api-head\")[\n")
            print(io, "== `", s, "` <", label(s), ">\n\n")
            print(io, "#el(\"span\", class: \"se-api-kind se-api-kind-", kind(s), "\")[", uppercasefirst(kind(s)), "]\n")
            print(io, "]\n\n")
            print(io, docstring(s), "\n")
            print(io, "]\n\n")
        end
    end
    mkpath(dirname(path))
    write(path, String(take!(io)))
    return path
end

println("Written ", write_page(OUT))
