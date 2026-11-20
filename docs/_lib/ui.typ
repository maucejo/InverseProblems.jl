// UI components of the SourceEstimation.jl documentation (HTML target).
// Styles live in theme/css/sourceestimation.css.

#let bm(x) = $bold(upright(#x))$

// Generic HTML element with an optional class. Outside of the HTML target (e.g. the query pass of Calepin, which collects the labels and headings of the page), the body is returned
// as is, so that the labels and headings it contains are not lost.
#let el(tag, class: none, attrs: (:), body) = context {
  if target() == "html" {
    let a = attrs
    if class != none { a.insert("class", class) }
    html.elem(tag, attrs: a, body)
  } else {
    body
  }
}

// <div class="…"> with an optional inline style
#let _div(class, body, style: none) = el("div", class: class,
  attrs: if style != none { (style: style) } else { (:) }, body)

// Julia colors, used for icons and badges
#let julia = (
  blue: "#4063D8",
  green: "#389826",
  purple: "#9558B2",
  red: "#CB3C33",
)

// VitePress-like hero: name, text, tagline, action buttons and image
#let hero(name: none, text: none, tagline: none, actions: (), image: none, install: none) = {
  html.elem("section", attrs: (class: "se-hero"))[
    #_div("se-hero-main")[
      #html.elem("h1", attrs: (class: "se-hero-name"))[#name]
      #html.elem("p", attrs: (class: "se-hero-text"))[#text]
      #html.elem("p", attrs: (class: "se-hero-tagline"))[#tagline]
      #_div("se-hero-actions")[
        #for a in actions {
          html.elem("a", attrs: (
            class: "landing-button " + if a.at("primary", default: false) { "landing-button-primary" } else { "landing-button-secondary" },
            href: a.href,
          ))[#a.label]
        }
      ]
      #if install != none {
        html.elem("div", attrs: (class: "se-install"))[
          #html.elem("span", attrs: (class: "se-prompt"))[pkg>]
          #html.elem("span")[#install]
        ]
      }
    ]
    #_div("se-hero-image")[
      #html.elem("img", attrs: (src: image, alt: "SourceEstimation.jl logo"))
    ]
  ]
}

// Feature card: icon (emoji or text), title, details, optional link and icon color
#let card(icon: none, title: none, color: julia.blue, href: none, body) = {
  let inner = _div("se-card", style: "--se-icon: " + color + ";")[
    #if icon != none { _div("se-card-icon")[#icon] }
    #html.elem("p", attrs: (class: "se-card-title"))[#title]
    #html.elem("p", attrs: (class: "se-card-details"))[#body]
  ]
  if href != none {
    html.elem("a", attrs: (class: "se-card-link", href: href), inner)
  } else {
    inner
  }
}

// Grid of feature cards (3 columns on the landing page)
#let features(..cards) = _div("se-features")[#cards.pos().join()]

// Grid of "next steps" cards, inside regular pages
#let next-steps(..cards) = _div("se-next")[#cards.pos().join()]

// Two-column layout
#let split(..cols) = _div("se-split")[#cols.pos().map(c => _div("se-split-col")[#c]).join()]

// Small colored pills
#let badge(body, color: julia.blue) = html.elem("span", attrs: (class: "se-badge", style: "--se-icon: " + color + ";"))[#body]
#let badges(..items) = _div("se-badges")[#items.pos().join()]

// Lead paragraph of a section
#let lead(body) = html.elem("p", attrs: (class: "se-section-lead"))[#body]

// ---------------------------------------------------------------------------
// Equations in HTML pages. Use in a page: #show: math-html
//  - block equations are rendered by Typst as SVG (exact layout, numbering
//    included), wrapped in <div class="mx-equation"> so the theme can centre them;
//  - inline equations stay MathML. Typst emits combining marks (U+0307,
//    U+0308…) for dot(), dot.double()…, which browsers draw misplaced, so they
//    are swapped for the spacing variants, which MathML positions correctly.
// ---------------------------------------------------------------------------
#let _spacing-accents = (
  "\u{300}": $\u{60}$,   // grave
  "\u{301}": $\u{B4}$,   // acute
  "\u{303}": $\u{2DC}$,  // tilde
  "\u{304}": $\u{AF}$,   // macron
  "\u{306}": $\u{2D8}$,  // breve
  "\u{307}": $\u{2D9}$,  // dot
  "\u{308}": $\u{A8}$,   // dot.double
  "\u{30A}": $\u{2DA}$,  // circle
  "\u{30C}": $\u{2C7}$,  // caron
)

#let math-font = "STIX Two Math"

// URL prefix from the current page back to the site root (same logic as
// Calepin's runtime, which cannot be imported from a helper module).
#let _site-root() = {
  let override = sys.inputs.at("calepin-site-root", default: "")
  if override != "" { return override }
  let href = sys.inputs.at("calepin-current-href", default: "")
  let depth = href.split("/").filter(part => part != "").len() - 1
  if depth <= 0 { "" } else { "../" * depth }
}

// Web font of the inline (MathML) equations. The theme CSS is inlined in every
// page, so a relative URL in it would depend on the depth of the page: the
// @font-face is emitted here, with the path resolved by Calepin for this page.
#let _math-font-face() = context {
  if target() == "html" {
    let href = _site-root() + "assets/fonts/STIXTwoMath-Regular.woff2"
    el("style", "@font-face { font-family: \"STIX Two Math\"; src: local(\"STIX Two Math\"), local(\"STIXTwoMath-Regular\"), url(\"" + href + "\") format(\"woff2\"); font-display: swap; }")
  }
}

// Fences of the inline (MathML) equations. Browsers stretch the fences (parentheses, brackets, bars) to the height of their content, and with the large ascent/descent of the math italic letters they pick an enlarged glyph with wide side bearings, which looks like spaces around the content: (μ, Σ) → ( μ, Σ ). Typst cannot emit `stretchy="false"`, so this script sets it, except when the fence encloses tall content (fraction, limits, root, matrix), where stretching is needed.
#let _fence-script = "document.addEventListener('DOMContentLoaded', function () {
  var fences = '()[]{}|‖⟨⟩⌊⌋⌈⌉';
  var tall = 'mfrac, munder, mover, munderover, msqrt, mroot, mtable';
  document.querySelectorAll('math mo').forEach(function (mo) {
    if (fences.indexOf(mo.textContent.trim()) < 0) return;
    var row = mo.parentElement;
    if (row && !row.querySelector(tall)) mo.setAttribute('stretchy', 'false');
  });
});"

#let _math-fences() = context {
  if target() == "html" { el("script", _fence-script) }
}

#let math-html(body) = {
  // Font of the SVG equations (Typst draws the glyphs, so it only has to be installed where `calepin compile` runs). Inline MathML uses the web font declared by _math-font-face.
  show math.equation: set text(font: math-font)
  show math.equation.where(block: true): it => context {
    if target() != "html" { return it }
    if it.numbering == none {
      el("div", class: "mx-equation", html.frame(it))
    } else {
      // The frame is as wide as the equation: render the number outside of it,
      // so that the theme can align it on the right of the column.
      let num = numbering(it.numbering, ..counter(math.equation).at(it.location()))
      el("div", class: "mx-equation mx-numbered")[
        #el("span", class: "mx-eq-body",
          html.frame(math.equation(block: true, numbering: none, it.body)))
        #el("span", class: "mx-eq-num")[#num]
      ]
    }
  }
  show math.accent: it => context {
    if target() == "html" and it.accent in _spacing-accents {
      html.elem("mover", attrs: (accent: "true"),
        { it.base; _spacing-accents.at(it.accent) })
    } else { it }
  }
  _math-font-face()
  _math-fences()
  body
}
