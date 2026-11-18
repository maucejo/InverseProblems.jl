/* Function calls in the Julia code blocks.
 *
 * The Julia grammar used by Typst does not scope the names of called functions
 * (`modefreq(beam, 5000.)`): they are left as plain text, while the themes in
 * theme/syntax/ color everything else. This script wraps the identifiers followed by
 * `(` in <span class="se-call">, only in the unscoped text of the code (never in strings,
 * comments, or other highlighted tokens). Colors: .se-call in theme/css/sourceestimation.css.
 */
(() => {
  "use strict";
  const IDENT = "[\\p{L}_][\\p{L}\\p{N}_!₀-₉ᵢ-ᵥ′]*";
  // A call inside a text node: `name(` ...
  const CALL = new RegExp("(" + IDENT + ")(?=\\.?\\()", "gu");
  // ... or a name ending a text node whose next sibling (the highlighted `(`) opens the call
  const TRAILING = new RegExp("(" + IDENT + ")$", "u");
  const KEYWORDS = new Set(["if", "elseif", "for", "while", "function", "return", "macro", "begin", "let", "do", "quote", "in", "isa"]);

  // `name(` or the broadcast form `name.(`, the dot being a separate highlighted token
  function opensCall(node) {
    const next = node.nextSibling;
    if (!next) return false;
    const t = next.textContent;
    if (t.startsWith("(") || t.startsWith(".(")) return true;
    return t === "." && Boolean(next.nextSibling && next.nextSibling.textContent.startsWith("("));
  }

  function wrapCalls(code) {
    for (const node of Array.from(code.childNodes)) {
      if (node.nodeType !== Node.TEXT_NODE) continue;
      const text = node.nodeValue;
      const spans = [];                              // [start, end] of the call names
      for (const m of text.matchAll(CALL)) {
        if (!KEYWORDS.has(m[1])) spans.push([m.index, m.index + m[1].length]);
      }
      const t = opensCall(node) ? text.match(TRAILING) : null;
      if (t && !KEYWORDS.has(t[1])) spans.push([t.index, t.index + t[1].length]);
      if (!spans.length) continue;

      const fragment = document.createDocumentFragment();
      let last = 0;
      for (const [a, b] of spans) {
        fragment.append(text.slice(last, a));
        const span = document.createElement("span");
        span.className = "se-call";
        span.textContent = text.slice(a, b);
        fragment.append(span);
        last = b;
      }
      fragment.append(text.slice(last));
      node.replaceWith(fragment);
    }
  }

  function run() {
    document.querySelectorAll('pre code[data-lang="julia"]').forEach(wrapCalls);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", run);
  } else {
    run();
  }
})();
