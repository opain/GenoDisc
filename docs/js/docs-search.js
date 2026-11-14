/* GenoDisc docs search.
   Dependency-free client-side search over section headings across the docs
   pages. The index (search-index.json) is generated from the rendered HTML by
   build_search_index.py and lists { page, url, heading, id } per heading.
   A small compass-styled search box is injected top-right, beside the theme
   toggle. Substring + token matching is plenty for a handful of pages. */
(function () {
  "use strict";

  function debounce(fn, ms) {
    var t;
    return function () {
      var args = arguments, self = this;
      clearTimeout(t);
      t = setTimeout(function () { fn.apply(self, args); }, ms);
    };
  }

  function build() {
    var box = document.createElement("div");
    box.className = "docs-search";
    var input = document.createElement("input");
    input.type = "search";
    input.placeholder = "Search docs…";
    input.setAttribute("aria-label", "Search documentation");
    var results = document.createElement("div");
    results.className = "docs-search-results";
    box.appendChild(input);
    box.appendChild(results);
    document.body.appendChild(box);
    return { input: input, results: results, box: box };
  }

  function score(entry, terms) {
    var h = entry.heading.toLowerCase();
    var hit = 0;
    for (var i = 0; i < terms.length; i++) {
      var idx = h.indexOf(terms[i]);
      if (idx === -1) return -1;            // every term must appear
      hit += (idx === 0 ? 3 : 1);           // prefix matches rank higher
    }
    return hit;
  }

  function render(results, matches, query) {
    results.innerHTML = "";
    if (!matches.length) {
      var e = document.createElement("div");
      e.className = "ds-empty";
      e.textContent = 'No matches for "' + query + '"';
      results.appendChild(e);
      results.classList.add("open");
      return;
    }
    matches.slice(0, 20).forEach(function (m) {
      var a = document.createElement("a");
      a.href = m.url + (m.id ? "#" + m.id : "");
      var page = document.createElement("span");
      page.className = "ds-page";
      page.textContent = m.page;
      a.appendChild(page);
      a.appendChild(document.createTextNode(m.heading));
      results.appendChild(a);
    });
    results.classList.add("open");
  }

  function wire(ui, index) {
    var onType = debounce(function () {
      var q = ui.input.value.trim().toLowerCase();
      if (q.length < 2) { ui.results.classList.remove("open"); ui.results.innerHTML = ""; return; }
      var terms = q.split(/\s+/);
      var matches = [];
      for (var i = 0; i < index.length; i++) {
        var s = score(index[i], terms);
        if (s >= 0) matches.push({ e: index[i], s: s });
      }
      matches.sort(function (a, b) { return b.s - a.s; });
      render(ui.results, matches.map(function (m) { return m.e; }), ui.input.value.trim());
    }, 120);

    ui.input.addEventListener("input", onType);
    ui.input.addEventListener("keydown", function (ev) {
      if (ev.key === "Escape") { ui.results.classList.remove("open"); ui.input.blur(); }
    });
    document.addEventListener("click", function (ev) {
      if (!ui.box.contains(ev.target)) ui.results.classList.remove("open");
    });
  }

  function init() {
    var ui = build();
    fetch("search-index.json")
      .then(function (r) { return r.ok ? r.json() : []; })
      .then(function (index) {
        if (!Array.isArray(index) || !index.length) { ui.box.style.display = "none"; return; }
        wire(ui, index);
      })
      .catch(function () { ui.box.style.display = "none"; });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
