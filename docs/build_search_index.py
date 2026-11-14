#!/usr/bin/env python3
"""Generate search-index.json for the docs search box.

Parses the rendered documentation pages and emits one entry per section
heading: {"page": ..., "url": ..., "heading": ..., "id": ...}. Run this after
rmarkdown::render_site() whenever page headings change:

    python3 build_search_index.py

Stdlib only (html.parser) so it runs anywhere R/pandoc already runs.
"""
import json
import os
import re
from html.parser import HTMLParser

# Published pages only (the internal source-only notebooks are not on the site).
PAGES = [
    ("index.html", "Home"),
    ("pipeline_guide.html", "User Guide"),
    ("pipeline_technical.html", "Technical Details"),
]

HERE = os.path.dirname(os.path.abspath(__file__))


class HeadingExtractor(HTMLParser):
    """rmarkdown wraps each header in <div class="section ..." id="X"><hN>...</hN>,
    so the anchor id lives on the enclosing section div, not the heading. Track a
    stack of div ids (None for non-section divs) and use the nearest section id."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.headings = []            # (id, text)
        self._divs = []               # stack: section id, or None per open <div>
        self._in_h = False
        self._buf = []

    def handle_starttag(self, tag, attrs):
        if tag == "div":
            a = dict(attrs)
            cls = a.get("class", "")
            self._divs.append(a.get("id") if "section" in cls else None)
        elif tag in ("h1", "h2", "h3", "h4"):
            self._in_h = True
            self._buf = []

    def handle_startendtag(self, tag, attrs):
        # self-closing (e.g. <a ... />) inside a heading contributes no text
        pass

    def handle_endtag(self, tag):
        if tag == "div" and self._divs:
            self._divs.pop()
        elif tag in ("h1", "h2", "h3", "h4") and self._in_h:
            self._in_h = False
            text = re.sub(r"\s+", " ", "".join(self._buf)).strip().replace("¶", "").strip()
            sec_id = next((d for d in reversed(self._divs) if d), None)
            if text and sec_id:
                self.headings.append((sec_id, text))
            self._buf = []

    def handle_data(self, data):
        if self._in_h:
            self._buf.append(data)


def main():
    index = []
    for fname, page in PAGES:
        path = os.path.join(HERE, fname)
        if not os.path.exists(path):
            print("skip (missing):", fname)
            continue
        with open(path, encoding="utf-8") as fh:
            parser = HeadingExtractor()
            parser.feed(fh.read())
        for hid, text in parser.headings:
            # skip the hidden document-title h1 and the auto "Contents" TOC header
            if text.lower() in ("contents",):
                continue
            index.append({"page": page, "url": fname, "heading": text, "id": hid})
        print("%-26s %d headings" % (fname, len(parser.headings)))

    out = os.path.join(HERE, "search-index.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(index, fh, ensure_ascii=False, indent=0, separators=(",", ":"))
    print("wrote", out, "(%d entries)" % len(index))


if __name__ == "__main__":
    main()
