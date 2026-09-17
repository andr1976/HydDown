# Configuration file for the HydDown CO2 Release Technical Reference (Sphinx).
# Mirrors the ThermoCalc techref build (mathjax + sphinxcontrib.bibtex).

project = "HydDown CO2 Release - Technical Reference"
author = "Anders Andreasen"
copyright = "2025, Anders Andreasen"
release = "2.0"

extensions = [
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"

templates_path = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# --- HTML ---
html_theme = "alabaster"
html_static_path = []
html_title = "HydDown CO2 Release Technical Reference"

# --- LaTeX / PDF ---
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    "preamble": r"""
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{tikz}
\usetikzlibrary{shapes.geometric,arrows.meta,positioning,fit,backgrounds,calc}
% Decision-tree / phase-map node styles (used in the Model Map chapter).
\tikzset{
  phase/.style   = {rectangle, rounded corners, draw=black!70, very thick,
                     fill=blue!6, align=center, inner sep=5pt, minimum height=1.0cm,
                     text width=3.4cm, font=\small},
  model/.style   = {rectangle, draw=black!35, thin, fill=black!3, align=left,
                     inner sep=4pt, text width=4.6cm, font=\scriptsize},
  branch/.style  = {diamond, aspect=2, draw=black!70, very thick, fill=orange!12,
                     align=center, inner sep=1pt, text width=2.6cm, font=\scriptsize},
  flowline/.style= {-{Latex[length=2.2mm]}, thick, black!75},
}
% NB: keep :cite: roles out of figure/table captions. sphinxcontrib-bibtex renders a
% citation as a fragile bracket-protected \hyperlink, which breaks hyperref/nameref
% caption-title extraction (\Hy@tempa "extra }") inside a \caption moving argument.
% Attribution is placed in a legend paragraph or the introducing sentence instead.
""",
}
latex_documents = [
    ("index", "hyddown_co2_techref.tex",
     "HydDown CO2 Release --- Technical Reference",
     "Anders Andreasen", "manual"),
]

numfig = True
math_number_all = False
