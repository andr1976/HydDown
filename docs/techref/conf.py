# Configuration file for the HydDown CO2 Release Technical Reference (Sphinx).
# Mirrors the ThermoCalc techref build (mathjax + sphinxcontrib.bibtex).

project = "HydDown CO2 Release - Technical Reference"
author = "Anders Andreasen"
copyright = "2025, Anders Andreasen"
release = "1.0"

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
""",
}
latex_documents = [
    ("index", "hyddown_co2_techref.tex",
     "HydDown CO2 Release --- Technical Reference",
     "Anders Andreasen", "manual"),
]

numfig = True
math_number_all = False
