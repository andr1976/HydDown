# HydDown CO₂ Release — Technical Reference

Sphinx/reStructuredText technical reference for the thermopack-based CO₂
release / dry-ice add-on developed on the `co2-release-hem` branch. It documents
the physical models, thermodynamics, heat/mass transfer, discharge (HEM/HNE),
dry-ice estimation, and the CARDICE validation.

## Contents

| File | Chapter |
|------|---------|
| `index.rst` | Title + table of contents |
| `introduction.rst` | Motivation, scope, architecture |
| `experimental.rst` | CARDICE / Ineris setup, test matrix, external-HT derivation |
| `thermodynamics.rst` | thermopack backend, triple point, CoolProp↔thermopack hand-overs |
| `discharge.rst` | HEM discharge, non-equilibrium (HNE) liquid boost, orifice/Cd |
| `dry_ice.rst` | Atmospheric and in-vessel dry-ice estimation |
| `heat_transfer.rst` | Two-node wall, natural convection, initialisation, zone balances |
| `validation.rst` | CARDICE tests 5–10 vs 1 Hz data |
| `references.bib` | Bibliography |
| `figures/` | Paper illustrations + validation figures |

## Building

Requires `sphinx` and `sphinxcontrib-bibtex`:

```bash
pip install sphinx sphinxcontrib-bibtex
cd docs/techref
make html        # -> _build/html/index.html
make latexpdf    # -> _build/latex/hyddown_co2_techref.pdf  (needs a LaTeX toolchain)
```

## Regenerating the figures

Paper illustrations in `figures/` (`vaillant_*.png`, `jamois_*.png`) are cropped from
the two source papers (Vaillant et al. GHGT-15 / SSRN 3823308; Jamois et al. IJGGC
129:103974) and included with attribution. The `cardice_*.png` validation figures are
rendered from `validation/CARDICE_*_full.pdf` (produced by the model-vs-data overlay
scripts).
