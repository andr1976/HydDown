# Supplementary material

Input files and the two-dimensional wall-conduction solver that reproduce the
results in the manuscript. They run with HydDown (development branch
`co2-release-hem`): <https://github.com/andr1976/HydDown>.

## Contents

```
cardice/    CARDICE (Ineris 2 m3 sphere) input files, tests 5-10
sintef/     ECCSEL/SINTEF dense-phase input files (9 tests) + two 2-D-wall variants
wall2d.py   Standalone 2-D axisymmetric conjugate wall-conduction solver
```

All input files use the two-zone non-equilibrium (partial-phase-equilibrium)
description with the solid-CO2 (dry-ice) extension. The vapour-contact wall uses
the Churchill--Chu correlation and the wetted wall the Cooper boiling correlation;
the CARDICE liquid tests (6, 8, 10) use a fixed wetted-wall coefficient.

### CARDICE (`cardice/`)

Ineris 2 m3 horizontal sphere, low- and medium-pressure blowdown.

| File | Phase |
|------|-------|
| `CARDICE_test5.yml`  | gas |
| `CARDICE_test6.yml`  | liquid |
| `CARDICE_test7.yml`  | gas |
| `CARDICE_test8.yml`  | liquid |
| `CARDICE_test9.yml`  | gas |
| `CARDICE_test10.yml` | liquid |

### ECCSEL/SINTEF (`sintef/`)

Vertical 273 mm ID / 1 m cylinder, dense-phase (~120 bar) blowdown. `False`
denotes a no-riser (gas-space) release, `True` a riser (liquid-space) release.

| File | P0 (bar) | T0 (C) | Nozzle (mm) | Release |
|------|----------|--------|-------------|---------|
| `Exp71.yml` | 122.6 | 25.2 | 8.0 | gas |
| `Exp72.yml` | 119.0 | 24.9 | 6.5 | gas |
| `Exp75.yml` | 119.0 | 25.0 | 4.5 | gas |
| `Exp52.yml` | 119.9 | 15.4 | 8.0 | liquid |
| `Exp53.yml` | 119.5 | 24.4 | 8.0 | liquid |
| `Exp56.yml` | 119.1 | 15.2 | 6.5 | liquid |
| `Exp57.yml` | 116.8 | 24.5 | 6.5 | liquid |
| `Exp45.yml` | 119.9 | 14.5 | 4.5 | liquid |
| `Exp46.yml` | 116.7 | 24.4 | 4.5 | liquid |

`Exp53_2dwall.yml` and `Exp72_2dwall.yml` repeat those two tests with the
two-dimensional conjugate wall enabled (`vessel.wall_model: 2d`, with the bottom
plate, flange and lid dimensions), reproducing the wall-temperature comparison in
Section 5.5. Every other sub-model is unchanged, so the pair isolates the effect
of the wall representation.

## Running

Install HydDown (branch `co2-release-hem`) and run any file:

```
python scripts/hyddown_main.py sintef/Exp72.yml
```

If a near-triple flash step fails at `time_step: 0.1`, reduce it to `0.02`.

## `wall2d.py`

The standalone 2-D axisymmetric conjugate wall-conduction solver used by the
`wall_model: 2d` option (masked structured finite-volume r-z grid, adiabatic
outer surface, Robin inner condition split at the moving liquid/condensate level:
pool boiling below, free convection above). The model and its relation to the
ECCSEL/SINTEF reference model are described in Appendix A of the manuscript. In
HydDown it lives at `src/hyddown/wall2d.py`; the copy here is for reference.
