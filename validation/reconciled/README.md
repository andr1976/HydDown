# CARDICE reconciled discharge — sensitivity case

A reconciliation of the pure-CO₂ CARDICE discharge model against the **true**
restriction-orifice sizes, kept **alongside** the baseline single-`Cd` calibration in
`validation/` (which stays unchanged as the fallback).

## What changed vs the baseline

The baseline inferred each orifice from the gas rate at one `Cd = 0.68`. SSRN 4292065
(GHGT-16, Drescher et al., **Table 1**) reports the actual restriction orifices —
**3, 4 and 6 mm** — so the discharge is now re-parameterised with those fixed sizes and
a physically-split discharge coefficient:

| Quantity | Value | Basis |
|---|---|---|
| **Orifice** | T5,T6 = 3 mm; T7,T9,T10 = 4 mm; **T8 = 5 mm effective** | SSRN 4292065 Table 1 (T8 corrected — see below) |
| **`Cd_gas`** | **0.84** | canonical sharp-edged gas-orifice coefficient. The raw data gives 0.87–0.92 (mean 0.89), *constant along each pressure decline*; **0.84 adopted** (within the physical range). Applied to gas releases **and** the gas tail of a liquid release. |
| **`Cd_liquid`** | **0.62** (fixed) | two-phase / flashing sharp-orifice coefficient (Darby / API 520 — cited by the paper). Applied to the liquid discharge only. |
| **`N` (HNE)** | T6 0.174, T8 0.175, T10 0.020, **pressure-scaled** | delayed/metastable flashing boost on top of `Cd_liquid`, calibrated to the measured steady liquid-drain rate; faded with (P−Ptr) as the vessel blows down (`liquid_ne_pressure_scaled`). |

All values are **calibrated to the raw 1 Hz Ineris data** (`background/CO2_blowdown-DJa.zip`,
folders `cardice-05..10`). Rates are the −d*M*/d*t* plateau (rolling 61 s slope), skipping
the initial valve-opening dead period so the estimate is not overcompensated by the
start-up spike.

## Key findings

- The raw data gives a **single gas coefficient ≈ 0.89** (0.87–0.92) that stays constant
  across each blowdown's whole pressure decline — an independent confirmation of the HEM +
  single-`Cd` gas model. The reconciled set adopts the **canonical `Cd_gas = 0.84`**;
  the gas rate then sits ~3–9 % below the (noisy) measured band while the gas-test
  durations land close to measured (T5 10.4 vs 10.9 h, T7 5.8 vs 5.9 h, T9 7.5 vs 6.7 h).
- With a **physical `Cd_liquid = 0.62`**, the liquid discharge needs a non-equilibrium
  boost `N`. Across the saturated-CO₂ regime (CARDICE + Ineris) `N` scales **linearly with
  the distance above the triple point**, `N ≈ 0.013·(P₀−Ptr[bar])` (R²≈0.88; a free power
  exponent returns 1.02, i.e. linear), rising from ~0 at the triple point to ~0.31 at
  28 bar. This is applied as `liquid_ne_pressure_scaled: true`, which fades each test's
  calibrated `N` linearly with (P−Ptr) during the blowdown — so the metastable boost
  vanishes near the triple point (why T10 at 10 bar already reads ~equilibrium).
- **T8 orifice correction:** SSRN Table 1 lists 4 mm, but T8's measured 0.314 kg/s exceeds
  Pursell's controlled 49-bar lab test (0.200 kg/s) through the *same* 4 mm hole —
  physically impossible (flow must rise with upstream pressure). An effective **~5 mm**
  orifice reconciles it and drops `N` from 0.525 to **0.175**, in line with T6.
- Every observable (pressure, inventory, discharge rate, internal fluid temperatures,
  inner/outer wall temperatures) tracks the data, and the central physics is preserved:
  gas releases retain 260–380 kg of in-vessel dry ice over multi-hour triple-point
  plateaus; liquid releases retain none.

## Engine support

Two new **optional** release keys (both default off / backward-compatible):
- `release.discharge_coef_gas` — a separate gas Cd; a liquid release uses `discharge_coef`
  (+ `liquid_nonequilibrium`) for the liquid and `discharge_coef_gas` for its gas tail.
- `release.liquid_ne_pressure_scaled` — when `true`, the non-equilibrium boost is faded
  linearly with (P−Ptr): `N_eff(P) = N · clip((P−Ptr)/(P₀−Ptr), 0, 1)`, so `N = N_user`
  at the initial pressure and → 0 at the triple point.

See `src/hyddown/co2_release.py`, `src/hyddown/hdclass.py` and `src/hyddown/validator.py`.

## Files

| File | Contents |
|---|---|
| `CARDICE_test{5..10}.yml` | reconciled inputs (true orifice + split `Cd` + calibrated `N`) |
| `CARDICE_test{5..10}_reconciled.pdf` | model-vs-data comparison figures |
| `../../scripts/cardice_reconcile_coefficients.py` | derives `Cd_gas`, `Cd_liquid`, `N` from the raw 1 Hz data |
| `../../scripts/cardice_run_reconciled.py` | runs the six cases and overlays the data (6-panel comparison) |

## Reproduce

```bash
# needs the raw data zip at background/CO2_blowdown-DJa.zip (~160 MB, not committed);
# both scripts extract what they need into validation/reconciled/_data on first run.
python scripts/cardice_reconcile_coefficients.py   # prints the reconciliation table
python scripts/cardice_run_reconciled.py           # regenerates the six comparison PDFs
```
