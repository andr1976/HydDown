# Munkejord / SINTEF Test 71 - dense-phase CO2 vessel depressurisation

HydDown validation against the SINTEF/NTNU ECCSEL dense-phase CO2 blowdown experiment
**Test 71** (8 mm nozzle, no riser, 122.6 bar / 25.2 C).

**Reference:** Hoydalsvik, Austegard, Deng, Blakseth, Hafner, Munkejord, *"Experiments and
modelling of CO2 vessel depressurization: flash boiling, heat transfer and critical flow through
nozzles"*, Applied Thermal Engineering 306 (2026) 133032 (open access, CC BY 4.0).
**Data:** Zenodo record 19589510 (CC BY 4.0). `exp71_blowdown.csv` is a downsampled extract of the
Test 71 1 Hz channels (pressure PT163/PT162, load-cell weight, fluid/wall thermocouples), time-
gridded to 0.5 s over the blowdown; see the paper/dataset for the full record.

## Files
- `test71_nem.yml` - HydDown NEM input. Because `hdclass` cannot yet initialise a dense single
  phase, the case starts at the **two-phase onset** (~52 bar, the paper's 25 C flash point) where
  NEM applies, with the ~47 kg inventory remaining after the ~1 s dense prelude. A dense-phase-start
  extension (single-phase -> flash handoff -> NEM) is planned to run it from the true 122.6 bar.
- `run_test71_nem.py` - run the NEM case (dt 0.02 s) and save arrays to `test71_nem_out.npz`.
- `plot_test71_nem.py` - overlay the NEM run vs the 1 Hz data -> `test71_nem_overlay.pdf`.
- `trial_test71.py` - standalone 0-D dense-phase trial driver (reuses `CO2ReleaseModelCP`;
  CoolProp above triple + thermopack-free solid below). Env vars: `EOS` (CP/tcPR/GERG2008),
  `HBOIL`, `MSTEEL`, `ROHSENOW`, `RCOND`.
- `test71_full_plot.py`, `mass_sweep_plot.py` - trial-based diagnostics (full overlay; retained
  dry ice vs lumped wall mass).

Note: some scripts contain absolute repo paths - adjust `BASE`/`sys.path` for your checkout.
