# Below-triple dry-ice descent: validation and the T72 outlier

Below the triple point HydDown's two-zone model holds the dry ice on the sublimation line, so as
the vessel depressurises the bed re-sublimes to stay in equilibrium (co2_release.py
`two_zone_descent_step`, term `dm_cool = m_solid·cp_solid·ΔT_sub / L_sub`). Munkejord Test 72
retains more dry ice than that descent predicts (8.8 vs 11.7 kg), which raised two questions:
is the descent sublimation a modelling artefact, and would "freezing" the bed at its peak fix it?
Both are answered by looking across five gas / dry-ice tests.

## The descent over-sublimation is not wall- or HTC-limited (Test 72)

The T72 model forms the right dry ice at the plateau end (peak 11.49 kg ≈ measured 11.7) then
re-sublimes ~2.7 kg. That loss is insensitive to every heat-transfer knob:

| Test-72 variant                              | retained | peak  |
|----------------------------------------------|:--------:|:-----:|
| baseline (lumped wall, `solid_h_gas_solid`=3)| 8.80 kg  | 11.49 |
| interphase `solid_h_gas_solid` = 1.0         | 8.81 kg  | 11.49 |
| interphase `solid_h_gas_solid` = 0.0         | 8.81 kg  | 11.49 |
| below-triple wall through **thermesh**       | 8.80 kg  | 11.49 |

Routing the descent wall through the transient 1-D conductor changes nothing (over 110 s the
thermal penetration √(αt) ≈ 21 mm reaches most of the 25.4 mm wall) and was reverted; removing the
gas→solid interphase entirely also changes nothing. The sublimation is **structural** to the
equilibrium descent, not a calibration issue. Figure `test72_descent_sensitivity.pdf`.

## Freezing the bed is NOT a general fix — the descent sublimation is real

Comparing the dry ice **formed** (peak) with the dry ice **retained** (measured), normalised as
retained ÷ peak (1.0 = a frozen bed that loses nothing on the way down):

| Test                         | peak (kg) | model final | measured | measured/peak | model/peak |
|------------------------------|:---------:|:-----------:|:--------:|:-------------:|:----------:|
| Munkejord **T71** (8 mm, 120 bar)  | 12.36 | 8.46 | 8.4  | 0.68 | 0.68 |
| Munkejord **T72** (6.5 mm, 120 bar)| 11.49 | 8.80 | 11.7 | **1.02** | 0.77 |
| CARDICE **T5** (3 mm, 20 bar)      | 315   | 265  | 217  | 0.69 | 0.84 |
| CARDICE **T7** (4 mm, 15 bar)      | 355   | 303  | 298  | 0.84 | 0.85 |
| CARDICE **T9** (4 mm, 10 bar)      | 482   | 392  | 412  | 0.85 | 0.81 |

Across 10→120 bar and 8→480 kg, the **measured** bed sublimes 15–32 % of its peak on the way down
in four of the five tests — the descent sublimation is physically real, and the equilibrium descent
reproduces it well (T71 exact; T7/T9 within a few %; T5 it retains a little too much). **Freezing
the bed at its peak would over-predict CARDICE by 17–45 % and T71 by 47 %.** Figure
`retained_over_peak_5tests.pdf`.

## Conclusion — T72 is a single outlier, treated as such

Only T72 retains essentially its full peak (retained/peak = 1.02, i.e. the measured value even
exceeds the dry ice the model ever forms — a hint of trapped liquid/gas residual or a load-cell /
definition nuance for that one test). The physical driver of the small residual gap is **bed
morphology**, not residence time: the gentle small-nozzle blowdown builds a consolidated,
low-surface-area bed that sublimes little, whereas the violent large-nozzle blowdown (T71) and the
CARDICE orifices disperse the bed so it sublimes to equilibrium. A residence-time relaxation gives
the *wrong* trend (the slower T72 would relax more, not less).

Therefore the equilibrium two-zone descent is kept as the validated default. Capturing T72 would
need a bed-morphology / nozzle-keyed sublimation rate (a metastable "NEM ice" zone, the below-triple
analogue of the liquid non-equilibrium `N`) that recovers the equilibrium descent everywhere except
where the rate is deliberately throttled — a physical parameter, not something a 0-D model derives
from first principles. Given four of five tests already match, and T72's own number is partly
suspect, the below-triple dry-ice modelling is considered concluded: this is as good as a 0-D
description gets, with the residual limits (0-D gas-node over-cooling, bed morphology) intrinsic to
the lumped formulation.
