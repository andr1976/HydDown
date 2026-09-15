# Test 72 residual dry-ice gap: a structural limit of the equilibrium descent

**Question (task 3):** thermesh closes Test 71 (8.5 vs 8.4 kg) but only partly closes Test 72
(8.8 vs measured 11.7 kg). The model *forms* the right amount of dry ice at the plateau end
(peak 11.49 kg ≈ measured 11.7) but then re-sublimes ~2.7 kg during the long (~110 s) descent.
Is the descent over-sublimation limited by the wall conduction or by the gas→dry-ice interphase?

**Answer: neither.** The retained dry ice is insensitive to both:

| Model variant (Test 72)                     | retained dry ice | peak formed |
|---------------------------------------------|:----------------:|:-----------:|
| baseline (lumped wall, `solid_h_gas_solid`=3) | 8.80 kg          | 11.49 kg    |
| gas→solid interphase `solid_h_gas_solid`=1.0  | 8.81 kg          | 11.49 kg    |
| gas→solid interphase `solid_h_gas_solid`=0.0  | 8.81 kg          | 11.49 kg    |
| below-triple wall routed through **thermesh** | 8.80 kg          | 11.49 kg    |
| **measured**                                  | **11.7 kg**      | –           |

* **Wall model:** routing the below-triple descent wall through the transient 1-D conductor
  (`thermesh`) instead of the lumped node changes nothing (8.80 → 8.80), because over the 110 s
  descent the thermal penetration √(αt) ≈ 21 mm reaches most of the 25.4 mm wall, so the two wall
  models deliver essentially the same heat. That change was implemented, verified to be a no-op,
  and reverted (it doubled runtime for zero accuracy gain).
* **Interphase HTC:** removing the gas→dry-ice heat path entirely (`solid_h_gas_solid`=0) *also*
  changes nothing (8.81 kg). The dry ice is adiabatic to the gas yet still loses the same ~2.7 kg.

The sublimation is therefore **structural**: the two-zone descent holds the dry ice on the
sublimation line, so as the vessel depressurises (5.18 → ~1 bar) the solid re-sublimes to stay in
equilibrium — independent of any heat-transfer coefficient. The real Test 72 bed retains that mass,
i.e. the physical dry ice departs from sublimation-line equilibrium (kinetically / mass-transfer
"frozen") in the slower small-nozzle blowdown. This is a limitation of the *equilibrium* 0-D
descent, not a wall- or HTC-calibration issue, and is consistent with the nozzle dependence the
reference model also does not fully capture. Test 71's fast descent (~40 s) re-sublimes little, so
it matches (8.5 vs 8.4).

Figure: `test72_descent_sensitivity.pdf`.
