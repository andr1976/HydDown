# Cleanup Plan

## Summary

During NEM implementation and debugging, we created many temporary diagnostic scripts, plots, and intermediate documentation files. This plan organizes what to keep vs. remove.

## Files to KEEP (Production)

### Core Documentation
- ✅ `NEM_IMPLEMENTATION_METHODOLOGY.md` - Complete NEM methodology reference
- ✅ `NEM_HEAT_TRANSFER_COMPLETE_FIX.md` - Summary of all heat transfer fixes
- ✅ `H_GL_CORRELATION_FINAL.md` - Final h_gl correlation documentation

### Original Files (unchanged)
- ✅ `CLAUDE.md` - Project instructions for Claude
- ✅ `CODE_OF_CONDUCT.md`
- ✅ `CONTRIBUTING.md`
- ✅ `Manual.md`
- ✅ `README.md`
- ✅ `run_validation_suite.py` - Official validation script
- ✅ `run_validation_with_plots.py` - Validation with plotting

## Files to REMOVE (Temporary/Diagnostic)

### Diagnostic Scripts (26 files)
```
analyze_nem_pressure.py
check_gas_liquid_heat_transfer.py
compare_enthalpy_vs_average.py
compare_h_vs_u.py
compare_nem.py
compare_nem_high_hgl.py
debug_conservation.py
debug_energy_detailed.py
debug_nem_energy.py
debug_nem_step.py
detailed_array_comparison.py
diagnostic_nem.py
gas_wall_htc_diagnostic.py
h_gl_sensitivity_study.py
investigate_ra_negative.py
nem_boiling_properties_diagnostic.py
psv_opening_diagnostic.py
test_enthalpy_transfer.py
test_h_gl_correlation.py
test_h_gl_detailed.py
test_h_gl_sensitivity.py
test_hgl_sensitivity.py
test_nem.py
test_pressure_residual.py
trace_nem_state.py
trace_phase_transfer.py
wetted_htc_comparison.py
```

### Diagnostic Plots (22 files)
```
equilibrium_only.pdf
equilibrium_only.png
h_gl_correlation_test.pdf
h_gl_correlation_test.png
h_gl_sensitivity.pdf
h_gl_sensitivity_study.pdf
hgl_sensitivity.pdf
nem_comparison.pdf
nem_comparison.png
nem_comparison_hgl5000.pdf
nem_comparison_high_hgl.pdf
nem_comparison_psv_fire.pdf
nem_pressure_analysis.pdf
nem_test_results.pdf
nem_test_results.png
nem_vs_equilibrium_comparison.pdf
nem_vs_equilibrium_comparison.png
psv_opening_diagnostic.pdf
psv_opening_diagnostic.png
wetted_htc_equilibrium_vs_nem_h500.pdf
wetted_htc_equilibrium_vs_nem_h500.png
```

### Intermediate Documentation (22 files)
```
ENTHALPY_FIX_RESULTS.md
ENTHALPY_VS_INTERNAL_ENERGY_COMPARISON.md
FINAL_NEM_SOLVER.md
FINAL_PHASE_TRANSFER_COMPARISON.md
GAS_WALL_HEAT_TRANSFER_FIX.md
H_GL_CORRELATION_IMPLEMENTATION.md  (superseded)
H_GL_SENSITIVITY_RESULTS.md
NEM_BOILING_FIX_RESULTS.md
NEM_BOILING_HEAT_TRANSFER_ISSUES.md
NEM_COMPARISON_RESULTS.md
NEM_CORRECTED_RESULTS.md
NEM_DISCHARGE_FIX.md
NEM_PSV_FIRE_RESULTS.md
NEM_STATUS.md
PHASE_TRANSFER_ANALYSIS.md
PHASE_TRANSFER_COMPARISON.md
PHASE_TRANSFER_ENERGY_ANALYSIS.md
PHASE_TRANSFER_LIMITS_ANALYSIS.md
PHASE_TRANSFER_VERIFICATION.md
RELAXATION_FACTOR_SENSITIVITY.md
solver_comparison.md
```

## Cleanup Options

### Option 1: ARCHIVE (Recommended)
Create `docs/nem_development/` directory and move all temporary files there for historical reference.

### Option 2: DELETE
Remove all temporary files completely (can recover from git if needed).

### Option 3: SELECTIVE
- Keep only the 3 core documentation files
- Archive diagnostic scripts (may be useful for future debugging)
- Delete plots (can regenerate if needed)
- Delete intermediate documentation

## Recommended Structure

```
HydDown/
├── docs/
│   ├── NEM_IMPLEMENTATION_METHODOLOGY.md  ← Keep
│   ├── NEM_HEAT_TRANSFER_COMPLETE_FIX.md  ← Keep
│   ├── H_GL_CORRELATION_FINAL.md          ← Keep
│   └── nem_development/                    ← Archive folder
│       ├── diagnostic_scripts/
│       └── intermediate_docs/
├── src/
│   └── hyddown/
│       ├── hdclass.py                      ← Modified (NEM implementation)
│       ├── transport.py                    ← Modified (h_gl correlation)
│       ├── validator.py                    ← Modified (h_gl validation)
│       └── examples/
│           └── nem_propane_psv_fire.yml   ← Modified (test case)
└── [original files unchanged]
```

## Summary

**Total files created during development:** 70
**Files to keep in production:** 3 core docs
**Files to clean up:** 67 temporary/diagnostic files

---

## ✅ CLEANUP EXECUTED (2026-03-26)

The following cleanup was performed:

### Actions Taken:

1. **Created archive directory structure:**
   ```
   docs/nem_development/
   ├── diagnostic_scripts/    (27 files)
   └── intermediate_docs/     (21 files)
   ```

2. **Moved diagnostic scripts to archive:**
   - 27 Python diagnostic/test scripts → `docs/nem_development/diagnostic_scripts/`
   - All scripts preserved for future debugging reference

3. **Moved intermediate documentation to archive:**
   - 21 intermediate markdown docs → `docs/nem_development/intermediate_docs/`
   - Historical documentation preserved

4. **Deleted diagnostic plots:**
   - ~22 PDF and PNG plot files removed
   - Can be regenerated from scripts if needed

5. **Core documentation kept in root:**
   - ✅ `NEM_IMPLEMENTATION_METHODOLOGY.md` (root directory)
   - ✅ `NEM_HEAT_TRANSFER_COMPLETE_FIX.md` (root directory)
   - ✅ `H_GL_CORRELATION_FINAL.md` (root directory)

### Result:

The repository root is now clean with only production-ready documentation. All development files are safely archived in `docs/nem_development/` and can be accessed if needed for future debugging or historical reference.

**Cleanup method used:** Selective archive (combination of Option 1 and Option 3)
