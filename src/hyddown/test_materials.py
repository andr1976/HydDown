# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

"""
Test suite for materials.py module.

Tests cover:
- von Mises stress calculations
- Ultimate Tensile Strength (UTS) lookups for all materials
- Allowable Tensile Strength (ATS) with safety factors
- Heat capacity (Cp) temperature dependencies
- Material property interpolation
- Edge cases and boundary conditions
"""

from hyddown import materials
import pytest
import numpy as np


class TestVonMisesStress:
    """Test von Mises stress calculations for pressure vessels."""

    def test_von_mises_basic_calculation(self):
        """Test basic von Mises stress calculation."""
        p = 100e5  # 100 bar
        d = 0.5  # 0.5 m inner diameter
        wt = 0.01  # 10 mm wall thickness
        stress = materials.von_mises(p, d, wt)

        # Verify stress is positive and reasonable
        assert stress > 0
        assert stress < 1e9  # Should be less than 1 GPa for typical cases

    def test_von_mises_higher_pressure(self):
        """Test von Mises stress increases with pressure."""
        d = 0.5
        wt = 0.01
        p_low = 50e5
        p_high = 150e5

        stress_low = materials.von_mises(p_low, d, wt)
        stress_high = materials.von_mises(p_high, d, wt)

        assert stress_high > stress_low

    def test_von_mises_thicker_wall(self):
        """Test von Mises stress decreases with wall thickness."""
        p = 100e5
        d = 0.5
        wt_thin = 0.005
        wt_thick = 0.020

        stress_thin = materials.von_mises(p, d, wt_thin)
        stress_thick = materials.von_mises(p, d, wt_thick)

        assert stress_thin > stress_thick

    def test_von_mises_with_custom_sigma_a(self):
        """Test von Mises calculation with custom axial stress."""
        p = 100e5
        d = 0.5
        wt = 0.01
        sigma_a_default = 30e6
        sigma_a_custom = 50e6

        stress_default = materials.von_mises(p, d, wt, sigma_a=sigma_a_default)
        stress_custom = materials.von_mises(p, d, wt, sigma_a=sigma_a_custom)

        # Higher axial stress should give higher equivalent stress
        assert stress_custom > stress_default


class TestUTS:
    """Test Ultimate Tensile Strength lookups."""

    def test_UTS_SS316_room_temperature(self):
        """Test UTS for SS316 at room temperature."""
        # Room temperature (20°C = 293.15 K)
        uts_room = materials.UTS(293.15, "SS316")
        assert uts_room == pytest.approx(575e6, rel=0.01)

    def test_UTS_SS316_elevated_temperature(self):
        """Test UTS for SS316 at elevated temperature."""
        # Test at 500°C (773.15 K)
        uts_high = materials.UTS(773.15, "SS316")
        uts_room = materials.UTS(293.15, "SS316")

        # Strength decreases with temperature
        assert uts_high < uts_room
        assert uts_high == pytest.approx(449e6, rel=0.05)

    def test_UTS_all_materials_positive(self):
        """Test UTS function returns positive values for all materials."""
        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_test = 373.15  # 100°C

        for mat in materials_list:
            uts = materials.UTS(T_test, mat)
            assert uts > 0
            assert uts < 1e9

    def test_UTS_temperature_degradation(self):
        """Test that UTS decreases with temperature for all materials."""
        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_low = 293.15  # 20°C
        T_high = 773.15  # 500°C

        for mat in materials_list:
            uts_low = materials.UTS(T_low, mat)
            uts_high = materials.UTS(T_high, mat)
            assert uts_high < uts_low

    def test_UTS_duplex_vs_carbon_steel(self):
        """Test that Duplex has higher UTS than carbon steel at room temp."""
        T = 293.15
        uts_duplex = materials.UTS(T, "Duplex")
        uts_cs235 = materials.UTS(T, "CS_235LT")
        uts_cs360 = materials.UTS(T, "CS_360LT")

        assert uts_duplex > uts_cs235
        assert uts_duplex > uts_cs360

    def test_UTS_CS360_higher_than_CS235(self):
        """Test that CS360LT has higher UTS than CS235LT."""
        T = 373.15  # 100°C
        uts_cs360 = materials.UTS(T, "CS_360LT")
        uts_cs235 = materials.UTS(T, "CS_235LT")

        assert uts_cs360 > uts_cs235

    def test_UTS_interpolation_consistency(self):
        """Test that interpolation gives smooth values."""
        mat = "SS316"
        # Test at intermediate temperature
        T_mid = 350  # Between tabulated values
        uts_mid = materials.UTS(T_mid, mat)

        # Should be between neighboring values
        T_low = 323.15  # 50°C
        T_high = 423.15  # 150°C
        uts_low = materials.UTS(T_low, mat)
        uts_high = materials.UTS(T_high, mat)

        assert uts_high < uts_mid < uts_low

    def test_UTS_extreme_temperatures(self):
        """Test UTS at extreme temperatures."""
        mat = "Duplex"

        # At minimum temperature (20°C)
        uts_min = materials.UTS(293.17, mat)
        assert uts_min == pytest.approx(730e6, rel=0.01)

        # At maximum temperature (1100°C)
        uts_max = materials.UTS(1373.17, mat)
        assert uts_max > 0
        assert uts_max < uts_min


class TestATS:
    """Test Allowable Tensile Strength calculations."""

    def test_ATS_default_safety_factor(self):
        """Test ATS with default safety factor k_s=0.85."""
        T = 373.15
        mat = "SS316"

        ats = materials.ATS(T, mat, k_s=0.85, k_y=1.0)
        uts = materials.UTS(T, mat)

        assert ats == pytest.approx(uts * 0.85, rel=0.01)

    def test_ATS_guaranteed_minimum(self):
        """Test ATS with guaranteed minimum (k_s=1.0)."""
        T = 373.15
        mat = "SS316"

        ats_default = materials.ATS(T, mat, k_s=0.85, k_y=1.0)
        ats_guaranteed = materials.ATS(T, mat, k_s=1.0, k_y=1.0)

        assert ats_guaranteed > ats_default
        assert ats_guaranteed == pytest.approx(materials.UTS(T, mat), rel=0.01)

    def test_ATS_uncertain_material_factor(self):
        """Test ATS with uncertain material data (k_y < 1.0)."""
        T = 373.15
        mat = "SS316"

        ats_normal = materials.ATS(T, mat, k_s=0.85, k_y=1.0)
        ats_uncertain = materials.ATS(T, mat, k_s=0.85, k_y=0.9)

        # Uncertain material should have lower allowable stress
        assert ats_uncertain < ats_normal

    def test_ATS_all_materials(self):
        """Test ATS for all material types."""
        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T = 473.15  # 200°C

        for mat in materials_list:
            ats = materials.ATS(T, mat)
            uts = materials.UTS(T, mat)
            # ATS should be less than UTS (with default k_s=0.85)
            assert ats < uts
            assert ats == pytest.approx(uts * 0.85, rel=0.01)


class TestSteelCp:
    """Test heat capacity lookups."""

    def test_steel_Cp_SS316_room_temperature(self):
        """Test heat capacity for SS316 at room temperature."""
        mat = "SS316"
        Cp = materials.steel_Cp(293.15, mat)  # 20°C
        assert Cp == pytest.approx(472, rel=0.01)

    def test_steel_Cp_temperature_dependence(self):
        """Test that Cp increases with temperature."""
        mat = "SS316"
        Cp_low = materials.steel_Cp(293.15, mat)  # 20°C
        Cp_high = materials.steel_Cp(773.15, mat)  # 500°C

        assert Cp_high > Cp_low

    def test_steel_Cp_all_materials(self):
        """Test heat capacity for all material types."""
        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_test = 473.15  # 200°C

        for mat in materials_list:
            cp = materials.steel_Cp(T_test, mat)
            assert cp > 400  # Reasonable lower bound for steel
            assert cp < 2000  # Reasonable upper bound

    def test_steel_Cp_carbon_steel_phase_transformation(self):
        """Test CS heat capacity discontinuity at phase transformation."""
        # Carbon steel has discontinuity at ~750°C due to austenite formation
        mat = "CS_235LT"

        # Test at temperatures around the phase transformation (750°C)
        Cp_700 = materials.steel_Cp(700 + 273.15, mat)  # 700°C
        Cp_750 = materials.steel_Cp(750 + 273.15, mat)  # 750°C

        # At 750°C there's a phase transformation, Cp jumps to 1450
        # At 700°C it's 900, so we should see a significant increase
        assert Cp_750 > Cp_700 * 1.4  # At least 40% increase
        assert Cp_750 > 1400  # Should be near the high value

    def test_steel_Cp_interpolation(self):
        """Test heat capacity interpolation between data points."""
        mat = "SS316"

        # Test at exact data points
        Cp_100 = materials.steel_Cp(373.15, mat)  # 100°C
        Cp_200 = materials.steel_Cp(473.15, mat)  # 200°C

        # Test at midpoint (should be between the two)
        Cp_150 = materials.steel_Cp(423.15, mat)  # 150°C

        assert Cp_100 < Cp_150 < Cp_200

    def test_steel_Cp_duplex_vs_SS316(self):
        """Test different materials have different heat capacities."""
        T = 473.15  # 200°C

        Cp_duplex = materials.steel_Cp(T, "Duplex")
        Cp_ss316 = materials.steel_Cp(T, "SS316")
        Cp_6mo = materials.steel_Cp(T, "6Mo")

        # All should be positive but different
        assert Cp_duplex > 0
        assert Cp_ss316 > 0
        assert Cp_6mo > 0
        # Duplex has higher Cp than SS316 at this temperature
        assert Cp_duplex > Cp_ss316

    def test_steel_Cp_extreme_temperatures(self):
        """Test heat capacity at extreme temperatures."""
        mat = "SS316"

        # At minimum temperature (20°C)
        Cp_min = materials.steel_Cp(293.15, mat)
        assert Cp_min > 0

        # At maximum temperature (1100°C)
        Cp_max = materials.steel_Cp(1373.15, mat)
        assert Cp_max > Cp_min

    def test_steel_Cp_6Mo_plateau(self):
        """Test that 6Mo Cp plateaus at high temperature."""
        mat = "6Mo"

        # 6Mo data shows plateau at 610 J/(kg·K) above ~750°C
        Cp_800 = materials.steel_Cp(1073.15, mat)  # 800°C
        Cp_1000 = materials.steel_Cp(1273.15, mat)  # 1000°C

        # Should be approximately the same (plateau)
        assert Cp_800 == pytest.approx(610, rel=0.05)
        assert Cp_1000 == pytest.approx(610, rel=0.05)


class TestDataConsistency:
    """Test consistency and validity of material property data."""

    def test_temperature_arrays_monotonic(self):
        """Test that temperature arrays are monotonically increasing."""
        # Check T_Cp array
        assert np.all(np.diff(materials.T_Cp) > 0)

        # Check T array
        assert np.all(np.diff(materials.T) > 0)

    def test_UTS_arrays_monotonic_decreasing(self):
        """Test that UTS decreases monotonically with temperature."""
        # UTS should decrease with temperature for all materials
        assert np.all(np.diff(materials.Duplex_UTS) <= 0)
        assert np.all(np.diff(materials.SS_UTS) <= 0)
        assert np.all(np.diff(materials.SMo_UTS) <= 0)
        assert np.all(np.diff(materials.CS_235LT_UTS) <= 0)
        assert np.all(np.diff(materials.CS_360LT_UTS) <= 0)

    def test_Cp_arrays_physically_reasonable(self):
        """Test that all Cp values are physically reasonable."""
        # All heat capacities should be between 400-2000 J/(kg·K) for steel
        all_cp_values = np.concatenate(
            [
                materials.SS316_Cp,
                materials.Duplex_Cp,
                materials.SMo_Cp,
                materials.CS_LT_Cp,
            ]
        )

        # Check all values are within reasonable range
        # Allow up to 1500 for CS phase transformation
        assert np.all(all_cp_values > 400)
        assert np.all(all_cp_values < 1500)

    def test_UTS_arrays_physically_reasonable(self):
        """Test that all UTS values are physically reasonable."""
        all_uts_values = np.concatenate(
            [
                materials.Duplex_UTS,
                materials.SS_UTS,
                materials.SMo_UTS,
                materials.CS_235LT_UTS,
                materials.CS_360LT_UTS,
            ]
        )

        # All UTS values should be positive and less than 1 GPa
        assert np.all(all_uts_values > 0)
        assert np.all(all_uts_values < 1e9)

    def test_temperature_array_lengths_match(self):
        """Test that temperature and property arrays have matching lengths."""
        # Cp data
        assert len(materials.T_Cp) == len(materials.SS316_Cp)
        assert len(materials.T_Cp) == len(materials.Duplex_Cp)
        assert len(materials.T_Cp) == len(materials.SMo_Cp)
        assert len(materials.T_Cp) == len(materials.CS_LT_Cp)

        # UTS data
        assert len(materials.T) == len(materials.Duplex_UTS)
        assert len(materials.T) == len(materials.SS_UTS)
        assert len(materials.T) == len(materials.SMo_UTS)
        assert len(materials.T) == len(materials.CS_235LT_UTS)
        assert len(materials.T) == len(materials.CS_360LT_UTS)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
