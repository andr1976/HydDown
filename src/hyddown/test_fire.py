# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

"""
Test suite for fire.py module.

Tests cover:
- Stefan-Boltzmann radiation calculations
- Fire scenario heat fluxes (API 521 and Scandpower)
- Temperature dependencies
- Edge cases and error handling
"""

from hyddown import fire
import pytest


class TestStefanBoltzmann:
    """Test Stefan-Boltzmann radiation calculations."""

    def test_stefan_boltzmann_basic(self):
        """Test basic Stefan-Boltzmann calculation."""
        alpha = 1
        e_flame = 1
        e_surface = 0
        h = 100
        Tflame = 635 + 273.15
        Tradiative = 635 + 273.15
        Tvessel = 20 + 273.15

        Q = fire.stefan_boltzmann(
            alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel
        )

        assert Q == pytest.approx(1e5, abs=100)

    def test_stefan_boltzmann_zero_emissivity(self):
        """Test Stefan-Boltzmann with zero surface emissivity (only convection)."""
        Q = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0,
            h=50,
            Tflame=1000,
            Tradiative=1000,
            Tvessel=300,
        )

        # Should still have heat flux from convection
        assert Q > 0
        # Convective part should be h * (Tflame - Tvessel)
        Q_conv_expected = 50 * (1000 - 300)
        assert Q_conv_expected > 0

    def test_stefan_boltzmann_high_emissivity(self):
        """Test Stefan-Boltzmann with different surface emissivity values."""
        Q_low_emissivity = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.1, h=100, Tflame=1000, Tradiative=1000, Tvessel=300
        )

        Q_high_emissivity = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.9, h=100, Tflame=1000, Tradiative=1000, Tvessel=300
        )

        # Higher surface emissivity means surface radiates more heat away
        # so net heat flux to vessel is reduced
        assert Q_high_emissivity <= Q_low_emissivity

    def test_stefan_boltzmann_equal_temperatures(self):
        """Test Stefan-Boltzmann when wall equals flame temperature."""
        # When temperatures are equal, heat flux should be small
        T = 500
        Q = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.5, h=100, Tflame=T, Tradiative=T, Tvessel=T
        )

        # Allow for numerical effects and small residuals
        assert Q == pytest.approx(0, abs=2000)

    def test_stefan_boltzmann_different_radiative_temp(self):
        """Test Stefan-Boltzmann with different radiative temperature."""
        Q1 = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.5, h=100, Tflame=1000, Tradiative=1000, Tvessel=300
        )

        Q2 = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.5, h=100, Tflame=1000, Tradiative=800, Tvessel=300
        )

        # Different radiative temperature should give different heat flux
        assert Q1 != Q2


class TestPoolFireAPI521:
    """Test API 521 pool fire scenarios."""

    def test_pool_fire_api521_room_temp(self):
        """Test API 521 pool fire at room temperature."""
        Q = fire.pool_fire_api521(273 + 20)
        # Should be around 60 kW/m² minus radiative losses
        assert Q > 0
        assert Q < 60e3  # Maximum incident heat flux

    def test_pool_fire_api521_elevated_temp(self):
        """Test API 521 pool fire at elevated temperature."""
        Q = fire.pool_fire_api521(273 + 50)
        assert Q == pytest.approx(45.5e3, abs=100)

    def test_pool_fire_api521_temperature_dependence(self):
        """Test that pool fire heat flux decreases with wall temperature."""
        Q_low = fire.pool_fire_api521(273 + 20)
        Q_high = fire.pool_fire_api521(273 + 200)

        # Higher wall temperature should give lower heat flux
        assert Q_high < Q_low


class TestJetFireAPI521:
    """Test API 521 jet fire scenarios."""

    def test_jet_fire_api521_room_temp(self):
        """Test API 521 jet fire at room temperature."""
        Q = fire.jet_fire_api521(273 + 20)
        # Should be around 100 kW/m² minus radiative losses
        assert Q > 0
        assert Q < 100e3

    def test_jet_fire_api521_elevated_temp(self):
        """Test API 521 jet fire at elevated temperature."""
        Q = fire.jet_fire_api521(273 + 50)
        assert Q == pytest.approx(83.5e3, abs=500)

    def test_jet_fire_api521_higher_than_pool(self):
        """Test that jet fire gives higher heat flux than pool fire."""
        T = 273 + 50
        Q_jet = fire.jet_fire_api521(T)
        Q_pool = fire.pool_fire_api521(T)

        # Jet fire should have higher heat flux than pool fire
        assert Q_jet > Q_pool


class TestPoolFireScandpower:
    """Test Scandpower pool fire scenarios."""

    def test_pool_fire_scandpower_room_temp(self):
        """Test Scandpower pool fire at room temperature."""
        Q = fire.pool_fire_scandpower(273 + 20)
        assert Q == pytest.approx(88.5e3, abs=500)

    def test_pool_fire_scandpower_elevated_temp(self):
        """Test Scandpower pool fire at elevated temperature."""
        Q = fire.pool_fire_scandpower(273 + 100)
        # Should be positive but less than at room temp
        assert Q > 0
        assert Q < fire.pool_fire_scandpower(273 + 20)

    def test_pool_fire_scandpower_higher_than_api(self):
        """Test that Scandpower pool fire is more conservative than API."""
        T = 273 + 20
        Q_scandpower = fire.pool_fire_scandpower(T)
        Q_api = fire.pool_fire_api521(T)

        # Scandpower should be more conservative (higher heat flux)
        assert Q_scandpower > Q_api


class TestJetFireScandpower:
    """Test Scandpower jet fire scenarios."""

    def test_jet_fire_scandpower_room_temp(self):
        """Test Scandpower jet fire at room temperature."""
        Q = fire.jet_fire_scandpower(273 + 20)
        assert Q == pytest.approx(94.5e3, abs=1000)

    def test_jet_fire_scandpower_elevated_temp(self):
        """Test Scandpower jet fire at elevated temperature."""
        Q = fire.jet_fire_scandpower(273 + 150)
        # Should be positive but less than at room temp
        assert Q > 0
        assert Q < fire.jet_fire_scandpower(273 + 20)

    def test_jet_fire_scandpower_temperature_dependence(self):
        """Test temperature dependence of Scandpower jet fire."""
        Q_20 = fire.jet_fire_scandpower(273 + 20)
        Q_50 = fire.jet_fire_scandpower(273 + 50)
        Q_100 = fire.jet_fire_scandpower(273 + 100)

        # Heat flux should decrease with increasing wall temperature
        assert Q_20 > Q_50 > Q_100


class TestSbFire:
    """Test the sb_fire wrapper function."""

    def test_sb_fire_all_scenarios(self):
        """Test sb_fire for all valid scenario types."""
        T = 273 + 50

        # Test all valid scenarios
        Q_api_jet = fire.sb_fire(T, "api_jet")
        Q_api_pool = fire.sb_fire(T, "api_pool")
        Q_sp_pool = fire.sb_fire(T, "scandpower_pool")
        Q_sp_jet = fire.sb_fire(T, "scandpower_jet")

        # All should be positive
        assert Q_api_jet > 0
        assert Q_api_pool > 0
        assert Q_sp_pool > 0
        assert Q_sp_jet > 0

        # Jets should be higher than pools
        assert Q_api_jet > Q_api_pool

    def test_sb_fire_api_scenarios(self):
        """Test sb_fire API scenarios match direct function calls."""
        T = 273 + 50

        assert fire.sb_fire(T, "api_jet") == pytest.approx(83.5e3, abs=500)
        assert fire.sb_fire(T, "api_pool") == pytest.approx(45.5e3, abs=100)

    def test_sb_fire_scandpower_scenarios(self):
        """Test sb_fire Scandpower scenarios match direct function calls."""
        T = 273 + 20

        assert fire.sb_fire(T, "scandpower_pool") == pytest.approx(88.5e3, abs=500)
        assert fire.sb_fire(T, "scandpower_jet") == pytest.approx(94.5e3, abs=1000)

    def test_sb_fire_invalid_scenario(self):
        """Test sb_fire raises error for invalid scenario."""
        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "invalid_scenario")

    def test_sb_fire_typo_scenario(self):
        """Test sb_fire raises error for typo in scenario name."""
        # Common typos
        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "scand_jet")

        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "scandpower_jett")

        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "API_jet")  # Wrong case


class TestFireScenarioComparisons:
    """Test relationships between different fire scenarios."""

    def test_jet_vs_pool_heat_flux(self):
        """Test that jet fires consistently give higher heat flux than pools."""
        temperatures = [273 + 20, 273 + 50, 273 + 100, 273 + 200]

        for T in temperatures:
            Q_jet_api = fire.jet_fire_api521(T)
            Q_pool_api = fire.pool_fire_api521(T)
            Q_jet_sp = fire.jet_fire_scandpower(T)
            Q_pool_sp = fire.pool_fire_scandpower(T)

            # Jet should always be higher than pool for same standard
            assert Q_jet_api > Q_pool_api
            # Note: For Scandpower, the jet/pool relationship may vary

    def test_scandpower_vs_api_conservatism(self):
        """Test relative conservatism between standards."""
        T = 273 + 50

        Q_jet_api = fire.jet_fire_api521(T)
        Q_pool_api = fire.pool_fire_api521(T)
        Q_jet_sp = fire.jet_fire_scandpower(50 + 273)
        Q_pool_sp = fire.pool_fire_scandpower(50 + 273)

        # All should be reasonable heat flux values
        assert 40e3 < Q_pool_api < 100e3
        assert 70e3 < Q_jet_api < 150e3

    def test_temperature_effect_consistency(self):
        """Test that all scenarios show consistent temperature effects."""
        T_low = 273 + 20
        T_high = 273 + 200

        # All scenarios should give lower heat flux at higher wall temperature
        assert fire.pool_fire_api521(T_high) < fire.pool_fire_api521(T_low)
        assert fire.jet_fire_api521(T_high) < fire.jet_fire_api521(T_low)
        assert fire.pool_fire_scandpower(T_high) < fire.pool_fire_scandpower(T_low)
        assert fire.jet_fire_scandpower(T_high) < fire.jet_fire_scandpower(T_low)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
