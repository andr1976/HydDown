# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Tests for the thermopack CO2 release model (co2_release.py) and its HydDown integration.

Covers:
- the atmospheric dry-ice split reproducing the DRY_ICE_HANDOVER.md table (CI drift check),
- the isentrope lever carrying dry ice below the CO2 triple point (the "dry ice at the
  throat" mechanism CoolProp cannot evaluate),
- the HEM release rate (positive, choked, throat below the stagnation pressure), and
- an end-to-end HydDown run populating the atmospheric-state arrays with mass leaving the
  correct (liquid) inventory.
"""

import math
import os

import pytest

# The whole module needs thermopack; skip cleanly if it is not installed.
pytest.importorskip("thermopack")

from hyddown.co2_release import CO2ReleaseModel, P_TRIPLE, _scal


def get_example_input(fname):
    import yaml

    fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "examples", fname)
    with open(fname) as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


@pytest.fixture(scope="module")
def model():
    return CO2ReleaseModel()


def test_frost_point(model):
    # Gas<->solid Gibbs-equality frost point at 1 atm (DRY_ICE_HANDOVER.md sec. 5.3)
    assert model.T_frost == pytest.approx(194.14, abs=0.5)


# DRY_ICE_HANDOVER.md sec. 3 - live-tcPR dry-ice mass fractions at 1 atm.
# (release phase, T [degC], isenthalpic, isentropic)
HANDOVER_TABLE = [
    ("liquid", 17.0, 0.29, 0.41),
    ("liquid", -30.0, 0.50, 0.56),
    ("gas", 17.0, 0.01, 0.22),
    ("gas", -30.0, 0.00, 0.13),
]


@pytest.mark.parametrize("phase,TC,isenth,isentr", HANDOVER_TABLE)
def test_atm_split_matches_handover(model, phase, TC, isenth, isentr):
    T = 273.15 + TC
    P0 = _scal(model.eos.bubble_pressure(T, model.z))
    h0, s0, _rho0, _T0 = model.stagnation(P0, phase)
    assert model.atm_split(h0)["solid_frac"] == pytest.approx(isenth, abs=0.03)
    assert model.atm_split_isentropic(s0)["solid_frac"] == pytest.approx(isentr, abs=0.03)


def test_atm_split_fractions_sum_to_one(model):
    h0, _s0, _r, _T = model.stagnation(20e5, "liquid")
    split = model.atm_split(h0)
    assert split["vapour_frac"] + split["solid_frac"] == pytest.approx(1.0, abs=1e-9)
    assert 0.0 <= split["solid_frac"] <= 1.0


def test_isentrope_lever_carries_dry_ice_below_triple_point(model):
    # This is the state CoolProp cannot evaluate. Along the sat-liquid isentrope the
    # sub-triple-point states must be vapour + solid with a physical solid fraction and
    # finite enthalpy/density (no crash).
    _h0, s0, _rho0, _T0 = model.stagnation(20e5, "liquid")
    for P in (5.0e5, 3.0e5, 2.0e5, 1.01325e5):
        assert P < P_TRIPLE  # all below the CO2 triple point (5.18 bar)
        h_kg, rho, T, solid = model._iso_props(P, s0)
        assert 0.0 < solid < 1.0
        assert rho > 0.0
        assert math.isfinite(h_kg)
        assert T < 217.0  # below the triple-point temperature


@pytest.mark.parametrize("P0bar", [6.0, 10.0, 20.0, 50.0])
def test_hem_rate_positive_and_choked(model, P0bar):
    area = math.pi * 0.02 ** 2 / 4
    r = model.hem_rate(P0bar * 1e5, "liquid", 0.62, area)
    assert r["mdot"] > 0.0
    assert r["choked"]
    assert model.p_back < r["P_throat"] < P0bar * 1e5
    assert 0.0 <= r["solid_frac_throat"] <= 1.0


def test_hem_rate_monotonic_in_pressure(model):
    # Higher tank pressure -> higher release rate (no solver artefacts across the kink).
    area = math.pi * 0.02 ** 2 / 4
    rates = [model.hem_rate(P * 1e5, "liquid", 0.62, area)["mdot"] for P in (6, 8, 12, 20, 40)]
    assert all(b > a for a, b in zip(rates, rates[1:]))


def test_hem_rate_no_flow_when_below_backpressure(model):
    area = math.pi * 0.02 ** 2 / 4
    r = model.hem_rate(0.9 * model.p_back, "gas", 0.62, area)
    assert r["mdot"] == 0.0


def test_release_run_smoke():
    from hyddown import HydDown

    inp = get_example_input("co2_release.yml")
    inp["calculation"]["end_time"] = 60
    inp["calculation"]["time_step"] = 5
    hd = HydDown(inp)
    hd.run(disable_pbar=True)

    n = int(inp["calculation"]["end_time"] / inp["calculation"]["time_step"])
    # Release rate is positive and drives the mass balance
    assert hd.release_rate[0] > 0.0
    # A liquid release draws from the liquid inventory (gas roughly unchanged)
    assert hd.m_liquid[n - 1] < hd.m_liquid[0]
    # Atmospheric state populated: dry ice forms at the frost point
    assert hd.T_atm[0] == pytest.approx(model_frost(), abs=0.5)
    assert hd.x_solid_atm[0] > 0.0
    assert hd.x_vap_atm[0] + hd.x_solid_atm[0] == pytest.approx(1.0, abs=1e-6)
    # Cumulative dry-ice mass is accumulating
    assert hd.m_dryice_cum[n - 1] > 0.0

    # Output plumbing
    df = hd.get_dataframe()
    assert "Cumulative dry-ice mass (kg)" in df.columns
    assert "Atmospheric dry-ice mass fraction (-)" in df.columns
    hd.generate_report()
    assert hd.report["total_dryice_mass"] > 0.0
    assert hd.report["max_release_rate"] > 0.0


def test_release_liquid_to_gas_switch_and_triple_point_floor():
    # The shipped example drains the liquid, switches the release to vapour, and then
    # approaches the CO2 triple point - which must freeze the release gracefully rather
    # than push the vessel into the (unsupported) solid-in-vessel regime and crash.
    from hyddown import HydDown

    inp = get_example_input("co2_release.yml")
    inp["calculation"]["end_time"] = 1500  # far enough to reach the floor (~800 s)
    hd = HydDown(inp)
    hd.run(disable_pbar=True)  # must not raise

    # Liquid was exhausted and the release switched to vapour
    assert hd.m_liquid.min() <= 1e-3
    assert hd.release_phase == "gas"
    # The triple-point floor guard fired and the tank never crossed below the triple point
    assert hd.release_frozen is True
    assert hd.P[hd.P > 0].min() > P_TRIPLE
    # Dry ice was produced both at the throat (near the floor) and at the atmosphere
    assert hd.solid_frac_throat.max() > 0.0
    assert hd.x_solid_atm.max() > 0.0


def model_frost():
    return CO2ReleaseModel().T_frost
