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

import numpy as np
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


# --------------------------------------------------------------------------------------
# Solid-in-vessel fallback (below the triple point)
# --------------------------------------------------------------------------------------

def test_triple_point_self_consistent(model):
    # tcPR triple point near the literature 216.6 K / 5.18 bar
    assert model.T_TRIPLE_EOS == pytest.approx(216.6, abs=0.5)
    assert model.P_TRIPLE_EOS == pytest.approx(5.2e5, rel=0.05)


def test_triple_lever_roundtrip(model):
    # build (M, U, V) from a chosen phase split, recover it
    ms, ml, mg = 300.0, 5000.0, 400.0
    M = ms + ml + mg
    U = ms * model.u_s + ml * model.u_l + mg * model.u_g
    V = ms * model.v_s + ml * model.v_l + mg * model.v_g
    r = model.triple_lever(M, U, V)
    assert r["m_s"] == pytest.approx(ms, abs=1e-6)
    assert r["m_l"] == pytest.approx(ml, abs=1e-6)
    assert r["m_g"] == pytest.approx(mg, abs=1e-6)


def test_sublimation_state_roundtrip(model):
    P, vg, ug, vs, us = model._sublimation_props(205.0)
    ms, mg = 200.0, 150.0
    M, U, V = ms + mg, ms * us + mg * ug, ms * vs + mg * vg
    st = model.sublimation_state(M, U, V)
    assert st["T"] == pytest.approx(205.0, abs=0.2)
    assert st["m_s"] == pytest.approx(ms, rel=0.02)
    assert st["m_g"] == pytest.approx(mg, rel=0.02)
    assert st["P"] < model.P_TRIPLE_EOS  # below the triple point


def test_triple_LG_from_MV_volume(model):
    M, V = 300.0, 20.0
    m_l, m_g, U = model.triple_LG_from_MV(M, V)
    assert m_l + m_g == pytest.approx(M, abs=1e-9)
    assert m_l * model.v_l + m_g * model.v_g == pytest.approx(V, rel=1e-6)
    assert U == pytest.approx(m_l * model.u_l + m_g * model.u_g)


def test_dispatcher_regimes(model):
    V = 20.0
    # three-phase: an interior (v,u) point -> "triple"
    ms, ml, mg = 100.0, 100.0, 100.0
    M = ms + ml + mg
    st = model.vessel_state_below_triple(
        M, ms * model.u_s + ml * model.u_l + mg * model.u_g,
        ms * model.v_s + ml * model.v_l + mg * model.v_g)
    assert st["regime"] == "triple"
    # solid+gas (colder, no liquid) -> "sublimation"
    P, vg, ug, vs, us = model._sublimation_props(200.0)
    ms2, mg2 = 150.0, 120.0
    st2 = model.vessel_state_below_triple(ms2 + mg2, ms2 * us + mg2 * ug, ms2 * vs + mg2 * vg)
    assert st2["regime"] == "sublimation"


def test_solid_in_vessel_run_accumulates_dry_ice():
    from hyddown import HydDown

    inp = get_example_input("co2_solid_in_vessel.yml")
    inp["calculation"]["end_time"] = 1500  # reaches the triple point (~700 s) and forms ice
    hd = HydDown(inp)
    hd.run(disable_pbar=True)  # must not raise

    assert hd.solid_regime is True
    # dry ice accumulated inside the vessel
    assert hd.m_solid.max() > 100.0
    # the tank descended below the (frozen) triple point instead of freezing there
    assert hd.P[hd.P > 0].min() < model_triple_P() * 1.02
    # solid-regime mass conservation: M[handoff] - leaked == M[end]
    solid = np.where(hd.m_solid > 0)[0]
    h0 = solid[0] - 1
    last = np.where(hd.mass_fluid > 0)[0][-1]
    leaked = float(np.sum(hd.release_rate[h0 + 1:last + 1]) * hd.tstep)
    assert hd.mass_fluid[h0] - leaked == pytest.approx(hd.mass_fluid[last], rel=1e-3)


def test_default_freezes_without_solid_in_vessel():
    from hyddown import HydDown

    inp = get_example_input("co2_gas_release.yml")
    inp["calculation"]["end_time"] = 1000  # long enough to hit the floor
    hd = HydDown(inp)
    hd.run(disable_pbar=True)
    # default behaviour: freeze at the floor, never enter the solid regime
    assert getattr(hd, "solid_regime", False) is False
    assert hd.release_frozen is True
    assert hd.m_solid.max() == 0.0
    assert hd.P[hd.P > 0].min() > P_TRIPLE  # never crossed below the triple point


def model_frost():
    return CO2ReleaseModel().T_frost


def model_triple_P():
    return CO2ReleaseModel().P_TRIPLE_EOS
