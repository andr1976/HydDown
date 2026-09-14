# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license
"""
Tests for the thermopack-free CoolProp CO2 release model (co2_release_cp.py).

Unlike test_co2_release.py these do NOT require thermopack - that is the whole point: the
CoolProp backend must run on a CoolProp + numpy install alone. The one cross-check against
thermopack GERG-2008 is skipped when thermopack is absent.
"""
import math
import os

import numpy as np
import pytest

from hyddown.co2_release_cp import CO2ReleaseModelCP


@pytest.fixture(scope="module")
def model():
    return CO2ReleaseModelCP()


def test_constructs_without_thermopack():
    """The CoolProp backend must import and build with thermopack blocked."""
    import builtins
    real_import = builtins.__import__

    def guard(name, *a, **k):
        if name.startswith("thermopack"):
            raise ImportError("thermopack must not be imported by the CoolProp backend")
        return real_import(name, *a, **k)

    builtins.__import__ = guard
    try:
        m = CO2ReleaseModelCP()
        assert m.eos is None
    finally:
        builtins.__import__ = real_import


def test_frost_point(model):
    # CoolProp / Span-Wagner sublimation temperature at 1 atm = -78.5 C
    assert model.T_frost == pytest.approx(194.66, abs=0.3)


def test_triple_point(model):
    assert model.T_TRIPLE_EOS == pytest.approx(216.59, abs=0.3)
    assert model.P_TRIPLE_EOS == pytest.approx(5.18e5, rel=0.02)


@pytest.mark.parametrize("phase,TC", [("liquid", 17.0), ("liquid", -30.0),
                                      ("gas", 17.0), ("gas", -30.0)])
def test_atm_split_sums_to_one(model, phase, TC):
    from CoolProp.CoolProp import PropsSI
    P0 = PropsSI("P", "T", 273.15 + TC, "Q", 0, "CO2")
    h0, s0, _r, _T = model.stagnation(P0, phase)
    a = model.atm_split(h0)
    assert a["vapour_frac"] + a["solid_frac"] == pytest.approx(1.0, abs=1e-9)
    ai = model.atm_split_isentropic(s0)
    assert ai["vapour_frac"] + ai["solid_frac"] == pytest.approx(1.0, abs=1e-9)


def test_liquid_release_makes_dry_ice(model):
    # a saturated-liquid release to atmosphere forms dry ice (isenthalpic)
    from CoolProp.CoolProp import PropsSI
    P0 = PropsSI("P", "T", 273.15 - 30.0, "Q", 0, "CO2")
    h0 = model.stagnation(P0, "liquid")[0]
    assert model.atm_split(h0)["solid_frac"] > 0.3


@pytest.mark.parametrize("P0bar", [50.0, 20.0, 10.0, 6.0])
def test_hem_rate_positive_and_choked(model, P0bar):
    area = math.pi * 0.02 ** 2 / 4
    r = model.hem_rate(P0bar * 1e5, "gas", 0.85, area)
    assert r["mdot"] > 0.0
    assert r["choked"] is True


def test_hem_rate_monotonic_in_pressure(model):
    area = math.pi * 0.02 ** 2 / 4
    rates = [model.hem_rate(P * 1e5, "gas", 0.85, area)["mdot"] for P in (8, 15, 30, 55)]
    assert all(x < y for x, y in zip(rates, rates[1:]))


def test_no_flow_below_backpressure(model):
    area = math.pi * 0.02 ** 2 / 4
    r = model.hem_rate(model.p_back * 0.5, "gas", 0.85, area)
    assert r["mdot"] == 0.0


def test_throat_carries_dry_ice_below_triple(model):
    # a low-pressure gas stagnation chokes below the triple point -> vapour+solid throat
    area = math.pi * 0.02 ** 2 / 4
    r = model.hem_rate(6.0e5, "gas", 0.85, area)
    assert r["solid_frac_throat"] > 0.0


def test_triple_lever_roundtrip(model):
    ms, ml, mg = 3.0, 10.0, 2.0
    V = ms * model.v_s + ml * model.v_l + mg * model.v_g
    U = ms * model.u_s + ml * model.u_l + mg * model.u_g
    r = model.triple_lever(ms + ml + mg, U, V)
    assert r["m_s"] == pytest.approx(ms, abs=1e-6)
    assert r["m_l"] == pytest.approx(ml, abs=1e-6)
    assert r["m_g"] == pytest.approx(mg, abs=1e-6)


def test_sublimation_state_roundtrip(model):
    T = 205.0
    P, vg, ug, vs, us = model._sublimation_props(T)
    ms, mg = 4.0, 1.5
    V = ms * vs + mg * vg
    U = ms * us + mg * ug
    st = model.sublimation_state(ms + mg, U, V)
    assert st["T"] == pytest.approx(T, abs=0.3)
    assert st["m_s"] == pytest.approx(ms, rel=0.02)


def test_matches_thermopack_gerg():
    """Cross-check: CoolProp-only reproduces thermopack GERG rates + dry-ice fractions."""
    pytest.importorskip("thermopack")
    from thermopack.multiparameter import multiparam
    from hyddown.co2_release import CO2ReleaseModel

    mg = CO2ReleaseModel.__new__(CO2ReleaseModel)
    mg.eos = multiparam("CO2", "GERG2008"); mg.eos.init_solid("CO2")
    mg.z = np.array([1.0]); mg.LIQ = mg.eos.LIQPH; mg.VAP = mg.eos.VAPPH
    mg.M = mg.eos.compmoleweight(1) / 1000.0
    mg.p_back = 1e5; mg.p_atm = 1e5
    mg.liquid_nonequilibrium = 0.0; mg.liquid_ne_pressure_scaled = False; mg.liquid_ne_pref = None
    mg.P_TRIPLE = 5.18e5
    mg._init_atm_endpoints(); mg._init_triple_point(); mg._init_gas_table()
    mc = CO2ReleaseModelCP(back_pressure=1e5, atm_pressure=1e5)

    area = math.pi * 0.008 ** 2 / 4
    for P0bar in (60, 20, 8, 6):
        for ph in ("gas", "liquid"):
            rg = mg.hem_rate(P0bar * 1e5, ph, 0.8, area)
            rc = mc.hem_rate(P0bar * 1e5, ph, 0.8, area)
            assert rc["mdot"] == pytest.approx(rg["mdot"], rel=0.01)
            assert rc["solid_frac_throat"] == pytest.approx(rg["solid_frac_throat"], abs=0.01)


def test_release_run_smoke_no_thermopack():
    """A full HydDown CO2 release run on the CoolProp backend, thermopack import blocked."""
    import builtins
    import yaml
    from hyddown import HydDown

    path = os.path.join(os.path.dirname(__file__), "examples", "co2_release_coolprop.yml")
    with open(path) as f:
        inp = yaml.safe_load(f)
    inp["calculation"]["end_time"] = 200  # keep the smoke test short

    real_import = builtins.__import__

    def guard(name, *a, **k):
        if name.startswith("thermopack"):
            raise ImportError("thermopack blocked")
        return real_import(name, *a, **k)

    builtins.__import__ = guard
    try:
        hd = HydDown(inp)
        hd.run()
    finally:
        builtins.__import__ = real_import
    assert type(hd.release_model).__name__ == "CO2ReleaseModelCP"
    assert hd.x_vap_atm[0] + hd.x_solid_atm[0] == pytest.approx(1.0, abs=1e-6)
    assert hd.T_atm[0] == pytest.approx(194.66, abs=0.6)
