# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license
"""Unit tests for the 2-D axisymmetric conjugate wall solver (:mod:`hyddown.wall2d`)."""
import numpy as np
import pytest

from hyddown import wall2d as w2

RHO, CP, T0 = 7950.0, 500.0, 298.15


def _steel_volume_analytic(g):
    r_in, r_out, r_lid, H = g["r_in"], g["r_out"], g["r_lid"], g["H"]
    v_bot = np.pi * r_out ** 2 * g["t_bot"]
    v_shell = np.pi * (r_out ** 2 - r_in ** 2) * H
    v_flange = np.pi * (r_lid ** 2 - r_in ** 2) * g["t_flange"]
    v_lid = np.pi * r_lid ** 2 * g["t_lid"]
    return v_bot + v_shell + v_flange + v_lid


def test_geometry_volume_matches_analytic():
    """The masked FV grid must reproduce the analytic steel volume of the stepped domain."""
    g = w2.default_sintef_geometry()
    W = w2.WallConduction2D(g, RHO, CP, T0)
    assert W.vol_u.sum() == pytest.approx(_steel_volume_analytic(g), rel=1e-6)


def test_adiabatic_energy_conserved():
    """With both HTCs zero (fully adiabatic) the field is frozen and energy is conserved."""
    W = w2.WallConduction2D(w2.default_sintef_geometry(), RHO, CP, T0)
    E0 = W.energy()
    for _ in range(50):
        W.step(1.0, liquid_level=0.5, liquid_present=True,
               h_gas=0.0, T_gas=200.0, h_liq=0.0, T_liq=150.0)
    assert W.T.max() - W.T.min() == pytest.approx(0.0, abs=1e-9)
    assert W.energy() == pytest.approx(E0, rel=1e-12)


def test_steady_state_approaches_fluid_temperature():
    """A uniform inner Robin BC to one fluid temperature, adiabatic outside, drives the whole
    wall to that fluid temperature at steady state (Laplace with adiabatic outer + Robin inner)."""
    W = w2.WallConduction2D(w2.default_sintef_geometry(), RHO, CP, T0)
    Tf = 250.0
    for _ in range(20000):
        W.step(5.0, liquid_level=2.0, liquid_present=True,
               h_gas=500.0, T_gas=Tf, h_liq=500.0, T_liq=Tf)
    assert W.T.mean() == pytest.approx(Tf, abs=0.2)
    assert W.T.max() == pytest.approx(Tf, abs=0.5)


def test_split_level_gradient_signs():
    """A cold boiling liquid below and a mild warm gas above must give a wetted wall colder than
    the dry wall, and a through-wall gradient (outer warmer than inner) on the wetted side."""
    W = w2.WallConduction2D(w2.default_sintef_geometry(), RHO, CP, T0)
    for _ in range(120):
        s = W.step(1.0, liquid_level=0.30, liquid_present=True,
                   h_gas=15.0, T_gas=273.0, h_liq=800.0, T_liq=233.0)
    assert s["T_inner_wet"] < s["T_inner_dry"]          # wetted wall colder
    assert s["T_outer_wet"] > s["T_inner_wet"]          # radial gradient (outer lags)
    assert s["Q_liq"] > s["Q_gas"] > 0                  # boiling extracts more than gas
    assert 233.0 < s["T_inner_wet"] < 273.0             # bracketed by the two fluid temps


def test_cylinder_report_keys():
    """The cylinder-only report (thermocouple-comparable) is ordered and colder on the wetted
    side, and lies within the full-domain dry/wetted averages' range."""
    W = w2.WallConduction2D(w2.default_sintef_geometry(), RHO, CP, T0)
    for _ in range(120):
        s = W.step(1.0, liquid_level=0.30, liquid_present=True,
                   h_gas=15.0, T_gas=273.0, h_liq=800.0, T_liq=233.0)
    assert s["T_cyl_in_min"] <= s["T_cyl_in_med"] <= s["T_cyl_in_max"]
    assert s["T_cyl_out_min"] <= s["T_cyl_out_med"] <= s["T_cyl_out_max"]
    # wetted-region cylinder wall colder than the dry-region cylinder wall
    assert s["T_cyl_in_wet"] < s["T_cyl_in_dry"]
    # cylinder inner report excludes the warm lid/flange, so it is colder than the full-domain
    # dry average (which the lid/flange bias upward)
    assert s["T_cyl_in_dry"] <= s["T_inner_dry"] + 1e-6


def test_flat_end_geometry_reduces_to_cylinder():
    """With zero bottom/flange/lid the domain is a plain cylindrical shell."""
    g = w2.build_geometry(r_in=0.1365, thickness=0.0254, H=1.0,
                          t_bot=0.0, t_flange=0.0, t_lid=0.0, r_lid=0.0)
    W = w2.WallConduction2D(g, RHO, CP, T0)
    v_shell = np.pi * ((0.1365 + 0.0254) ** 2 - 0.1365 ** 2) * 1.0
    assert W.vol_u.sum() == pytest.approx(v_shell, rel=1e-6)
