# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

from hyddown import transport as tp
from CoolProp.CoolProp import PropsSI
import pytest


def test_orifice():
    P1 = 10.e5
    P2 = 5.5e5
    D = PropsSI('D', 'P', P1, 'T', 298.15, 'HEOS::N2')
    cpcv = PropsSI('CP0MOLAR', 'T', 298.15, 'P', P1, 'HEOS::N2') / PropsSI('CVMOLAR', 'T', 298.15, 'P', P1, 'HEOS::N2')
    assert tp.gas_release_rate(P1, P2, D, cpcv, 0.85, 0.01**2/4 * 3.1415) == pytest.approx(9.2 / 60, 0.001)


def test_controlvalve():
    P1 = 10.e5
    P2 = 5.5e5
    T1 = 20. + 273.15
    MW = PropsSI('M', 'P', P1, 'T', T1, 'HEOS::N2')
    Z1 = PropsSI('Z', 'P', P1, 'T', T1, 'HEOS::N2')
    gamma = PropsSI('CP0MOLAR', 'T', T1, 'P', P1, 'HEOS::N2') / PropsSI('CVMOLAR', 'T', T1, 'P', P1, 'HEOS::N2')
    assert tp.control_valve(P1, P2, T1, Z1, MW, gamma, 500) == pytest.approx(21.92,0.07)


def test_psv():
    P1 = 100e5
    Pback = 1e5
    Pset = 90e5
    blowdown = 0.1
    T1 = 27. + 273.15
    rho = PropsSI('D', 'P', P1, 'T', T1, 'HEOS::N2')
    gamma = PropsSI('CP0MOLAR', 'T', T1, 'P', P1, 'HEOS::N2') / PropsSI('CVMOLAR', 'T', T1, 'P', P1, 'HEOS::N2')
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(P1, Pback, Pset, blowdown, rho, gamma, CD, area) == pytest.approx(1.57, 0.02)


def test_hinner():
    h = tp.h_inner(0.305, 311, 505.4, 1e5, 'HEOS::air')
    assert h == pytest.approx(7,0.1)


def test_NNu():
    NGr = tp.Gr(0.305, 311, 505.4, 1e5, 'HEOS::air')
    NPr = tp.Pr((311 + 505.4) / 2, 1e5, 'HEOS::air')
    NRa = NGr * NPr
    NNu = tp.Nu(NRa, NPr)
    print(NNu)
    assert NNu == pytest.approx(7.5e7, 0.1e7)


def test_NRa():
    NGr = tp.Gr(0.305, 311, 505.4, 1e5, 'HEOS::air')
    NPr = tp.Pr((311 + 505.4) / 2, 1e5, 'HEOS::air')
    assert (NGr * NPr) == pytest.approx(1.3e8, 0.1e8)


def test_NGr():
    assert tp.Gr(0.305, 311, 505.4, 1e5, 'HEOS::air') == pytest.approx(1.8e8, 0.1e8)

if __name__ == "__main__":
    test_NGr()