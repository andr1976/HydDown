# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

from hyddown import transport as tp
from hyddown import fire
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


def test_stefan_boltzmann():
    alpha = 1
    e_flame = 1
    e_surface = 0
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    assert (fire.stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, 20+273.15)) == pytest.approx(1e5,100)


def test_pool_fire_api521():
    assert fire.pool_fire_api521(273+50) == pytest.approx(45.5e3, 100)


def test_jet_fire_api521():
    assert fire.jet_fire_api521(273+50) == pytest.approx(83.5e3, 500)


def test_jet_fire_scandpower():
    assert fire.jet_fire_scandpower(273+20) == pytest.approx(94.5e3,100)


def test_pool_fire_scandpower():
    assert fire.pool_fire_scandpower(273+20) == pytest.approx(88.5e3,50)



if __name__ == "__main__":
    test_NGr()