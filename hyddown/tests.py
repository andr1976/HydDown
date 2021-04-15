# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

from hyddown import transport as tp
from hyddown import fire
from hyddown import validator
from CoolProp.CoolProp import PropsSI
import pytest


def get_example_input(fname):
    import os
    import yaml

    if "C:\\Users\\ANRA" in os.getcwd():
        fname = r"C:\\Users\\ANRA\\Documents\\GitHub\\HydDown\\examples\\" + fname
    else:
        fname = r"//home//runner//work//HydDown//HydDown//examples//" + fname

    with open(fname) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)

    return input


def test_orifice():
    P1 = 10.0e5
    P2 = 5.5e5
    D = PropsSI("D", "P", P1, "T", 298.15, "HEOS::N2")
    cpcv = PropsSI("CP0MOLAR", "T", 298.15, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", 298.15, "P", P1, "HEOS::N2"
    )
    assert tp.gas_release_rate(
        P1, P2, D, cpcv, 0.85, 0.01 ** 2 / 4 * 3.1415
    ) == pytest.approx(9.2 / 60, rel=0.001)


def test_orifice1():
    P1 = 10.0e5
    P2 = 6.5e5
    D = PropsSI("D", "P", P1, "T", 298.15, "HEOS::N2")
    cpcv = PropsSI("CP0MOLAR", "T", 298.15, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", 298.15, "P", P1, "HEOS::N2"
    )
    assert tp.gas_release_rate(
        P1, P2, D, cpcv, 0.85, 0.01 ** 2 / 4 * 3.1415
    ) == pytest.approx(9.2 / 60, rel=0.2)


def test_controlvalve():
    P1 = 10.0e5
    P2 = 5.5e5
    T1 = 20.0 + 273.15
    MW = PropsSI("M", "P", P1, "T", T1, "HEOS::N2")
    Z1 = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    assert tp.control_valve(P1, P2, T1, Z1, MW, gamma, 500) == pytest.approx(
        21.92, rel=0.05
    )

def test_psv3():
    Pback = 1e5
    Pset = 18.2e5
    blowdown = 0.1
    P1 = 0.99 * Pset * (1 - blowdown)
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M",  "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == 0
    psv_state = "closed"
    assert tp.relief_valve(
        P1*1.01, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == 0


def test_psv2():
    P1 = 21.0e5
    Pback = 1e5
    Pset = 20.99e5
    blowdown = 0.1
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M",  "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1046 / 3600, rel=0.03)
    psv_state = "open"
    assert tp.relief_valve(
        Pset*0.99, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1046 / 3600, rel=0.03)


def test_psv():
    P1 = 100e5
    Pback = 1e5
    Pset = 99.2e5
    blowdown = 0.1
    T1 = 25.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M",  "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1.57, rel=0.02)


def test_api_psv_relief():
    assert tp.api_psv_release_rate(121.9e5, 71e5, 1.39, 0.975, 298.15, 1.01, 2/1e3, 71e-6) == pytest.approx(1846/3600, rel=0.01)
    assert tp.api_psv_release_rate(121.9e5, 1e5, 1.39, 0.975, 298.15, 1.01, 2/1e3, 71e-6) == pytest.approx(1860/3600, rel=0.01)


def test_hinner():
    h = tp.h_inner(0.305, 311, 505.4, 1e5, "HEOS::air")
    assert h == pytest.approx(7, abs=0.1)


def test_hinner_mixed():
    mdot = 1e-10
    D = 0.010
    h = tp.h_inner_mixed(0.305, 311, 505.4, 1e5, "HEOS::air", mdot, D)
    assert h == pytest.approx(7, abs=0.1)


def test_NNu():
    NGr = tp.Gr(0.305, 311, 505.4, 1e5, "HEOS::air")
    NPr = tp.Pr((311 + 505.4) / 2, 1e5, "HEOS::air")
    NRa = NGr * NPr
    NNu = tp.Nu(NRa, NPr)
    print(NNu)
    assert NNu == pytest.approx(62.2, abs=1.0)


def test_NRa():
    NGr = tp.Gr(0.305, 311, 505.4, 1e5, "HEOS::air")
    NPr = tp.Pr((311 + 505.4) / 2, 1e5, "HEOS::air")
    assert (NGr * NPr) == pytest.approx(1.3e8, abs=0.1e8)


def test_NGr():
    assert tp.Gr(0.305, 311, 505.4, 1e5, "HEOS::air") == pytest.approx(1.8e8, abs=0.1e8)


def test_NPr():
    assert tp.Pr((311 + 505.4) / 2, 1e5, "HEOS::air") == pytest.approx(0.7, rel=0.01)


def test_stefan_boltzmann():
    alpha = 1
    e_flame = 1
    e_surface = 0
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    assert (
        fire.stefan_boltzmann(
            alpha, e_flame, e_surface, h, Tflame, Tradiative, 20 + 273.15
        )
    ) == pytest.approx(1e5, abs=100)


def test_pool_fire_api521():
    assert fire.pool_fire_api521(273 + 50) == pytest.approx(45.5e3, abs=100)


def test_jet_fire_api521():
    assert fire.jet_fire_api521(273 + 50) == pytest.approx(83.5e3, abs=500)


def test_jet_fire_scandpower():
    assert fire.jet_fire_scandpower(273 + 20) == pytest.approx(94.5e3, abs=1000)


def test_pool_fire_scandpower():
    assert fire.pool_fire_scandpower(273 + 20) == pytest.approx(88.5e3, abs=500)


def test_validator():
    import os 
    import yaml

    for fname in os.listdir("examples/"):
        with open("examples//" + fname) as infile:
            input = yaml.load(infile, Loader=yaml.FullLoader)
        assert validator.validate_mandatory_ruleset(input) == True


def test_validator2():
    import os 
    import yaml

    for fname in os.listdir("examples/"):
        with open("examples//" + fname) as infile:
            input = yaml.load(infile, Loader=yaml.FullLoader)
        assert validator.heat_transfer_validation(input) == True
        assert validator.valve_validation(input) == True


def test_sim_orifice_full():
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()
    hdown.plot()
    hdown.generate_report()
    assert hdown.report['final_mass'] == pytest.approx(0.14, rel=0.01)
    assert hdown.report['min_wall_temp'] == pytest.approx(284.9, rel=0.01)
    assert hdown.report['min_fluid_temp'] == pytest.approx(193.7, rel=0.01)
    assert hdown.report['time_min_fluid_temp'] == pytest.approx(38.2, rel=0.01)
    assert hdown.report['max_mass_rate'] == pytest.approx(0.869, rel=0.01)


def test_sim_controlvalve():
    from hyddown import HydDown

    input = get_example_input("controlvalve.yml")
    hdown = HydDown(input)
    hdown.run()


def test_sim_psv():
    from hyddown import HydDown

    input = get_example_input("psv.yml")
    hdown = HydDown(input)
    hdown.run()


def test_sim_cv_filling():
    from hyddown import HydDown

    input = get_example_input("cv_filling.yml")
    hdown = HydDown(input)
    hdown.run()

def test_sim_stefan_boltzmann():
    from hyddown import HydDown

    input = get_example_input("psv_sb.yml")
    hdown = HydDown(input)
    hdown.run()

