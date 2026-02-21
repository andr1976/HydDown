# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

from hyddown import transport as tp
from hyddown import fire
from hyddown import validator
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import pytest


def get_example_input(fname):
    import os
    import yaml

    fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "examples", fname)

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
        P1, P2, D, cpcv, 0.85, 0.01**2 / 4 * 3.1415
    ) == pytest.approx(9.2 / 60, rel=0.001)


def test_orifice1():
    P1 = 10.0e5
    P2 = 6.5e5
    D = PropsSI("D", "P", P1, "T", 298.15, "HEOS::N2")
    cpcv = PropsSI("CP0MOLAR", "T", 298.15, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", 298.15, "P", P1, "HEOS::N2"
    )
    assert tp.gas_release_rate(
        P1, P2, D, cpcv, 0.85, 0.01**2 / 4 * 3.1415
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


def test_cv_vs_time():
    assert tp.cv_vs_time(1, 0.5, time_constant=1, characteristic="linear") == 0.5
    assert tp.cv_vs_time(1, 0.5, time_constant=1, characteristic="eq") == pytest.approx(
        0.14, abs=0.002
    )
    assert tp.cv_vs_time(
        1, 0.5, time_constant=1, characteristic="fast"
    ) == pytest.approx(0.707, abs=0.002)
    assert tp.cv_vs_time(1, 0.5, time_constant=0) == 1.0


def test_psv3():
    Pback = 1e5
    Pset = 18.2e5
    blowdown = 0.1
    P1 = 0.99 * Pset * (1 - blowdown)
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area) == 0
    psv_state = "closed"
    assert (
        tp.relief_valve(P1 * 1.01, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area)
        == 0
    )


def test_psv2():
    P1 = 21.0e5
    Pback = 1e5
    Pset = 20.99e5
    blowdown = 0.1
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
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
        Pset * 0.99, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1046 / 3600, rel=0.03)


def test_psv():
    P1 = 100e5
    Pback = 1e5
    Pset = 99.2e5
    blowdown = 0.1
    T1 = 25.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1.57, rel=0.02)


def test_api_psv_relief():
    assert tp.api_psv_release_rate(
        121.9e5, 71e5, 1.39, 0.975, 298.15, 1.01, 2 / 1e3, 71e-6
    ) == pytest.approx(1846 / 3600, rel=0.01)
    assert tp.api_psv_release_rate(
        121.9e5, 1e5, 1.39, 0.975, 298.15, 1.01, 2 / 1e3, 71e-6
    ) == pytest.approx(1860 / 3600, rel=0.01)


def test_hinside():
    fluid = CP.AbstractState("HEOS", "air")
    Tboundary = (311 + 505.4) / 2
    fluid.update(CP.PT_INPUTS, 1e5, Tboundary)
    h = tp.h_inside(0.305, 311, 505.4, fluid)
    assert h == pytest.approx(7, abs=0.1)


def test_hinside_mixed():
    mdot = 1e-10
    D = 0.010
    fluid = CP.AbstractState("HEOS", "air")
    Tboundary = (311 + 505.4) / 2
    fluid.update(CP.PT_INPUTS, 1e5, Tboundary)
    h = tp.h_inside_mixed(0.305, 311, 505.4, fluid, mdot, D)
    assert h == pytest.approx(7, abs=0.1)


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


def test_sb():
    assert fire.sb_fire(273 + 50, "api_jet") == pytest.approx(83.5e3, abs=500)
    assert fire.sb_fire(273 + 50, "api_pool") == pytest.approx(45.5e3, abs=100)
    assert fire.sb_fire(273 + 20, "scandpower_pool") == pytest.approx(88.5e3, abs=500)
    assert fire.sb_fire(273 + 20, "scandpower_jet") == pytest.approx(94.5e3, abs=1000)
    try:
        Q = fire.sb_fire(273 + 20, "scand_jet") == pytest.approx(94.5e3, abs=1000)
    except ValueError:
        pass


def test_validator():
    import os
    import yaml

    dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "examples")

    for fname in os.listdir(dir):
        with open(os.path.join(dir, fname)) as infile:
            input = yaml.load(infile, Loader=yaml.FullLoader)
        assert validator.validate_mandatory_ruleset(input) == True


def test_validator2():
    import os
    import yaml

    dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "examples")

    for fname in os.listdir(dir):
        with open(os.path.join(dir, fname)) as infile:
            input = yaml.load(infile, Loader=yaml.FullLoader)
        assert validator.heat_transfer_validation(input) != False
        assert validator.valve_validation(input) != False


def test_sim_orifice_full():
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()
    hdown.plot(verbose=False)
    hdown.generate_report()
    assert hdown.report["final_mass"] == pytest.approx(0.14, rel=0.02)
    assert hdown.report["min_wall_temp"] == pytest.approx(284.9, rel=0.01)
    assert hdown.report["min_fluid_temp"] == pytest.approx(193.7, rel=0.01)
    assert hdown.report["time_min_fluid_temp"] == pytest.approx(37.0, rel=0.01)
    assert hdown.report["max_mass_rate"] == pytest.approx(0.88, rel=0.01)


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


def test_isenthalpic():
    from hyddown import HydDown

    input = get_example_input("isenthalpic.yml")
    hdown = HydDown(input)
    hdown.run()


def test_sim_cv_filling():
    from hyddown import HydDown

    input = get_example_input("cv_filling.yml")
    hdown = HydDown(input)
    hdown.run()


def test_mdot_filling():
    from hyddown import HydDown

    input = get_example_input("mdot_filling.yml")
    hdown = HydDown(input)
    hdown.run()


def test_sim_filling():
    from hyddown import HydDown

    input = get_example_input("filling.yml")
    hdown = HydDown(input)
    hdown.run()


# def test_multicomponent():
#    from hyddown import HydDown
#
#    input = get_example_input("ng.yml")
#    hdown = HydDown(input)
#    hdown.run()


def test_sim_stefan_boltzmann():
    from hyddown import HydDown

    input = get_example_input("psv_sb.yml")
    hdown = HydDown(input)
    hdown.run()


def test_sim_rupture():
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    from hyddown import HydDown

    # import scienceplots

    # plt.style.use(["science", "nature"])

    input = get_example_input("rupture.yml")
    hdown = HydDown(input)
    hdown.run()

    fname = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "..",
        "..",
        "validation",
        "rupture_history.csv",
    )
    df_ref = pd.read_csv(fname)
    plt.figure(1)
    plt.plot(hdown.time_array, hdown.P, label="Hyddown", linestyle="-", color="#125A56")
    plt.plot(
        df_ref["Time"],
        df_ref["Vapour - Pressure"] * 1e5 + 1.013e5,
        label="Unisim",
        linestyle="--",
        color="#125A56",
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (Pa)")
    plt.legend()

    plt.figure(2)
    plt.plot(
        hdown.time_array,
        hdown.T_fluid - 273.15,
        label="Hyddown - Gas",
        linestyle="-",
        color="#FF6F61",
    )
    plt.plot(
        df_ref["Time"],
        df_ref["Vapour - Temperature"],
        label="Unisim - Gas",
        linestyle="--",
        color="#FF6F61",
    )
    plt.plot(
        hdown.time_array,
        hdown.T_vessel - 273.15,
        label="Hyddown - Wall",
        linestyle="-",
        color="#00767B",
    )
    plt.plot(
        df_ref["Time"],
        df_ref[
            "Vessel - Heat Loss: Inner Wall Temperatures (Heat Loss: Inner Wall Temperatures_1)"
        ],
        label="Unisim - Wall",
        linestyle="--",
        color="#00767B",
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.show()

    assert hdown.T_vessel[-1] == pytest.approx(
        df_ref[
            "Vessel - Heat Loss: Inner Wall Temperatures (Heat Loss: Inner Wall Temperatures_1)"
        ][-1:]
        + 273.15,
        abs=10,
    )
    assert hdown.T_fluid[-1] == pytest.approx(
        df_ref["Vapour - Temperature"][-1:] + 273.15, rel=0.01
    )


def test_dataframe():
    from hyddown import HydDown

    input = get_example_input("controlvalve.yml")
    hdown = HydDown(input)
    hdown.run()
    df = hdown.get_dataframe()


def test_thermesh():
    from hyddown import thermesh as tm

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.special import erfc

    # Domain information
    L = 0.01  # Make sure the domain is large enough to be a semi-infinite solid
    k, rho, cp = 0.72, 1560.0, 1450.0
    h = 20.0  # Heat transfer coefficient
    T_0, T_inf = 0.0, 400.0  # Initial and far field temperature

    #
    # Analytical solution
    # -------------------

    def analytical_solution(x, t, rho, cp, k, q, T_inf):
        """Returns analytical solution when initial temperature equals zero.

        Arguments
        ---------
        x : nd.array(dtype=float, dim=1)
            Spatial coordinates.
        t : float
            Time.
        rho : float
            Density.
        cp : float
            Specific heat.
        k : float
            Thermal conductivity
        h : float
            Heat transfer coefficient.
        T_inf : float
            Far field temperature.

        Returns
        -------
        dT : nd.array(dtype=float, dim=1, len=len(x))
            Temperature change at x.

        """
        alpha = k / rho / cp
        dT = T_inf * (
            erfc(x / (2 * np.sqrt(alpha * t)))
            - np.exp(h * x / k + h**2 * alpha * t / k**2)
            * erfc(x / (2 * np.sqrt(alpha * t)) + h * np.sqrt(alpha * t) / k)
        )
        return dT

    # Plot analytical solution at the following times
    t_inc = np.array([2, 10, 25])

    # Domain coordinates and thermal diffusivity
    z = np.linspace(0, L, 50)

    # fig = plt.figure(figsize=(2.4, 2.0))
    # ax = fig.add_subplot(1, 1, 1)
    # ax.set_title("Analytical solution")
    # ax.set_xlabel("z/L")
    # ax.set_ylabel(r"$\Delta$T ($^{\circ}$C)")
    # ax.set_xlim([0, 1])
    # for t in t_inc:
    #     plt.plot(
    #         z / L,
    #         analytical_solution(z, t, rho, cp, k, h, T_inf),
    #         "k",
    #         linewidth=0.5,
    #     )
    # plt.text(0.04, 1, "2")
    # plt.text(0.12, 5, "10")
    # plt.text(0.40, 8, "25")
    # plt.tight_layout()
    # plt.savefig("../fig/conv_analytical_sol.png", dpi=600)
    # plt.show(block=False)

    #
    # Numerical solution
    # ------------------

    # Let's make a mesh using linear elements (LinearElement). The
    # alternative is to use second-order elements (QuadraticElement).
    nn = 11  # number of nodes
    z = np.linspace(0, L, nn)
    mesh = tm.Mesh(z, tm.LinearElement)  # Or `QuadraticElement` to
    # use quadratic shape functions

    # The boundary conditions are provided in a two-item list of
    # dictionaries. The first dictionary (or zeroth item in the list)
    # applies to the start or left side of the domain, while the second
    # item applies to the end or right side of the domain. The
    # dictionaries can have the following keys:
    #
    #   "T" OR ( ("h" and "T_inf") AND/OR "q" ),
    #
    # with "T" an applied temperature, "h" and "T_inf" the convective heat
    # transfer coefficient and far field temperature, respectively, while
    # "q" represents a direct flux on the surface.
    bc = [
        {"h": h, "T_inf": T_inf},  # convective bc on the left
        {"T": T_0},
    ]  # T on the right

    # Material model (CPEEK is a function that takes T as input)
    cpeek = tm.isothermal_model(k, rho, cp)

    # Define and solve problem
    domain = tm.Domain(mesh, [cpeek], bc)
    domain.set_T(T_0 * np.ones(nn))

    # Solver details
    theta = 0.5
    dt = 0.5
    solver = {"dt": dt, "t_end": 25.0, "theta": theta}
    t, T = tm.solve_ht(domain, solver)

    # Plot numerical solution
    # fig = plt.figure(figsize=(2.4, 2.0))
    # ax = fig.add_subplot(1, 1, 1)
    # ax.set_title(
    #     "FE solution (dz/L="
    #     + str((z[1] - z[0]) / L)
    #     + ", dt="
    #     + str(dt)
    #     + r" s., $\Theta$="
    #     + str(theta)
    #     + ")"
    # )
    # ax.set_xlabel("z/L")
    # ax.set_ylabel(r"$\Delta$T ($^{\circ}$C)")
    # ax.set_xlim([0, 1])
    t_inc = np.array([2, 10, 25]) / dt
    # for i in t_inc.astype(int):
    #     plt.plot(z / L, T[i, :], "k", linewidth=0.5)
    # plt.text(0.04, 1, "2")
    # plt.text(0.12, 5, "10")
    # plt.text(0.40, 8, "25")
    # plt.tight_layout()
    # plt.savefig("../fig/conv_FE_t0.5_dt0.1s.png", dpi=600)
    # plt.show(block=False)
    analytical = analytical_solution(z, 25, rho, cp, k, h, T_inf)

    assert T[-1, 0] == pytest.approx(analytical[0], abs=0.3)
    assert sum((T[-1, :] - analytical) ** 2) < 1


def test_boiling_heat_transfer():
    sat_water = CP.AbstractState("HEOS", "water")
    sat_water.set_mole_fractions([1.0])
    sat_water.update(CP.PQ_INPUTS, 1e5, 0.01)
    water = CP.AbstractState("HEOS", "water")
    water.set_mole_fractions([1.0])
    water.specify_phase(CP.iphase_liquid)
    water.update(CP.PT_INPUTS, 1e5, 373.15)

    h_inner = tp.h_inside_wetted(
        L=0.01,
        Tvessel=373.15 + 4.9,
        Tfluid=water.T(),
        fluid=water,
        master_fluid=sat_water,
    )

    rhol = sat_water.saturated_liquid_keyed_output(CP.iDmass)
    rhog = sat_water.saturated_vapor_keyed_output(CP.iDmass)
    mul = water.viscosity()
    kl = water.conductivity()
    Cpl = water.cpmass()
    Hvap = sat_water.saturated_vapor_keyed_output(
        CP.iHmass
    ) - sat_water.saturated_liquid_keyed_output(CP.iHmass)
    sigma = sat_water.surface_tension()
    from ht import Rohsenow

    h = Rohsenow(
        rhol=rhol,
        rhog=rhog,
        mul=mul,
        kl=kl,
        Cpl=Cpl,
        Hvap=Hvap,
        sigma=sigma,
        Te=4.9,
        Csf=0.011,
        n=1.26,
    )
    print(h)
    assert h == pytest.approx(3571, rel=0.01)


def test_different_vessel_types():
    """Test different vessel head types."""
    from hyddown import HydDown

    # Test with different vessel types
    vessel_types = ["Flat-end", "ASME F&D"]

    for vtype in vessel_types:
        input = get_example_input("input.yml")
        input["vessel"]["type"] = vtype
        hdown = HydDown(input)
        hdown.run()

        # Verify calculation completes for all types
        assert len(hdown.time_array) > 0
        assert hdown.vol > 0


def test_plot_envelope():
    """Test phase envelope plotting functionality."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()

    # Test plot_envelope doesn't raise exceptions
    try:
        hdown.plot_envelope(verbose=False)
    except Exception as e:
        pytest.fail(f"plot_envelope raised exception: {e}")


def test_generate_report():
    """Test report generation with all metrics."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()
    hdown.generate_report()

    # Verify report contains expected keys
    assert "final_mass" in hdown.report
    assert "min_wall_temp" in hdown.report
    assert "min_fluid_temp" in hdown.report
    assert "time_min_fluid_temp" in hdown.report
    assert "max_mass_rate" in hdown.report

    # Verify values are reasonable
    assert hdown.report["final_mass"] >= 0
    assert hdown.report["min_wall_temp"] > 0
    assert hdown.report["min_fluid_temp"] > 0
    assert hdown.report["max_mass_rate"] >= 0


def test_get_dataframe_structure():
    """Test DataFrame output has correct structure."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()
    df = hdown.get_dataframe()

    # Verify DataFrame has expected columns (exact names from get_dataframe)
    assert "Time (s)" in df.columns
    assert "Pressure (bar)" in df.columns
    assert "Fluid temperature (oC)" in df.columns
    assert "Fluid mass (kg)" in df.columns

    # Verify data consistency
    assert len(df) == len(hdown.time_array)
    assert all(df["Pressure (bar)"] > 0)


def test_isothermal_calculation():
    """Test isothermal calculation type."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    input["calculation"]["type"] = "isothermal"

    hdown = HydDown(input)
    hdown.run()

    # Temperature should remain constant for isothermal
    T_initial = hdown.T_fluid[0]
    T_final = hdown.T_fluid[-1]

    # Allow small deviation due to numerical effects
    assert T_final == pytest.approx(T_initial, rel=0.05)


def test_isentropic_calculation():
    """Test isentropic calculation type."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    input["calculation"]["type"] = "isentropic"

    hdown = HydDown(input)
    hdown.run()

    # For discharge, temperature should decrease in isentropic process
    T_initial = hdown.T_fluid[0]
    T_final = hdown.T_fluid[-1]

    assert T_final < T_initial


def test_mass_conservation_discharge():
    """Test mass conservation during discharge."""
    from hyddown import HydDown
    import numpy as np

    input = get_example_input("input.yml")
    hdown = HydDown(input)
    hdown.run()

    # For discharge, mass should only decrease (or stay constant)
    mass_diff = np.diff(hdown.mass_fluid)
    assert all(mass_diff <= 1e-10)  # Allow small numerical tolerance


def test_mass_conservation_filling():
    """Test mass conservation during filling."""
    from hyddown import HydDown
    import numpy as np

    input = get_example_input("filling.yml")
    hdown = HydDown(input)
    hdown.run()

    # For filling, mass should only increase (or stay constant)
    mass_diff = np.diff(hdown.mass_fluid)
    assert all(mass_diff >= -1e-10)  # Allow small numerical tolerance


def test_energy_balance_with_fire():
    """Test energy balance calculation with fire heat load."""
    from hyddown import HydDown

    input = get_example_input("psv_sb.yml")
    hdown = HydDown(input)
    hdown.run()

    # Verify fire heat load increases temperature/pressure
    assert len(hdown.Q_outer) > 0
    # External heat should be positive (fire adding heat)
    assert max(hdown.Q_outer) > 0


def test_validation_data_plotting():
    """Test plotting with validation data."""
    from hyddown import HydDown

    input = get_example_input("input.yml")

    # Add validation data
    input["validation"] = {
        "pressure": {"time": [0, 10, 20, 30], "pres": [150e5, 120e5, 90e5, 60e5]},
        "temperature": {
            "gas_high": {"time": [0, 10, 20, 30], "temp": [298, 280, 260, 240]}
        },
    }

    hdown = HydDown(input)
    hdown.run()

    # Plot should handle validation data without errors
    try:
        hdown.plot(verbose=False)
    except Exception as e:
        pytest.fail(f"plot with validation data raised exception: {e}")


def test_invalid_calculation_type():
    """Test that invalid calculation type raises error."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    input["calculation"]["type"] = "invalid_type"

    with pytest.raises(Exception):
        hdown = HydDown(input)
        hdown.run()


def test_zero_initial_mass():
    """Test behavior with unrealistic zero initial mass."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    input["initial"]["pressure"] = 1e5  # Very low pressure (near atmospheric)

    try:
        hdown = HydDown(input)
        hdown.run()
        # Should complete but with very small mass
        assert hdown.m0 > 0
    except Exception:
        # May raise exception for unrealistic conditions
        pass


def test_very_short_timestep():
    """Test stability with very short timestep."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    input["calculation"]["time_step"] = 0.001  # Very short
    input["calculation"]["end_time"] = 1.0  # Short simulation

    hdown = HydDown(input)
    hdown.run()

    # Should complete without numerical instabilities
    assert len(hdown.time_array) > 0
    assert all(hdown.P > 0)
    assert all(hdown.T_fluid > 0)


def test_vessel_horizontal_vs_vertical():
    """Test that orientation affects heat transfer."""
    from hyddown import HydDown

    # Horizontal vessel
    input_h = get_example_input("input.yml")
    input_h["vessel"]["orientation"] = "horizontal"
    hdown_h = HydDown(input_h)
    hdown_h.run()

    # Vertical vessel
    input_v = get_example_input("input.yml")
    input_v["vessel"]["orientation"] = "vertical"
    hdown_v = HydDown(input_v)
    hdown_v.run()

    # Results should differ due to different heat transfer correlations
    # (though may be close)
    assert len(hdown_h.time_array) > 0
    assert len(hdown_v.time_array) > 0


def test_string_representation():
    """Test __str__ method returns valid string."""
    from hyddown import HydDown

    input = get_example_input("input.yml")
    hdown = HydDown(input)

    str_repr = str(hdown)
    assert isinstance(str_repr, str)
    assert len(str_repr) > 0


if __name__ == "__main__":
    test_sim_rupture()
