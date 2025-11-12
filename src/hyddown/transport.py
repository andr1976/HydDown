# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

import math
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from ht import Rohsenow


def Gr(L, Tfluid, Tvessel, P, species):
    """
    Calculation of Grasshof number. See eq. 4.7-4 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    P : float
        Pressure of fluid inventory

    Returns
    ----------
    Gr : float
        Grasshof number
    """
    # Estimating the temperature at the fluid film interface
    T = (Tfluid + Tvessel) / 2
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T|gas", T, "P", P, species)
    nu = PropsSI("V", "T|gas", T, "P", P, "HEOS::" + species.split("::")[1]) / PropsSI(
        "D", "T|gas", T, "P", P, species
    )
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / nu**2
    return Gr


def Pr(T, P, species):
    """
    Calculation of Prandtl number, eq. 4.5-6 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993

    Parameters
    ----------
    T : float
        Temperature of the fluid film interface
    P : float
        Pressure of fluid inventory

    Returns
    ----------
    Pr : float
        Prantdl number
    """
    C = PropsSI("C", "T|gas", T, "P", P, species)
    V = PropsSI("V", "T|gas", T, "P", P, "HEOS::" + species.split("::")[1])
    L = PropsSI("L", "T|gas", T, "P", P, "HEOS::" + species.split("::")[1])
    Pr = C * V / L

    return Pr


def Nu(Ra, Pr):
    """
    Calculation of Nusselt number for natural convection. See eq. 4.7-4  and Table 4.7-1 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993

    Parameters
    ----------
    Ra : float
        Raleigh number
    Pr : float
        Prandtl number

    Returns
    ----------
    Nu : float
        Nusselt numebr
    """
    if Ra >= 1e9:
        NNu = 0.13 * Ra**0.333
    elif Ra < 1e9 and Ra > 1e4:
        NNu = 0.59 * Ra**0.25
    else:
        NNu = 1.36 * Ra**0.20
    return NNu


def h_inside(L, Tvessel, Tfluid, fluid):
    """
    Calculation of internal natural convective heat transfer coefficient from Nusselt number
    and using the coolprop low level interface.

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
        Coolprop fluid object

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient
    """
    cond = fluid.conductivity()
    visc = fluid.viscosity()
    cp = fluid.cpmass()
    Pr = cp * visc / cond
    beta = fluid.isobaric_expansion_coefficient()
    nu = visc / fluid.rhomass()
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / nu**2
    Ra = Pr * Gr
    NNu = Nu(Ra, Pr)
    h_inner = NNu * cond / L
    return h_inner


def h_inner(L, Tfluid, Tvessel, P, species):
    """
    Calculation of internal natural convective heat transfer coefficient from Nusselt number
    and using the coolprop high level interface. Not currently in use.

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    species : str
        Fluid definition string

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient
    """
    NPr = Pr((Tfluid + Tvessel) / 2, P, species)
    NGr = Gr(L, Tfluid, Tvessel, P, species)
    NRa = NPr * NGr
    NNu = Nu(NRa, NPr)
    h_inner = (
        NNu
        * PropsSI(
            "L",
            "T|gas",
            (Tfluid + Tvessel) / 2,
            "P",
            P,
            "HEOS::" + species.split("::")[1],
        )
        / L
    )
    return h_inner


def h_inside_mixed(L, Tvessel, Tfluid, fluid, mdot, D):
    """
    Calculation of internal mixed natural/forced convective heat transfer coefficient from Nusselt number
    and using the coolprop low level interface.

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
        Coolprop fluid object
    mdot : float
        Mass flow
    D : float
        Characteristic diameter for Reynolds number estimation

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient
    """
    cond = fluid.conductivity()
    visc = fluid.viscosity()
    cp = fluid.cpmass()
    Pr = cp * visc / cond

    T = (Tfluid + Tvessel) / 2
    beta = fluid.isobaric_expansion_coefficient()
    nu = visc / fluid.rhomass()
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / nu**2
    Ra = Pr * Gr

    NNu_free = Nu(Ra, Pr)  # 0.13 * NRa**0.333
    Re = 4 * abs(mdot) / (visc * math.pi * D)
    NNu_forced = 0.56 * Re**0.67
    return (NNu_free + NNu_forced) * cond / L


def h_inner_mixed(L, Tfluid, Tvessel, P, species, mdot, D):
    """
    Calculation of internal mixed (nutural/forced convective) heat transfer coefficient from Nusselt number
    and using the coolprop high level interface. Not currently in use.

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    P : float
        Fluid pressure
    species : str
        Fluid definition string
    mdot : float
        Mass flow
    D : float
        Characteristic diameter for Reynolds number estimation


    Returns
    ----------
    h_inner : float
        Heat transfer coefficient
    """
    NPr = Pr((Tfluid + Tvessel) / 2, P, species)
    NGr = Gr(L, Tfluid, Tvessel, P, species)
    NRa = NPr * NGr
    NNu_free = Nu(NRa, NPr)  # 0.13 * NRa**0.333
    Re = (
        4
        * abs(mdot)
        / (
            PropsSI("V", "T|gas", Tfluid, "P", P, "HEOS::" + species.split("::")[1])
            * math.pi
            * D
        )
    )
    NNu_forced = 0.56 * Re**0.67
    return (
        (NNu_free + NNu_forced)
        * PropsSI(
            "L",
            "T|gas",
            (Tfluid + Tvessel) / 2,
            "P",
            P,
            "HEOS::" + species.split("::")[1],
        )
        / L
    )


def h_inside_liquid(L, Tvessel, Tfluid, fluid, master_fluid):
    """
    Calculation of internal natural convective heat transfer coefficient from Nusselt number

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
            Liquid Fluid object equilibrated at film temperature
    master_fluid : obj
            Master Fluid object (used for surface tension etc.)

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient (W/m2 K)
    """
    cond = fluid.conductivity()
    visc = fluid.viscosity()
    cp = fluid.cpmass()
    nu = visc / fluid.rhomass()
    Pr = cp * visc / cond
    beta = fluid.isobaric_expansion_coefficient()

    Pr = cp * visc / cond
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / nu**2
    Ra = Pr * Gr
    NNu = Nu(Ra, Pr)
    h_inner = NNu * cond / L
    return h_inner


def h_inside_wetted(L, Tvessel, Tfluid, fluid, master_fluid):
    """
    Calculation of internal heat transfer coefficient for boiling liquid

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
            Gas object equilibrated at film temperature

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient (W/m2 K)
    """
    kl = fluid.conductivity()
    mul = fluid.viscosity()
    sigma = master_fluid.surface_tension()

    h_boil = Rohsenow(
        rhol=master_fluid.saturated_liquid_keyed_output(CP.iDmass),
        rhog=master_fluid.saturated_vapor_keyed_output(CP.iDmass),
        mul=mul,
        kl=kl,
        Cpl=fluid.cpmass(),
        Hvap=(
            master_fluid.saturated_vapor_keyed_output(CP.iHmass)
            - master_fluid.saturated_liquid_keyed_output(CP.iHmass)
        ),
        sigma=sigma,
        Te=max((Tvessel - Tfluid), 0),
        Csf=0.013,
        # n=1.3,
        # Csf=0.018,
        n=1.7,
    )

    if math.isnan(h_boil):
        h_boil = 0

    h_conv = h_inside_liquid(L, Tvessel, Tfluid, fluid, master_fluid)
    # return 3000
    # return h_conv
    return max(h_boil, h_conv)
    # return min(max(h_boil, h_conv, 1000), 3000)


def gas_release_rate(P1, P2, rho, k, CD, area):
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream pressure
    rho : float
        Fluid density
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    are : float
        Orifice area

    Returns
    ----------
        : float
        Gas release rate / mass flow of discharge
    """
    if P1 > P2:
        if P1 / P2 > ((k + 1) / 2) ** ((k) / (k - 1)):
            flow_coef = 1
        else:
            flow_coef = (
                2
                / (k - 1)
                * (((k + 1) / 2) ** ((k + 1) / (k - 1)))
                * ((P2 / P1) ** (2 / k))
                * (1 - (P2 / P1) ** ((k - 1) / k))
            )

        return (
            math.sqrt(flow_coef)
            * CD
            * area
            * math.sqrt(rho * P1 * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        )
    else:
        return 0


def relief_valve(P1, Pback, Pset, blowdown, k, CD, T1, Z, MW, area):
    """
    Pop action relief valve model including hysteresis.
    The pressure shall rise above P_set to open and
    decrease below P_reseat (P_set*(1-blowdown)) to close

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pback : float
        Downstream / backpressure
    Pset : float
        Set pressure of the PSV / relief valve
    blowdown : float
        The percentage of the set pressure at which the valve reseats
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    T1 : float
        Upstream temperature
    Z : float
        Compressibility
    MW : float
        Molecular weight of the gas relieved
    area : float
        PSV orifice area

    Returns
    ----------
        : float
        Relief rate / mass flow
    """

    global psv_state
    if P1 > Pset:
        eff_area = area
        psv_state = "open"
    elif P1 < Pset * (1 - blowdown):
        eff_area = 0
        psv_state = "closed"
    else:
        if psv_state == "open":
            eff_area = area
        elif psv_state == "closed":
            eff_area = 0
        else:
            raise ValueError("Unknown PSV open/close state.")

    if eff_area > 0:
        return api_psv_release_rate(P1, Pback, k, CD, T1, Z, MW, area)
    else:
        return 0.0


def api_psv_release_rate(P1, Pback, k, CD, T1, Z, MW, area):
    """
    PSV vapour relief rate calculated according to API 520 Part I 2014
    Eq. 5, 9, 15, 18

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pback : float
        Downstream / backpressure
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    T1 : float
        Upstream temperature
    Z : float
        Compressibility
    MW : float
        Molecular weight of the gas relieved
    area : float
        PSV orifice area

    Returns
    ----------
        : float
        Relief rate / mass flow
    """

    P1 = P1 / 1000
    Pback = Pback / 1000
    area = area * 1e6
    MW = MW * 1000
    C = 0.03948 * (k * (2 / (k + 1)) ** ((k + 1) / (k - 1))) ** 0.5
    if P1 / Pback > ((k + 1) / 2) ** ((k) / (k - 1)):
        w = CD * area * C * P1 / math.sqrt(T1 * Z / MW)
    else:
        r = Pback / P1
        f2 = ((k / (k - 1)) * r ** (2 / k) * (1 - r ** ((k - 1) / k)) / (1 - r)) ** 0.5
        print(f2)
        w = CD * area * f2 / (T1 * Z / (MW * P1 * (P1 - Pback))) ** 0.5 / 17.9
    return w / 3600


def cv_vs_time(Cv_max, t, time_constant=0, characteristic="linear"):
    """
    Control valve flow coefficient vs time / actuator postion
    assuming a linear rate of actuator for the three archetypes of
    characteristics: linear, equal percentage and fast/quick opening.

    Parameters
    ----------
    Cv_max : float
        Valve flow coefficient at full open position
    t : float
        Time
    time_constant : float (optional)
        The time required for the actuator to fully open.
        Default to instant open
    characteristic : string (optional)
        Valve characteristic
        Default to linear.
    """

    if time_constant == 0:
        return Cv_max
    else:
        if characteristic == "linear":
            return Cv_max * min(t / time_constant, 1)
        elif characteristic == "eq":
            # https://www.spiraxsarco.com/learn-about-steam/control-hardware-electric-pneumatic-actuation/control-valve-characteristics
            tau = 50
            travel = min(t / time_constant, 1)
            return Cv_max * math.exp(travel * math.log(tau)) / tau
        elif characteristic == "fast":
            # square root function used
            return Cv_max * min(t / time_constant, 1) ** (0.5)
        else:
            return Cv_max


def control_valve(P1, P2, T, Z, MW, gamma, Cv, xT=0.75, FP=1):
    """
    Flow calculated from ANSI/ISA control valve equations for single phase gas flow.
    Equation 19 pp. 132 in
    Control Valves / Guy Borden, editor; Paul Friedmann, style editor

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream / backpressure
    T : float
        Upstream temperature
    Z : float
        Upstream compressibility
    MW : float
        Molecular weight of the gas relieved
    gamma : float
        Upstream Ideal gas k (Cp/Cv)
    Cv : float
        Valve coefficient
    xT : float
        Value of xT for valve fitting assembly, default value
    FP : float
        Piping geometry factor

    Returns
    ----------
        : float
        Mass flow
    """

    P1 = P1 / 1e5
    P2 = P2 / 1e5
    MW = MW * 1000
    N8 = 94.8
    Fk = gamma / 1.4
    x = (P1 - P2) / P1
    if x < 0:
        x = 0
    Y = 1.0 - min(x, Fk * xT) / (3.0 * Fk * xT)
    mass_flow = N8 * FP * Cv * P1 * Y * (MW * min(x, xT * Fk) / T / Z) ** 0.5
    return mass_flow / 3600  # kg/s
