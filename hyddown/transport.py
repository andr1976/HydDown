# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import math
from CoolProp.CoolProp import PropsSI


def Gr(L, Tfluid, Tvessel, P, species):
    """
    Calculation of Grasshof number. See eq. 4.7-4 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993
    """
    T = (Tfluid + Tvessel) / 2
    beta = PropsSI("ISOBARIC_EXPANSION_COEFFICIENT", "T|gas", T, "P", P, species)
    nu = PropsSI("V", "T|gas", T, "P", P, 'HEOS::'+species.split('::')[1]) / PropsSI("D", "T|gas", T, "P", P, species)
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L ** 3 / nu ** 2
    return Gr


def Pr(T, P, species):
    """
    Calculation of Prandtl number, eq. 4.5-6 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993
    """
    C = PropsSI("C", "T|gas", T, "P", P, species)
    V = PropsSI("V", "T|gas", T, "P", P, 'HEOS::'+species.split('::')[1])
    L = PropsSI("L", "T|gas", T, "P", P, 'HEOS::'+species.split('::')[1])
    Pr = C * V / L

    return Pr


def Nu(Ra, Pr):
    """
    Calculation of Nusselt number for natural convection. See eq. 4.7-4  and Table 4.7-1 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993
    """
    if Ra >= 1e9:
        NNu = 0.13 * Ra ** 0.333
    elif Ra < 1e9 and Ra > 1e4:
        NNu = 0.59 * Ra ** 0.25
    else:
        NNu = 1.36 * Ra ** 0.20
    return NNu

def h_inside(L, Tvessel, Tfluid, fluid):
    cond = fluid.conductivity()
    visc = fluid.viscosity()
    cp = fluid.cpmass()
    Pr = cp * visc / cond

    T = (Tfluid + Tvessel) / 2
    beta = fluid.isobaric_expansion_coefficient()
    nu = visc / fluid.rhomass()
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L ** 3 / nu ** 2
    Ra = Pr * Gr
    NNu = Nu(Ra, Pr)
    h_inner = NNu * cond / L
    return h_inner


def h_inner(L, Tfluid, Tvessel, P, species):
    """
    Calculation of heat transfer coefficient from Nusselt number
    """
    NPr = Pr((Tfluid + Tvessel) / 2, P, species)
    NGr = Gr(L, Tfluid, Tvessel, P, species)
    NRa = NPr * NGr
    NNu = Nu(NRa, NPr)
    h_inner = NNu * PropsSI("L", "T|gas", (Tfluid + Tvessel) / 2, "P", P, 'HEOS::'+species.split('::')[1]) / L
    return h_inner


def h_inner_mixed(L, Tfluid, Tvessel, P, species, mdot, D):
    NPr = Pr((Tfluid + Tvessel) / 2, P, species)
    NGr = Gr(L, Tfluid, Tvessel, P, species)
    NRa = NPr * NGr
    NNu_free = Nu(NRa,NPr)  # 0.13 * NRa**0.333
    Re = 4 * abs(mdot) / (PropsSI("V", "T|gas", Tfluid, "P", P, 'HEOS::'+species.split('::')[1]) * math.pi * D)
    NNu_forced = 0.56 * Re ** 0.67
    return (
        (NNu_free + NNu_forced)
        * PropsSI("L", "T|gas", (Tfluid + Tvessel) / 2, "P", P, 'HEOS::'+species.split('::')[1])
        / L
    )


def gas_release_rate(P1, P2, rho, k, CD, area):
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996
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
        f2 = ((k / (k - 1)) * r ** (2 / k) * (1 - r**((k - 1) / k)) / (1 - r))**0.5
        print(f2)
        w = CD * area * f2 / (T1 * Z / (MW * P1 * (P1 - Pback)))**0.5 / 17.9
    return w/3600


def control_valve(P1, P2, T, Z, MW, gamma, Cv, xT=0.75, FP=1):
    """
    Flow calculated from ANSI/ISA control valve equations for single phase gas flow.
    Equation 19 pp. 132 in
    Control Valves / Guy Borden, editor; Paul Friedmann, style editor
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