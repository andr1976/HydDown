# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Heat and mass transfer correlations for vessel depressurization/pressurization.

This module provides functions for calculating:
- Dimensionless numbers (Grashof, Prandtl, Nusselt, Rayleigh)
- Heat transfer coefficients for natural and forced convection
- Mass flow rates through orifices, control valves, and relief valves
- Boiling heat transfer (pool boiling, film boiling)

The correlations are based on established literature including:
- Geankoplis, Transport Processes and Unit Operations
- API Standard 520/521 for relief valve sizing
- Rohsenow pool boiling correlation

All functions use SI units unless otherwise specified.
CoolProp is used as the thermodynamic backend for fluid property calculations.
"""

import math
from scipy import optimize
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


def h_inside_liquid(L, Tvessel, Tfluid, fluid):
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

    # Check if sigma is very small (near critical point) to avoid division by zero
    if sigma < 1e-6:
        h_boil = 0
    else:
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

    h_conv = h_inside_liquid(L, Tvessel, Tfluid, fluid)
    # return min(h_boil, h_conv)
    # return h_conv
    return min(max(h_boil, h_conv), 3000)


def h_gas_liquid_interface(T_gas, T_liquid, P, L_char, fluid_gas, fluid_liquid):
    """
    Calculate gas-liquid interfacial heat transfer coefficient for stratified two-phase flow.

    Uses natural convection correlations for horizontal interface based on temperature
    configuration:
    - If T_gas > T_liquid (typical): Hot plate lower surface, liquid properties
    - If T_liquid > T_gas (rare): Hot plate upper surface, gas properties

    Physical situation:
    - Horizontal interface area (vessel cross-section at liquid level)
    - Natural convection driven by buoyancy from temperature difference
    - Characteristic length typically vessel diameter

    Parameters
    ----------
    T_gas : float
        Gas phase bulk temperature [K]
    T_liquid : float
        Liquid phase bulk temperature [K]
    P : float
        System pressure [Pa]
    L_char : float
        Characteristic length for interface [m]
        Typically vessel diameter for horizontal or vertical vessels
    fluid_gas : obj
        CoolProp fluid object for gas phase (already updated at current state)
    fluid_liquid : obj
        CoolProp fluid object for liquid phase (already updated at current state)

    Returns
    -------
    h_gl : float
        Gas-liquid interfacial heat transfer coefficient [W/m²K]

    Notes
    -----
    - Correlation automatically handles both stable and unstable stratification
    - Stable (T_gas > T_liquid): Modest h_gl, heat flows downward
    - Unstable (T_liquid > T_gas): Higher h_gl, heat flows upward
    - Can be extended to include forced/mixed convection effects
    - Typical h_gl range: 20-500 W/m²K for natural convection

    References
    ----------
    - Churchill and Chu (1975): Natural convection from horizontal surfaces
    - Incropera & DeWitt: Fundamentals of Heat and Mass Transfer
    - McAdams: Heat Transmission
    """
    # Temperature difference (driving force for heat transfer)
    dT = abs(T_gas - T_liquid)

    # Safety check: if temperatures are equal, return minimum value
    if dT < 0.01:
        return 100.0  # Minimum baseline (accounts for residual mixing/turbulence)

    # =========================================================================
    # SELECT CORRELATION BASED ON TEMPERATURE CONFIGURATION
    # =========================================================================

    if T_gas > T_liquid:
        # TYPICAL CASE: Hot gas above cold liquid (stable stratification)
        # Heat flows downward from gas to liquid
        # Use correlation for LOWER surface of hot plate with LIQUID properties
        # This is stable - reduced natural convection

        try:
            # Get liquid properties
            k = fluid_liquid.conductivity()
            mu = fluid_liquid.viscosity()
            cp = fluid_liquid.cpmass()
            rho = fluid_liquid.rhomass()
            beta = fluid_liquid.isobaric_expansion_coefficient()

            # Safety check: beta can be negative near saturation for liquids
            # Use absolute value since we're interested in buoyancy magnitude
            beta = abs(beta)

            # Dimensionless numbers
            Pr = cp * mu / k
            nu = mu / rho
            Gr = 9.81 * beta * dT * L_char**3 / nu**2
            Ra = Gr * Pr

            # Detailed check for invalid values
            if math.isnan(Pr) or math.isnan(nu) or math.isnan(Gr) or math.isnan(Ra):
                # NaN detected - property calculation issue
                h_gl = 100.0
            elif math.isinf(Ra) or Ra > 1e20:
                # Extremely high Ra - cap at maximum
                Nu = 0.15 * (1e15)**0.333  # Use very high but finite Ra
                h_gl = Nu * k / L_char
            elif Ra <= 0:
                # Non-positive Ra (shouldn't happen with abs(beta), but safety check)
                h_gl = 100.0
            else:
                # Nusselt number for LOWER surface of hot plate (stable)
                # From Churchill & Chu / Incropera & DeWitt
                if Ra < 1e7:
                    Nu = 0.27 * Ra**0.25  # Laminar, stable
                else:
                    Nu = 0.15 * Ra**0.333  # Turbulent, stable

                # Heat transfer coefficient
                h_gl = Nu * k / L_char

        except Exception as e:
            # Fallback if calculation fails
            h_gl = 50.0  # Conservative default for stable stratification

    else:
        # RARE CASE: Hot liquid below cold gas (unstable stratification)
        # Heat flows upward from liquid to gas
        # Use correlation for UPPER surface of hot plate with GAS properties
        # This is unstable - enhanced natural convection

        try:
            # Get gas properties
            k = fluid_gas.conductivity()
            mu = fluid_gas.viscosity()
            cp = fluid_gas.cpmass()
            rho = fluid_gas.rhomass()
            beta = fluid_gas.isobaric_expansion_coefficient()

            # Safety check: use absolute value of beta
            beta = abs(beta)

            # Dimensionless numbers
            Pr = cp * mu / k
            nu = mu / rho
            Gr = 9.81 * beta * dT * L_char**3 / nu**2
            Ra = Gr * Pr

            # Detailed check for invalid values
            if math.isnan(Pr) or math.isnan(nu) or math.isnan(Gr) or math.isnan(Ra):
                # NaN detected - property calculation issue
                h_gl = 200.0
            elif math.isinf(Ra) or Ra > 1e20:
                # Extremely high Ra - cap at maximum
                Nu = 0.15 * (1e15)**0.333  # Use very high but finite Ra
                h_gl = Nu * k / L_char
            elif Ra <= 0:
                # Non-positive Ra (shouldn't happen with abs(beta), but safety check)
                h_gl = 200.0
            else:
                # Nusselt number for UPPER surface of hot plate (unstable)
                # From Churchill & Chu / Incropera & DeWitt
                if Ra < 1e7:
                    Nu = 0.54 * Ra**0.25  # Laminar, unstable
                else:
                    Nu = 0.15 * Ra**0.333  # Turbulent, unstable

                # Heat transfer coefficient
                h_gl = Nu * k / L_char

        except Exception as e:
            # Fallback if calculation fails
            h_gl = 200.0  # Higher default for unstable stratification

    # Final safety checks
    if h_gl < 100.0:
        h_gl = 100.0  # Minimum realistic value (accounts for residual mixing, turbulence)
    if h_gl > 2000.0:
        h_gl = 2000.0  # Maximum for natural convection

    return h_gl


def h_gas_liquid_interface_two_sided(T_gas, T_liquid, P, L_char, fluid_gas, fluid_liquid):
    """
    Calculate gas-liquid interfacial heat transfer coefficient using a two-sided
    (series resistance) model.

    Computes natural convection HTCs for both the gas side and liquid side of
    the interface independently, then combines them as resistances in series:

        1/h_overall = 1/h_gas + 1/h_liquid

    This captures the thermal resistance of both boundary layers at the interface.

    Parameters
    ----------
    T_gas : float
        Gas phase bulk temperature [K]
    T_liquid : float
        Liquid phase bulk temperature [K]
    P : float
        System pressure [Pa]
    L_char : float
        Characteristic length for interface [m]
    fluid_gas : obj
        CoolProp fluid object for gas phase
    fluid_liquid : obj
        CoolProp fluid object for liquid phase

    Returns
    -------
    h_gl : float
        Overall gas-liquid interfacial heat transfer coefficient [W/m²K]

    Notes
    -----
    For stable stratification (T_gas > T_liquid, typical):
    - Gas side: hot gas above interface, stable (lower surface of hot plate)
    - Liquid side: cold liquid below interface, stable (lower surface of hot plate)

    For unstable stratification (T_liquid > T_gas, rare):
    - Both sides use upper surface of hot plate correlations (enhanced convection)

    References
    ----------
    - Churchill and Chu (1975): Natural convection from horizontal surfaces
    - Incropera & DeWitt: Fundamentals of Heat and Mass Transfer
    """
    dT = abs(T_gas - T_liquid)

    if dT < 0.01:
        return 100.0

    stable = T_gas > T_liquid

    # --- Gas side HTC ---
    h_gas = _natconv_htc(dT, L_char, fluid_gas, stable)

    # --- Liquid side HTC ---
    h_liquid = _natconv_htc(dT, L_char, fluid_liquid, stable)

    # Series resistance: 1/h_overall = 1/h_gas + 1/h_liquid
    if h_gas > 0 and h_liquid > 0:
        h_gl = 1.0 / (1.0 / h_gas + 1.0 / h_liquid)
    else:
        h_gl = min(h_gas, h_liquid) if max(h_gas, h_liquid) > 0 else 50.0

    # Safety bounds
    if h_gl < 5.0:
        h_gl = 5.0
    if h_gl > 2000.0:
        h_gl = 2000.0

    return h_gl


def _natconv_htc(dT, L_char, fluid, stable):
    """
    Calculate natural convection HTC for one side of a horizontal interface.

    Parameters
    ----------
    dT : float
        Temperature difference across interface [K]
    L_char : float
        Characteristic length [m]
    fluid : obj
        CoolProp fluid object (already updated at current state)
    stable : bool
        True for stable stratification (lower surface of hot plate correlations),
        False for unstable (upper surface of hot plate correlations)

    Returns
    -------
    h : float
        Heat transfer coefficient [W/m²K]
    """
    try:
        k = fluid.conductivity()
        mu = fluid.viscosity()
        cp = fluid.cpmass()
        rho = fluid.rhomass()
        beta = abs(fluid.isobaric_expansion_coefficient())

        Pr = cp * mu / k
        nu = mu / rho
        Gr = 9.81 * beta * dT * L_char**3 / nu**2
        Ra = Gr * Pr

        if math.isnan(Ra) or math.isinf(Ra) or Ra <= 0:
            return 50.0

        if stable:
            # Lower surface of hot plate (stable, reduced convection)
            if Ra < 1e7:
                Nu = 0.27 * Ra**0.25
            else:
                Nu = 0.15 * Ra**0.333
        else:
            # Upper surface of hot plate (unstable, enhanced convection)
            if Ra < 1e7:
                Nu = 0.54 * Ra**0.25
            else:
                Nu = 0.15 * Ra**0.333

        h = Nu * k / L_char
        return h

    except Exception:
        return 50.0


def hem_release_rate(P1, Pback, Cd, area, fluid):
    """
    Fluid mass flow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions calculated applying the HEM (Homogenous Equilibrium Model)
    assumption.

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pback : float
        Back/downstream pressure
    Cd : float
        Coefficient of discharge
    area : float
        Orifice area
    fluid : obj
        Fluid object

    Returns
    ----------
        : float
        Gas release rate / mass flow of discharge
    """
    P0 = fluid.p()
    h0 = fluid.hmass()
    s0 = fluid.smass()

    def negflux(P):
        fluid.update(CP.PSmass_INPUTS, P, s0)
        return -fluid.rhomass() * math.sqrt(2 * (h0 - fluid.hmass()))

    P = optimize.minimize_scalar(negflux, bounds=(Pback, P1), method="bounded")["x"]

    G = -negflux(P)

    mass_flow = Cd * area * G
    fluid.update(CP.PSmass_INPUTS, P0, s0)
    return mass_flow


def gas_release_rate(P1, P2, rho, k, CD, area):
    """
    Gas mass flow (kg/s) trough a hole at critical (sonic) or subcritical
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
    area : float
        Orifice area

    Returns
    ----------
        : float
        Gas release rate / mass flow of discharge
    """
    # Only calculate if there is positive pressure difference (P1 > P2)
    if P1 > P2:
        # Determine if flow is critical (sonic/choked) or subcritical
        # Critical pressure ratio from isentropic flow theory
        critical_pressure_ratio = ((k + 1) / 2) ** ((k) / (k - 1))

        if P1 / P2 > critical_pressure_ratio:
            # Critical (choked) flow - velocity at throat equals speed of sound
            # Flow coefficient = 1 for sonic conditions
            flow_coef = 1
        else:
            # Subcritical flow - expansion is not complete to sonic velocity
            # Apply Yellow Book correction factor for subcritical flow
            flow_coef = (
                2
                / (k - 1)
                * (((k + 1) / 2) ** ((k + 1) / (k - 1)))
                * ((P2 / P1) ** (2 / k))
                * (1 - (P2 / P1) ** ((k - 1) / k))
            )

        # Calculate mass flow rate using Yellow Book equation 2.22
        return (
            math.sqrt(flow_coef)
            * CD
            * area
            * math.sqrt(rho * P1 * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        )
    else:
        # No flow if downstream pressure exceeds upstream pressure
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

    # Convert units for API 520 equations
    P1 = P1 / 1000  # Pa to kPa
    Pback = Pback / 1000  # Pa to kPa
    area = area * 1e6  # m² to mm²
    MW = MW * 1000  # kg/mol to g/mol

    # API 520 critical flow coefficient
    C = 0.03948 * (k * (2 / (k + 1)) ** ((k + 1) / (k - 1))) ** 0.5

    # Check if flow is critical (choked) or subcritical
    if P1 / Pback > ((k + 1) / 2) ** ((k) / (k - 1)):
        # Critical flow (choked at throat)
        w = CD * area * C * P1 / math.sqrt(T1 * Z / MW)
    else:
        # Subcritical flow (not choked)
        r = Pback / P1
        # Subcritical flow correction factor f2
        f2 = ((k / (k - 1)) * r ** (2 / k) * (1 - r ** ((k - 1) / k)) / (1 - r)) ** 0.5
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
