# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Mass flow calculator for valve discharge and filling operations.

This module provides the MassFlowCalculator class which encapsulates all
mass flow rate calculations for different valve types used in HydDown
simulations. It handles compressible flow through orifices, relief valves,
control valves, and constant mass flow devices.

The MassFlowCalculator handles:
- Orifice flow (Yellow Book equation 2.22)
- Relief valve flow with hysteresis (pop-action)
- PSV flow (API 520/521)
- Control valve flow with Cv characteristics
- Constant mass flow rate (mdot)
- Time-varying mass flow profiles

Key Features:
- Manages PSV state (open/closed) with hysteresis
- Handles both critical (choked) and subcritical flow regimes
- Supports filling and discharge modes
- Provides control valve characteristics (linear, equal percentage, fast opening)

Typical usage:
    calc = MassFlowCalculator(
        valve_type="orifice",
        flow_direction="discharge",
        diameter=0.01,
        discharge_coef=0.8,
        back_pressure=101325
    )
    mdot = calc.calculate(P=1e6, T=300, rho=10, k=1.4)
"""

import math
import numpy as np


class MassFlowCalculator:
    """
    Mass flow rate calculator for various valve types.

    This class encapsulates all mass flow calculations for discharge and filling
    operations through different valve types. It manages valve-specific state
    (e.g., PSV open/closed status) and provides a unified interface for calculating
    mass flow rates.

    Parameters
    ----------
    valve_type : str
        Type of valve: "orifice", "psv", "relief", "controlvalve", "mdot"
    flow_direction : str
        Flow direction: "discharge" or "filling"
    diameter : float, optional
        Valve/orifice diameter [m]. Required for orifice, psv.
    discharge_coef : float, optional
        Discharge coefficient [-]. Required for orifice, psv. Default 0.62.
    back_pressure : float
        Downstream/back pressure [Pa]
    set_pressure : float, optional
        Set pressure for relief/PSV [Pa]. Required for psv, relief.
    blowdown : float, optional
        Blowdown fraction for PSV [-]. Required for psv. Typically 0.1-0.2.
    Cv : float, optional
        Control valve flow coefficient. Required for controlvalve.
    xT : float, optional
        Pressure drop ratio factor for control valve. Default 0.75.
    FP : float, optional
        Piping geometry factor for control valve. Default 1.0.
    mdot : float or array-like, optional
        Constant mass flow rate [kg/s] or time-varying profile. Required for mdot.
    time : array-like, optional
        Time points for mdot profile [s]. Required if mdot is array.
    time_constant : float, optional
        Control valve actuation time constant [s]. Default 0.
    characteristic : str, optional
        Control valve characteristic: "linear", "eq", "fast". Default "linear".

    Attributes
    ----------
    psv_state : str
        Current PSV state: "open" or "closed". Maintains hysteresis.
    """

    def __init__(
        self,
        valve_type,
        flow_direction,
        diameter=None,
        discharge_coef=0.62,
        back_pressure=101325,
        set_pressure=None,
        blowdown=None,
        Cv=None,
        xT=0.75,
        FP=1.0,
        mdot=None,
        time=None,
        time_constant=0,
        characteristic="linear",
    ):
        """
        Initialize mass flow calculator with valve parameters.

        Raises
        ------
        ValueError
            If required parameters for the valve type are missing
        """
        self.valve_type = valve_type
        self.flow_direction = flow_direction
        self.diameter = diameter
        self.discharge_coef = discharge_coef
        self.back_pressure = back_pressure
        self.set_pressure = set_pressure
        self.blowdown = blowdown
        self.Cv = Cv
        self.xT = xT
        self.FP = FP
        self.mdot_value = mdot
        self.time_profile = time
        self.time_constant = time_constant
        self.characteristic = characteristic

        # PSV state management (for hysteresis)
        self.psv_state = "closed"

        # Calculate orifice area if diameter provided
        self.area = math.pi * diameter**2 / 4 if diameter is not None else None

        # Validate required parameters
        self._validate_parameters()

    def _validate_parameters(self):
        """Validate that required parameters are provided for valve type."""
        if self.valve_type == "orifice":
            if self.diameter is None or self.discharge_coef is None:
                raise ValueError("Orifice valve requires 'diameter' and 'discharge_coef'")

        elif self.valve_type == "psv":
            if any(
                p is None
                for p in [self.diameter, self.discharge_coef, self.set_pressure, self.blowdown]
            ):
                raise ValueError(
                    "PSV requires 'diameter', 'discharge_coef', 'set_pressure', 'blowdown'"
                )

        elif self.valve_type == "relief":
            if self.set_pressure is None:
                raise ValueError("Relief valve requires 'set_pressure'")

        elif self.valve_type == "controlvalve":
            if self.Cv is None:
                raise ValueError("Control valve requires 'Cv'")

        elif self.valve_type == "mdot":
            if self.mdot_value is None:
                raise ValueError("Mdot valve requires 'mdot' parameter")

    def calculate(self, P, T=None, rho=None, k=None, Z=None, MW=None, t=0):
        """
        Calculate mass flow rate for current valve type and conditions.

        Parameters
        ----------
        P : float
            Upstream pressure [Pa]
        T : float, optional
            Upstream temperature [K]. Required for psv, relief, controlvalve.
        rho : float, optional
            Upstream density [kg/m³]. Required for orifice.
        k : float, optional
            Heat capacity ratio Cp/Cv [-]. Required for orifice, psv, relief, controlvalve.
        Z : float, optional
            Compressibility factor [-]. Required for psv, relief, controlvalve.
        MW : float, optional
            Molecular weight [kg/mol]. Required for psv, relief, controlvalve.
        t : float, optional
            Current time [s]. Required for controlvalve, mdot with time profile.

        Returns
        -------
        float
            Mass flow rate [kg/s]. Positive for discharge, negative for filling.
        """
        if self.valve_type == "orifice":
            mdot = self._orifice_flow(P, rho, k)

        elif self.valve_type == "psv":
            mdot = self._psv_flow(P, T, k, Z, MW)

        elif self.valve_type == "relief":
            mdot = self._relief_valve_flow(P, T, k, Z, MW)

        elif self.valve_type == "controlvalve":
            mdot = self._control_valve_flow(P, T, k, Z, MW, t)

        elif self.valve_type == "mdot":
            mdot = self._constant_mdot(t)

        else:
            raise ValueError(f"Unknown valve type: {self.valve_type}")

        # Apply flow direction (negative for filling)
        if self.flow_direction == "filling":
            mdot = -mdot

        return mdot

    def _orifice_flow(self, P1, rho, k):
        """
        Calculate orifice flow using Yellow Book equation 2.22.

        Handles both critical (choked) and subcritical flow regimes.

        Parameters
        ----------
        P1 : float
            Upstream pressure [Pa]
        rho : float
            Upstream density [kg/m³]
        k : float
            Heat capacity ratio Cp/Cv [-]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        P2 = self.back_pressure

        # Only calculate if there is positive pressure difference
        if P1 <= P2:
            return 0.0

        # Critical pressure ratio from isentropic flow theory
        critical_pressure_ratio = ((k + 1) / 2) ** (k / (k - 1))

        if P1 / P2 > critical_pressure_ratio:
            # Critical (choked) flow - velocity at throat equals speed of sound
            flow_coef = 1.0
        else:
            # Subcritical flow - apply Yellow Book correction factor
            flow_coef = (
                2
                / (k - 1)
                * (((k + 1) / 2) ** ((k + 1) / (k - 1)))
                * ((P2 / P1) ** (2 / k))
                * (1 - (P2 / P1) ** ((k - 1) / k))
            )

        # Yellow Book equation 2.22
        mdot = (
            math.sqrt(flow_coef)
            * self.discharge_coef
            * self.area
            * math.sqrt(rho * P1 * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        )

        return mdot

    def _psv_flow(self, P1, T1, k, Z, MW):
        """
        Calculate PSV flow with hysteresis using API 520/521.

        Implements pop-action relief valve model where:
        - Valve opens when P1 > set_pressure
        - Valve closes when P1 < set_pressure * (1 - blowdown)
        - Maintains state during hysteresis region

        Parameters
        ----------
        P1 : float
            Upstream pressure [Pa]
        T1 : float
            Upstream temperature [K]
        k : float
            Heat capacity ratio Cp/Cv [-]
        Z : float
            Compressibility factor [-]
        MW : float
            Molecular weight [kg/mol]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        # Determine effective area based on PSV state and hysteresis
        if P1 > self.set_pressure:
            # Pressure above set point - open valve
            eff_area = self.area
            self.psv_state = "open"
        elif P1 < self.set_pressure * (1 - self.blowdown):
            # Pressure below reseat point - close valve
            eff_area = 0.0
            self.psv_state = "closed"
        else:
            # In hysteresis region - maintain current state
            if self.psv_state == "open":
                eff_area = self.area
            else:
                eff_area = 0.0

        if eff_area > 0:
            return self._api_psv_release_rate(P1, T1, k, Z, MW)
        else:
            return 0.0

    def _api_psv_release_rate(self, P1, T1, k, Z, MW):
        """
        Calculate PSV vapor relief rate per API 520 Part I 2014.

        Implements API 520 equations 5, 9, 15, 18 for compressible flow
        through pressure safety valves.

        Parameters
        ----------
        P1 : float
            Upstream pressure [Pa]
        T1 : float
            Upstream temperature [K]
        k : float
            Heat capacity ratio Cp/Cv [-]
        Z : float
            Compressibility factor [-]
        MW : float
            Molecular weight [kg/mol]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        Pback = self.back_pressure

        # Convert units for API 520 equations
        P1_kPa = P1 / 1000  # Pa to kPa
        Pback_kPa = Pback / 1000  # Pa to kPa
        area_mm2 = self.area * 1e6  # m² to mm²
        MW_gmol = MW * 1000  # kg/mol to g/mol

        # API 520 critical flow coefficient
        C = 0.03948 * (k * (2 / (k + 1)) ** ((k + 1) / (k - 1))) ** 0.5

        # Check if flow is critical (choked) or subcritical
        critical_ratio = ((k + 1) / 2) ** (k / (k - 1))

        if P1_kPa / Pback_kPa > critical_ratio:
            # Critical flow (choked at throat)
            w = self.discharge_coef * area_mm2 * C * P1_kPa / math.sqrt(T1 * Z / MW_gmol)
        else:
            # Subcritical flow (not choked)
            r = Pback_kPa / P1_kPa
            # Subcritical flow correction factor f2
            f2 = (
                (k / (k - 1)) * r ** (2 / k) * (1 - r ** ((k - 1) / k)) / (1 - r)
            ) ** 0.5
            w = (
                self.discharge_coef
                * area_mm2
                * f2
                / (T1 * Z / (MW_gmol * P1_kPa * (P1_kPa - Pback_kPa))) ** 0.5
                / 17.9
            )

        # Convert from kg/h to kg/s
        return w / 3600

    def _relief_valve_flow(self, P1, T1, k, Z, MW):
        """
        Calculate relief valve flow (simple model without PSV hysteresis).

        For relief valves, flow is calculated using API methods but without
        the hysteresis behavior of PSVs.

        Parameters
        ----------
        P1 : float
            Upstream pressure [Pa]
        T1 : float
            Upstream temperature [K]
        k : float
            Heat capacity ratio Cp/Cv [-]
        Z : float
            Compressibility factor [-]
        MW : float
            Molecular weight [kg/mol]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        # For relief valves, calculate flow if above set pressure
        if P1 > self.set_pressure:
            # Use simplified API calculation (assuming area from set pressure)
            # This would need to be enhanced based on specific relief valve sizing
            return self._api_psv_release_rate(P1, T1, k, Z, MW)
        else:
            return 0.0

    def _control_valve_flow(self, P1, T1, k, Z, MW, t):
        """
        Calculate control valve flow using Cv characteristic.

        Implements control valve sizing equations with time-varying Cv based
        on valve characteristic (linear, equal percentage, or fast opening).

        Parameters
        ----------
        P1 : float
            Upstream pressure [Pa]
        T1 : float
            Upstream temperature [K]
        k : float
            Heat capacity ratio Cp/Cv [-]
        Z : float
            Compressibility factor [-]
        MW : float
            Molecular weight [kg/mol]
        t : float
            Current time [s]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        P2 = self.back_pressure

        # Get time-varying Cv based on valve characteristic
        Cv = self._cv_vs_time(t)

        # Control valve sizing equation for compressible flow
        # Based on ISA/IEC standards
        dP = P1 - P2

        if dP <= 0:
            return 0.0

        # Pressure drop ratio
        x = dP / P1

        # Check if flow is choked
        if x > self.xT:
            x = self.xT

        # Mass flow rate (adapted from ISA 75.01.01)
        # Simplified equation for gas flow
        Y = 1 - x / (3 * self.xT)  # Expansion factor

        # Convert Cv to metric units and calculate mass flow
        # Cv in US units, need to convert
        # Using simplified correlation
        mdot = 94.8 * Cv * self.FP * Y * math.sqrt(P1 * dP / (T1 * Z * MW))

        return mdot

    def _cv_vs_time(self, t):
        """
        Calculate control valve Cv vs time based on valve characteristic.

        Models three common valve characteristics:
        - Linear: Cv proportional to valve position
        - Equal percentage: Logarithmic Cv vs position
        - Fast/quick opening: Front-loaded Cv vs position

        Parameters
        ----------
        t : float
            Current time [s]

        Returns
        -------
        float
            Flow coefficient Cv at time t
        """
        Cv_max = self.Cv

        # Handle time constant (valve actuation time)
        if self.time_constant > 0:
            # Valve position increases from 0 to 1 over time_constant
            position = min(t / self.time_constant, 1.0)
        else:
            # Instantaneous opening
            position = 1.0

        # Apply valve characteristic
        if self.characteristic == "linear":
            Cv = Cv_max * position

        elif self.characteristic == "eq":
            # Equal percentage: exponential relationship
            # Cv = Cv_min * R^position, where R is rangeability (typically 50)
            R = 50  # Rangeability
            Cv_min = Cv_max / R
            Cv = Cv_min * R**position

        elif self.characteristic == "fast":
            # Fast/quick opening: square root relationship
            Cv = Cv_max * math.sqrt(position)

        else:
            # Default to linear
            Cv = Cv_max * position

        return Cv

    def _constant_mdot(self, t):
        """
        Return constant or time-varying mass flow rate.

        Parameters
        ----------
        t : float
            Current time [s]

        Returns
        -------
        float
            Mass flow rate [kg/s]
        """
        if self.time_profile is not None:
            # Time-varying mdot profile
            mdot_array = np.asarray(self.mdot_value)
            time_array = np.asarray(self.time_profile)
            return np.interp(t, time_array, mdot_array)
        else:
            # Constant mdot
            return self.mdot_value

    def reset_psv_state(self):
        """Reset PSV state to closed. Useful for initialization or resets."""
        self.psv_state = "closed"
