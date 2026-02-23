# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Thermodynamic solver for CoolProp fluid property calculations.

This module provides the ThermodynamicSolver class which encapsulates all
thermodynamic state calculations for HydDown simulations. It manages CoolProp
AbstractState objects and provides methods for solving thermodynamic problems
(PH-problems, UD-problems) for both single-component and multicomponent fluids.

The ThermodynamicSolver handles:
- CoolProp AbstractState initialization and management
- PH-problems: Finding temperature at constant pressure and enthalpy
- UD-problems: Finding pressure and temperature at constant internal energy and density
- Liquid level calculations for two-phase systems
- Single-component (fast, direct CoolProp calls) and multicomponent (slower, numerical optimization) fluids

Key Methods:
- PHproblem(): Solve for temperature at given P, H (isenthalpic expansion)
- UDproblem(): Solve for P, T at given U, density (constant internal energy)
- calc_liquid_level(): Calculate liquid height in two-phase systems

Typical usage:
    solver = ThermodynamicSolver(
        species="HEOS::Hydrogen",
        mole_fractions=[1.0],
        vessel_geometry=inner_vol
    )
    T = solver.PHproblem(H=1e6, P=1e5, Tguess=300)
"""

import numpy as np
from scipy.optimize import minimize, root_scalar
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from hyddown import safety_checks as sc


class ThermodynamicSolver:
    """
    Thermodynamic solver for CoolProp-based fluid property calculations.

    This class encapsulates all thermodynamic calculations required for HydDown
    simulations, managing CoolProp AbstractState objects and providing methods
    to solve for thermodynamic states given various property pairs.

    Parameters
    ----------
    species : str
        CoolProp fluid name (e.g., "HEOS::Hydrogen") or mixture
        (e.g., "HEOS::Methane[0.9]&Ethane[0.1]")
    mole_fractions : list of float
        Mole fractions for each component (must sum to 1.0)
    vessel_geometry : fluids.TANK, optional
        Vessel geometry object for liquid level calculations.
        Required if calc_liquid_level() will be used.
    species_SRK : str, optional
        CoolProp species string using SRK backend for transport properties.
        If not provided, defaults to species.

    Attributes
    ----------
    fluid : CoolProp.AbstractState
        Main fluid state for thermodynamic calculations
    vent_fluid : CoolProp.AbstractState
        Fluid state for vent/discharge calculations (gas phase)
    res_fluid : CoolProp.AbstractState
        Fluid state for reservoir/filling calculations
    transport_fluid : CoolProp.AbstractState
        Fluid state for transport property calculations (gas phase)
    transport_fluid_wet : CoolProp.AbstractState
        Fluid state for wetted region transport properties (liquid phase)
    is_multicomponent : bool
        True if fluid is multicomponent mixture
    """

    def __init__(self, species, mole_fractions, vessel_geometry=None, species_SRK=None):
        """
        Initialize the thermodynamic solver with fluid composition and geometry.

        Parameters
        ----------
        species : str
            CoolProp fluid name
        mole_fractions : list of float
            Mole fractions for each component
        vessel_geometry : fluids.TANK, optional
            Vessel geometry for liquid level calculations
        species_SRK : str, optional
            Species string for SRK backend (transport properties)
        """
        self.species = species
        self.mole_fractions = mole_fractions
        self.vessel_geometry = vessel_geometry
        self.is_multicomponent = "&" in species

        # Set up SRK species for transport properties
        if species_SRK is None:
            self.species_SRK = species
        else:
            self.species_SRK = species_SRK

        # Initialize CoolProp AbstractState objects
        self._initialize_fluid_states()

    def _initialize_fluid_states(self):
        """
        Create CoolProp AbstractState objects for different calculation needs.

        Creates:
        - fluid: Main thermodynamic state
        - vent_fluid: Vent/discharge calculations (gas phase)
        - res_fluid: Reservoir/filling calculations
        - transport_fluid: Transport properties (gas phase, SRK backend)
        - transport_fluid_wet: Transport properties (liquid phase, SRK backend)
        """
        # Main fluid state
        self.fluid = CP.AbstractState("HEOS", self.species)
        if self.is_multicomponent:
            self.fluid.specify_phase(CP.iphase_gas)
        self.fluid.set_mole_fractions(self.mole_fractions)

        # Vent fluid (for discharge calculations, gas phase)
        self.vent_fluid = CP.AbstractState("HEOS", self.species)
        self.vent_fluid.specify_phase(CP.iphase_gas)
        self.vent_fluid.set_mole_fractions(self.mole_fractions)

        # Reservoir fluid (for filling calculations)
        self.res_fluid = CP.AbstractState("HEOS", self.species)
        self.res_fluid.set_mole_fractions(self.mole_fractions)

        # Transport properties (gas phase, SRK backend for robustness)
        self.transport_fluid = CP.AbstractState("HEOS", self.species_SRK)
        self.transport_fluid.specify_phase(CP.iphase_gas)
        self.transport_fluid.set_mole_fractions(self.mole_fractions)

        # Transport properties (liquid phase for wetted regions)
        self.transport_fluid_wet = CP.AbstractState("HEOS", self.species_SRK)
        self.transport_fluid_wet.specify_phase(CP.iphase_liquid)
        self.transport_fluid_wet.set_mole_fractions(self.mole_fractions)

    def calc_liquid_level(self):
        """
        Calculate liquid level height based on current two-phase fluid state.

        For two-phase systems (0 ≤ quality ≤ 1), calculates the height of liquid
        phase in the vessel based on vapor quality, phase densities, and vessel geometry.
        Uses vessel geometry from fluids.TANK to convert liquid volume to height.

        Returns
        -------
        float
            Liquid level height from vessel bottom [m].
            Returns 0.0 for single-phase gas (quality > 1) or subcooled liquid (quality < 0).

        Raises
        ------
        ValueError
            If vessel_geometry was not provided during initialization
        """
        if self.vessel_geometry is None:
            raise ValueError(
                "vessel_geometry must be provided to ThermodynamicSolver "
                "for liquid level calculations"
            )

        # Check if in two-phase region
        quality = self.fluid.Q()
        if quality >= 0 and quality <= 1:
            # Get saturated phase densities
            rho_liq = self.fluid.saturated_liquid_keyed_output(CP.iDmass)
            rho_vap = self.fluid.saturated_vapor_keyed_output(CP.iDmass)

            # Calculate liquid mass and volume
            m_liq = self.fluid.rhomass() * self.vessel_geometry.V_total * (1 - quality)
            V_liq = m_liq / rho_liq

            # Convert volume to height using vessel geometry
            h_liq = self.vessel_geometry.h_from_V(V_liq)
            return h_liq
        else:
            return 0.0

    def PHres(self, T, P, H):
        """
        Residual enthalpy function to be minimized during PH-problem.

        Used by numerical optimizer (scipy.optimize.minimize) to find temperature
        that satisfies constant pressure-enthalpy constraints for multicomponent
        fluids. The optimizer adjusts T until residual approaches zero.

        Parameters
        ----------
        T : float or array-like
            Temperature estimate [K]. Optimizer may pass as array.
        P : float
            Pressure [Pa]
        H : float
            Target enthalpy [J/kg]

        Returns
        -------
        float
            Squared normalized enthalpy residual (dimensionless).
            Zero when correct temperature is found.
        """
        # Extract scalar from array (scipy optimizers pass arrays, CoolProp needs scalars)
        T_scalar = float(T.item()) if hasattr(T, "item") else float(T)
        self.vent_fluid.update(CP.PT_INPUTS, P, T_scalar)
        return ((H - self.vent_fluid.hmass()) / H) ** 2

    def PHres_relief(self, T, P, H):
        """
        Residual enthalpy function for PH-problem during relief valve calculations.

        Used by numerical optimizer (scipy.optimize.root_scalar) to find temperature
        for relief valve discharge calculations. Similar to PHres() but uses the
        main fluid state instead of vent_fluid state.

        Parameters
        ----------
        T : float or array-like
            Temperature estimate [K]
        P : float
            Pressure [Pa]
        H : float
            Target enthalpy [J/kg]

        Returns
        -------
        float
            Normalized enthalpy residual (dimensionless).
            Zero when correct temperature is found.
        """
        # Extract scalar from array (scipy optimizers pass arrays, CoolProp needs scalars)
        T_scalar = float(T.item()) if hasattr(T, "item") else float(T)
        self.fluid.update(CP.PT_INPUTS, P, T_scalar)
        return (H - self.fluid.hmass()) / H

    def PHproblem(self, H, P, Tguess, relief=False):
        """
        Solve constant pressure, constant enthalpy problem (isenthalpic expansion).

        Finds temperature at specified pressure and enthalpy. Typical use case is
        modeling adiabatic throttling/expansion through valves without work extraction.
        For multicomponent mixtures, uses numerical optimization to find temperature.
        For single component fluids, uses direct CoolProp methods for speed.

        Parameters
        ----------
        H : float
            Enthalpy [J/kg]
        P : float
            Pressure [Pa]
        Tguess : float
            Initial guess for temperature [K]
        relief : bool, optional
            If True, use relief valve solver (root_scalar with Newton method).
            If False, use general solver (minimize with Nelder-Mead).
            Default is False.

        Returns
        -------
        float
            Temperature at (P, H) [K]

        Raises
        ------
        ThermodynamicConvergenceError
            If numerical optimization fails to converge for multicomponent fluid
        """
        # Multicomponent case: requires numerical optimization
        if self.is_multicomponent:
            x0 = Tguess

            if not relief:
                # Use Nelder-Mead minimization
                res = minimize(
                    self.PHres,
                    x0,
                    args=(P, H),
                    method="Nelder-Mead",
                    options={"xatol": 0.1, "fatol": 0.001},
                )
                # Check convergence
                sc.check_optimization_convergence(
                    res, solver_name="PHproblem", state_vars={"P": P, "H": H, "T_guess": x0}
                )
                T1 = res.x[0]
            else:
                # Use Newton root finding for relief valve
                res = root_scalar(
                    self.PHres_relief,
                    args=(P, H),
                    x0=x0,
                    method="newton",
                )
                # Check convergence
                sc.check_optimization_convergence(
                    res, solver_name="PHproblem_relief", state_vars={"P": P, "H": H}
                )
                T1 = res.root

        # Single component case: direct CoolProp calculation
        else:
            T1 = PropsSI("T", "P", P, "H", H, self.species)

        return T1

    def UDres(self, x, U, rho):
        """
        Residual function for UD-problem (constant internal energy and density).

        Used by numerical optimizer to find pressure and temperature that satisfy
        constant internal energy and density constraints. Minimizes sum of squared
        normalized residuals for both U and rho.

        Parameters
        ----------
        x : array-like of float
            [Pressure [Pa], Temperature [K]]
        U : float
            Target internal energy [J/kg]
        rho : float
            Target density [kg/m³]

        Returns
        -------
        float
            Sum of squared normalized residuals (dimensionless)
        """
        self.fluid.update(CP.PT_INPUTS, x[0], x[1])
        return ((U - self.fluid.umass()) / U) ** 2 + ((rho - self.fluid.rhomass()) / rho) ** 2

    def UDproblem(self, U, rho, Pguess, Tguess):
        """
        Solve constant internal energy, constant density problem.

        Finds pressure and temperature at specified internal energy and density.
        Relevant for 1st law of thermodynamics with constant volume. For multicomponent
        mixtures, uses numerical optimization to find P and T. For single component
        fluids, uses direct CoolProp methods for speed.

        Parameters
        ----------
        U : float
            Internal energy [J/kg]
        rho : float
            Density [kg/m³]
        Pguess : float
            Initial guess for pressure [Pa]
        Tguess : float
            Initial guess for temperature [K]

        Returns
        -------
        P1 : float
            Pressure at (U, rho) [Pa]
        T1 : float
            Temperature at (U, rho) [K]
        Ures : float
            Internal energy residual [J/kg]. Zero for single component fluids.

        Raises
        ------
        ThermodynamicConvergenceError
            If numerical optimization fails to converge for multicomponent fluid
        """
        # Multicomponent case: requires numerical optimization
        if self.is_multicomponent:
            x0 = [Pguess, Tguess]
            res = minimize(
                self.UDres,
                x0,
                args=(U, rho),
                method="Nelder-Mead",
                options={"xatol": 0.1, "fatol": 0.001},
            )
            # Check convergence
            sc.check_optimization_convergence(
                res,
                solver_name="UDproblem",
                state_vars={"U": U, "rho": rho, "P_guess": Pguess, "T_guess": Tguess},
            )
            P1 = res.x[0]
            T1 = res.x[1]
            Ures = U - self.fluid.umass()

        # Single component case: direct CoolProp calculation
        else:
            P1 = PropsSI("P", "D", rho, "U", U, self.species)
            T1 = PropsSI("T", "D", rho, "U", U, self.species)
            Ures = 0

        return P1, T1, Ures
