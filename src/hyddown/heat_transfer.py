# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Heat transfer models for vessel wall calculations.

This module provides heat transfer calculation classes for vessel walls in
HydDown simulations. It handles convective heat transfer, fire scenarios,
and optional detailed wall conduction modeling.

The module implements a strategy pattern where different external heat sources
(convection, fire, etc.) share common internal calculations but differ in how
external heat flux is calculated.

Key Features:
- Unified interface for all heat transfer types
- Two-phase system support (wetted vs unwetted regions)
- Optional detailed 1-D transient wall conduction
- Internal convection correlations (natural/forced)
- External boundary conditions (convection, fire)

Typical usage:
    model = ConvectiveHeatTransfer(
        vessel_geometry=inner_vol,
        h_out=10.0,
        Tamb=300,
        vessel_cp=500,
        vessel_density=7800,
        vol_solid=0.05,
        surf_area_inner=10.0,
        surf_area_outer=11.0
    )
    model.calculate(i, time[i], P, T_fluid, ...)
"""

import math
import numpy as np
import CoolProp.CoolProp as CP
from hyddown import transport as tp
from hyddown import fire
from hyddown import thermesh as tm
from hyddown import safety_checks as sc


class WallHeatTransfer:
    """
    Base class for vessel wall heat transfer calculations.

    This class provides common functionality for calculating heat transfer
    between vessel walls and fluid, including:
    - Internal convection (wall to fluid)
    - External heat sources (convection, fire, etc.)
    - Wall temperature evolution (lumped capacitance or detailed conduction)
    - Two-phase system handling (wetted vs unwetted regions)

    Subclasses must implement _calculate_external_heat_flux() to define
    the specific external boundary condition (convection, fire, etc.).

    Parameters
    ----------
    vessel_geometry : TANK
        Vessel geometry object from fluids library (provides SA_from_h, etc.)
    vessel_orientation : str
        Vessel orientation: "horizontal" or "vertical"
    diameter : float
        Vessel diameter [m]
    length : float
        Vessel length [m]
    vessel_cp : float
        Vessel material specific heat capacity [J/kg/K]
    vessel_density : float
        Vessel material density [kg/m³]
    vol_solid : float
        Volume of vessel wall material [m³]
    surf_area_inner : float
        Inner surface area [m²]
    surf_area_outer : float
        Outer surface area [m²]
    tstep : float
        Time step [s]
    species : str
        CoolProp fluid name
    flow_direction : str
        Flow direction: "discharge" or "filling"
    h_in : float or str
        Internal heat transfer coefficient [W/m²/K] or "calc" to calculate
    use_detailed_wall_model : bool, optional
        Enable detailed 1-D wall conduction (requires thermal_conductivity). Default False.
    thermal_conductivity : float, optional
        Wall thermal conductivity [W/m/K]. Required if use_detailed_wall_model=True.
    thickness : float, optional
        Wall thickness [m]. Required if use_detailed_wall_model=True.
    T0 : float, optional
        Initial temperature [K]. Required if use_detailed_wall_model=True.
    Tamb : float, optional
        Ambient temperature [K]. Required if use_detailed_wall_model=True.
    liner_thermal_conductivity : float, optional
        Liner thermal conductivity [W/m/K]. For multi-layer walls.
    liner_density : float, optional
        Liner density [kg/m³]. For multi-layer walls.
    liner_heat_capacity : float, optional
        Liner specific heat capacity [J/kg/K]. For multi-layer walls.
    liner_thickness : float, optional
        Liner thickness [m]. For multi-layer walls.

    Attributes
    ----------
    psv_state : str
        PSV state tracking (not used in heat transfer, but maintained for compatibility)
    T_profile : int or ndarray
        Temperature profile through wall (unwetted region). Initialized to 0.
    T_profile_w : int or ndarray
        Temperature profile through wall (wetted region). Initialized to 0.
    T_profile2 : int or ndarray
        Temperature profile for multi-layer wall (unwetted). Initialized to 0.
    T_profile2_w : int or ndarray
        Temperature profile for multi-layer wall (wetted). Initialized to 0.
    """

    def __init__(
        self,
        vessel_geometry,
        vessel_orientation,
        diameter,
        length,
        vessel_cp,
        vessel_density,
        vol_solid,
        surf_area_inner,
        surf_area_outer,
        tstep,
        species,
        flow_direction,
        h_in="calc",
        use_detailed_wall_model=False,
        thermal_conductivity=None,
        thickness=None,
        T0=None,
        Tamb=None,
        liner_thermal_conductivity=None,
        liner_density=None,
        liner_heat_capacity=None,
        liner_thickness=None,
    ):
        """Initialize wall heat transfer model with vessel parameters."""
        self.inner_vol = vessel_geometry
        self.vessel_orientation = vessel_orientation
        self.diameter = diameter
        self.length = length
        self.vessel_cp = vessel_cp
        self.vessel_density = vessel_density
        self.vol_solid = vol_solid
        self.surf_area_inner = surf_area_inner
        self.surf_area_outer = surf_area_outer
        self.tstep = tstep
        self.species = species
        self.flow_direction = flow_direction
        self.h_in = h_in

        # Detailed wall model parameters
        self.use_detailed_wall_model = use_detailed_wall_model
        self.thermal_conductivity = thermal_conductivity
        self.thickness = thickness
        self.T0 = T0
        self.Tamb = Tamb

        # Multi-layer wall parameters (liner + shell)
        self.liner_thermal_conductivity = liner_thermal_conductivity
        self.liner_density = liner_density
        self.liner_heat_capacity = liner_heat_capacity
        self.liner_thickness = liner_thickness
        self.has_liner = liner_thermal_conductivity is not None

        # Validate detailed wall model parameters
        if use_detailed_wall_model:
            if any(
                p is None
                for p in [thermal_conductivity, thickness, T0, Tamb, vessel_cp, vessel_density]
            ):
                raise ValueError(
                    "Detailed wall model requires: thermal_conductivity, thickness, "
                    "T0, Tamb, vessel_cp, vessel_density"
                )

            # Validate multi-layer parameters if liner is specified
            if self.has_liner:
                if any(
                    p is None
                    for p in [liner_density, liner_heat_capacity, liner_thickness]
                ):
                    raise ValueError(
                        "Multi-layer wall model requires: liner_thermal_conductivity, "
                        "liner_density, liner_heat_capacity, liner_thickness"
                    )

        # Initialize thermesh variables for detailed model
        if use_detailed_wall_model:
            if self.has_liner:
                # Multi-layer wall
                self.T_profile2 = 0  # Will be initialized in first call
                self.T_profile2_w = 0
            else:
                # Single-layer wall
                self.T_profile = 0  # Will be initialized in first call
                self.T_profile_w = 0
            self.temp_profile = []
            self.z = None  # Node locations through thickness

    def calculate(
        self,
        i,
        t,
        liquid_level,
        T_fluid,
        P,
        T_vessel,
        T_vessel_wetted,
        T_inner_wall,
        T_inner_wall_wetted,
        T_outer_wall,
        T_outer_wall_wetted,
        fluid,
        transport_fluid,
        transport_fluid_wet,
        mass_rate,
    ):
        """
        Calculate heat transfer for current time step.

        This is the main template method that orchestrates the heat transfer
        calculation sequence. It calls various helper methods in order:
        1. Calculate characteristic length
        2. Calculate wetted area
        3. Calculate internal heat transfer coefficients
        4. Calculate internal heat flux
        5. Calculate external heat flux (strategy - varies by subclass)
        6. Update wall temperatures
        7. Optionally solve detailed wall conduction

        Parameters
        ----------
        i : int
            Current time step index
        t : float
            Current time [s]
        liquid_level : array
            Liquid level history [m]
        T_fluid : array
            Fluid temperature history [K]
        P : array
            Pressure history [Pa]
        T_vessel : array
            Vessel wall temperature (unwetted) history [K]
        T_vessel_wetted : array
            Vessel wall temperature (wetted) history [K]
        T_inner_wall : array
            Inner wall temperature (unwetted) history [K]
        T_inner_wall_wetted : array
            Inner wall temperature (wetted) history [K]
        T_outer_wall : array
            Outer wall temperature (unwetted) history [K]
        T_outer_wall_wetted : array
            Outer wall temperature (wetted) history [K]
        fluid : AbstractState
            CoolProp AbstractState for fluid
        transport_fluid : AbstractState
            CoolProp AbstractState for transport property calculations
        transport_fluid_wet : AbstractState
            CoolProp AbstractState for wetted region transport properties
        mass_rate : array
            Mass flow rate history [kg/s]

        Returns
        -------
        dict
            Dictionary containing calculated values:
            - 'h_inside': Internal heat transfer coefficient [W/m²/K]
            - 'h_inside_wetted': Internal heat transfer coefficient wetted [W/m²/K]
            - 'Q_inner': Internal heat flux unwetted [W]
            - 'Q_inner_wetted': Internal heat flux wetted [W]
            - 'q_inner': Internal heat flux density unwetted [W/m²]
            - 'q_inner_wetted': Internal heat flux density wetted [W/m²]
            - 'Q_outer': External heat flux unwetted [W]
            - 'Q_outer_wetted': External heat flux wetted [W]
            - 'q_outer': External heat flux density unwetted [W/m²]
            - 'q_outer_wetted': External heat flux density wetted [W/m²]
            - 'T_vessel': Updated vessel temperature unwetted [K]
            - 'T_vessel_wetted': Updated vessel temperature wetted [K]
            - 'T_inner_wall': Updated inner wall temperature unwetted [K]
            - 'T_inner_wall_wetted': Updated inner wall temperature wetted [K]
            - 'T_outer_wall': Updated outer wall temperature unwetted [K]
            - 'T_outer_wall_wetted': Updated outer wall temperature wetted [K]
        """
        # Step 1: Calculate characteristic length
        L = self._calculate_characteristic_length()

        # Step 2: Calculate wetted area for two-phase systems
        wetted_area = self._calculate_wetted_area(liquid_level, i)

        # Step 3: Calculate internal heat transfer coefficients
        hi, hiw = self._calculate_internal_convection(
            i, L, T_fluid, P, T_inner_wall, T_inner_wall_wetted,
            fluid, transport_fluid, transport_fluid_wet, mass_rate
        )

        # Step 4: Calculate internal heat flux (wall to fluid)
        Q_inner, Q_inner_wetted, q_inner, q_inner_wetted = self._calculate_internal_heat_flux(
            i, wetted_area, hi, hiw, T_inner_wall, T_inner_wall_wetted, T_fluid
        )

        # Step 5: Calculate external heat flux (STRATEGY - varies by subclass)
        Q_outer, Q_outer_wetted, q_outer, q_outer_wetted = self._calculate_external_heat_flux(
            i, wetted_area, T_vessel, T_vessel_wetted, T_outer_wall, T_outer_wall_wetted
        )

        # Step 6: Update wall temperatures (energy balance)
        T_vessel_new, T_vessel_wetted_new = self._update_wall_temperature(
            i, wetted_area, Q_inner, Q_outer, Q_inner_wetted, Q_outer_wetted,
            T_vessel, T_vessel_wetted, liquid_level
        )

        # Step 7: Detailed wall conduction (if enabled)
        if self.use_detailed_wall_model:
            result = self._solve_wall_conduction(
                i, wetted_area, Q_inner, Q_outer, Q_inner_wetted, Q_outer_wetted
            )
            T_inner_wall_new = result[0]
            T_outer_wall_new = result[1]
            T_inner_wall_wetted_new = result[2]
            T_outer_wall_wetted_new = result[3]
            T_bonded_wall_new = result[4]
            T_bonded_wall_wetted_new = result[5]
        else:
            # Lumped capacitance: inner/outer = vessel temperature
            T_inner_wall_new = T_vessel_new
            T_outer_wall_new = T_vessel_new
            T_inner_wall_wetted_new = T_vessel_wetted_new
            T_outer_wall_wetted_new = T_vessel_wetted_new
            T_bonded_wall_new = None
            T_bonded_wall_wetted_new = None

        # Return all calculated values
        return {
            'h_inside': hi,
            'h_inside_wetted': hiw,
            'Q_inner': Q_inner,
            'Q_inner_wetted': Q_inner_wetted,
            'q_inner': q_inner,
            'q_inner_wetted': q_inner_wetted,
            'Q_outer': Q_outer,
            'Q_outer_wetted': Q_outer_wetted,
            'q_outer': q_outer,
            'q_outer_wetted': q_outer_wetted,
            'T_vessel': T_vessel_new,
            'T_vessel_wetted': T_vessel_wetted_new,
            'T_inner_wall': T_inner_wall_new,
            'T_inner_wall_wetted': T_inner_wall_wetted_new,
            'T_outer_wall': T_outer_wall_new,
            'T_outer_wall_wetted': T_outer_wall_wetted_new,
            'T_bonded_wall': T_bonded_wall_new,
            'T_bonded_wall_wetted': T_bonded_wall_wetted_new,
        }

    def _calculate_characteristic_length(self):
        """
        Calculate characteristic length for convection correlations.

        For horizontal vessels, use diameter.
        For vertical vessels, use length.

        Returns
        -------
        float
            Characteristic length [m]
        """
        if self.vessel_orientation == "horizontal":
            return self.diameter
        else:
            return self.length

    def _calculate_wetted_area(self, liquid_level, i):
        """
        Calculate wall surface area in contact with liquid phase.

        For two-phase systems, the vessel wall is divided into:
        - Wetted region: In contact with liquid
        - Unwetted region: In contact with vapor

        Parameters
        ----------
        liquid_level : array
            Liquid level history [m]
        i : int
            Current time step index

        Returns
        -------
        float
            Wetted surface area [m²]
        """
        wetted_area = self.inner_vol.SA_from_h(liquid_level[i - 1])
        if np.isnan(wetted_area):
            wetted_area = 0
        return wetted_area

    def _calculate_internal_convection(
        self, i, L, T_fluid, P, T_inner_wall, T_inner_wall_wetted,
        fluid, transport_fluid, transport_fluid_wet, mass_rate
    ):
        """
        Calculate internal heat transfer coefficients (wall to fluid).

        Handles both discharge and filling modes, as well as two-phase systems.

        Parameters
        ----------
        i : int
            Current time step index
        L : float
            Characteristic length [m]
        T_fluid : array
            Fluid temperature history [K]
        P : array
            Pressure history [Pa]
        T_inner_wall : array
            Inner wall temperature (unwetted) history [K]
        T_inner_wall_wetted : array
            Inner wall temperature (wetted) history [K]
        fluid : AbstractState
            CoolProp AbstractState for fluid
        transport_fluid : AbstractState
            CoolProp AbstractState for transport property calculations
        transport_fluid_wet : AbstractState
            CoolProp AbstractState for wetted region transport properties
        mass_rate : array
            Mass flow rate history [kg/s]

        Returns
        -------
        tuple of float
            (hi, hiw) - Internal heat transfer coefficients [W/m²/K]
            hi: Gas/vapor region coefficient
            hiw: Liquid (wetted) region coefficient
        """
        if self.h_in != "calc":
            # User-specified heat transfer coefficient
            hi = self.h_in
            hiw = self.h_in
        else:
            # Calculate heat transfer coefficient
            if self.flow_direction == "filling":
                # Filling mode: forced convection
                T_film = (T_fluid[i - 1] + T_inner_wall[i - 1]) / 2
                transport_fluid.update(CP.PT_INPUTS, P[i - 1], T_film)

                hi = tp.h_inside_mixed(
                    L,
                    T_inner_wall[i - 1],
                    T_fluid[i - 1],
                    transport_fluid,
                    mass_rate,
                    self.diameter,
                )

                # Wetted region
                if fluid.Q() >= 0 and fluid.Q() <= 1:
                    transport_fluid_wet.update(CP.PT_INPUTS, P[i - 1], T_film)
                    hiw = tp.h_inside_wetted(
                        L,
                        T_inner_wall_wetted[i - 1],
                        T_fluid[i - 1],
                        transport_fluid_wet,
                        fluid,
                    )
                else:
                    hiw = hi

            else:
                # Discharge mode: natural convection
                T_film = (T_fluid[i - 1] + T_inner_wall[i - 1]) / 2
                try:
                    transport_fluid.update(CP.PT_INPUTS, P[i - 1], T_film)
                    if fluid.Q() >= 0 and fluid.Q() <= 1:
                        transport_fluid_wet.update(
                            CP.PT_INPUTS, P[i - 1], T_fluid[i - 1]
                        )
                except:
                    transport_fluid.update(CP.PQ_INPUTS, P[i - 1], 1.0)

                hi = tp.h_inside(
                    L,
                    T_inner_wall[i - 1],
                    T_fluid[i - 1],
                    transport_fluid,
                )

                # Wetted region
                if fluid.Q() >= 0 and fluid.Q() <= 1:
                    hiw = tp.h_inside_wetted(
                        L,
                        T_inner_wall_wetted[i - 1],
                        T_fluid[i - 1],
                        transport_fluid_wet,
                        fluid,
                    )
                else:
                    hiw = hi

        return hi, hiw

    def _calculate_internal_heat_flux(
        self, i, wetted_area, hi, hiw, T_inner_wall, T_inner_wall_wetted, T_fluid
    ):
        """
        Calculate internal heat flux from wall to fluid.

        Uses Newton's law of cooling: Q = h * A * ΔT

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        hi : float
            Internal heat transfer coefficient unwetted [W/m²/K]
        hiw : float
            Internal heat transfer coefficient wetted [W/m²/K]
        T_inner_wall : array
            Inner wall temperature (unwetted) history [K]
        T_inner_wall_wetted : array
            Inner wall temperature (wetted) history [K]
        T_fluid : array
            Fluid temperature history [K]

        Returns
        -------
        tuple of float
            (Q_inner, Q_inner_wetted, q_inner, q_inner_wetted)
            Q_*: Total heat flux [W]
            q_*: Heat flux density [W/m²]
        """
        # Unwetted (gas/vapor) region
        Q_inner = (
            (self.surf_area_inner - wetted_area)
            * hi
            * (T_inner_wall[i - 1] - T_fluid[i - 1])
        )
        q_inner = hi * (T_inner_wall[i - 1] - T_fluid[i - 1])

        # Wetted (liquid) region
        Q_inner_wetted = (
            wetted_area
            * hiw
            * (T_inner_wall_wetted[i - 1] - T_fluid[i - 1])
        )

        # Avoid division by zero when fully vapor (wetted_area = 0)
        if wetted_area > 0:
            q_inner_wetted = Q_inner_wetted / wetted_area
        else:
            q_inner_wetted = 0.0

        if np.isnan(Q_inner_wetted):
            Q_inner_wetted = 0

        return Q_inner, Q_inner_wetted, q_inner, q_inner_wetted

    def _calculate_external_heat_flux(
        self, i, wetted_area, T_vessel, T_vessel_wetted, T_outer_wall, T_outer_wall_wetted
    ):
        """
        Calculate external heat flux from environment to vessel.

        This is an ABSTRACT method that must be implemented by subclasses.
        Different subclasses provide different external boundary conditions:
        - ConvectiveHeatTransfer: h_out * (T_amb - T_wall)
        - FireHeatTransfer: fire.sb_fire(T_wall, fire_type)

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        T_vessel : array
            Vessel wall temperature (unwetted) history [K]
        T_vessel_wetted : array
            Vessel wall temperature (wetted) history [K]
        T_outer_wall : array
            Outer wall temperature (unwetted) history [K]
        T_outer_wall_wetted : array
            Outer wall temperature (wetted) history [K]

        Returns
        -------
        tuple of float
            (Q_outer, Q_outer_wetted, q_outer, q_outer_wetted)
            Q_*: Total heat flux [W]
            q_*: Heat flux density [W/m²]

        Raises
        ------
        NotImplementedError
            If not overridden by subclass
        """
        raise NotImplementedError(
            "Subclass must implement _calculate_external_heat_flux()"
        )

    def _update_wall_temperature(
        self, i, wetted_area, Q_inner, Q_outer, Q_inner_wetted, Q_outer_wetted,
        T_vessel, T_vessel_wetted, liquid_level
    ):
        """
        Update vessel wall temperature using energy balance.

        Uses lumped capacitance model: dT/dt = (Q_in - Q_out) / (m * cp)

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        Q_inner : float
            Internal heat flux unwetted [W]
        Q_outer : float
            External heat flux unwetted [W]
        Q_inner_wetted : float
            Internal heat flux wetted [W]
        Q_outer_wetted : float
            External heat flux wetted [W]
        T_vessel : array
            Vessel wall temperature (unwetted) history [K]
        T_vessel_wetted : array
            Vessel wall temperature (wetted) history [K]
        liquid_level : array
            Liquid level history [m]

        Returns
        -------
        tuple of float
            (T_vessel_new, T_vessel_wetted_new) - Updated temperatures [K]
        """
        # Unwetted region energy balance
        T_vessel_new = T_vessel[i - 1] + (Q_outer - Q_inner) * self.tstep / (
            self.vessel_cp
            * self.vessel_density
            * self.vol_solid
            * (self.inner_vol.A - wetted_area)
            / self.inner_vol.A
        )

        # Wetted region energy balance
        if Q_outer_wetted == 0 and Q_inner_wetted == 0:
            T_vessel_wetted_new = T_vessel_wetted[i - 1]
        else:
            if wetted_area > 0:
                T_vessel_wetted_new = T_vessel_wetted[i - 1] + (
                    Q_outer_wetted - Q_inner_wetted
                ) * self.tstep / (
                    self.vessel_cp
                    * self.vessel_density
                    * self.vol_solid
                    * wetted_area
                    / self.inner_vol.A
                )
            else:
                T_vessel_wetted_new = T_vessel_wetted[i - 1]

        if np.isnan(T_vessel_wetted_new):
            T_vessel_wetted_new = T_vessel_new

        return T_vessel_new, T_vessel_wetted_new

    def _solve_wall_conduction(self, i, wetted_area, Q_inner, Q_outer, Q_inner_wetted, Q_outer_wetted):
        """
        Solve 1-D transient heat conduction through vessel wall.

        Uses thermesh finite element solver for detailed temperature distribution
        through wall thickness. Supports both single-layer and multi-layer
        (liner + shell) walls.

        This method is only called if use_detailed_wall_model=True.

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        Q_inner : float
            Internal heat flux unwetted [W]
        Q_outer : float
            External heat flux unwetted [W]
        Q_inner_wetted : float
            Internal heat flux wetted [W]
        Q_outer_wetted : float
            External heat flux wetted [W]

        Returns
        -------
        tuple of float
            (T_inner_wall, T_outer_wall, T_inner_wall_wetted, T_outer_wall_wetted)
            Wall surface temperatures [K]
        """
        # Crank-Nicolson scheme parameters
        theta = 0.5  # Unconditionally stable, 2nd order accurate
        dt = self.tstep / 10  # Sub-step for thermal solver (finer time resolution)

        k = self.thermal_conductivity
        rho = self.vessel_density
        cp = self.vessel_cp

        if not self.has_liner:
            # ================================================================
            # SINGLE-LAYER WALL
            # ================================================================
            nn = 11  # Number of nodes through wall thickness
            z = np.linspace(0, self.thickness, nn)
            self.z = z

            # Create meshes for unwetted and wetted regions
            mesh = tm.Mesh(z, tm.LinearElement)
            mesh_w = tm.Mesh(z, tm.LinearElement)

            # Material model with constant properties
            material = tm.isothermal_model(k, rho, cp)
            material_w = tm.isothermal_model(k, rho, cp)

            # Initialize temperature profile on first time step
            if type(self.T_profile) == type(int()) and self.T_profile == 0:
                # Initial steady-state solution
                bc = [
                    {"T": self.T0},
                    {"T": self.Tamb},
                ]
                domain = tm.Domain(mesh, [material], bc)
                domain.set_T((self.Tamb + self.T0) / 2 * np.ones(len(mesh.nodes)))
                solver = {
                    "dt": 100,
                    "t_end": 10000,
                    "theta": theta,
                }
                t_bonded, self.T_profile = tm.solve_ht(domain, solver)

                # Wetted region initial solution
                bc_w = [
                    {"T": self.T0},
                    {"T": self.Tamb},
                ]
                domain_w = tm.Domain(mesh_w, [material_w], bc_w)
                domain_w.set_T((self.Tamb + self.T0) / 2 * np.ones(len(mesh_w.nodes)))
                solver_w = {
                    "dt": 100,
                    "t_end": 10000,
                    "theta": theta,
                }
                t_bonded_w, self.T_profile_w = tm.solve_ht(domain_w, solver_w)
            else:
                # Transient solution with heat flux boundary conditions
                # Unwetted region
                bc = [
                    {
                        "q": Q_outer / (self.surf_area_inner - wetted_area) if (self.surf_area_inner - wetted_area) > 0 else 0.0
                    },
                    {
                        "q": -Q_inner / (self.surf_area_inner - wetted_area) if (self.surf_area_inner - wetted_area) > 0 else 0.0
                    },
                ]
                domain = tm.Domain(mesh, [material], bc)
                domain.set_T(self.T_profile[-1, :])
                solver = {
                    "dt": dt,
                    "t_end": self.tstep,
                    "theta": theta,
                }
                t_bonded, self.T_profile = tm.solve_ht(domain, solver)

                # Wetted region
                bc_w = [
                    {"q": Q_outer_wetted / wetted_area if wetted_area > 0 else 0.0},
                    {"q": -Q_inner_wetted / wetted_area if wetted_area > 0 else 0.0},
                ]
                domain_w = tm.Domain(mesh_w, [material_w], bc_w)
                domain_w.set_T(self.T_profile_w[-1, :])
                solver_w = {
                    "dt": dt,
                    "t_end": self.tstep,
                    "theta": theta,
                }
                t_bonded_w, self.T_profile_w = tm.solve_ht(domain_w, solver_w)

            # Final solve for current timestep
            solver = {"dt": dt, "t_end": self.tstep, "theta": theta}
            t, self.T_profile = tm.solve_ht(domain, solver)
            solver_w = {"dt": dt, "t_end": self.tstep, "theta": theta}
            t_w, self.T_profile_w = tm.solve_ht(domain_w, solver_w)

            # Store temperature profile
            self.temp_profile.append(self.T_profile[-1, :])

            # Extract wall surface temperatures
            T_outer_wall = self.T_profile[-1, 0]
            T_inner_wall = self.T_profile[-1, -1]
            T_outer_wall_wetted = self.T_profile_w[-1, 0]
            T_inner_wall_wetted = self.T_profile_w[-1, -1]
            T_bonded_wall = None  # No bonded interface for single-layer
            T_bonded_wall_wetted = None

        else:
            # ================================================================
            # MULTI-LAYER WALL (Liner + Shell)
            # ================================================================
            k_liner = self.liner_thermal_conductivity
            rho_liner = self.liner_density
            cp_liner = self.liner_heat_capacity

            # Material models
            liner = tm.isothermal_model(k_liner, rho_liner, cp_liner)
            shell = tm.isothermal_model(k, rho, cp)
            liner_w = tm.isothermal_model(k_liner, rho_liner, cp_liner)
            shell_w = tm.isothermal_model(k, rho, cp)

            # Create composite mesh
            nn = 11  # Number of nodes per layer
            z_shell = np.linspace(0, self.thickness, nn)
            z_liner = np.linspace(-self.liner_thickness, 0, nn)
            z2 = np.hstack((z_liner, z_shell[1:]))  # Combine meshes
            self.z = z2

            # Create meshes for unwetted and wetted regions
            mesh2 = tm.Mesh(z2, tm.LinearElement)
            mesh2_w = tm.Mesh(z2, tm.LinearElement)

            # Mark shell subdomain (positive z)
            for j, elem in enumerate(mesh2.elem):
                if elem.nodes.mean() > 0.0:
                    mesh2.subdomain[j] = 1
                    mesh2_w.subdomain[j] = 1

            # Initialize temperature profile on first time step
            if type(self.T_profile2) == type(int()) and self.T_profile2 == 0:
                # Initial steady-state solution
                bc = [
                    {"T": self.T0},
                    {"T": self.Tamb},
                ]
                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                domain2.set_T((self.Tamb + self.T0) / 2 * np.ones(len(mesh2.nodes)))
                solver2 = {
                    "dt": 100,
                    "t_end": 10000,
                    "theta": theta,
                }
                t_bonded, self.T_profile2 = tm.solve_ht(domain2, solver2)

                # Wetted region initial solution
                bc_w = [
                    {"T": self.T0},
                    {"T": self.Tamb},
                ]
                domain2_w = tm.Domain(mesh2_w, [liner_w, shell_w], bc_w)
                domain2_w.set_T((self.Tamb + self.T0) / 2 * np.ones(len(mesh2.nodes)))
                solver2_w = {
                    "dt": 100,
                    "t_end": 10000,
                    "theta": theta,
                }
                t_bonded_w, self.T_profile2_w = tm.solve_ht(domain2_w, solver2_w)
            else:
                # Transient solution with heat flux boundary conditions
                # Unwetted region - use safe_divide to handle fully liquid case
                bc = [
                    {
                        "q": sc.safe_divide(
                            -Q_inner,
                            self.surf_area_inner - wetted_area,
                            name="q_inner_unwetted",
                            default=0.0
                        )
                    },
                    {
                        "q": sc.safe_divide(
                            Q_outer,
                            self.surf_area_outer - wetted_area,
                            name="q_outer_unwetted",
                            default=0.0
                        )
                    },
                ]
                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                domain2.set_T(self.T_profile2[-1, :])
                solver2 = {
                    "dt": dt,
                    "t_end": self.tstep,
                    "theta": theta,
                }
                t_bonded, self.T_profile2 = tm.solve_ht(domain2, solver2)

                # Wetted region - use safe_divide to handle fully vapor case
                bc_w = [
                    {
                        "q": sc.safe_divide(
                            -Q_inner_wetted,
                            wetted_area,
                            name="q_inner_wetted",
                            default=0.0
                        )
                    },
                    {
                        "q": sc.safe_divide(
                            Q_outer_wetted,
                            wetted_area,
                            name="q_outer_wetted",
                            default=0.0
                        )
                    },
                ]
                domain2_w = tm.Domain(mesh2_w, [liner_w, shell_w], bc_w)
                domain2_w.set_T(self.T_profile2_w[-1, :])
                solver2_w = {
                    "dt": dt,
                    "t_end": self.tstep,
                    "theta": theta,
                }
                t_bonded_w, self.T_profile2_w = tm.solve_ht(domain2_w, solver2_w)

            # Extract wall surface temperatures
            T_outer_wall = self.T_profile2[-1, -1]  # Shell outer surface
            T_inner_wall = self.T_profile2[-1, 0]    # Liner inner surface
            T_bonded_wall = self.T_profile2[-1, nn - 1]  # Liner/shell interface
            T_outer_wall_wetted = self.T_profile2_w[-1, -1]
            T_inner_wall_wetted = self.T_profile2_w[-1, 0]
            T_bonded_wall_wetted = self.T_profile2_w[-1, nn - 1]

            # Store temperature profile
            self.temp_profile.append(self.T_profile2[-1, :])

        return (
            T_inner_wall, T_outer_wall, T_inner_wall_wetted, T_outer_wall_wetted,
            T_bonded_wall, T_bonded_wall_wetted
        )


class ConvectiveHeatTransfer(WallHeatTransfer):
    """
    Heat transfer model with external convection boundary condition.

    Used for "specified_h" and "detailed" heat transfer types where external
    heat transfer is via natural/forced convection to ambient conditions.

    External heat flux: Q = h_out * A * (T_amb - T_wall)

    Parameters
    ----------
    h_out : float
        External heat transfer coefficient [W/m²/K]
    Tamb : float
        Ambient temperature [K]
    **kwargs
        Additional parameters passed to WallHeatTransfer base class

    Examples
    --------
    >>> model = ConvectiveHeatTransfer(
    ...     vessel_geometry=inner_vol,
    ...     vessel_orientation="vertical",
    ...     diameter=0.273,
    ...     length=1.524,
    ...     vessel_cp=500,
    ...     vessel_density=7800,
    ...     vol_solid=0.05,
    ...     surf_area_inner=1.5,
    ...     surf_area_outer=1.7,
    ...     tstep=0.1,
    ...     species="HEOS::Hydrogen",
    ...     flow_direction="discharge",
    ...     h_out=10.0,
    ...     Tamb=300
    ... )
    """

    def __init__(self, h_out, Tamb, **kwargs):
        """Initialize convective heat transfer model."""
        # Pass Tamb to parent for validation if detailed wall model is used
        kwargs['Tamb'] = Tamb
        super().__init__(**kwargs)
        self.h_out = h_out
        # Tamb already set by parent, but ensure it's set here too for clarity
        self.Tamb = Tamb

    def _calculate_external_heat_flux(
        self, i, wetted_area, T_vessel, T_vessel_wetted, T_outer_wall, T_outer_wall_wetted
    ):
        """
        Calculate external convective heat flux from ambient.

        Q = h_out * A * (T_amb - T_wall)

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        T_vessel : array
            Vessel wall temperature (unwetted) history [K]
        T_vessel_wetted : array
            Vessel wall temperature (wetted) history [K]
        T_outer_wall : array
            Outer wall temperature (unwetted) history [K]
        T_outer_wall_wetted : array
            Outer wall temperature (wetted) history [K]

        Returns
        -------
        tuple of float
            (Q_outer, Q_outer_wetted, q_outer, q_outer_wetted)
        """
        # Unwetted region - convection from ambient
        Q_outer = (
            (self.surf_area_inner - wetted_area)
            * self.surf_area_outer
            / self.surf_area_inner
            * self.h_out
            * (self.Tamb - T_outer_wall[i - 1])
        )
        q_outer = self.h_out * (self.Tamb - T_outer_wall[i - 1])

        # Wetted region - convection from ambient
        Q_outer_wetted = (
            wetted_area
            * self.surf_area_outer
            / self.surf_area_inner
            * self.h_out
            * (self.Tamb - T_outer_wall_wetted[i - 1])
        )
        q_outer_wetted = self.h_out * (self.Tamb - T_outer_wall_wetted[i - 1])

        if np.isnan(Q_outer_wetted):
            Q_outer_wetted = 0

        return Q_outer, Q_outer_wetted, q_outer, q_outer_wetted


class FireHeatTransfer(WallHeatTransfer):
    """
    Heat transfer model with fire boundary condition (Stefan-Boltzmann).

    Used for "s-b" (Stefan-Boltzmann) heat transfer type where external
    heat load is from fire scenarios using radiative + convective heat transfer.

    Fire types supported:
    - API 521 pool fire (60 kW/m²)
    - API 521 jet fire (100 kW/m²)
    - Scandpower pool fire (100 kW/m²)
    - Scandpower jet fire (250 kW/m²)
    - Scandpower peak fire (350 kW/m²)

    External heat flux: Q = fire.sb_fire(T_wall, fire_type) * A

    Parameters
    ----------
    fire_type : str
        Type of fire scenario (e.g., "api521_pool", "scandpower_jet")
    **kwargs
        Additional parameters passed to WallHeatTransfer base class

    Examples
    --------
    >>> model = FireHeatTransfer(
    ...     vessel_geometry=inner_vol,
    ...     vessel_orientation="vertical",
    ...     diameter=0.273,
    ...     length=1.524,
    ...     vessel_cp=500,
    ...     vessel_density=7800,
    ...     vol_solid=0.05,
    ...     surf_area_inner=1.5,
    ...     surf_area_outer=1.7,
    ...     tstep=0.1,
    ...     species="HEOS::Hydrogen",
    ...     flow_direction="discharge",
    ...     fire_type="api521_pool"
    ... )
    """

    def __init__(self, fire_type, **kwargs):
        """Initialize fire heat transfer model."""
        super().__init__(**kwargs)
        self.fire_type = fire_type

    def _calculate_external_heat_flux(
        self, i, wetted_area, T_vessel, T_vessel_wetted, T_outer_wall, T_outer_wall_wetted
    ):
        """
        Calculate external fire heat flux using Stefan-Boltzmann equation.

        Q = fire.sb_fire(T_wall, fire_type) * A

        The fire.sb_fire() function accounts for:
        - Incident radiative heat flux from fire
        - Re-radiation from hot vessel surface
        - Convective heat transfer

        Parameters
        ----------
        i : int
            Current time step index
        wetted_area : float
            Wetted surface area [m²]
        T_vessel : array
            Vessel wall temperature (unwetted) history [K]
        T_vessel_wetted : array
            Vessel wall temperature (wetted) history [K]
        T_outer_wall : array
            Outer wall temperature (unwetted) history [K]
        T_outer_wall_wetted : array
            Outer wall temperature (wetted) history [K]

        Returns
        -------
        tuple of float
            (Q_outer, Q_outer_wetted, q_outer, q_outer_wetted)
        """
        # Unwetted region - fire heat load
        Q_outer = (
            fire.sb_fire(T_vessel[i - 1], self.fire_type)
            * (self.surf_area_inner - wetted_area)
            * self.surf_area_outer
            / self.surf_area_inner
        )
        q_outer = fire.sb_fire(T_vessel[i - 1], self.fire_type)

        # Wetted region - fire heat load
        Q_outer_wetted = (
            fire.sb_fire(T_vessel_wetted[i - 1], self.fire_type)
            * wetted_area
            * self.surf_area_outer
            / self.surf_area_inner
        )
        q_outer_wetted = fire.sb_fire(T_vessel_wetted[i - 1], self.fire_type)

        if np.isnan(Q_outer_wetted):
            Q_outer_wetted = 0

        return Q_outer, Q_outer_wetted, q_outer, q_outer_wetted
