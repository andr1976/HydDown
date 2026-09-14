# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Main calculation engine for pressure vessel filling and discharge simulations.

This module contains the HydDown class, which is the core of the HydDown package.
It integrates mass and energy balances over time to simulate pressure vessel
depressurization (discharge) and pressurization (filling) with heat transfer effects.

The HydDown class:
- Reads and validates YAML input defining vessel geometry, initial conditions,
  calculation type, valve parameters, and heat transfer settings
- Initializes thermodynamic state using CoolProp for fluid properties
- Integrates mass and energy balances using explicit Euler time stepping
- Calculates heat transfer between fluid and vessel wall
- Handles various thermodynamic paths: isothermal, isenthalpic, isentropic,
  constant internal energy, and full energy balance
- Supports multiple valve types: orifice, control valve, relief valve, constant mass flow
- Models fire heat loads using Stefan-Boltzmann approach
- Tracks two-phase systems with separate gas/liquid temperatures
- Can model 1-D transient heat conduction through vessel walls (composite materials)
- Stores time-series results and provides plotting capabilities

Calculation types:
- isothermal: Constant temperature (very slow process with large heat reservoir)
- isenthalpic: Constant enthalpy (adiabatic, no work)
- isentropic: Constant entropy (adiabatic with PV work)
- specified_U: Constant internal energy
- energybalance: Full energy balance with heat transfer and work

Heat transfer modes (for energybalance):
- fixed_U: Fixed overall heat transfer coefficient
- fixed_Q: Fixed heat input rate
- specified_h: Specified internal/external heat transfer coefficients
- detailed: 1-D transient conduction through vessel wall
- fire: External fire heat load using Stefan-Boltzmann equation

Valve types:
- orifice: Compressible flow through orifice (Yellow Book equation)
- control_valve: Control valve with Cv characteristic
- relief_valve: API 520/521 relief valve sizing
- mdot: Constant mass flow rate

The integration scheme uses explicit Euler method with user-specified time step.
Results are stored in numpy arrays for time, pressure, temperature, mass flow, etc.

Typical usage:
    import yaml
    from hyddown import HydDown

    with open('input.yml') as f:
        input_data = yaml.load(f, Loader=yaml.FullLoader)

    hdown = HydDown(input_data)
    hdown.run()
    hdown.plot()
"""

import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.optimize import minimize
from scipy.optimize import root_scalar
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from hyddown import transport as tp
from hyddown import validator
from hyddown import fire
from hyddown import thermesh as tm
import fluids


class HydDown:
    """
    Main class to to hold problem definition, running problem, storing results etc.
    """

    def __init__(self, input):
        """
        Parameters
        ----------
        input : dict
            Dict holding problem definition
        """
        self.input = input
        self.verbose = 0
        self.isrun = False
        self.validate_input()
        self.read_input()
        self.initialize()
        # sb_heat_load = np.loadtxt(
        #    "C:\\Users\\AndersAndreasen\\Documents\\GitHub\\HydDown\\src\\hyddown\\examples\\LPG-heat_load.txt"
        # )
        # self.sb_heat_load = lambda x: np.interp(
        #    x, sb_heat_load[:, 0], sb_heat_load[:, 1]
        # )

    def validate_input(self):
        """
        Validating the provided problem definition dict

        Raises
        ------
        ValueError
            If missing input is detected.
        """
        valid = validator.validation(self.input)
        if valid is False:
            raise ValueError("Error in input file")

    def read_input(self):
        """
        Reading in input/ problem definition dict and assigning to classs
        attributes.
        """
        self.length = self.input["vessel"]["length"]
        self.diameter = self.input["vessel"]["diameter"]
        if "type" in self.input["vessel"]:
            self.vessel_type = self.input["vessel"]["type"]
        else:
            self.vessel_type = "Flat-end"

        if "orientation" in self.input["vessel"]:
            if self.input["vessel"]["orientation"] == "horizontal":
                horizontal = True
            else:
                horizontal = False
        else:
            horizontal = True
            # Orientation

        if self.vessel_type == "Flat-end":
            self.inner_vol = fluids.TANK(
                D=self.diameter, L=self.length, horizontal=horizontal
            )
        elif self.vessel_type == "ASME F&D":
            self.inner_vol = fluids.TANK(
                D=self.diameter,
                L=self.length,
                sideA="torispherical",
                sideB="torispherical",
                horizontal=horizontal,
            )
        elif self.vessel_type == "DIN":
            self.inner_vol = fluids.TANK(
                D=self.diameter,
                L=self.length,
                sideA="torispherical",
                sideB="torispherical",
                sideA_f=1,
                sideA_k=0.1,
                sideB_f=1,
                sideB_k=0.1,
                horizontal=horizontal,
            )
        elif self.vessel_type == "Hemispherical":
            self.inner_vol = fluids.TANK(
                D=self.diameter,
                L=self.length,
                sideA="spherical",
                sideB="spherical",
                sideA_a=0.5 * self.diameter,
                sideB_a=0.5 * self.diameter,
                horizontal=horizontal,
            )

        t_wall = self.input["vessel"].get("thickness", 0.0)
        if t_wall <= 0 or self.vessel_type == "Flat-end":
            # Flat ends: add_thickness (outer length = L + 2t) is the correct concentric shell.
            self.outer_vol = self.inner_vol.add_thickness(t_wall)
        else:
            # Head-type tanks (Hemispherical / ASME F&D / DIN): fluids.TANK.add_thickness
            # spuriously adds 2*thickness to the cylinder length, inserting a phantom
            # cylindrical band (for a sphere this over-states the wall volume by ~50%).
            # Rebuild the outer shell concentrically with the SAME cylinder length.
            Do = self.diameter + 2 * t_wall
            if self.vessel_type == "Hemispherical":
                self.outer_vol = fluids.TANK(
                    D=Do, L=self.length, sideA="spherical", sideB="spherical",
                    sideA_a=0.5 * Do, sideB_a=0.5 * Do, horizontal=horizontal,
                )
            elif self.vessel_type == "ASME F&D":
                self.outer_vol = fluids.TANK(
                    D=Do, L=self.length, sideA="torispherical", sideB="torispherical",
                    horizontal=horizontal,
                )
            else:  # DIN
                self.outer_vol = fluids.TANK(
                    D=Do, L=self.length, sideA="torispherical", sideB="torispherical",
                    sideA_f=1, sideA_k=0.1, sideB_f=1, sideB_k=0.1, horizontal=horizontal,
                )

        self.p0 = self.input["initial"]["pressure"]
        self.T0 = self.input["initial"]["temperature"]

        self.species = "HEOS::" + self.input["initial"]["fluid"]

        # Detects if a multi component fluid is specified using & for separation of components
        if "&" in self.input["initial"]["fluid"]:
            comp_frac_pair = [
                str.replace("[", " ").replace("]", "").split(" ")
                for str in self.input["initial"]["fluid"].split("&")
            ]
            comp = [pair[0] for pair in comp_frac_pair]
            compSRK = [pair[0] + "-SRK" for pair in comp_frac_pair]
            molefracs = np.asarray([float(pair[1]) for pair in comp_frac_pair])
            molefracs = molefracs / sum(molefracs)
            self.molefracs = molefracs
            sep = "&"
            self.comp = sep.join(comp)
            self.compSRK = sep.join(compSRK)
        # Normally single component fluid is specified
        else:
            self.comp = self.input["initial"]["fluid"]
            self.molefracs = [1.0]
            self.compSRK = self.input["initial"]["fluid"]

        self.tstep = self.input["calculation"]["time_step"]
        self.time_tot = self.input["calculation"]["end_time"]
        self.method = self.input["calculation"]["type"]

        # Check for non-equilibrium model flag
        if "non_equilibrium" in self.input["calculation"]:
            self.non_equilibrium = self.input["calculation"]["non_equilibrium"]
        else:
            self.non_equilibrium = False

        # Non-equilibrium model only works with single component and energybalance
        if self.non_equilibrium:
            if "&" in self.input["initial"]["fluid"]:
                raise ValueError("Non-equilibrium model only supports single component fluids")
            if self.method != "energybalance":
                raise ValueError("Non-equilibrium model requires calculation.type = 'energybalance'")

        # Reading valve specific data
        if (
            self.input["valve"]["type"] == "orifice"
            or self.input["valve"]["type"] == "psv"
            or self.input["valve"]["type"] == "hem_release"
        ):
            self.p_back = self.input["valve"]["back_pressure"]
            self.D_orifice = self.input["valve"]["diameter"]
            self.CD = self.input["valve"]["discharge_coef"]
            if self.input["valve"]["type"] == "psv":
                self.Pset = self.input["valve"]["set_pressure"]
                self.blowdown = self.input["valve"]["blowdown"]
                self.psv_state = "closed"
        elif self.input["valve"]["type"] == "relief":
            self.p_back = self.input["valve"]["back_pressure"]
            self.Pset = self.input["valve"]["set_pressure"]
        elif self.input["valve"]["type"] == "controlvalve":
            self.p_back = self.input["valve"]["back_pressure"]
            self.Cv = self.input["valve"]["Cv"]
            if "xT" in self.input["valve"]:
                self.xT = self.input["valve"]["xT"]
            if "Fp" in self.input["valve"]:
                self.Fp = self.input["valve"]["Fp"]
            if (
                "characteristic" in self.input["valve"]
                and "time_constant" in self.input["valve"]
            ):
                self.valve_characteristic = self.input["valve"]["characteristic"]
                self.valve_time_constant = self.input["valve"]["time_constant"]
            else:
                self.valve_characteristic = "linear"
                self.valve_time_constant = 0

        elif (
            self.input["valve"]["type"]
            == "mdot"
            # and self.input["valve"]["flow"] == "filling"
        ):
            self.p_back = self.input["valve"]["back_pressure"]
        elif self.input["valve"]["type"] == "none":
            # No throttling device: the outflow is driven by a top-level "release" block.
            self.p_back = self.input["valve"].get("back_pressure", 101325.0)

        # Reading release specific data (thermopack CO2 HEM release + dry-ice state).
        # All thermopack work is delegated to hyddown.co2_release.CO2ReleaseModel; only
        # the plain input values are read here (no thermopack import in the HydDown class).
        self.has_release = "release" in self.input
        if self.has_release:
            rel = self.input["release"]
            self.release_type = rel["type"]  # 'liquid' (liquid space) or 'gas' (vapour space)
            self.D_release = rel["diameter"]
            self.CD_release = rel["discharge_coef"]
            # Separate discharge coefficient for a GAS discharge - a gas-space release, or the
            # gas tail of a liquid release once the liquid is exhausted. Defaults to the single
            # discharge_coef so existing inputs are unchanged. This lets a liquid release use a
            # two-phase/flashing Cd for the liquid (with liquid_nonequilibrium) while its gas
            # tail keeps the reliable single-phase (HEM) gas Cd. See docs/techref discharge.
            self.CD_release_gas = rel.get("discharge_coef_gas", self.CD_release)
            self.release_back_pressure = rel.get("back_pressure", 101325.0)
            self.release_atm_pressure = rel.get(
                "atm_pressure", self.release_back_pressure
            )
            # Thermodynamic backend for the CO2 release model. "CoolProp" (default) is the
            # thermopack-free implementation (CoolProp + solid table); "tcPR" uses the
            # thermopack backend (requires a thermopack install).
            self.release_eos = rel.get("eos", "CoolProp")
            # Opt-in solid-in-vessel fallback: once the tank reaches the triple point,
            # continue with a thermopack three-phase / sublimation model instead of
            # freezing. solid_h_inner is the (simplified) internal HTC used there.
            self.solid_in_vessel = rel.get("solid_in_vessel", False)
            # solid_h_inner: a fixed wall->cold-phase HTC [W/m2 K], OR a boiling correlation
            # evaluated at the triple-point saturated liquid with the wall superheat -
            # "calc" -> Rohsenow (h_inside_wetted); "cooper" -> Cooper reduced-pressure
            # correlation (carries CO2's high reduced pressure natively; the SINTEF choice).
            # In correlation mode the numeric value below is the fallback used when the
            # correlation cannot be evaluated.
            _shi = rel.get("solid_h_inner", 20.0)
            self.solid_h_inner_mode = _shi.lower() if isinstance(_shi, str) else "fixed"
            self.solid_h_inner = 150.0 if isinstance(_shi, str) else _shi
            # Two-zone plateau HTCs: wall->gas keeps the gas warm; gas->liquid/solid
            # interphase is kept ~0 so heat does not melt the freezing dry ice.
            self.solid_h_gas_wall = rel.get("solid_h_gas_wall", 15.0)
            self.solid_h_gas_liquid = rel.get("solid_h_gas_liquid", 0.0)  # plateau (gas<->liquid)
            self.solid_h_gas_solid = rel.get("solid_h_gas_solid", 0.0)  # descent (gas<->dry ice)
            self.solid_gas_wall_frac = rel.get("solid_gas_wall_frac", 0.5)
            # Non-equilibrium factor for a LIQUID discharge (delayed/metastable flashing
            # through a short orifice): 0 = equilibrium HEM, 1 = frozen all-liquid.
            self.liquid_nonequilibrium = rel.get("liquid_nonequilibrium", 0.0)
            # Optionally fade that boost linearly with (P - P_triple) as the vessel
            # depressurises (metastable flashing is a low-pressure saturated effect that
            # vanishes near the triple point); N = liquid_nonequilibrium at the initial
            # pressure. Default off, so existing inputs keep a constant N.
            self.liquid_ne_pressure_scaled = rel.get("liquid_ne_pressure_scaled", False)

        # valve type
        # - constant_mass
        # - functional mass flow
        self.thickness = 0

        if "rupture" in self.input:
            self.rupture_material = self.input["rupture"]["material"]
            if "fire" in self.input["rupture"]:
                self.rupture_fire = self.input["rupture"]["fire"]
            else:
                self.rupture_fire = "api_jet"
            self.rupture_k_s = self.input["rupture"].get("k_s", 0.85)

        # Reading heat transfer related data/information
        if "heat_transfer" in self.input:
            self.heat_method = self.input["heat_transfer"]["type"]
            if self.heat_method == "specified_h" or self.heat_method == "specified_U":
                self.Tamb = self.input["heat_transfer"]["temp_ambient"]
            if self.heat_method == "specified_U":
                self.Ufix = self.input["heat_transfer"]["U_fix"]
            if self.heat_method == "specified_Q":
                self.Qfix = self.input["heat_transfer"]["Q_fix"]
            if self.heat_method == "specified_h":
                self.vessel_cp = self.input["vessel"]["heat_capacity"]
                self.vessel_density = self.input["vessel"]["density"]
                self.vessel_orientation = self.input["vessel"]["orientation"]
                self.thickness = self.input["vessel"]["thickness"]
                self.h_out = self.input["heat_transfer"]["h_outer"]
                self.h_in = self.input["heat_transfer"]["h_inner"]
                if self.input["valve"]["flow"] == "filling":
                    if "D_throat" in self.input["heat_transfer"]:
                        self.D_throat = self.input["heat_transfer"]["D_throat"]
                    else:
                        self.D_throat = self.input["vessel"]["diameter"]
            if self.heat_method == "s-b":
                self.fire_type = self.input["heat_transfer"]["fire"]
                self.h_in = "calc"
                self.vessel_cp = self.input["vessel"]["heat_capacity"]
                self.vessel_density = self.input["vessel"]["density"]
                self.vessel_orientation = self.input["vessel"]["orientation"]
                self.thickness = self.input["vessel"]["thickness"]
                if "scaling" in self.input["heat_transfer"]:
                    self.scaling = self.input["heat_transfer"]["scaling"]
                else:
                    self.scaling = 1.0
                if self.input["valve"]["flow"] == "filling":
                    raise ValueError("Filling and Fire heat load not implemented")
            if self.heat_method == "specified_q":
                self.vessel_cp = self.input["vessel"]["heat_capacity"]
                self.vessel_density = self.input["vessel"]["density"]
                self.vessel_orientation = self.input["vessel"]["orientation"]
                self.thickness = self.input["vessel"]["thickness"]
                # Handle q_outer: can be a number (fixed) or dict (time-dependent)
                q_outer_input = self.input["heat_transfer"]["q_outer"]
                if isinstance(q_outer_input, (int, float)):
                    # Fixed heat flux
                    self.q_outer_func = lambda t: q_outer_input
                elif isinstance(q_outer_input, dict):
                    # Time-dependent heat flux from dict
                    time_data = np.array(q_outer_input["time"])
                    heat_flux_data = np.array(q_outer_input["heat_flux"])
                    self.q_outer_func = lambda t: np.interp(t, time_data, heat_flux_data)
                else:
                    raise ValueError("q_outer must be a number or dict with time/heat_flux")
                # Handle h_inner
                if "h_inner" in self.input["heat_transfer"]:
                    self.h_in = self.input["heat_transfer"]["h_inner"]
                else:
                    self.h_in = "calc"
                if self.input["valve"]["flow"] == "filling":
                    if "D_throat" in self.input["heat_transfer"]:
                        self.D_throat = self.input["heat_transfer"]["D_throat"]
                    else:
                        self.D_throat = self.input["vessel"]["diameter"]

    def initialize(self):
        """
        Preparing for running problem by creating the fluid objects required
        instantiating arrays for storing time-dependent results, setting additional
        required class attributes.
        """
        self.vol = self.inner_vol.V_total
        self.vol_tot = self.outer_vol.V_total
        self.vol_solid = self.vol_tot - self.vol
        self.surf_area_outer = self.outer_vol.A
        self.surf_area_inner = self.inner_vol.A

        self.fluid = CP.AbstractState("HEOS", self.comp)
        if "&" in self.comp:
            self.fluid.specify_phase(CP.iphase_gas)
        self.fluid.set_mole_fractions(self.molefracs)

        self.transport_fluid = CP.AbstractState("HEOS", self.compSRK)
        self.transport_fluid.specify_phase(CP.iphase_gas)
        self.transport_fluid.set_mole_fractions(self.molefracs)

        self.transport_fluid_wet = CP.AbstractState("HEOS", self.compSRK)
        self.transport_fluid_wet.specify_phase(CP.iphase_liquid)
        self.transport_fluid_wet.set_mole_fractions(self.molefracs)

        self.vent_fluid = CP.AbstractState("HEOS", self.comp)
        self.vent_fluid.specify_phase(CP.iphase_gas)
        self.vent_fluid.set_mole_fractions(self.molefracs)

        if "liquid_level" in self.input["vessel"]:
            ll = self.input["vessel"]["liquid_level"]
            V_liquid = self.inner_vol.V_from_h(ll)
            self.fluid.update(CP.PQ_INPUTS, self.p0, 0)
            liq_rho = self.fluid.rhomass()
            m_liq = V_liquid * liq_rho

            self.fluid.update(CP.PQ_INPUTS, self.p0, 1)
            gas_rho = self.fluid.rhomass()
            V_vapour = self.inner_vol.V_total - V_liquid
            m_vap = V_vapour * gas_rho

            m_tot = m_liq + m_vap
            rho0 = m_tot / self.inner_vol.V_total
            self.Q0 = m_vap / m_tot
            self.fluid.update(CP.PQ_INPUTS, self.p0, self.Q0)
            self.T0 = self.fluid.T()
            self.liquid_level0 = ll
            self.vent_fluid.update(CP.PQ_INPUTS, self.p0, 1.0)

        else:
            self.fluid.update(CP.PT_INPUTS, self.p0, self.T0)
            self.liquid_level0 = 0.0
            self.Q0 = self.fluid.Q()
            self.vent_fluid.update(CP.PT_INPUTS, self.p0, self.T0)

        self.res_fluid = CP.AbstractState("HEOS", self.comp)
        self.res_fluid.set_mole_fractions(self.molefracs)
        if self.input["valve"]["flow"] == "filling":
            self.res_fluid.update(CP.PT_INPUTS, self.p_back, self.T0)

        # Non-equilibrium model: Create separate AbstractState objects for gas and liquid
        if self.non_equilibrium:
            # Create parent phase objects - these CAN enter two-phase region
            # Do NOT use specify_phase - we want to detect phase change via quality
            self.fluid_gas = CP.AbstractState("HEOS", self.comp)
            self.fluid_gas.set_mole_fractions(self.molefracs)

            self.fluid_liquid = CP.AbstractState("HEOS", self.comp)
            self.fluid_liquid.set_mole_fractions(self.molefracs)

            # Initialize at saturation conditions (equilibrium at t=0)
            if "liquid_level" in self.input["vessel"]:
                # Two-phase initial condition. The liquid zone is saturated liquid at p0.
                # The gas zone is saturated vapour by default, or - if
                # initial.gas_temperature is given - a warmer (superheated) gas at p0: a
                # lower gas density and mass, so at a matched total mass the fill (liquid
                # mass) takes up the difference. Pressure and mass balance still close.
                self.fluid_liquid.update(CP.PQ_INPUTS, self.p0, 0.0)
                self.T_liquid0 = self.fluid_liquid.T()  # Saturated liquid
                if "gas_temperature" in self.input["initial"]:
                    self.T_gas0 = self.input["initial"]["gas_temperature"]
                    self.fluid_gas.update(CP.PT_INPUTS, self.p0, self.T_gas0)
                else:
                    self.fluid_gas.update(CP.PQ_INPUTS, self.p0, 1.0)
                    self.T_gas0 = self.fluid_gas.T()  # Saturated vapour

                # Calculate initial masses
                ll = self.input["vessel"]["liquid_level"]
                V_liquid = self.inner_vol.V_from_h(ll)
                self.m_liquid0 = V_liquid * self.fluid_liquid.rhomass()
                V_vapour = self.inner_vol.V_total - V_liquid
                self.m_gas0 = V_vapour * self.fluid_gas.rhomass()
            else:
                # Single-phase gas initial condition
                self.fluid_gas.update(CP.PT_INPUTS, self.p0, self.T0)
                self.T_gas0 = self.T0
                self.m_gas0 = self.fluid_gas.rhomass() * self.vol
                self.m_liquid0 = 0.0
                self.T_liquid0 = self.T0  # Not used if no liquid

        # data storage
        data_len = int(self.time_tot / self.tstep)
        self.rho = np.zeros(data_len)
        self.T_fluid = np.zeros(data_len)
        self.T_vent = np.zeros(data_len)
        self.T_vessel = np.zeros(data_len)
        self.T_vessel_wetted = np.zeros(data_len)
        self.T_inner_wall = np.zeros(data_len)
        self.T_inner_wall_wetted = np.zeros(data_len)
        self.T_outer_wall = np.zeros(data_len)
        self.T_outer_wall_wetted = np.zeros(data_len)
        self.T_bonded_wall = np.zeros(data_len)
        self.T_bonded_wall_wetted = np.zeros(data_len)
        self.Q_outer = np.zeros(data_len)
        self.Q_inner = np.zeros(data_len)
        self.Q_outer_wetted = np.zeros(data_len)
        self.Q_inner_wetted = np.zeros(data_len)
        self.q_outer = np.zeros(data_len)
        self.q_inner = np.zeros(data_len)
        self.q_outer_wetted = np.zeros(data_len)
        self.q_inner_wetted = np.zeros(data_len)
        self.h_inside = np.zeros(data_len)
        self.h_inside_wetted = np.zeros(data_len)
        self.h_gas_liquid = np.zeros(data_len)  # Gas-liquid interfacial HTC for NEM
        # Initialize Biot arrays as NaN (will be calculated if thermal_conductivity_biot specified)
        self.Biot = np.full(
            data_len, np.nan
        )  # Biot number for lumped capacitance validation
        self.Biot_wetted = np.full(data_len, np.nan)  # Biot number for wetted region
        self.T_vent = np.zeros(data_len)
        self.H_mass = np.zeros(data_len)
        self.S_mass = np.zeros(data_len)
        self.U_mass = np.zeros(data_len)
        self.U_tot = np.zeros(data_len)
        self.U_res = np.zeros(data_len)
        self.P = np.zeros(data_len)
        self.mass_fluid = np.zeros(data_len)
        self.mass_rate = np.zeros(data_len)
        self.time_array = np.zeros(data_len)
        self.relief_area = np.zeros(data_len)
        self.temp_profile = []
        self.rho0 = self.fluid.rhomass()
        self.m0 = self.rho0 * self.vol
        self.MW = self.fluid.molar_mass()
        self.vapour_mole_fraction = np.zeros(data_len)
        self.vapour_mass_fraction = np.zeros(data_len)
        self.vapour_volume_fraction = np.zeros(data_len)
        self.liquid_level = np.zeros(data_len)

        # Non-equilibrium model data storage
        if self.non_equilibrium:
            self.T_gas = np.zeros(data_len)  # Gas phase temperature
            self.T_liquid = np.zeros(data_len)  # Liquid phase temperature
            self.m_gas = np.zeros(data_len)  # Gas phase mass
            self.m_liquid = np.zeros(data_len)  # Liquid phase mass
            self.U_gas = np.zeros(data_len)  # Gas phase internal energy (specific)
            self.U_liquid = np.zeros(data_len)  # Liquid phase internal energy (specific)
            self.rho_gas = np.zeros(data_len)  # Gas phase density
            self.rho_liquid = np.zeros(data_len)  # Liquid phase density
            self.mdot_phase_transfer = np.zeros(data_len)  # Mass transfer rate (condensation > 0, evaporation < 0)

        # Release model (thermopack CO2 HEM release rate + dry-ice atmospheric state).
        # Imported lazily so non-release runs never import thermopack.
        if self.has_release:
            if str(self.release_eos).lower().replace("-", "") in ("coolprop", "cp"):
                # thermopack-free backend (CoolProp + solid table)
                from hyddown.co2_release_cp import CO2ReleaseModelCP as _ReleaseModel
            else:
                # thermopack backend (tcPR / GERG2008 / MEOS)
                from hyddown.co2_release import CO2ReleaseModel as _ReleaseModel

            self.release_model = _ReleaseModel(
                back_pressure=self.release_back_pressure,
                atm_pressure=self.release_atm_pressure,
                eos=self.release_eos,
                liquid_nonequilibrium=self.liquid_nonequilibrium,
                liquid_ne_pressure_scaled=self.liquid_ne_pressure_scaled,
                liquid_ne_pref=self.p0,
            )
            # A 'liquid' release switches to 'gas' once the liquid inventory is exhausted.
            self.release_phase = self.release_type
            # Set once the tank approaches the CO2 triple point: the release is frozen so
            # the vessel state never crosses into the solid-in-vessel regime (out of scope,
            # and where the CoolProp vessel solver fails). See DRY_ICE_HANDOVER.md sec. 5.5.
            self.release_frozen = False
            # Atmospheric (1 atm) end-state time series
            self.T_atm = np.zeros(data_len)  # atmospheric temperature [K]
            self.x_vap_atm = np.zeros(data_len)  # atmospheric vapour mass fraction [-]
            self.x_solid_atm = np.zeros(data_len)  # atmospheric dry-ice (solid) mass fraction [-]
            self.solid_frac_throat = np.zeros(data_len)  # dry-ice fraction at choked throat [-]
            self.m_dryice_cum = np.zeros(data_len)  # cumulative dry-ice mass released [kg]
            self.release_choked = np.zeros(data_len)  # 1.0 if choked flow, else 0.0
            self.release_rate = np.zeros(data_len)  # release-hole mass flow only [kg/s]
            # Solid-in-vessel fallback (below the triple point)
            self.m_solid = np.zeros(data_len)  # in-vessel dry-ice (solid CO2) mass [kg]
            self.solid_regime = False  # True once handed off to the solid-in-vessel model
            self.M_vessel = 0.0  # total vessel mass [kg] (thermopack-basis state)
            self.U_vessel = 0.0  # total vessel internal energy [J] (thermopack basis)
            # Two-zone plateau/descent state (gas zone + liquid/solid zone)
            self.tz_plateau = False
            self.tz_descent = False
            self.tz_m_gas = 0.0
            self.tz_U_gas = 0.0
            self.tz_M_ls = 0.0
            self.tz_U_ls = 0.0
            self.tz_m_solid = 0.0
            # Below-triple wetted (liquid/solid-contact) wall node. Seeded at handoff from
            # the cold above-triple wetted wall and then relaxed toward the liquid/solid
            # temperature, so the wetted wall keeps tracking the cold phase instead of
            # jumping to the (warm) gas-contact wall. None until the handoff.
            self.tz_T_wall_wet = None

    def calc_liquid_level(self):
        """
        Calculate liquid level height based on current two-phase fluid state.

        For two-phase systems (0 ≤ quality ≤ 1), calculates the height of liquid
        phase in the vessel based on vapor quality, phase densities, and vessel geometry.
        Uses vessel geometry from fluids.TANK to convert liquid volume to height.

        Parameters
        ----------
        fluid : CoolProp AbstractState
            Current fluid state

        Returns
        -------
        float
            Liquid level height from vessel bottom [m].
            Returns 0.0 for single-phase gas (quality > 1 or subcooled liquid).
        """
        if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
            rho_liq = self.fluid.saturated_liquid_keyed_output(CP.iDmass)
            rho_vap = self.fluid.saturated_vapor_keyed_output(CP.iDmass)
            m_liq = self.fluid.rhomass() * self.inner_vol.V_total * (1 - self.fluid.Q())
            V_liq = m_liq / rho_liq
            h_liq = self.inner_vol.h_from_V(V_liq)
            return h_liq
        else:
            return 0.0

    def compute_release(self, P, i):
        """
        Compute the thermopack HEM release mass rate at tank pressure ``P`` and store the
        atmospheric (1 atm) dry-ice end state at time index ``i``.

        The stagnation branch (saturated liquid / vapour) follows ``self.release_phase``,
        which is set from ``release.type`` and flips ``liquid`` -> ``gas`` once the liquid
        inventory is exhausted (handled in the mass balance).

        Parameters
        ----------
        P : float
            Current tank pressure [Pa].
        i : int
            Time index at which to store the atmospheric state.

        Returns
        -------
        float
            Release mass flow [kg/s] (>= 0).
        """
        # Validity floor: once the tank approaches/drops to the triple point, dry ice starts
        # forming inside the vessel (solid-in-vessel regime, out of scope) - stop the
        # release rather than extrapolate. Also stop if there is no driving pressure.
        if (
            self.release_frozen
            or P <= self.release_model.P_TRIPLE
            or P <= self.release_back_pressure
        ):
            return 0.0

        rate, atm = self.release_model.release_state(
            P,
            self.release_phase,
            self.CD_release_gas if self.release_phase == "gas" else self.CD_release,
            self.D_release ** 2 / 4 * math.pi,
        )
        self.T_atm[i] = atm["T"]
        self.x_vap_atm[i] = atm["vapour_frac"]
        self.x_solid_atm[i] = atm["solid_frac"]
        self.solid_frac_throat[i] = rate["solid_frac_throat"]
        self.release_choked[i] = 1.0 if rate["choked"] else 0.0
        self.release_rate[i] = rate["mdot"]
        return rate["mdot"]

    def _h_gas_wall(self, T_gas, T_wall, P):
        """Gas->wall heat-transfer coefficient below the triple point [W/m2K].

        ``release.solid_h_gas_wall`` is either a fixed number (default 15) or the string
        ``"calc"`` for a natural-convection estimate. The correlation (transport.h_inner:
        Pr, Gr -> Ra -> Nu, h = Nu*k/L) uses CoolProp's ``T|gas`` phase spec, which is valid
        for CO2 vapour below the triple point (verified down to ~1 bar). Falls back to 15 on
        any CoolProp failure or a negligible gas/wall temperature difference.
        """
        hgw = self.solid_h_gas_wall
        if isinstance(hgw, str) and hgw.lower() == "calc":
            if abs(T_wall - T_gas) < 0.1:
                return 15.0
            try:
                return tp.h_inner(self.inner_vol.D, T_gas, T_wall, P, "HEOS::CO2")
            except Exception:
                return 15.0
        return float(hgw)

    def _wetted_wall_htc(self, Tw, T_cold):
        """Below-triple wetted-wall HTC [W/m2 K]: the fixed ``solid_h_inner``, or a boiling
        correlation selected by ``solid_h_inner``:

          * ``"calc"``   -> Rohsenow (``h_inside_wetted``);
          * ``"cooper"`` -> Cooper reduced-pressure correlation, which carries CO2's high
                            reduced pressure natively (the correlation the SINTEF reference
                            model uses; needs only P_r, M, q and roughness).

        No boiling liquid exists below the triple point, so the correlation is evaluated at the
        triple-point saturated liquid and only the driving superheat ``Tw - T_cold`` varies with
        the descent. Both are capped at 3000 W/m2 K (as ``h_inside_wetted`` is): above that the
        wall is conduction-limited and its temperature is insensitive to the exact coefficient,
        and the cap keeps the explicit wall step stable. Falls back to the fixed value when the
        correlation cannot be evaluated (or when there is no superheat to drive boiling).
        """
        mode = getattr(self, "solid_h_inner_mode", "fixed")
        if mode == "fixed":
            return self.solid_h_inner
        Te = Tw - T_cold
        if Te <= 0.0:
            return self.solid_h_inner
        Ptp = self.release_model.P_TRIPLE_EOS
        try:
            if mode == "cooper":
                import ht
                from CoolProp.CoolProp import PropsSI
                if not hasattr(self, "_co2_Pc"):
                    self._co2_Pc = PropsSI("Pcrit", "CO2")
                    self._co2_MW = PropsSI("molar_mass", "CO2") * 1000.0  # g/mol
                h = ht.Cooper(P=Ptp, Pc=self._co2_Pc, MW=self._co2_MW, Te=Te, Rp=1e-6)
            else:  # "calc" -> Rohsenow
                self.fluid_liquid.update(CP.PQ_INPUTS, Ptp, 0.0)       # sat liquid at triple
                self.transport_fluid_wet.update(CP.PQ_INPUTS, Ptp, 0.0)
                L = self.diameter if self.vessel_orientation == "horizontal" else self.length
                h = tp.h_inside_wetted(L, Tw, T_cold, self.transport_fluid_wet, self.fluid_liquid)
            if h and h > 0:
                return min(h, 3000.0)
        except Exception:
            pass
        return self.solid_h_inner

    def _wetted_wall_step(self, i, T_cold, has_cold):
        """Evolve the below-triple wetted (liquid/solid-contact) wall node and store it in
        ``T_vessel_wetted[i]``.

        The wetted wall carries its own energy balance: it exchanges with the cold phase
        (liquid at the triple point on the plateau, dry ice on the sublimation line during
        the descent) through ``release.solid_h_inner`` over the wetted (1 - gas) fraction of
        the inner area, plus ambient over the same fraction of the outer area. It uses the
        wetted fraction of the wall mass. This keeps the reported wetted wall tracking the
        cold phase down toward the measured bottom-of-vessel wall temperatures, instead of
        collapsing to the (warm) gas-contact wall as it did before.

        The reciprocal wall->solid heat is intentionally NOT removed from the liquid/solid
        zone: that zone is kept adiabatic to the wall, matching the existing two-zone
        simplification and leaving the calibrated retained-dry-ice mass unchanged.
        """
        dt = self.tstep
        m_wall = getattr(self, "vessel_density", 0.0) * self.vol_solid
        cp_wall = getattr(self, "vessel_cp", 500.0)
        h_out = getattr(self, "h_out", 0.0)
        frac_wet = max(1.0 - self.solid_gas_wall_frac, 0.0)
        if self.tz_T_wall_wet is None:
            self.tz_T_wall_wet = self.T_vessel_wetted[i - 1]
        Tw = self.tz_T_wall_wet
        m_ww = m_wall * frac_wet
        if has_cold and m_ww > 0:
            Tamb = getattr(self, "Tamb", Tw)
            A_wet = self.surf_area_inner * frac_wet
            Q_ws = self._wetted_wall_htc(Tw, T_cold) * A_wet * (Tw - T_cold)  # wall -> cold phase
            Q_out = h_out * self.surf_area_outer * frac_wet * (Tamb - Tw)
            Tw = Tw + dt * (Q_out - Q_ws) / (m_ww * cp_wall)
        self.tz_T_wall_wet = Tw
        self.T_vessel_wetted[i] = Tw

    def _two_zone_plateau_step(self, i):
        """Two-zone triple-point plateau step: warm gas zone + adiabatic liquid/solid lever.

        The gas keeps its own (superheated) temperature - heated by the wall, decoupled
        from the cold liquid/solid by a near-zero interphase HTC - so the vessel reproduces
        the measured gas superheat while the liquid freezes to dry ice at the pinned triple
        point. Transitions to the sublimation descent once the liquid is exhausted.
        """
        rm = self.release_model
        dt = self.tstep
        V = self.vol
        area = self.D_release ** 2 / 4 * math.pi

        T_g_prev = rm.gas_T_from_u(self.tz_U_gas / self.tz_m_gas)
        Twall_prev = self.T_vessel[i - 1]
        A_g = self.surf_area_inner * self.solid_gas_wall_frac
        Q_wg = self._h_gas_wall(T_g_prev, Twall_prev, self.P[i - 1]) * A_g * (Twall_prev - T_g_prev)  # wall -> gas
        Q_gl = self.solid_h_gas_liquid * A_g * (T_g_prev - rm.T_TRIPLE_EOS)  # gas -> L/S

        r = rm.two_zone_plateau_step(
            self.tz_m_gas, self.tz_U_gas, self.tz_M_ls, self.tz_U_ls,
            Q_wg, Q_gl, dt, self.CD_release_gas, area, V,
        )
        self.tz_m_gas, self.tz_U_gas = r["m_g"], r["U_g"]
        self.tz_M_ls, self.tz_U_ls = r["M_ls"], r["U_ls"]

        # lumped wall: exchanges with the gas and ambient (liquid/solid zone is adiabatic)
        m_wall = getattr(self, "vessel_density", 0.0) * self.vol_solid
        cp_wall = getattr(self, "vessel_cp", 500.0)
        h_out = getattr(self, "h_out", 0.0)
        Tamb = getattr(self, "Tamb", T_g_prev)
        Q_out = h_out * self.surf_area_outer * (Tamb - Twall_prev)
        self.T_vessel[i] = (
            Twall_prev + dt * (Q_out - Q_wg) / (m_wall * cp_wall) if m_wall > 0 else Twall_prev
        )

        m_s = max(r["m_s"], 0.0)
        m_l = max(r["m_l"], 0.0)
        m_g = max(r["m_g"], 0.0)
        atm = rm.atm_split(rm.gas_h_at(r["T_g"]))  # released vapour flashes to 1 atm

        self.P[i] = r["P"]
        self.T_gas[i] = r["T_g"]
        self.T_liquid[i] = rm.T_TRIPLE_EOS
        self.T_fluid[i] = r["T_g"]
        # wetted wall tracks the cold liquid/solid at the triple point (not the gas wall)
        self._wetted_wall_step(i, rm.T_TRIPLE_EOS, True)
        self.m_solid[i] = m_s
        self.m_liquid[i] = m_l
        self.m_gas[i] = m_g
        self.mass_fluid[i] = m_g + self.tz_M_ls
        self.rho[i] = self.mass_fluid[i] / V
        self.mass_rate[i] = r["mdot"]
        self.release_rate[i] = r["mdot"]
        self.solid_frac_throat[i] = 0.0
        self.release_choked[i] = 1.0
        self.liquid_level[i] = self.inner_vol.h_from_V(m_l * rm.v_l) if m_l > 1e-6 else 0.0
        self.T_atm[i] = atm["T"]
        self.x_vap_atm[i] = atm["vapour_frac"]
        self.x_solid_atm[i] = atm["solid_frac"]

        # transition to the two-zone sublimation descent once the liquid is exhausted
        if m_l <= 1e-3:
            self.tz_plateau = False
            self.tz_descent = True
            self.release_phase = "gas"
            self.tz_m_solid = m_s  # dry ice carried into the descent (warm gas kept)

    def _two_zone_descent_step(self, i):
        """Two-zone sublimation descent: warm gas leaks/depressurises; the near-adiabatic
        dry ice sublimes only enough to cool itself down the sublimation line, so most of
        it is retained."""
        rm = self.release_model
        dt = self.tstep
        V = self.vol
        area = self.D_release ** 2 / 4 * math.pi

        T_g_prev = rm._gas_T_from_u_P(self.tz_U_gas / self.tz_m_gas, self.P[i - 1])
        Twall_prev = self.T_vessel[i - 1]
        A_g = self.surf_area_inner * self.solid_gas_wall_frac
        Q_wg = self._h_gas_wall(T_g_prev, Twall_prev, self.P[i - 1]) * A_g * (Twall_prev - T_g_prev)
        UA_gs = self.solid_h_gas_solid * A_g  # gas->dry-ice interphase conductance [W/K]

        if self.P[i - 1] <= self.release_back_pressure * 1.002:
            r = {"m_g": self.tz_m_gas, "U_g": self.tz_U_gas, "m_solid": self.tz_m_solid,
                 "T_g": T_g_prev, "T_s": rm._T_sub_of_P(self.P[i - 1]),
                 "P": self.P[i - 1], "mdot": 0.0}
        else:
            r = rm.two_zone_descent_step(self.tz_m_gas, self.tz_U_gas, self.tz_m_solid,
                                         self.P[i - 1], Q_wg, UA_gs, dt, self.CD_release_gas, area, V)
        self.tz_m_gas, self.tz_U_gas, self.tz_m_solid = r["m_g"], r["U_g"], r["m_solid"]

        m_wall = getattr(self, "vessel_density", 0.0) * self.vol_solid
        cp_wall = getattr(self, "vessel_cp", 500.0)
        h_out = getattr(self, "h_out", 0.0)
        Tamb = getattr(self, "Tamb", T_g_prev)
        Q_out = h_out * self.surf_area_outer * (Tamb - Twall_prev)
        self.T_vessel[i] = (
            Twall_prev + dt * (Q_out - Q_wg) / (m_wall * cp_wall) if m_wall > 0 else Twall_prev
        )

        atm = rm.atm_split(rm._gas2d(rm._g2_h, r["T_g"], r["P"]))
        self.P[i] = r["P"]
        self.T_gas[i] = r["T_g"]
        # T_s is pinned to the sublimation line, so it only means something while dry ice is
        # actually present. For a liquid drain that empties to residual gas before the triple
        # point (m_solid -> 0), reporting T_s would be a phantom cold "solid" for a zero-mass
        # zone; report the real residual (gas) temperature instead.
        self.T_liquid[i] = r["T_s"] if r["m_solid"] > 1e-2 else r["T_g"]
        self.T_fluid[i] = r["T_g"]
        # wetted wall tracks the dry ice down the sublimation line while solid is present;
        # once the solid is gone (liquid drain -> residual gas) the wetted wall is meaningless,
        # so relax it toward the gas-contact wall instead.
        if r["m_solid"] > 1e-2:
            self._wetted_wall_step(i, r["T_s"], True)
        else:
            self.T_vessel_wetted[i] = self.T_vessel[i]
            self.tz_T_wall_wet = self.T_vessel[i]
        self.m_solid[i] = r["m_solid"]
        self.m_liquid[i] = 0.0
        self.m_gas[i] = r["m_g"]
        self.mass_fluid[i] = r["m_g"] + r["m_solid"]
        self.rho[i] = self.mass_fluid[i] / V
        self.mass_rate[i] = r["mdot"]
        self.release_rate[i] = r["mdot"]
        self.solid_frac_throat[i] = 0.0
        self.release_choked[i] = 1.0 if r["mdot"] > 0 else 0.0
        self.liquid_level[i] = 0.0
        self.T_atm[i] = atm["T"]
        self.x_vap_atm[i] = atm["vapour_frac"]
        self.x_solid_atm[i] = atm["solid_frac"]

    def _solid_regime_step(self, i):
        """One timestep of the opt-in solid-in-vessel model (below the triple point).

        Integrates the vessel total mass and internal energy (thermopack basis) under the
        continuing leak and a simplified lumped-wall heat transfer, then resolves the
        state analytically: a three-phase invariant point (T, P pinned at the triple
        point, dry ice accumulating) or a solid+gas point riding the sublimation line
        once the liquid is exhausted. All thermodynamics come from the release model.

        Two regimes:
        * while liquid remains - a two-zone plateau: a warm (superheated) gas zone plus
          an adiabatic liquid/solid zone that does the freezing lever, coupled so the gas
          fills the vapour volume at the triple-point pressure. This reproduces the
          measured gas superheat and the liquid-freezes-to-dry-ice plateau.
        * once the liquid is exhausted - the single-zone sublimation descent (solid+gas
          riding the sublimation line down to the back pressure).

        Simplifications: lumped wall; the liquid/solid zone is adiabatic to the wall; the
        gas-wall and interphase HTCs are ``release.solid_h_gas_wall`` /
        ``release.solid_h_gas_liquid``; phase equilibrium is instantaneous.
        """
        rm = self.release_model
        dt = self.tstep
        V = self.vol
        area = self.D_release ** 2 / 4 * math.pi

        if self.tz_plateau:
            self._two_zone_plateau_step(i)
            return
        if self.tz_descent:
            self._two_zone_descent_step(i)
            return

        T_prev = self.T_fluid[i - 1]
        P_prev = self.P[i - 1]

        # --- leak from the previous state (explicit Euler) ---
        if self.m_liquid[i - 1] <= 1e-3 and self.release_phase == "liquid":
            self.release_phase = "gas"  # liquid exhausted -> vapour leak (regime C)
        if P_prev <= self.release_back_pressure * 1.002:
            rate = {"mdot": 0.0, "h0": rm.h_g, "solid_frac_throat": 0.0, "choked": False}
            atm = {"T": 0.0, "vapour_frac": 0.0, "solid_frac": 0.0}
        elif self.release_phase == "liquid" and self.m_liquid[i - 1] > 1e-3:
            # regime B: saturated-liquid leak at the triple point
            rate = rm.hem_rate(rm.P_TRIPLE_EOS, "liquid", self.CD_release, area)
            atm = rm.atm_split(rate["h0"])
        else:
            # vapour leak from the current vessel state (triple point or sublimation line)
            rate = rm.gas_leak_rate(T_prev, P_prev, self.CD_release_gas, area)
            atm = rm.atm_split(rate["h0"])
        mdot = rate["mdot"]
        h_leak = rate["h0"]

        # --- simplified lumped-wall heat transfer ---
        m_wall = getattr(self, "vessel_density", 0.0) * self.vol_solid
        cp_wall = getattr(self, "vessel_cp", 500.0)
        h_out = getattr(self, "h_out", 0.0)
        Tamb = getattr(self, "Tamb", T_prev)
        Twall_prev = self.T_vessel[i - 1]
        Q_out = h_out * self.surf_area_outer * (Tamb - Twall_prev)  # W into the wall
        Q_in = self.solid_h_inner * self.surf_area_inner * (Twall_prev - T_prev)  # W into fluid
        if m_wall > 0:
            self.T_vessel[i] = Twall_prev + dt * (Q_out - Q_in) / (m_wall * cp_wall)
        else:
            self.T_vessel[i] = Twall_prev

        # --- overall mass & energy balance (thermopack basis) ---
        self.M_vessel = self.M_vessel - mdot * dt
        self.U_vessel = self.U_vessel + dt * (Q_in - mdot * h_leak)

        # --- resolve the phase state ---
        st = rm.vessel_state_below_triple(self.M_vessel, self.U_vessel, V)
        if st["regime"] == "above":
            # Warmed back onto the L+G saturation at the triple point: clamp there with
            # the volume-consistent split (mass-conserving).
            m_l, m_g, self.U_vessel = rm.triple_LG_from_MV(self.M_vessel, V)
            st = {"m_s": 0.0, "m_l": m_l, "m_g": m_g,
                  "T": rm.T_TRIPLE_EOS, "P": rm.P_TRIPLE_EOS}

        m_s = max(st["m_s"], 0.0)
        m_l = max(st["m_l"], 0.0)
        m_g = max(st["m_g"], 0.0)

        # --- store results ---
        self.P[i] = st["P"]
        self.T_fluid[i] = st["T"]
        self.T_gas[i] = st["T"]
        self.T_liquid[i] = st["T"]
        self.T_vessel_wetted[i] = self.T_vessel[i]
        self.m_solid[i] = m_s
        self.m_liquid[i] = m_l
        self.m_gas[i] = m_g
        self.mass_fluid[i] = self.M_vessel
        self.rho[i] = self.M_vessel / V
        self.mass_rate[i] = mdot
        self.release_rate[i] = mdot
        self.solid_frac_throat[i] = rate.get("solid_frac_throat", 0.0)
        self.release_choked[i] = 1.0 if rate.get("choked", False) else 0.0
        self.liquid_level[i] = (
            self.inner_vol.h_from_V(m_l * rm.v_l) if m_l > 1e-6 else 0.0
        )
        # atmospheric (released-stream) dry ice
        self.T_atm[i] = atm["T"]
        self.x_vap_atm[i] = atm["vapour_frac"]
        self.x_solid_atm[i] = atm["solid_frac"]

    def PHres(self, T, P, H):
        """
        Residual enthalpy function to be minimized during PH-problem.

        Used by numerical optimizer (scipy.optimize.minimize) to find temperature
        that satisfies constant pressure-enthalpy constraints for multicomponent
        fluids. The optimizer adjusts T until residual approaches zero.

        Parameters
        ----------
        H : float
            Enthalpy at initial/final conditions [J/kg]
        P : float
            Pressure at final conditions [Pa]
        T : float
            Updated estimate for the final temperature at (P,H) [K]

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
        H : float
            Enthalpy at initial/final conditions [J/kg]
        P : float
            Pressure at final conditions [Pa]
        T : float
            Updated estimate for the final temperature at (P,H) [K]

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
        Defining a constant pressure, constant enthalpy problem i.e. typical adiabatic
        problem like e.g. valve flow for the vented flow (during discharge).
        For multicomponent mixture the final temperature is changed/optimised until the residual
        enthalpy is near zero in an optimisation step. For single component fluid the coolprop
        built in methods are used for speed.

        Parameters
        ----------
        H : float
            Enthalpy at initial/final conditions
        P : float
            Pressure at final conditions.
        Tguess : float
            Initial guess for the final temperature at P,H
        """

        # Multicomponent case
        if "&" in self.species:
            x0 = Tguess
            if relief == False:
                res = minimize(
                    self.PHres,
                    x0,
                    args=(P, H),
                    method="Nelder-Mead",
                    options={"xatol": 0.1, "fatol": 0.001},
                )
                T1 = res.x[0]
            else:
                res = root_scalar(
                    self.PHres_relief,
                    args=(P, H),
                    x0=x0,
                    method="newton",
                )
                T1 = res.root
        # single component fluid case
        else:
            T1 = PropsSI("T", "P", P, "H", H, self.species)
        return T1

    def UDres(self, x, U, rho):
        """
        Residual U-rho to be minimised during a U-rho/UV-problem

        Parameters
        ----------
        U : float
            Internal energy at final conditions
        rho : float
            Density at final conditions
        """
        self.fluid.update(CP.PT_INPUTS, x[0], x[1])
        return ((U - self.fluid.umass()) / U) ** 2 + (
            (rho - self.fluid.rhomass()) / rho
        ) ** 2

    def UDproblem(self, U, rho, Pguess, Tguess):
        """
        Defining a constant UV problem i.e. constant internal energy and density/volume
        problem relevant for the 1. law of thermodynamics.
        For multicomponent mixture the final temperature/pressure is changed/optimised until the
        residual U/rho is near zero. For single component fluid the coolprop
        built in methods are used for speed.

        Parameters
        ----------
        U : float
            Internal energy at final conditions
        rho : float
            Density at final conditions.
        Pguess : float
            Initial guess for the final pressure at U, rho
        Tguess : float
            Initial guess for the final temperature at U, rho
        """
        if "&" in self.species:
            x0 = [Pguess, Tguess]
            res = minimize(
                self.UDres,
                x0,
                args=(U, rho),
                method="Nelder-Mead",
                options={"xatol": 0.1, "fatol": 0.001},
            )
            P1 = res.x[0]
            T1 = res.x[1]
            Ures = U - self.fluid.umass()
        else:
            P1 = PropsSI("P", "D", rho, "U", U, self.species)
            T1 = PropsSI("T", "D", rho, "U", U, self.species)
            Ures = 0
        return P1, T1, Ures

    def nem_objective_function(self, x, m_gas, m_liquid, U_gas_tot, U_liquid_tot, T_gas_prev, T_liquid_prev):
        """
        Objective function for NEM solver (inspired by ORS-openthermo).

        Solves for (T_gas, T_liquid, P) simultaneously by minimizing
        normalized residuals for volume and energy constraints.

        Parameters
        ----------
        x : array
            [T_gas, T_liquid, P]
        m_gas : float
            Gas mass [kg]
        m_liquid : float
            Liquid mass [kg]
        U_gas_tot : float
            Total gas internal energy [J]
        U_liquid_tot : float
            Total liquid internal energy [J]
        T_gas_prev : float
            Previous gas temperature (for bounds checking) [K]
        T_liquid_prev : float
            Previous liquid temperature (for bounds checking) [K]

        Returns
        -------
        float
            Objective value (sum of normalized squared residuals)
        """
        T_gas, T_liquid, P = x

        # Sanity checks
        if P <= 0 or T_gas <= 0 or T_liquid <= 0:
            return 1e10

        # Sanity checks on inputs - return high penalty if out of reasonable range
        if P <= 1e4 or P > 1e8:  # 0.1 to 1000 bar
            return 1e10
        if T_gas <= 100 or T_gas > 1000:  # 100 K to 1000 K
            return 1e10
        if T_liquid <= 100 or T_liquid > 1000:
            return 1e10

        # Update gas state with T, P
        # Handle case where PT is at saturation (CoolProp will fail)
        if m_gas > 1e-6:
            try:
                self.fluid_gas.update(CP.PT_INPUTS, P, T_gas)
                rho_gas = self.fluid_gas.rhomass()
                U_gas_calc = self.fluid_gas.umass() * m_gas
                V_gas = m_gas / rho_gas
            except ValueError as e:
                # PT failed - likely at or near saturation
                # Penalize this state heavily to push optimizer away from saturation line
                return 1e10
        else:
            U_gas_calc = 0.0
            V_gas = 0.0

        # Update liquid state with T, P
        if m_liquid > 1e-6:
            try:
                self.fluid_liquid.update(CP.PT_INPUTS, P, T_liquid)
                rho_liquid = self.fluid_liquid.rhomass()
                U_liquid_calc = self.fluid_liquid.umass() * m_liquid
                V_liquid = m_liquid / rho_liquid
            except ValueError as e:
                # PT failed - likely at or near saturation
                # Penalize this state heavily
                return 1e10
        else:
            U_liquid_calc = 0.0
            V_liquid = 0.0

        # Calculate normalized squared residuals
        V_total = V_gas + V_liquid

        # Volume constraint (normalized by vessel volume)
        vol_res = ((V_total - self.vol) / self.vol) ** 2

        # Energy constraints (normalized by target energy)
        if m_gas > 1e-6 and U_gas_tot > 0:
            energy_res_gas = ((U_gas_calc - U_gas_tot) / U_gas_tot) ** 2
        else:
            energy_res_gas = 0.0

        if m_liquid > 1e-6 and U_liquid_tot > 0:
            energy_res_liq = ((U_liquid_calc - U_liquid_tot) / U_liquid_tot) ** 2
        else:
            energy_res_liq = 0.0

        # Combined objective with equal weights
        # All residuals are normalized squared differences
        objective = vol_res + energy_res_gas + energy_res_liq

        # Debug: Store individual residuals for diagnosis (only if requested)
        if hasattr(self, '_debug_nem_residuals') and self._debug_nem_residuals:
            print(f"  P={P/1e5:.2f} bar, T_g={T_gas:.1f} K, T_l={T_liquid:.1f} K")
            print(f"    vol_res={vol_res:.2e}, U_gas_res={energy_res_gas:.2e}, U_liq_res={energy_res_liq:.2e}")
            print(f"    V_calc={V_total:.4f} m³ vs V_target={self.vol:.4f} m³")
            print(f"    U_gas_calc={U_gas_calc/1e6:.3f} MJ vs U_gas_target={U_gas_tot/1e6:.3f} MJ")
            print(f"    U_liq_calc={U_liquid_calc/1e6:.3f} MJ vs U_liq_target={U_liquid_tot/1e6:.3f} MJ")

        return objective

    def nem_residual(self, P, m_gas, m_liquid, U_gas_tot, U_liquid_tot):
        """
        Residual function for non-equilibrium model (NEM) solver.

        Solves for common pressure P such that total volume equals vessel volume.

        Key insight: Both phases are in mechanical equilibrium (same pressure P).
        Given P and U for each phase, CoolProp gives us rho and T via DmassUmass inputs.

        The constraint is: V_gas(P,U_gas) + V_liquid(P,U_liquid) = V_vessel

        Parameters
        ----------
        P : float
            Common pressure for both phases [Pa]
        m_gas : float
            Gas phase mass [kg]
        m_liquid : float
            Liquid phase mass [kg]
        U_gas_tot : float
            Total gas phase internal energy [J]
        U_liquid_tot : float
            Total liquid phase internal energy [J]

        Returns
        -------
        float
            Volume residual (V_total - V_vessel) / V_vessel
        """
        # Ensure positive pressure
        if P <= 0:
            return 1e10

        V_total = 0.0

        # Gas phase: given P and U, find rho (and T implicitly)
        # We need to iterate to find rho_gas such that fluid_gas(rho_gas, U_gas) gives pressure P
        if m_gas > 1e-6:
            U_gas_spec = U_gas_tot / m_gas
            try:
                # For given P and U, find the density
                # This requires iteration since we can't directly specify P,U to CoolProp
                # We use PUmass_INPUTS if available, otherwise iterate
                self.fluid_gas.update(CP.PUmass_INPUTS, P, U_gas_spec)
                rho_gas = self.fluid_gas.rhomass()
                V_gas = m_gas / rho_gas
                V_total += V_gas
            except:
                # PUmass might not work for all fluids, try iteration
                # Initial guess: use ideal gas for first estimate
                from scipy.optimize import brentq

                def rho_residual(rho):
                    try:
                        self.fluid_gas.update(CP.DmassUmass_INPUTS, rho, U_gas_spec)
                        return self.fluid_gas.p() - P
                    except:
                        return 1e10

                try:
                    # Search for density that gives the target pressure
                    rho_gas = brentq(rho_residual, 0.01, 1000.0, xtol=1e-6)
                    V_gas = m_gas / rho_gas
                    V_total += V_gas
                except:
                    return 1e10

        # Liquid phase: given P and U, find rho (and T implicitly)
        if m_liquid > 1e-6:
            U_liquid_spec = U_liquid_tot / m_liquid
            try:
                self.fluid_liquid.update(CP.PUmass_INPUTS, P, U_liquid_spec)
                rho_liquid = self.fluid_liquid.rhomass()
                V_liquid = m_liquid / rho_liquid
                V_total += V_liquid
            except:
                # PUmass might not work, try iteration
                from scipy.optimize import brentq

                def rho_residual(rho):
                    try:
                        self.fluid_liquid.update(CP.DmassUmass_INPUTS, rho, U_liquid_spec)
                        return self.fluid_liquid.p() - P
                    except:
                        return 1e10

                try:
                    # Search for density that gives the target pressure
                    rho_liquid = brentq(rho_residual, 100.0, 2000.0, xtol=1e-6)
                    V_liquid = m_liquid / rho_liquid
                    V_total += V_liquid
                except:
                    return 1e10

        # Return volume residual
        vol_residual = (V_total - self.vol) / self.vol
        return vol_residual

    def nem_solve_state(self, m_gas, m_liquid, U_gas_tot, U_liquid_tot, P_guess, T_gas_guess, T_liquid_guess):
        """
        Solve for (T_gas, T_liquid, P) in non-equilibrium model.

        Uses optimization approach inspired by ORS-openthermo: minimize combined
        residuals for volume and energy constraints.

        Parameters
        ----------
        m_gas : float
            Gas phase mass [kg]
        m_liquid : float
            Liquid phase mass [kg]
        U_gas_tot : float
            Total gas internal energy [J]
        U_liquid_tot : float
            Total liquid internal energy [J]
        P_guess : float
            Initial pressure guess [Pa]
        T_gas_guess : float
            Initial gas temperature guess [K]
        T_liquid_guess : float
            Initial liquid temperature guess [K]

        Returns
        -------
        tuple
            (P, T_gas, T_liquid) - solved state [Pa, K, K]
        """
        from scipy.optimize import minimize

        # Initial guess
        x0 = [T_gas_guess, T_liquid_guess, P_guess]

        # Bounds: allow reasonable deviations from previous state
        bounds = [
            (max(T_gas_guess - 50, 200), T_gas_guess + 100),     # T_gas bounds
            (max(T_liquid_guess - 50, 200), T_liquid_guess + 100),  # T_liquid bounds
            (max(P_guess * 0.5, 1e5), P_guess * 2.0)              # P bounds
        ]

        # Use Nelder-Mead (doesn't require derivatives, works well for non-smooth problems)
        result = minimize(
            fun=self.nem_objective_function,
            x0=x0,
            args=(m_gas, m_liquid, U_gas_tot, U_liquid_tot, T_gas_guess, T_liquid_guess),
            method='Nelder-Mead',
            bounds=bounds,
            options={'maxiter': 1000, 'xatol': 1e-6, 'fatol': 1e-9}
        )

        # Accept solution if objective is reasonable
        # Note: Objective = vol_res² + U_gas_res² + U_liq_res² (all normalized)
        # Objective < 1e-4 means ~0.01% error in each constraint
        # Objective < 1e-2 means ~1% error in each constraint
        if result.fun < 1e-2:
            T_gas, T_liquid, P = result.x
            return P, T_gas, T_liquid

        # If convergence is poor, print warning but still return result
        # This can happen when constraints are difficult to satisfy exactly
        if result.fun < 0.1:  # ~10% error - still usable
            print(f"Warning: NEM solver converged with objective={result.fun:.2e} at P={P_guess/1e5:.1f} bar")
            T_gas, T_liquid, P = result.x
            return P, T_gas, T_liquid

        # If Nelder-Mead didn't converge well enough, raise error
        raise RuntimeError(f"NEM solver failed. Objective={result.fun:.2e}, P_guess={P_guess/1e5:.2f} bar, T_gas={T_gas_guess:.1f} K, T_liq={T_liquid_guess:.1f} K")

    def nem_calc_phase_transfer(self, P, T_gas, T_liquid, m_gas, m_liquid, dt):
        """
        Calculate phase transfer based on vapor quality from CoolProp thermodynamic updates.

        After updating each phase with its internal energy, CoolProp returns a vapor quality:
        - For liquid phase: quality > 0 means some liquid has evaporated (transfer to gas)
        - For gas phase: quality < 1 means some gas has condensed (transfer to liquid)

        This approach uses thermodynamics directly rather than relaxation assumptions.

        Parameters
        ----------
        P : float
            Current pressure [Pa]
        T_gas : float
            Gas phase temperature [K]
        T_liquid : float
            Liquid phase temperature [K]
        m_gas : float
            Gas phase mass [kg]
        m_liquid : float
            Liquid phase mass [kg]
        dt : float
            Time step [s]

        Returns
        -------
        dm_transfer : float
            Net mass transfer [kg] (positive = condensation gas→liquid, negative = evaporation liquid→gas)

        Note: This function signature kept for compatibility, but actual phase transfer
        is now calculated directly in the main integration loop using vapor quality.
        """
        # This function is now a placeholder - actual phase transfer calculated
        # in main loop using vapor quality from CoolProp updates
        return 0.0

    def run(self, disable_pbar=True):
        """
        Routine for running the actual problem defined i.e. integrating the mass and energy balances
        """
        # Inititialise / setting initial values for t=0
        if self.isrun is True:
            self.initialize()
        input = self.input
        self.rho[0] = self.rho0
        self.T_fluid[0] = self.T0
        # Initial wall temperature. Default: uniform at T0. For a two-phase start with a
        # superheated gas zone (initial.gas_temperature), the wall has equilibrated with the
        # phase it touches before the blowdown, so initialise the gas-contact (unwetted) wall
        # at the gas temperature and the wetted (liquid-contact) wall at the liquid temperature.
        # An explicit initial.wall_temperature overrides both with a uniform value.
        T_wall_gas = self.T0
        T_wall_wet = self.T0
        if "wall_temperature" in self.input["initial"]:
            T_wall_gas = T_wall_wet = self.input["initial"]["wall_temperature"]
        elif ("gas_temperature" in self.input["initial"]
              and "liquid_level" in self.input["vessel"]):
            T_wall_gas = self.T_gas0
            T_wall_wet = self.T_liquid0
        self.T_vessel[0] = T_wall_gas
        self.T_inner_wall[0] = T_wall_gas
        self.T_outer_wall[0] = T_wall_gas
        self.T_bonded_wall[0] = T_wall_gas
        self.T_vessel_wetted[0] = T_wall_wet
        self.T_inner_wall_wetted[0] = T_wall_wet
        self.T_outer_wall_wetted[0] = T_wall_wet
        self.T_bonded_wall_wetted[0] = T_wall_wet
        self.liquid_level[0] = self.liquid_level0
        if self.input["valve"]["flow"] == "discharge":
            self.T_vent[0] = self.T0
        self.H_mass[0] = self.fluid.hmass()
        self.S_mass[0] = self.fluid.smass()
        self.U_mass[0] = self.fluid.umass()
        self.U_tot[0] = self.fluid.umass() * self.m0
        self.P[0] = self.p0
        self.mass_fluid[0] = self.m0
        if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
            self.vapour_mole_fraction[0] = self.fluid.Q()
        else:
            self.vapour_mole_fraction[0] = 1.0

        # Non-equilibrium model initialization
        if self.non_equilibrium:
            self.T_gas[0] = self.T_gas0
            self.T_liquid[0] = self.T_liquid0
            self.m_gas[0] = self.m_gas0
            self.m_liquid[0] = self.m_liquid0
            if self.m_liquid0 > 0:
                # Two-phase: try saturation, fall back to +1°C superheat
                self.fluid_liquid.update(CP.PQ_INPUTS, self.p0, 0.0)
                self.U_liquid[0] = self.fluid_liquid.umass()
                self.rho_liquid[0] = self.fluid_liquid.rhomass()
                if "gas_temperature" in self.input["initial"]:
                    # Specified superheated gas zone: set U/rho from the warm state so the
                    # first solver step keeps the gas warm (a saturated overwrite here would
                    # collapse T_gas back to saturation at step 1). Superheated gas is off
                    # the Q=1 boundary, so the DmassUmass round-trip below is stable.
                    self.fluid_gas.update(CP.PT_INPUTS, self.p0, self.T_gas0)
                else:
                    try:
                        self.fluid_gas.update(CP.PQ_INPUTS, self.p0, 1.0)
                        _rho = self.fluid_gas.rhomass()
                        _U = self.fluid_gas.umass()
                        # Verify DmassUmass round-trip at saturation
                        self.fluid_gas.update(CP.DmassUmass_INPUTS, _rho, _U)
                        # Verify quality is valid (Q=1 at saturation boundary)
                        _Q = self.fluid_gas.Q()
                        if _Q < 0 or _Q > 1:
                            raise ValueError("Invalid quality at saturation")
                    except Exception:
                        # Saturation boundary unstable - add slight superheat
                        self.T_gas0 += 1.0
                        self.T_gas[0] = self.T_gas0
                        self.fluid_gas.update(CP.PT_INPUTS, self.p0, self.T_gas0)
            else:
                # Single-phase gas
                self.fluid_gas.update(CP.PT_INPUTS, self.p0, self.T_gas0)
                self.U_liquid[0] = 0.0
                self.rho_liquid[0] = 0.0
            self.U_gas[0] = self.fluid_gas.umass()
            self.rho_gas[0] = self.fluid_gas.rhomass()
            self.mdot_phase_transfer[0] = 0.0

        if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
            cpcv = self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) / (
                self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) - 8.314
            )
            rho0 = self.fluid.saturated_vapor_keyed_output(CP.iDmass)
            Z = self.fluid.saturated_vapor_keyed_output(CP.iZ)
        else:
            cpcv = self.fluid.cp0molar() / (self.fluid.cp0molar() - 8.314)
            rho0 = self.rho0
            Z = self.fluid.compressibility_factor()

        massflow_stop_switch = 0

        # Calculating initial mass rate for t=0 depending on mass flow device
        # and filling/discharge mode
        if input["valve"]["type"] == "orifice":
            if input["valve"]["flow"] == "filling":
                k = self.res_fluid.cp0molar() / (self.res_fluid.cp0molar() - 8.314)
                self.mass_rate[0] = -tp.gas_release_rate(
                    self.p_back,
                    self.p0,
                    self.res_fluid.rhomass(),
                    k,
                    self.CD,
                    self.D_orifice**2 / 4 * math.pi,
                )
            else:
                self.mass_rate[0] = tp.gas_release_rate(
                    self.p0,
                    self.p_back,
                    rho0,
                    cpcv,
                    self.CD,
                    self.D_orifice**2 / 4 * math.pi,
                )

        elif input["valve"]["type"] == "hem_release":
            if input["valve"]["flow"] == "filling":
                raise ValueError(
                    "Unsupported valve: ",
                    input["valve"]["type"],
                    " for vessel filling.",
                )
            else:
                self.mass_rate[0] = tp.hem_release_rate(
                    self.p0,
                    self.p_back,
                    self.CD,
                    self.D_orifice**2 / 4 * math.pi,
                    self.fluid,
                )

        elif input["valve"]["type"] == "mdot":
            if "mdot" in input["valve"].keys() and "time" in input["valve"].keys():
                mdot = np.asarray(input["valve"]["mdot"])
                time = np.asarray(input["valve"]["time"])
                max_i = int(time[-1] / self.tstep)
                interp_time = np.linspace(
                    0,
                    self.tstep * len(self.time_array),
                    len(self.time_array),
                    endpoint=False,
                )[:max_i]
                self.mass_rate[:max_i] = np.interp(interp_time, time, mdot)
                if input["valve"]["flow"] == "filling":
                    self.mass_rate *= -1

            else:
                if input["valve"]["flow"] == "filling":
                    self.mass_rate[:] = -input["valve"]["mdot"]
                else:
                    self.mass_rate[:] = input["valve"]["mdot"]

        elif input["valve"]["type"] == "controlvalve":
            Cv = tp.cv_vs_time(
                self.Cv, 0, self.valve_time_constant, self.valve_characteristic
            )
            if input["valve"]["flow"] == "filling":
                Z = self.res_fluid.compressibility_factor()
                MW = self.MW
                k = self.res_fluid.cp0molar() / (self.res_fluid.cp0molar() - 8.314)
                self.mass_rate[0] = -tp.control_valve(
                    self.p_back, self.p0, self.T0, Z, MW, k, Cv
                )
            else:
                MW = self.MW
                k = cpcv
                self.mass_rate[0] = tp.control_valve(
                    self.p0, self.p_back, self.T0, Z, MW, k, Cv
                )
        elif input["valve"]["type"] == "psv":
            if input["valve"]["flow"] == "filling":
                raise ValueError(
                    "Unsupported valve: ",
                    input["valve"]["type"],
                    " for vessel filling.",
                )
            self.mass_rate[0] = tp.relief_valve(
                self.p0,
                self.p_back,
                self.Pset,
                self.blowdown,
                cpcv,
                self.CD,
                self.T0,
                Z,
                self.MW,
                self.D_orifice**2 / 4 * math.pi,
            )

        # Release outflow (thermopack CO2 HEM). Additive: for release-only the valve rate
        # is 0 (valve type "none"); a concurrent valve (future extension) would sum in here.
        if self.has_release:
            self.mass_rate[0] = self.mass_rate[0] + self.compute_release(self.p0, 0)

        self.time_array[0] = 0

        # ============================================================================
        # MAIN TIME INTEGRATION LOOP
        # ============================================================================
        # Integrate mass and energy balances forward in time using explicit Euler method
        # At each time step:
        #   1. Update mass inventory from mass flow rate
        #   2. Calculate new density (mass/volume)
        #   3. Determine thermodynamic state based on calculation method:
        #      - Simple methods (isenthalpic/isentropic/isothermal/constantU):
        #        Direct CoolProp update with (density, H/S/T/U)
        #      - Energy balance method: Solve for temperature from energy balance
        #        accounting for heat transfer, work, and enthalpy changes
        #   4. Calculate heat transfer (if energybalance method)
        #   5. Update vessel wall temperatures
        #   6. Calculate mass flow rate for next time step
        # ============================================================================

        # Initialize heat transfer variables
        T_profile, T_profile2 = 0, 0  # Temperature profiles for detailed wall model
        relief_area = []  # Track relief valve area changes

        for i in tqdm(
            range(1, len(self.time_array)),
            desc="hyddown",
            disable=disable_pbar,
            total=len(self.time_array),
        ):
            self.time_array[i] = self.time_array[i - 1] + self.tstep

            # ---- CO2 triple-point handling (release) ----
            # The CoolProp vessel solver fails once the tank drops below the triple point.
            # Two behaviours: the default freezes the release there; the opt-in
            # solid_in_vessel model hands off to a thermopack three-phase / sublimation
            # model and keeps going (dry ice forms in the vessel).
            if self.has_release and self.solid_in_vessel:
                if (
                    not self.solid_regime
                    and self.P[i - 1] <= self.release_model.P_TRIPLE_EOS * 1.05
                ):
                    # Hand off to the below-triple model. Keep the NEM's two zones: a warm
                    # (superheated) gas zone and a liquid/solid zone, both re-based into
                    # thermopack's energy basis (masses are basis-free). The single-zone
                    # (M_vessel, U_vessel) state is also seeded for the sublimation descent.
                    rm = self.release_model
                    self.solid_regime = True
                    self.tz_plateau = True
                    self.tz_m_gas = max(self.m_gas[i - 1], 1e-6)
                    T_gas_h = min(max(self.T_gas[i - 1], rm._gt_T[0]), rm._gt_T[-1])
                    self.tz_U_gas = self.tz_m_gas * rm.gas_u_at(T_gas_h)
                    self.tz_M_ls = self.m_liquid[i - 1]
                    self.tz_U_ls = self.tz_M_ls * rm.u_l
                    # Carry the cold liquid-contact wall temperature from the above-triple
                    # detailed wall model into the below-triple wetted-wall node, so it keeps
                    # cooling with the liquid/solid instead of resetting to the gas-side wall.
                    self.tz_T_wall_wet = self.T_vessel_wetted[i - 1]
                    self.M_vessel = self.m_liquid[i - 1] + self.m_gas[i - 1]
                    _ml, _mg, self.U_vessel = rm.triple_LG_from_MV(
                        self.M_vessel, self.vol
                    )
                if self.solid_regime:
                    self._solid_regime_step(i)
                    continue
            elif self.has_release and not self.release_frozen:
                # Default freeze: the 15% margin covers a single step's pressure drop so
                # the state never crosses below the triple point. The tank then holds.
                if self.P[i - 1] <= self.release_model.P_TRIPLE * 1.15:
                    self.release_frozen = True
                    self.mass_rate[i - 1] = 0.0

            self.mass_fluid[i] = (
                self.mass_fluid[i - 1] - self.mass_rate[i - 1] * self.tstep
            )

            self.rho[i] = self.mass_fluid[i] / self.vol

            # ------------------------------------------------------------------------
            # THERMODYNAMIC STATE UPDATE
            # ------------------------------------------------------------------------
            # For simple methods, density changes but one other property remains constant
            # CoolProp can directly calculate (T, P) from (density, H/S/U/T) for single components
            # For multicomponent fluids, numerical optimization is required (slower)
            if self.method == "isenthalpic":
                self.fluid.update(CP.DmassHmass_INPUTS, self.rho[i], self.H_mass[i - 1])
                self.T_fluid[i] = self.fluid.T()
                self.P[i] = self.fluid.p()

            elif self.method == "isentropic":
                self.fluid.update(CP.DmassSmass_INPUTS, self.rho[i], self.S_mass[i - 1])
                self.T_fluid[i] = self.fluid.T()
                self.P[i] = self.fluid.p()

            elif self.method == "isothermal":
                self.fluid.update(CP.DmassT_INPUTS, self.rho[i], self.T0)
                self.T_fluid[i] = self.T0
                self.P[i] = self.fluid.p()

            elif self.method == "constantU":
                self.fluid.update(CP.DmassUmass_INPUTS, self.rho[i], self.U_mass[i - 1])
                self.T_fluid[i] = self.fluid.T()
                self.P[i] = self.fluid.p()

            # ------------------------------------------------------------------------
            # ENERGY BALANCE METHOD (Most Complex Case)
            # ------------------------------------------------------------------------
            # Full energy balance accounting for:
            #   - Heat transfer between fluid and vessel wall (convection)
            #   - Heat transfer between vessel and environment (convection, radiation, fire)
            #   - Enthalpy change from mass flow (inlet/outlet streams)
            #   - Work done by expanding/compressing fluid
            #   - 1-D transient conduction through vessel wall (if detailed method)
            #
            # Energy balance on fluid:
            #   dU/dt = Q_transfer + H_inlet*dm_inlet/dt - H_outlet*dm_outlet/dt
            #
            # where U = internal energy, Q = heat transfer rate, H = specific enthalpy
            elif self.method == "energybalance":
                # ====================================================================
                # HEAT TRANSFER COEFFICIENT CALCULATIONS
                # ====================================================================
                # Calculate convective heat transfer coefficients for specified_h or detailed methods
                if self.heat_method == "specified_h" or self.heat_method == "detailed" or self.heat_method == "specified_q":
                    if self.h_in == "calc":
                        if self.vessel_orientation == "horizontal":
                            L = self.diameter
                        else:
                            L = self.length
                        if input["valve"]["flow"] == "filling":
                            # T_film = (self.T_fluid[i - 1] + self.T_vessel[i - 1]) / 2
                            # For NEM: use gas temperature for unwetted (gas-side) heat transfer
                            # For equilibrium: use bulk fluid temperature
                            if self.non_equilibrium and hasattr(self, 'T_gas'):
                                T_for_gas_side_htc = self.T_gas[i - 1]
                            else:
                                T_for_gas_side_htc = self.T_fluid[i - 1]

                            T_film = (
                                T_for_gas_side_htc + self.T_inner_wall[i - 1]
                            ) / 2
                            self.transport_fluid.update(
                                CP.PT_INPUTS, self.P[i - 1], T_film
                            )

                            hi = tp.h_inside_mixed(
                                L,
                                # self.T_vessel[i - 1],
                                self.T_inner_wall[i - i],
                                T_for_gas_side_htc,
                                self.transport_fluid,
                                self.mass_rate[i - 1],
                                self.diameter,
                            )
                            # NEM: Use liquid properties for wetted heat transfer
                            # Equilibrium: Use existing logic with equilibrium fluid
                            if self.non_equilibrium:
                                # For NEM, if liquid exists, use nucleate boiling correlation
                                # No need to check quality - liquid phase is always liquid
                                liquid_exists = self.m_liquid[i-1] > 1e-6
                                if liquid_exists:
                                    # Use liquid temperature for film temperature
                                    T_film_wet = (self.T_liquid[i-1] + self.T_inner_wall_wetted[i-1]) / 2
                                    try:
                                        self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], T_film_wet)
                                    except:
                                        # If film temperature fails, use liquid temperature
                                        try:
                                            self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], self.T_liquid[i-1])
                                        except:
                                            # Last resort: use saturation
                                            self.transport_fluid_wet.update(CP.PQ_INPUTS, self.P[i-1], 0.0)

                                    # Update fluid_liquid to saturated state for h_inside_wetted
                                    # (needed for surface tension and saturated properties)
                                    self.fluid_liquid.update(CP.PQ_INPUTS, self.P[i-1], 0.0)

                                    hiw = tp.h_inside_wetted(
                                        L,
                                        self.T_inner_wall_wetted[i - 1],
                                        self.T_liquid[i - 1],      # Use liquid temperature
                                        self.transport_fluid_wet,
                                        self.fluid_liquid,         # Use liquid phase object at saturation
                                    )
                                else:
                                    # No liquid - use gas-side coefficient
                                    hiw = hi
                            else:
                                # Equilibrium mode: use existing logic
                                if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                                    self.transport_fluid_wet.update(
                                        CP.PT_INPUTS, self.P[i - 1], T_film
                                    )
                                    hiw = tp.h_inside_wetted(
                                        L,
                                        self.T_inner_wall_wetted[i - 1],
                                        self.T_fluid[i - 1],
                                        self.transport_fluid_wet,
                                        self.fluid,
                                    )
                                else:
                                    hiw = hi
                        else:
                            # For NEM: use gas temperature for unwetted (gas-side) heat transfer
                            # For equilibrium: use bulk fluid temperature
                            if self.non_equilibrium and hasattr(self, 'T_gas'):
                                T_for_gas_side_htc = self.T_gas[i - 1]
                            else:
                                T_for_gas_side_htc = self.T_fluid[i - 1]

                            T_film = (
                                T_for_gas_side_htc + self.T_inner_wall[i - 1]
                            ) / 2
                            try:
                                self.transport_fluid.update(
                                    CP.PT_INPUTS, self.P[i - 1], T_film
                                )
                                # Update transport fluid for wetted surface
                                # For NEM, will be updated later with liquid temperature
                                if not self.non_equilibrium:
                                    if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                                        self.transport_fluid_wet.update(
                                            CP.PT_INPUTS, self.P[i - 1], self.T_fluid[i - 1]
                                        )
                            except:
                                self.transport_fluid.update(
                                    CP.PQ_INPUTS, self.P[i - 1], 1.0
                                )
                            hi = tp.h_inside(
                                L,
                                self.T_inner_wall[i - 1],
                                T_for_gas_side_htc,
                                self.transport_fluid,
                            )

                            # NEM: Check liquid phase for boiling, use liquid properties
                            # Equilibrium: Use existing logic with equilibrium fluid
                            if self.non_equilibrium:
                                # For NEM, if liquid exists, use nucleate boiling correlation
                                # No need to check quality - liquid phase is always liquid
                                liquid_exists = self.m_liquid[i-1] > 1e-6
                                if liquid_exists:
                                    # Use liquid temperature for film temperature
                                    T_film_wet = (self.T_liquid[i-1] + self.T_inner_wall_wetted[i-1]) / 2
                                    try:
                                        self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], T_film_wet)
                                    except:
                                        # If film temperature fails, use liquid temperature
                                        try:
                                            self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], self.T_liquid[i-1])
                                        except:
                                            # Last resort: use saturation
                                            self.transport_fluid_wet.update(CP.PQ_INPUTS, self.P[i-1], 0.0)

                                    # Update fluid_liquid to saturated state for h_inside_wetted
                                    # (needed for surface tension and saturated properties)
                                    self.fluid_liquid.update(CP.PQ_INPUTS, self.P[i-1], 0.0)

                                    hiw = tp.h_inside_wetted(
                                        L,
                                        self.T_inner_wall_wetted[i - 1],
                                        self.T_liquid[i - 1],      # Use liquid temperature
                                        self.transport_fluid_wet,
                                        self.fluid_liquid,         # Use liquid phase object at saturation
                                    )
                                else:
                                    # No liquid - use gas-side coefficient
                                    hiw = hi
                            else:
                                # Equilibrium mode: use existing logic
                                if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                                    hiw = tp.h_inside_wetted(
                                        L,
                                        self.T_inner_wall_wetted[i - 1],
                                        self.T_fluid[i - 1],
                                        self.transport_fluid_wet,
                                        self.fluid,
                                    )
                                else:
                                    hiw = hi
                    else:
                        hi = self.h_in
                        hiw = self.h_in  # Use same coefficient for wetted surface

                    self.h_inside[i] = hi
                    self.h_inside_wetted[i] = hiw

                    # ================================================================
                    # TWO-PHASE HEAT TRANSFER (Wetted vs Unwetted Areas)
                    # ================================================================
                    # For two-phase systems, split vessel into two regions:
                    #   1. Wetted area: Wall in contact with liquid phase
                    #   2. Unwetted area: Wall in contact with gas/vapor phase
                    # Different heat transfer coefficients and wall temperatures
                    # are calculated for each region based on liquid level
                    wetted_area = self.inner_vol.SA_from_h(self.liquid_level[i - 1])
                    if np.isnan(wetted_area):
                        wetted_area = 0

                    # Calculate outer wetted area for correct heat transfer on outer surface
                    # Add wall thickness to liquid level height since outer vessel reference
                    # is at bottom of wall, not bottom of inner surface
                    liquid_level_outer = self.liquid_level[i - 1] + self.thickness
                    wetted_area_outer = self.outer_vol.SA_from_h(liquid_level_outer)
                    if np.isnan(wetted_area_outer):
                        wetted_area_outer = 0

                    # Heat transfer from unwetted wall (gas side) to fluid
                    # For NEM: use gas temperature
                    # For equilibrium: use bulk fluid temperature
                    if self.non_equilibrium and hasattr(self, 'T_gas'):
                        T_for_gas_side_htc = self.T_gas[i - 1]
                    else:
                        T_for_gas_side_htc = self.T_fluid[i - 1]

                    self.Q_inner[i] = (
                        (self.surf_area_inner - wetted_area)
                        * hi
                        * (self.T_inner_wall[i - 1] - T_for_gas_side_htc)
                    )
                    self.q_inner[i] = hi * (
                        self.T_inner_wall[i - 1] - T_for_gas_side_htc
                    )

                    # Heat transfer from wetted wall (liquid side) to fluid
                    # Uses different heat transfer coefficient (hiw) for liquid contact
                    # For NEM: use liquid temperature
                    # For equilibrium: use bulk fluid temperature
                    if self.non_equilibrium and hasattr(self, 'T_liquid'):
                        T_fluid_wet = self.T_liquid[i - 1]
                    else:
                        T_fluid_wet = self.T_fluid[i - 1]

                    self.Q_inner_wetted[i] = (
                        wetted_area
                        * hiw
                        * (self.T_inner_wall_wetted[i - 1] - T_fluid_wet)
                    )

                    # Avoid division by zero when fully vapor (wetted_area = 0)
                    if wetted_area > 0:
                        self.q_inner_wetted[i] = self.Q_inner_wetted[i] / wetted_area
                    else:
                        self.q_inner_wetted[i] = 0.0

                    if np.isnan(self.Q_inner_wetted[i]):
                        self.Q_inner_wetted[i] = 0

                    # Heat transfer from environment to unwetted outer wall
                    # Use outer surface area directly (outer surface is exposed to environment)
                    if self.heat_method == "specified_q":
                        # User-defined heat flux (time-dependent or constant)
                        q_ext = self.q_outer_func(self.time_array[i])
                        self.Q_outer[i] = q_ext * (self.surf_area_outer - wetted_area_outer)
                        self.q_outer[i] = q_ext
                        self.Q_outer_wetted[i] = q_ext * wetted_area_outer
                        self.q_outer_wetted[i] = q_ext
                    else:
                        # specified_h or detailed: convective heat transfer
                        self.Q_outer[i] = (
                            (self.surf_area_outer - wetted_area_outer)
                            * self.h_out
                            * (self.Tamb - self.T_outer_wall[i - 1])
                        )
                        self.q_outer[i] = self.h_out * (
                            self.Tamb - self.T_outer_wall[i - 1]
                        )

                        self.Q_outer_wetted[i] = (
                            wetted_area_outer
                            * self.h_out
                            * (self.Tamb - self.T_outer_wall_wetted[i - 1])
                        )

                        self.q_outer_wetted[i] = self.h_out * (
                            self.Tamb - self.T_outer_wall_wetted[i - 1]
                        )
                    if np.isnan(self.Q_outer_wetted[i]):
                        self.Q_outer_wetted[i] = 0

                    self.T_vessel[i] = self.T_vessel[i - 1] + (
                        self.Q_outer[i] - self.Q_inner[i]
                    ) * self.tstep / (
                        self.vessel_cp
                        * self.vessel_density
                        * self.vol_solid
                        * (self.inner_vol.A - wetted_area)
                        / self.inner_vol.A
                    )

                    # Update wetted wall temperature only if wetted area exists
                    if wetted_area < 1e-10 or (self.Q_outer_wetted[i] == 0 and self.Q_inner_wetted[i] == 0):
                        self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1]
                    else:
                        self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1] + (
                            self.Q_outer_wetted[i] - self.Q_inner_wetted[i]
                        ) * self.tstep / (
                            self.vessel_cp
                            * self.vessel_density
                            * self.vol_solid
                            * wetted_area
                            / self.inner_vol.A
                        )

                    # ================================================================
                    # 1-D TRANSIENT HEAT CONDUCTION (Detailed Thermal Model)
                    # ================================================================
                    # Solve transient heat conduction through vessel wall using
                    # finite element method (thermesh module).
                    # Used for vessels with low thermal conductivity (Type III/IV)
                    # or composite materials where temperature gradient through
                    # wall thickness is significant.
                    #
                    # For single-layer walls: Uses isothermal material model
                    # For multi-layer walls: Uses piecewise linear model for liner+shell
                    #
                    # Time integration: Crank-Nicolson (theta=0.5) for stability
                    # Spatial discretization: Linear finite elements with 11 nodes
                    if "thermal_conductivity" in self.input["vessel"].keys():
                        theta = 0.5  # Crank-Nicolson scheme (unconditionally stable, 2nd order)
                        dt = (
                            self.tstep / 10
                        )  # Sub-step for thermal solver (finer time resolution)
                        k, rho, cp = (
                            self.input["vessel"]["thermal_conductivity"],
                            self.vessel_density,
                            self.vessel_cp,
                        )
                        # Check if single-layer or composite (liner + shell) construction
                        if (
                            "liner_thermal_conductivity"
                            not in self.input["vessel"].keys()
                        ):
                            # Single-layer wall construction
                            nn = 11  # number of nodes through wall thickness
                            z = np.linspace(0, self.thickness, nn)
                            self.z = z
                            # Create meshes for unwetted and wetted regions
                            mesh = tm.Mesh(z, tm.LinearElement)
                            mesh_w = tm.Mesh(z, tm.LinearElement)
                            # Material model with constant properties
                            cpeek = tm.isothermal_model(k, rho, cp)
                            cpeek_w = tm.isothermal_model(k, rho, cp)

                            # Initialize temperature profile on first time step
                            if type(T_profile) == type(int()) and T_profile == 0:
                                bc = [
                                    {"T": self.T0},
                                    {"T": self.Tamb},
                                ]
                                domain = tm.Domain(mesh, [cpeek], bc)
                                domain.set_T(
                                    (self.Tamb + self.T0) / 2 * np.ones(len(mesh.nodes))
                                )
                                solver = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded, T_profile = tm.solve_ht(domain, solver)

                                bc_w = [
                                    {"T": self.T0},
                                    {"T": self.Tamb},
                                ]
                                domain_w = tm.Domain(mesh_w, [cpeek_w], bc_w)
                                domain_w.set_T(
                                    (self.Tamb + self.T0)
                                    / 2
                                    * np.ones(len(mesh_w.nodes))
                                )
                                solver_w = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile_w = tm.solve_ht(
                                    domain_w, solver_w
                                )
                            else:
                                # Boundary conditions: z=0 is outer wall, z=L is inner wall
                                bc = [
                                    {
                                        "q": self.q_outer[i]
                                        # / (self.surf_area_outer - wetted_area_outer)
                                    },
                                    {
                                        "q": -self.q_inner[i]
                                        # / (self.surf_area_inner - wetted_area)
                                    },
                                ]
                                domain = tm.Domain(mesh, [cpeek], bc)
                                domain.set_T(T_profile[-1, :])
                                solver = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded, T_profile = tm.solve_ht(domain, solver)
                                # Wetted wall will be solved below if liquid is present

                            self.temp_profile.append(T_profile[-1, :])
                            self.T_outer_wall[i] = T_profile[-1, 0]
                            self.T_inner_wall[i] = T_profile[-1, -1]

                            # Only solve wetted wall if liquid is present
                            if self.liquid_level[i - 1] > 0 and wetted_area > 0:
                                bc_w = [
                                    {
                                        "q": self.q_outer_wetted[i]
                                    },  # / wetted_area_outer},
                                    {"q": -self.q_inner_wetted[i]},  # / wetted_area},
                                ]
                                domain_w = tm.Domain(mesh_w, [cpeek_w], bc_w)
                                domain_w.set_T(T_profile_w[-1, :])
                                solver_w = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_w, T_profile_w = tm.solve_ht(domain_w, solver_w)
                                self.T_outer_wall_wetted[i] = T_profile_w[-1, 0]
                                self.T_inner_wall_wetted[i] = T_profile_w[-1, -1]
                            else:
                                # No liquid - wetted wall same as unwetted
                                self.T_outer_wall_wetted[i] = T_profile[-1, 0]
                                self.T_inner_wall_wetted[i] = T_profile[-1, -1]
                        else:
                            k_liner = self.input["vessel"]["liner_thermal_conductivity"]
                            rho_liner = self.input["vessel"]["liner_density"]
                            cp_liner = self.input["vessel"]["liner_heat_capacity"]
                            liner = tm.isothermal_model(k_liner, rho_liner, cp_liner)
                            shell = tm.isothermal_model(k, rho, cp)
                            liner_w = tm.isothermal_model(k_liner, rho_liner, cp_liner)
                            shell_w = tm.isothermal_model(k, rho, cp)

                            thk = self.input["vessel"]["thickness"]  # thickness in m
                            nn = 11  # number of nodes
                            z_shell = np.linspace(0, thk, nn)  # node locations

                            thk = self.input["vessel"]["liner_thickness"]
                            z_liner = np.linspace(-thk, 0, nn)  # node locations
                            z2 = np.hstack((z_liner, z_shell[1:]))
                            self.z = z2
                            mesh2 = tm.Mesh(z2, tm.LinearElement)
                            mesh2_w = tm.Mesh(z2, tm.LinearElement)
                            for j, elem in enumerate(mesh2.elem):
                                if elem.nodes.mean() > 0.0:
                                    mesh2.subdomain[j] = 1
                                    mesh2_w.subdomain[j] = 1

                            if type(T_profile2) == type(int()) and T_profile2 == 0:
                                bc = [
                                    {"T": self.T0},
                                    {"T": self.Tamb},
                                ]
                                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                                domain2.set_T(
                                    (self.Tamb + self.T0)
                                    / 2
                                    * np.ones(len(mesh2.nodes))
                                )
                                solver2 = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded, T_profile2 = tm.solve_ht(domain2, solver2)
                                bc_w = [
                                    {"T": self.T0},
                                    {"T": self.Tamb},
                                ]
                                domain2_w = tm.Domain(mesh2_w, [liner_w, shell_w], bc_w)
                                domain2_w.set_T(
                                    (self.Tamb + self.T0)
                                    / 2
                                    * np.ones(len(mesh2.nodes))
                                )
                                solver2_w = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile2_w = tm.solve_ht(
                                    domain2_w, solver2_w
                                )
                            else:
                                # Boundary conditions: z=-liner_thickness is inner, z=thickness is outer
                                bc = [
                                    {
                                        "q": -self.q_inner[i]
                                        # / (self.surf_area_inner - wetted_area)
                                    },
                                    {
                                        "q": self.q_outer[i]
                                        # / (self.surf_area_outer - wetted_area_outer)
                                    },
                                ]
                                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                                domain2.set_T(T_profile2[-1, :])
                                solver2 = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded, T_profile2 = tm.solve_ht(domain2, solver2)
                                bc_w = [
                                    {"q": -self.q_inner_wetted[i]},  #  / wetted_area},
                                    {
                                        "q": self.q_outer_wetted[i]
                                    },  # / wetted_area_outer},
                                ]
                                domain2_w = tm.Domain(mesh2_w, [liner_w, shell_w], bc_w)
                                domain2_w.set_T(T_profile2_w[-1, :])
                                solver2_w = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile2_w = tm.solve_ht(
                                    domain2_w, solver2_w
                                )
                            self.T_outer_wall[i] = T_profile2[-1, -1]
                            self.T_inner_wall[i] = T_profile2[-1, 0]
                            self.T_bonded_wall[i] = T_profile2[-1, (nn - 1)]
                            self.T_outer_wall_wetted[i] = T_profile2_w[-1, -1]
                            self.T_inner_wall_wetted[i] = T_profile2_w[-1, 0]
                            self.T_bonded_wall_wetted[i] = T_profile2_w[-1, (nn - 1)]
                            self.temp_profile.append(T_profile2[-1, :])
                    else:
                        # Lumped capacitance model (no 1D heat transfer)
                        self.T_inner_wall[i] = self.T_vessel[i]
                        self.T_outer_wall[i] = self.T_vessel[i]
                        self.T_inner_wall_wetted[i] = self.T_vessel_wetted[i]
                        self.T_outer_wall_wetted[i] = self.T_vessel_wetted[i]

                        # Calculate Biot number to validate lumped capacitance assumption
                        # Bi = h*L/k where L is characteristic length (wall thickness)
                        # Bi << 0.1: lumped model valid (uniform wall temperature)
                        # Bi > 0.1: thermal gradient significant, 1D model recommended
                        if "thickness" in self.input["vessel"]:
                            L_char = self.thickness
                            # Check for thermal_conductivity_biot first (for Biot calc only)
                            # Otherwise check thermal_conductivity (which triggers 1D model)
                            k_wall = None
                            if "thermal_conductivity_biot" in self.input["vessel"]:
                                k_wall = self.input["vessel"][
                                    "thermal_conductivity_biot"
                                ]
                            elif "thermal_conductivity" in self.input["vessel"]:
                                k_wall = self.input["vessel"]["thermal_conductivity"]

                            if k_wall is None:
                                self.Biot[i] = np.nan
                                self.Biot_wetted[i] = np.nan
                            else:
                                self.Biot[i] = hi * L_char / k_wall
                                self.Biot_wetted[i] = (
                                    hiw * L_char / k_wall if hiw > 0 else 0.0
                                )

                                # Issue warning if Biot number exceeds threshold
                                if i == 1 or (i % 100 == 0 and self.Biot[i] > 0.1):
                                    if self.Biot[i] > 0.1:
                                        import warnings

                                        warnings.warn(
                                            f"t={self.time_array[i]:.1f}s: Biot number (unwetted) = {self.Biot[i]:.3f} > 0.1. "
                                            f"Lumped capacitance model may be inaccurate. Consider using 1D heat transfer "
                                            f"by specifying 'thermal_conductivity' in vessel properties.",
                                            UserWarning,
                                        )
                                if (
                                    self.liquid_level[i - 1] > 0
                                    and self.Biot_wetted[i] > 0.1
                                ):
                                    if i == 1 or i % 100 == 0:
                                        import warnings

                                        warnings.warn(
                                            f"t={self.time_array[i]:.1f}s: Biot number (wetted) = {self.Biot_wetted[i]:.3f} > 0.1. "
                                            f"Lumped capacitance model may be inaccurate.",
                                            UserWarning,
                                        )
                        else:
                            self.Biot[i] = np.nan
                            self.Biot_wetted[i] = np.nan

                elif self.heat_method == "s-b":
                    if self.vessel_orientation == "horizontal":
                        L = self.diameter
                    else:
                        L = self.length

                    wetted_area = self.inner_vol.SA_from_h(self.liquid_level[i - 1])
                    if np.isnan(wetted_area):
                        wetted_area = 0

                    # Calculate outer wetted area for correct heat transfer on outer surface
                    # Add wall thickness to liquid level height since outer vessel reference
                    # is at bottom of wall, not bottom of inner surface
                    liquid_level_outer = self.liquid_level[i - 1] + self.thickness
                    wetted_area_outer = self.outer_vol.SA_from_h(liquid_level_outer)
                    if np.isnan(wetted_area_outer):
                        wetted_area_outer = 0

                    # Determine which wall temperature to use for internal heat transfer
                    # For 1D heat transfer: Use actual inner wall surface temperature
                    # For 0D model: Use bulk vessel temperature
                    if "thermal_conductivity" in self.input["vessel"].keys():
                        # 1D model: Use inner wall surface temperature
                        T_wall_inner = self.T_inner_wall[i - 1]
                        T_wall_inner_wetted = self.T_inner_wall_wetted[i - 1]
                    else:
                        # 0D model: Use bulk vessel temperature
                        T_wall_inner = self.T_vessel[i - 1]
                        T_wall_inner_wetted = self.T_vessel_wetted[i - 1]

                    # For NEM: use gas temperature for unwetted (gas-side) heat transfer
                    # For equilibrium: use bulk fluid temperature
                    if self.non_equilibrium and hasattr(self, 'T_gas'):
                        T_for_gas_side = self.T_gas[i - 1]
                    else:
                        T_for_gas_side = self.T_fluid[i - 1]

                    hi = tp.h_inner(
                        L,
                        T_for_gas_side,
                        T_wall_inner,
                        self.P[i - 1],
                        self.species,
                    )
                    self.h_inside[i] = hi

                    # NEM: Check liquid phase for boiling, use liquid properties
                    # Equilibrium: Use existing logic with equilibrium fluid
                    if self.non_equilibrium:
                        # For NEM, if liquid exists, use nucleate boiling correlation
                        # No need to check quality - liquid phase is always liquid
                        liquid_exists = self.m_liquid[i-1] > 1e-6
                        if liquid_exists:
                            # Use liquid temperature and liquid phase object
                            try:
                                self.transport_fluid_wet.update(
                                    CP.PT_INPUTS, self.P[i - 1], self.T_liquid[i - 1]
                                )
                            except:
                                # If update fails, use saturation
                                self.transport_fluid_wet.update(CP.PQ_INPUTS, self.P[i - 1], 0.0)

                            # Update fluid_liquid to saturated state for h_inside_wetted
                            # (needed for surface tension and saturated properties)
                            self.fluid_liquid.update(CP.PQ_INPUTS, self.P[i - 1], 0.0)

                            hiw = tp.h_inside_wetted(
                                L,
                                T_wall_inner_wetted,
                                self.T_liquid[i - 1],      # Use liquid temperature
                                self.transport_fluid_wet,
                                self.fluid_liquid,         # Use liquid phase object at saturation
                            )
                        else:
                            # No liquid - use gas-side coefficient
                            hiw = hi
                    else:
                        # Equilibrium mode: use existing logic
                        if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                            self.transport_fluid_wet.update(
                                CP.PT_INPUTS, self.P[i - 1], self.T_fluid[i - 1]
                            )
                            hiw = tp.h_inside_wetted(
                                L,
                                T_wall_inner_wetted,
                                self.T_fluid[i - 1],
                                self.transport_fluid_wet,
                                self.fluid,
                            )
                        else:
                            hiw = hi

                    # For unwetted heat transfer, use gas temperature in NEM
                    if self.non_equilibrium and hasattr(self, 'T_gas'):
                        T_for_gas_side_htc = self.T_gas[i - 1]
                    else:
                        T_for_gas_side_htc = self.T_fluid[i - 1]

                    self.Q_inner[i] = self.scaling * (
                        (self.surf_area_inner - wetted_area)
                        * hi
                        * (T_wall_inner - T_for_gas_side_htc)
                    )

                    self.q_inner[i] = hi * (T_wall_inner - T_for_gas_side_htc)

                    # For wetted heat transfer, use liquid temperature in NEM
                    if self.non_equilibrium and hasattr(self, 'T_liquid'):
                        T_fluid_wet = self.T_liquid[i - 1]
                    else:
                        T_fluid_wet = self.T_fluid[i - 1]

                    self.Q_inner_wetted[i] = self.scaling * (
                        wetted_area * hiw * (T_wall_inner_wetted - T_fluid_wet)
                    )
                    self.q_inner_wetted[i] = hiw * (
                        T_wall_inner_wetted - T_fluid_wet
                    )
                    if np.isnan(self.Q_inner_wetted[i]):
                        self.Q_inner_wetted[i] = 0

                    # ================================================================
                    # FIRE HEAT LOAD CALCULATIONS (Stefan-Boltzmann Method)
                    # ================================================================
                    # Calculate external heat flux from fire using Stefan-Boltzmann
                    # equation accounting for radiative and convective heat transfer.
                    # Fire types: API 521 pool/jet fire, Scandpower pool/jet/peak fires
                    # Heat flux is temperature-dependent (higher vessel T → more re-radiation)
                    #
                    # For 1D heat transfer model: Use outer wall temperature (surface exposed to fire)
                    # For 0D model: Use bulk vessel temperature
                    # For two-phase systems: Calculate separate heat fluxes for
                    # wetted (liquid contact) and unwetted (gas contact) regions

                    # Determine which temperature to use for fire heat flux calculation
                    if "thermal_conductivity" in self.input["vessel"].keys():
                        # 1D model: Use actual outer wall surface temperature
                        T_fire_surface = self.T_outer_wall[i - 1]
                        T_fire_surface_wetted = self.T_outer_wall_wetted[i - 1]
                    else:
                        # 0D model: Use bulk vessel temperature
                        T_fire_surface = self.T_vessel[i - 1]
                        T_fire_surface_wetted = self.T_vessel_wetted[i - 1]

                    # Fire heats the outer surface - use outer surface area directly
                    self.Q_outer[i] = (
                        self.scaling
                        * fire.sb_fire(T_fire_surface, self.fire_type)
                        * (self.surf_area_outer - wetted_area_outer)
                    )

                    self.q_outer[i] = fire.sb_fire(T_fire_surface, self.fire_type)

                    self.Q_outer_wetted[i] = self.scaling * (
                        fire.sb_fire(T_fire_surface_wetted, self.fire_type)
                        * wetted_area_outer
                    )
                    self.q_outer_wetted[i] = fire.sb_fire(
                        T_fire_surface_wetted, self.fire_type
                    )

                    if np.isnan(self.Q_outer_wetted[i]):
                        self.Q_outer_wetted[i] = 0

                    # ================================================================
                    # 1-D TRANSIENT HEAT CONDUCTION (Fire Scenario with Detailed Thermal Model)
                    # ================================================================
                    # Solve transient heat conduction through vessel wall using
                    # finite element method (thermesh module) for fire scenarios.
                    # This section is activated when thermal_conductivity is specified.
                    # Otherwise, falls back to simple lumped capacitance model below.
                    if "thermal_conductivity" in self.input["vessel"].keys():
                        theta = 0.5  # Crank-Nicolson scheme (unconditionally stable, 2nd order)
                        dt = (
                            self.tstep / 10
                        )  # Sub-step for thermal solver (finer time resolution)
                        k, rho, cp = (
                            self.input["vessel"]["thermal_conductivity"],
                            self.vessel_density,
                            self.vessel_cp,
                        )
                        # Check if single-layer or composite (liner + shell) construction
                        if (
                            "liner_thermal_conductivity"
                            not in self.input["vessel"].keys()
                        ):
                            # Single-layer wall construction
                            nn = 11  # number of nodes through wall thickness
                            z = np.linspace(0, self.thickness, nn)
                            self.z = z
                            # Create meshes for unwetted and wetted regions
                            mesh = tm.Mesh(z, tm.LinearElement)
                            mesh_w = tm.Mesh(z, tm.LinearElement)
                            # Material model with constant properties
                            cpeek = tm.isothermal_model(k, rho, cp)
                            cpeek_w = tm.isothermal_model(k, rho, cp)

                            # Initialize temperature profile on first time step
                            if type(T_profile) == type(int()) and T_profile == 0:
                                bc = [
                                    {"T": self.T0},
                                    {"T": self.T0},
                                ]
                                domain = tm.Domain(mesh, [cpeek], bc)
                                domain.set_T(self.T0 * np.ones(len(mesh.nodes)))
                                solver = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded, T_profile = tm.solve_ht(domain, solver)

                                bc_w = [
                                    {"T": self.T0},
                                    {"T": self.T0},
                                ]
                                domain_w = tm.Domain(mesh_w, [cpeek_w], bc_w)
                                domain_w.set_T(self.T0 * np.ones(len(mesh_w.nodes)))
                                solver_w = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile_w = tm.solve_ht(
                                    domain_w, solver_w
                                )
                            else:
                                # Boundary conditions: z=0 is inner wall, z=L is outer wall
                                bc = [
                                    {
                                        "q": -self.q_inner[i]
                                        # / (self.surf_area_inner - wetted_area)
                                    },
                                    {
                                        "q": self.q_outer[i]
                                        # / (self.surf_area_outer - wetted_area_outer)
                                    },
                                ]
                                domain = tm.Domain(mesh, [cpeek], bc)
                                domain.set_T(T_profile[-1, :])
                                solver = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded, T_profile = tm.solve_ht(domain, solver)
                                # Wetted wall boundary conditions will be set below
                                # in the conditional block that checks liquid_level

                            # Only solve wetted wall if liquid is present
                            # Check liquid_level to handle both gas-only and liquid-depleted cases
                            if self.liquid_level[i - 1] > 0 and wetted_area > 0:
                                bc_w = [
                                    {"q": -self.q_inner_wetted[i]},  # / wetted_area},
                                    {
                                        "q": self.q_outer_wetted[i]
                                    },  # , / wetted_area_outer},
                                ]
                                domain_w = tm.Domain(mesh_w, [cpeek_w], bc_w)
                                domain_w.set_T(T_profile_w[-1, :])
                                solver_w = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_w, T_profile_w = tm.solve_ht(domain_w, solver_w)
                                self.T_inner_wall_wetted[i] = T_profile_w[-1, 0]
                                self.T_outer_wall_wetted[i] = T_profile_w[-1, -1]
                                self.T_vessel_wetted[i] = T_profile_w[-1, :].mean()
                            else:
                                # No liquid - wetted wall same as unwetted
                                self.T_inner_wall_wetted[i] = T_profile[-1, 0]
                                self.T_outer_wall_wetted[i] = T_profile[-1, -1]
                                self.T_vessel_wetted[i] = T_profile[-1, :].mean()

                            self.temp_profile.append(T_profile[-1, :])
                            self.T_inner_wall[i] = T_profile[-1, 0]
                            self.T_outer_wall[i] = T_profile[-1, -1]
                            # Update mean vessel temperature from wall temperatures
                            self.T_vessel[i] = T_profile[-1, :].mean()
                        else:
                            # Composite wall construction (liner + shell)
                            k_liner = self.input["vessel"]["liner_thermal_conductivity"]
                            rho_liner = self.input["vessel"]["liner_density"]
                            cp_liner = self.input["vessel"]["liner_heat_capacity"]
                            liner = tm.isothermal_model(k_liner, rho_liner, cp_liner)
                            shell = tm.isothermal_model(k, rho, cp)
                            liner_w = tm.isothermal_model(k_liner, rho_liner, cp_liner)
                            shell_w = tm.isothermal_model(k, rho, cp)

                            thk = self.input["vessel"]["thickness"]  # thickness in m
                            nn = 11  # number of nodes
                            z_shell = np.linspace(0, thk, nn)  # node locations

                            thk = self.input["vessel"]["liner_thickness"]
                            z_liner = np.linspace(-thk, 0, nn)  # node locations
                            z2 = np.hstack((z_liner, z_shell[1:]))
                            self.z = z2
                            mesh2 = tm.Mesh(z2, tm.LinearElement)
                            mesh2_w = tm.Mesh(z2, tm.LinearElement)
                            for j, elem in enumerate(mesh2.elem):
                                if elem.nodes.mean() > 0.0:
                                    mesh2.subdomain[j] = 1
                                    mesh2_w.subdomain[j] = 1

                            if type(T_profile2) == type(int()) and T_profile2 == 0:
                                bc = [
                                    {"T": self.T0},
                                    {"T": self.T0},
                                ]
                                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                                domain2.set_T(self.T0 * np.ones(len(mesh2.nodes)))
                                solver2 = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded, T_profile2 = tm.solve_ht(domain2, solver2)
                                bc_w = [
                                    {"T": self.T0},
                                    {"T": self.T0},
                                ]
                                domain2_w = tm.Domain(mesh2_w, [liner_w, shell_w], bc_w)
                                domain2_w.set_T(self.T0 * np.ones(len(mesh2.nodes)))
                                solver2_w = {
                                    "dt": 100,
                                    "t_end": 10000,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile2_w = tm.solve_ht(
                                    domain2_w, solver2_w
                                )
                            else:
                                # Boundary conditions: z=-liner_thickness is inner, z=thickness is outer
                                bc = [
                                    {
                                        "q": -self.q_inner[i]
                                        # / (self.surf_area_inner - wetted_area)
                                    },
                                    {
                                        "q": self.q_outer[i]
                                        # / (self.surf_area_outer - wetted_area_outer)
                                    },
                                ]
                                domain2 = tm.Domain(mesh2, [liner, shell], bc)
                                domain2.set_T(T_profile2[-1, :])
                                solver2 = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded, T_profile2 = tm.solve_ht(domain2, solver2)

                                # Only solve wetted wall if liquid is present
                                # Check liquid_level to handle both gas-only and liquid-depleted cases
                                if self.liquid_level[i - 1] > 0 and wetted_area > 0:
                                    bc_w = [
                                        {
                                            "q": -self.q_inner_wetted[i]
                                        },  # / wetted_area},
                                        {
                                            "q": self.q_outer_wetted[i]
                                            # / wetted_area_outer
                                        },
                                    ]
                                    domain2_w = tm.Domain(
                                        mesh2_w, [liner_w, shell_w], bc_w
                                    )
                                    domain2_w.set_T(T_profile2_w[-1, :])
                                    solver2_w = {
                                        "dt": dt,
                                        "t_end": self.tstep,
                                        "theta": theta,
                                    }
                                    t_bonded_w, T_profile2_w = tm.solve_ht(
                                        domain2_w, solver2_w
                                    )
                                    self.T_inner_wall_wetted[i] = T_profile2_w[-1, 0]
                                    self.T_outer_wall_wetted[i] = T_profile2_w[-1, -1]
                                    self.T_bonded_wall_wetted[i] = T_profile2_w[
                                        -1, (nn - 1)
                                    ]
                                    self.T_vessel_wetted[i] = T_profile2_w[-1, :].mean()
                                else:
                                    # No liquid - wetted wall same as unwetted
                                    self.T_inner_wall_wetted[i] = T_profile2[-1, 0]
                                    self.T_outer_wall_wetted[i] = T_profile2[-1, -1]
                                    self.T_bonded_wall_wetted[i] = T_profile2[
                                        -1, (nn - 1)
                                    ]
                                    self.T_vessel_wetted[i] = T_profile2[-1, :].mean()

                            self.T_inner_wall[i] = T_profile2[-1, 0]
                            self.T_outer_wall[i] = T_profile2[-1, -1]
                            self.T_bonded_wall[i] = T_profile2[-1, (nn - 1)]
                            self.temp_profile.append(T_profile2[-1, :])
                            # Update mean vessel temperature from wall temperatures
                            self.T_vessel[i] = T_profile2[-1, :].mean()
                    else:
                        # ================================================================
                        # SIMPLE LUMPED CAPACITANCE MODEL (No thermal gradient)
                        # ================================================================
                        # When thermal_conductivity is not specified, use simple 0D model
                        # with uniform vessel wall temperature (no spatial gradient)
                        self.T_vessel[i] = self.T_vessel[i - 1] + (
                            (self.Q_outer[i] - self.Q_inner[i]) / self.scaling
                        ) * self.tstep / (
                            self.vessel_cp
                            * self.vessel_density
                            * self.vol_solid
                            * (self.inner_vol.A - wetted_area)
                            / self.inner_vol.A
                        )
                        if self.liquid_level[i - 1] > 0:
                            self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1] + (
                                (self.Q_outer_wetted[i] - self.Q_inner_wetted[i])
                                / self.scaling
                            ) * self.tstep / (
                                self.vessel_cp
                                * self.vessel_density
                                * self.vol_solid
                                * wetted_area
                                / self.inner_vol.A
                            )
                        else:
                            # Hack to heat up previous liquid wetted surface
                            # Fire heats outer surface - use outer surface area directly
                            self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1] + (
                                fire.sb_fire(
                                    self.T_vessel_wetted[i - 1], self.fire_type
                                )
                                * (self.surf_area_outer - wetted_area_outer)
                                - (self.surf_area_inner - wetted_area)
                                * hi
                                * (self.T_vessel_wetted[i - 1] - self.T_fluid[i - 1])
                            ) * self.tstep / (
                                self.vessel_cp
                                * self.vessel_density
                                * self.vol_solid
                                * (self.inner_vol.A - wetted_area)
                                / self.inner_vol.A
                            )

                        if np.isnan(self.T_vessel_wetted[i]):
                            self.T_vessel_wetted[i] = self.T_vessel[i]

                        # Calculate Biot number to validate lumped capacitance assumption
                        # Bi = h*L/k where L is characteristic length (wall thickness)
                        # Bi << 0.1: lumped model valid (uniform wall temperature)
                        # Bi > 0.1: thermal gradient significant, 1D model recommended
                        if "thickness" in self.input["vessel"]:
                            L_char = (
                                self.thickness
                            )  # Characteristic length = wall thickness
                            # Use inner heat transfer coefficient (typically controls)
                            # and vessel thermal conductivity if available
                            # Check for thermal_conductivity_biot first (for Biot calc only)
                            # Otherwise check thermal_conductivity (which triggers 1D model)
                            k_wall = None
                            if "thermal_conductivity_biot" in self.input["vessel"]:
                                k_wall = self.input["vessel"][
                                    "thermal_conductivity_biot"
                                ]
                            elif "thermal_conductivity" in self.input["vessel"]:
                                k_wall = self.input["vessel"]["thermal_conductivity"]

                            if k_wall is None:
                                # Lumped model without k specified - can't calc Bi, set to NaN
                                self.Biot[i] = np.nan
                                self.Biot_wetted[i] = np.nan
                            else:
                                self.Biot[i] = hi * L_char / k_wall
                                self.Biot_wetted[i] = (
                                    hiw * L_char / k_wall if hiw > 0 else 0.0
                                )

                                # Issue warning if Biot number exceeds threshold
                                if i == 1 or (i % 100 == 0 and self.Biot[i] > 0.1):
                                    if self.Biot[i] > 0.1:
                                        import warnings

                                        warnings.warn(
                                            f"t={self.time_array[i]:.1f}s: Biot number (unwetted) = {self.Biot[i]:.3f} > 0.1. "
                                            f"Lumped capacitance model may be inaccurate. Consider using 1D heat transfer "
                                            f"by specifying 'thermal_conductivity' in vessel properties.",
                                            UserWarning,
                                        )
                                if (
                                    self.liquid_level[i - 1] > 0
                                    and self.Biot_wetted[i] > 0.1
                                ):
                                    if i == 1 or i % 100 == 0:
                                        import warnings

                                        warnings.warn(
                                            f"t={self.time_array[i]:.1f}s: Biot number (wetted) = {self.Biot_wetted[i]:.3f} > 0.1. "
                                            f"Lumped capacitance model may be inaccurate.",
                                            UserWarning,
                                        )
                        else:
                            self.Biot[i] = np.nan
                            self.Biot_wetted[i] = np.nan

                        self.T_inner_wall[i] = self.T_vessel[i]
                        self.T_outer_wall[i] = self.T_vessel[i]
                        self.T_inner_wall_wetted[i] = self.T_vessel_wetted[i]
                        self.T_outer_wall_wetted[i] = self.T_vessel_wetted[i]

                elif self.heat_method == "specified_U":
                    self.Q_inner[i] = (
                        self.surf_area_outer
                        * self.Ufix
                        * (self.Tamb - self.T_fluid[i - 1])
                    )
                    self.T_vessel[i] = self.T_vessel[0]
                elif self.heat_method == "specified_Q":
                    self.Q_inner[i] = self.Qfix
                    self.T_vessel[i] = self.T_vessel[0]
                else:
                    self.Q_inner[i] = 0.0
                    self.T_vessel[i] = self.T_vessel[0]

                # NMOL = self.mass_fluid[i - 1] / self.MW
                # NMOL_ADD = (self.mass_fluid[i] - self.mass_fluid[i - 1]) / self.MW
                # New
                U_start = self.U_mass[i - 1] * self.mass_fluid[i - 1]

                # Smooting finction for very early times /numerical trick
                # Might not be necessary.
                x = 1 - math.exp(-1 * self.time_array[i]) ** 0.66

                # Finding the inlet/outlet enthalpy rate for the energy balance
                if input["valve"]["flow"] == "filling":
                    # h_in = self.fluid.hmass()
                    h_in = x * self.res_fluid.hmass() + (1 - x) * self.res_fluid.umass()

                else:
                    # h_in = self.fluid.hmass()
                    if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                        h_in = self.fluid.saturated_vapor_keyed_output(CP.iHmass)
                    else:
                        h_in = self.fluid.hmass()

                if i > 1:
                    P1 = self.P[i - 2]
                else:
                    P1 = self.P[i - 1]

                U_end = (
                    U_start
                    - self.tstep * self.mass_rate[i - 1] * h_in
                    + self.tstep * self.Q_inner[i]
                    + self.tstep * self.Q_inner_wetted[i]
                )

                self.U_mass[i] = U_end / self.mass_fluid[i]

                # ====================================================================
                # NON-EQUILIBRIUM MODEL (NEM) ENERGY BALANCE
                # ====================================================================
                # For non-equilibrium model, solve separate energy balances for gas and liquid phases
                if self.non_equilibrium:
                    # First, update masses with valve/release flow only (no phase transfer yet)
                    if input["valve"]["flow"] == "discharge":
                        if self.has_release and self.release_phase == "liquid":
                            # Liquid-space release: the outflow leaves the liquid inventory.
                            # (release_phase flips to "gas" once the liquid is exhausted, so
                            # this branch stops draining liquid at that point.)
                            self.m_liquid[i] = self.m_liquid[i-1] - self.mass_rate[i-1] * self.tstep
                            self.m_gas[i] = self.m_gas[i-1]
                        else:
                            # Vapour-space release / standard blowdown: leaves the gas phase.
                            self.m_gas[i] = self.m_gas[i-1] - self.mass_rate[i-1] * self.tstep
                            self.m_liquid[i] = self.m_liquid[i-1]
                    else:  # filling
                        # Mass enters through gas phase
                        self.m_gas[i] = self.m_gas[i-1] - self.mass_rate[i-1] * self.tstep
                        self.m_liquid[i] = self.m_liquid[i-1]

                    # Prevent negative masses
                    if self.m_gas[i] < 0:
                        self.m_gas[i] = 0.0
                    if self.m_liquid[i] < 0:
                        self.m_liquid[i] = 0.0

                    # ====================================================================
                    # PHASE TRANSFER (Part of Mass Balance)
                    # ====================================================================
                    # Check if previous timestep had two-phase conditions in either phase
                    # Transfer mass and energy between phases based on vapor quality
                    # This happens BEFORE energy balance, as part of mass redistribution

                    dm_evap_mass = 0.0  # Mass evaporated (liquid → gas)
                    dm_cond_mass = 0.0  # Mass condensed (gas → liquid)
                    E_evap = 0.0  # Energy transferred with evaporation
                    E_cond = 0.0  # Energy transferred with condensation

                    # Check liquid phase from previous timestep
                    if self.m_liquid[i-1] > 1e-6:
                        try:
                            # Check if liquid was in two-phase region at previous state
                            self.fluid_liquid.update(CP.DmassUmass_INPUTS, self.rho_liquid[i-1], self.U_liquid[i-1])
                            quality_liquid = self.fluid_liquid.Q()

                            if quality_liquid > 0 and quality_liquid < 1:
                                # Liquid is in two-phase region - some should evaporate
                                relax_factor = 0.8  # Transfer 80% per timestep
                                dm_evap_mass = quality_liquid * self.m_liquid[i-1] * relax_factor

                                # Safety limit: evaporation rate limited by available heat
                                # Maximum dm/dt from heat input: Q / h_fg
                                h_liq_sat_check = self.fluid_liquid.saturated_liquid_keyed_output(CP.iHmass)
                                h_vap_sat_check = self.fluid_liquid.saturated_vapor_keyed_output(CP.iHmass)
                                h_fg = h_vap_sat_check - h_liq_sat_check

                                # Estimate available heat for evaporation (use PREVIOUS timestep)
                                Q_available = max(self.Q_inner_wetted[i-1], 0.0) * self.tstep  # J
                                dm_evap_max = Q_available / h_fg if h_fg > 0 else dm_evap_mass

                                # Limit to 50x heat-based rate to allow for stored energy
                                dm_evap_mass = min(dm_evap_mass, dm_evap_max * 50.0)

                                # Transfer mass
                                self.m_liquid[i] -= dm_evap_mass
                                self.m_gas[i] += dm_evap_mass

                                # Energy transferred with phase change
                                # Evaporated mass leaves liquid as VAPOR (phase change occurs)
                                # The evaporated liquid becomes saturated vapor
                                # Use (h+u)/2 compromise between enthalpy and internal energy
                                h_vap_sat = self.fluid_liquid.saturated_vapor_keyed_output(CP.iHmass)
                                u_vap_sat = self.fluid_liquid.saturated_vapor_keyed_output(CP.iUmass)
                                e_vap_sat = (h_vap_sat + u_vap_sat) / 2.0

                                # E_evap = energy carried by evaporated mass
                                # This energy is SUBTRACTED from liquid and ADDED to gas
                                E_evap = dm_evap_mass * e_vap_sat
                        except:
                            pass  # Not in two-phase, no transfer

                    # Check gas phase from previous timestep (only if no evaporation)
                    if self.m_gas[i-1] > 1e-6 and dm_evap_mass == 0:
                        try:
                            # Check if gas was in two-phase region at previous state
                            self.fluid_gas.update(CP.DmassUmass_INPUTS, self.rho_gas[i-1], self.U_gas[i-1])
                            quality_gas = self.fluid_gas.Q()

                            if quality_gas > 0 and quality_gas < 1:
                                # Gas is in two-phase region - some should condense
                                relax_factor = 0.8  # Transfer 80% per timestep
                                dm_cond_mass = (1.0 - quality_gas) * self.m_gas[i-1] * relax_factor

                                # Safety limit: condensation rate limited by heat removal
                                h_liq_sat_check = self.fluid_gas.saturated_liquid_keyed_output(CP.iHmass)
                                h_vap_sat_check = self.fluid_gas.saturated_vapor_keyed_output(CP.iHmass)
                                h_fg = h_vap_sat_check - h_liq_sat_check

                                # Estimate heat removal rate (use PREVIOUS timestep)
                                Q_removal = max(-self.Q_inner[i-1], 0.0) * self.tstep  # J (positive value)
                                dm_cond_max = Q_removal / h_fg if h_fg > 0 else dm_cond_mass

                                # Limit to 50x heat-based rate
                                dm_cond_mass = min(dm_cond_mass, dm_cond_max * 50.0)

                                # Transfer mass
                                self.m_gas[i] -= dm_cond_mass
                                self.m_liquid[i] += dm_cond_mass

                                # Energy transferred with phase change
                                # Condensed mass leaves gas as LIQUID (phase change occurs)
                                # The condensed vapor becomes saturated liquid
                                # Use (h+u)/2 compromise between enthalpy and internal energy
                                h_liq_sat = self.fluid_gas.saturated_liquid_keyed_output(CP.iHmass)
                                u_liq_sat = self.fluid_gas.saturated_liquid_keyed_output(CP.iUmass)
                                e_liq_sat = (h_liq_sat + u_liq_sat) / 2.0

                                # E_cond = energy carried by condensed mass
                                # This energy is SUBTRACTED from gas and ADDED to liquid
                                E_cond = dm_cond_mass * e_liq_sat
                        except:
                            pass  # Not in two-phase, no transfer

                    # Store phase transfer rate for output
                    net_transfer = dm_cond_mass - dm_evap_mass  # Positive = condensation
                    self.mdot_phase_transfer[i] = net_transfer / self.tstep

                    # Calculate gas-liquid interfacial heat transfer
                    # Q_gl = h_gl * A_interface * (T_gas - T_liquid)
                    # Positive = heat from gas to liquid
                    if self.m_liquid[i-1] > 1e-6 and self.m_gas[i-1] > 1e-6:
                        # Calculate interface area from liquid level using TANK geometry
                        V_liquid_prev = self.m_liquid[i-1] / self.rho_liquid[i-1]
                        ll_prev = self.inner_vol.h_from_V(V_liquid_prev)

                        # Use fluids.TANK.A_cross_sectional for exact interface area
                        # Handles horizontal/vertical, all head types (hemispherical, F&D, etc.)
                        A_interface = self.inner_vol.A_cross_sectional(ll_prev)

                        # Heat transfer coefficient (W/m²K)
                        # Can be specified in input file, calculated, or use default
                        if "h_gas_liquid" in input["calculation"]:
                            h_gl_input = input["calculation"]["h_gas_liquid"]

                            if isinstance(h_gl_input, str) and h_gl_input.lower() in ("calc", "calc_two_sided"):
                                # Calculate h_gl using natural convection correlation
                                L_char = self.inner_vol.D

                                try:
                                    if h_gl_input.lower() == "calc_two_sided":
                                        h_gl = tp.h_gas_liquid_interface_two_sided(
                                            self.T_gas[i-1],
                                            self.T_liquid[i-1],
                                            self.P[i-1],
                                            L_char,
                                            self.fluid_gas,
                                            self.fluid_liquid
                                        )
                                    else:
                                        h_gl = tp.h_gas_liquid_interface(
                                            self.T_gas[i-1],
                                            self.T_liquid[i-1],
                                            self.P[i-1],
                                            L_char,
                                            self.fluid_gas,
                                            self.fluid_liquid
                                        )
                                except Exception as e:
                                    if i == 1 or (i < 100 and i % 50 == 0):
                                        import warnings
                                        warnings.warn(f"h_gl correlation failed at t={self.time_array[i]:.1f}s: {str(e)[:80]}, using fallback h_gl=100 W/m²K")
                                    h_gl = 100.0
                            else:
                                # User-specified fixed value
                                h_gl = float(h_gl_input)
                        else:
                            h_gl = 50.0  # Default value

                        # Store h_gl in array for diagnostics
                        self.h_gas_liquid[i] = h_gl

                        # Heat transfer rate (W)
                        Q_gas_liquid = h_gl * A_interface * (self.T_gas[i-1] - self.T_liquid[i-1])
                    else:
                        Q_gas_liquid = 0.0
                        self.h_gas_liquid[i] = 0.0

                    # ====================================================================
                    # ENERGY BALANCE
                    # ====================================================================
                    # Apply heat transfer and flow work to each phase
                    # Phase transfer energy already accounted for in mass balance

                    # Gas phase energy balance:
                    # dU_gas = Q_wall_gas - Q_gas_liquid + h_valve*dm_valve + E_phase_transfer
                    U_gas_start = self.U_gas[i-1] * self.m_gas[i-1]

                    # Enthalpy for valve flow
                    if input["valve"]["flow"] == "filling":
                        h_valve = self.res_fluid.hmass()
                        U_gas_tentative = (U_gas_start
                                          - self.tstep * self.mass_rate[i-1] * h_valve
                                          + self.tstep * self.Q_inner[i]
                                          - self.tstep * Q_gas_liquid
                                          + E_evap - E_cond)  # Phase transfer energy
                    else:  # discharge
                        # Attribute the outflow enthalpy to the phase the mass actually
                        # leaves from. A liquid-space release draws liquid (release_phase
                        # == "liquid"); a vapour-space release / standard blowdown draws
                        # gas. Exactly one of h_gas_out / h_liq_out is non-zero.
                        liquid_release = self.has_release and self.release_phase == "liquid"
                        h_gas_out = 0.0
                        h_liq_out = 0.0
                        if liquid_release:
                            if self.m_liquid[i-1] > 1e-6:
                                self.fluid_liquid.update(
                                    CP.DmassUmass_INPUTS, self.rho_liquid[i-1], self.U_liquid[i-1]
                                )
                                h_liq_out = self.fluid_liquid.hmass()
                        else:
                            if self.m_gas[i-1] > 1e-6:
                                # Use DmassUmass to get enthalpy (avoids PT issues at saturation)
                                self.fluid_gas.update(CP.DmassUmass_INPUTS, self.rho_gas[i-1], self.U_gas[i-1])
                                h_gas_out = self.fluid_gas.hmass()
                        U_gas_tentative = (U_gas_start
                                          - self.tstep * self.mass_rate[i-1] * h_gas_out
                                          + self.tstep * self.Q_inner[i]
                                          - self.tstep * Q_gas_liquid
                                          + E_evap - E_cond)  # Phase transfer energy

                    # Liquid phase energy balance:
                    # dU_liquid = Q_wall_liquid + Q_gas_liquid - h_liq_out*dm_liquid - E_phase_transfer
                    U_liquid_start = self.U_liquid[i-1] * self.m_liquid[i-1]
                    # Outflow enthalpy leaving the liquid phase (0 unless a liquid release)
                    liquid_out_term = 0.0
                    if input["valve"]["flow"] == "discharge":
                        liquid_out_term = self.tstep * self.mass_rate[i-1] * h_liq_out
                    U_liquid_tentative = (U_liquid_start
                                         - liquid_out_term
                                         + self.tstep * self.Q_inner_wetted[i]
                                         + self.tstep * Q_gas_liquid
                                         - E_evap + E_cond)  # Phase transfer energy (opposite sign)

                    # Final internal energies for P-solver
                    U_gas_end = U_gas_tentative
                    U_liquid_end = U_liquid_tentative

                    # Safety check: if gas mass is very small, switch to single-phase liquid
                    if self.m_gas[i] < 1e-3:  # Less than 1 gram of gas
                        # Essentially single-phase liquid - add remaining gas energy to liquid
                        U_liquid_end += U_gas_end
                        U_gas_end = 0.0
                        self.m_gas[i] = 0.0

                    # Safety check: if liquid mass is very small, switch to single-phase gas
                    if self.m_liquid[i] < 1e-3:  # Less than 1 gram of liquid
                        # Essentially single-phase gas - add remaining liquid energy to gas
                        U_gas_end += U_liquid_end
                        U_liquid_end = 0.0
                        self.m_liquid[i] = 0.0

                    # Solve for pressure P such that volume constraint is satisfied
                    # Given P and U, CoolProp directly gives T and ρ (no nested iteration needed!)
                    from scipy.optimize import brentq, brenth, ridder, bisect, minimize, newton

                    def volume_residual(P):
                        """
                        For given pressure P, update phases with (P, U) and check volume constraint.
                        Returns (V_gas + V_liquid - V_vessel) / V_vessel
                        """
                        try:
                            U_gas_specific = U_gas_end / self.m_gas[i]
                            U_liquid_specific = U_liquid_end / self.m_liquid[i]

                            # Update gas with (P, U) → CoolProp gives ρ and T directly
                            self.fluid_gas.update(CP.PUmass_INPUTS, P, U_gas_specific)
                            rho_gas = self.fluid_gas.rhomass()

                            # Update liquid with (P, U) → CoolProp gives ρ and T directly
                            self.fluid_liquid.update(CP.PUmass_INPUTS, P, U_liquid_specific)
                            rho_liquid = self.fluid_liquid.rhomass()

                            # Calculate volumes
                            V_gas = self.m_gas[i] / rho_gas
                            V_liquid = self.m_liquid[i] / rho_liquid
                            V_total = V_gas + V_liquid

                            # Return normalized volume residual
                            return ((V_total - self.vol) / self.vol)

                        except Exception:
                            # CoolProp PUmass_INPUTS can fail near the critical point.
                            # Return a sign-consistent residual based on physical behavior:
                            #   low P → expanded gas → V_total > V_vessel → positive
                            #   high P → compressed gas → V_total < V_vessel → negative
                            # Using +1.0 for all failures creates false sign changes
                            # that mislead bracket-based solvers (brentq/bisect).
                            if P > self.P[i-1]:
                                return -10.0  # High P side: compressed → negative
                            else:
                                return 10.0   # Low P side: expanded → positive

                    try:
                        # Check if we're essentially single-phase
                        if self.m_gas[i] < 1e-6 and self.m_liquid[i] > 1e-6:
                            # Single-phase liquid - use DU update
                            U_total_specific = U_liquid_end / self.m_liquid[i]
                            rho_total = self.mass_fluid[i] / self.vol
                            self.fluid.update(CP.DmassUmass_INPUTS, rho_total, U_total_specific)
                            P_solution = self.fluid.p()
                        elif self.m_liquid[i] < 1e-6 and self.m_gas[i] > 1e-6:
                            # Single-phase gas - use DU update
                            U_total_specific = U_gas_end / self.m_gas[i]
                            rho_total = self.mass_fluid[i] / self.vol
                            self.fluid.update(CP.DmassUmass_INPUTS, rho_total, U_total_specific)
                            P_solution = self.fluid.p()
                        else:
                            # Two-phase - solve for pressure equilibrium
                            # Two-phase VLE physically cannot exist above P_crit,
                            # so capping bounds at P_crit*0.95 is correct here.
                            # (When liquid depletes, code switches to single-phase
                            # DU path above, which has no pressure cap.)
                            P_crit = self.fluid_gas.p_critical()
                            P_min = max(self.P[i-1] * 0.5, 1e5)  # At least 1 bar
                            P_max = min(self.P[i-1] * 2.0, P_crit * 0.95)

                            # Solve for pressure that gives correct volume
                            # Use fallback solver chain for robustness
                            P_solution = None
                            solver_method = None

                            try:
                                # Method 1: Try brentq (fastest)
                                P_solution = brentq(volume_residual, P_min, P_max, xtol=1e-5, maxiter=100)
                                solver_method = "brentq"
                            except ValueError:
                                # Root not bracketed - try ridder (more robust)
                                try:
                                    P_solution = ridder(volume_residual, P_min, P_max, xtol=1e-5, maxiter=100)
                                    solver_method = "ridder"
                                except ValueError:
                                    # Still not bracketed - expand bounds and try bisect
                                    P_min_expanded = max(P_min * 0.2, 1e5)  # Expand to 0.2x
                                    P_max_expanded = min(P_max * 5.0, P_crit * 0.95)

                                    try:
                                        P_solution = bisect(volume_residual, P_min_expanded, P_max_expanded, xtol=1e-5, maxiter=200)
                                        solver_method = "bisect_expanded"
                                    except ValueError:
                                        # Last resort: check if previous pressure gives acceptable residual
                                        res_prev = volume_residual(self.P[i-1])
                                        if abs(res_prev) < 0.05:  # Within 5% of target volume
                                            P_solution = self.P[i-1]
                                            solver_method = "prev_pressure"
                                        else:
                                            # Give up - cannot find solution
                                            raise ValueError(
                                                f"All solvers failed. Residual at P_prev={self.P[i-1]/1e5:.2f} bar: {res_prev:.4f}"
                                            )

                        # Store solved pressure
                        self.P[i] = P_solution

                        # Update gas state with solved pressure
                        if self.m_gas[i] > 1e-6:
                            U_gas_specific = U_gas_end / self.m_gas[i]
                            self.fluid_gas.update(CP.PUmass_INPUTS, P_solution, U_gas_specific)
                            self.rho_gas[i] = self.fluid_gas.rhomass()
                            self.U_gas[i] = self.fluid_gas.umass()
                            self.T_gas[i] = self.fluid_gas.T()
                        else:
                            self.rho_gas[i] = 0.0
                            self.U_gas[i] = 0.0
                            self.T_gas[i] = self.T_gas[i-1]

                        # Update liquid state with solved pressure
                        if self.m_liquid[i] > 1e-6:
                            U_liquid_specific = U_liquid_end / self.m_liquid[i]
                            self.fluid_liquid.update(CP.PUmass_INPUTS, P_solution, U_liquid_specific)
                            self.rho_liquid[i] = self.fluid_liquid.rhomass()
                            self.U_liquid[i] = self.fluid_liquid.umass()
                            self.T_liquid[i] = self.fluid_liquid.T()
                        else:
                            self.rho_liquid[i] = 0.0
                            self.U_liquid[i] = 0.0
                            self.T_liquid[i] = self.T_liquid[i-1]

                        # Update combined state for compatibility
                        self.T_fluid[i] = self.T_gas[i] if self.m_gas[i] > 1e-6 else self.T_liquid[i]
                        self.rho[i] = self.mass_fluid[i] / self.vol

                        # Update main fluid object for compatibility (best-effort). Below the
                        # CO2 triple point an equilibrium D,U flash of the (warm, non-equilibrium)
                        # inventory can land in the solid region and raise, even while the NEM
                        # zones legitimately hold a pressure above the triple point. The NEM
                        # tracks the two zones explicitly and the solid-regime handoff keys off
                        # the NEM pressure, so on failure keep the last valid self.fluid state.
                        total_U = (U_gas_end + U_liquid_end) / self.mass_fluid[i]
                        try:
                            self.fluid.update(CP.DmassUmass_INPUTS, self.rho[i], total_U)
                        except ValueError:
                            pass

                        # Calculate liquid level for NEM
                        if self.m_liquid[i] > 1e-6 and self.rho_liquid[i] > 1e-6:
                            V_liquid = self.m_liquid[i] / self.rho_liquid[i]
                            self.liquid_level[i] = self.inner_vol.h_from_V(V_liquid)
                        else:
                            self.liquid_level[i] = 0.0

                        # Phase transfer already handled in mass balance step
                        # No post-solver adjustments needed

                    except Exception as e:
                        # Fail hard - do not allow pressure equilibrium violations
                        raise RuntimeError(
                            f"NEM solver failed at t={self.time_array[i]:.2f}s: {e}\n"
                            f"  P_guess range: [{P_min:.2f}, {P_max:.2f}] Pa\n"
                            f"  m_gas: {self.m_gas[i]:.3f} kg, m_liquid: {self.m_liquid[i]:.3f} kg\n"
                            f"  U_gas_end: {U_gas_end:.1f} J, U_liquid_end: {U_liquid_end:.1f} J\n"
                            f"Solver must maintain strict pressure equilibrium. Check:\n"
                            f"  1. Timestep may be too large (current: {self.tstep}s)\n"
                            f"  2. Energy balance may have unphysical values\n"
                            f"  3. Phase transfer rate may be too aggressive"
                        ) from e

                # Not pretty if-statement and a hack for fire relief area estimation. Most cases go directly to the first ...else... clause
                elif input["valve"]["type"] == "relief":
                    if self.Pset <= self.P[i - 1]:
                        if self.liquid_level[i - 1] > 0:

                            self.fluid.update(
                                CP.HmassP_INPUTS,
                                self.fluid.hmass()
                                + self.tstep
                                * (self.Q_inner[i] + self.Q_inner_wetted[i])
                                / self.mass_fluid[i],
                                self.Pset,
                            )

                            total_mass = self.fluid.rhomass() * self.vol

                            Hvap = self.fluid.saturated_vapor_keyed_output(
                                CP.iHmass
                            ) - self.fluid.saturated_liquid_keyed_output(CP.iHmass)
                            BOG = (
                                (self.Q_inner_wetted[i] + self.Q_inner[i])
                                / (Hvap)
                                * 3600
                                * 24
                                / self.mass_fluid[i]
                                * 100
                            )
                            # print("BOG rate at relief (%%/day): ", BOG)

                            self.mass_rate[i] = (
                                self.Q_inner_wetted[i] + self.Q_inner[i]
                            ) / Hvap

                            P1 = self.Pset
                            self.P[i] = P1

                            self.T_fluid[i] = self.fluid.T()

                            Z = self.fluid.saturated_vapor_keyed_output(CP.iZ)
                            cp0molar = self.fluid.saturated_vapor_keyed_output(
                                CP.iCp0molar
                            )

                            self.relief_area[i] = fluids.API520_A_g(
                                self.mass_rate[i],
                                self.fluid.T(),
                                Z,
                                self.MW * 1000,
                                cp0molar / (cp0molar - 8.314),
                                P1,
                                self.p_back,
                                0.975,
                                1,
                                1,
                            )
                        else:
                            T1 = self.PHproblem(
                                h_in
                                + self.tstep * self.Q_inner[i] / self.mass_fluid[i],
                                self.Pset,
                                Tguess=self.T_fluid[i - 1] + 5,
                                relief=True,
                            )
                            self.T_fluid[i] = T1
                            P1 = self.Pset
                            self.P[i] = P1
                            self.T_fluid[i] = T1
                            self.fluid.update(CP.PT_INPUTS, self.P[i], T1)
                            self.mass_rate[i] = (
                                (
                                    (1 / self.fluid.rhomass() - 1 / self.rho[i])
                                    * self.mass_fluid[i]
                                )
                                * self.rho[i]
                                / self.tstep
                            )
                            self.relief_area[i] = fluids.API520_A_g(
                                self.mass_rate[i],
                                T1,
                                self.fluid.compressibility_factor(),
                                self.MW * 1000,
                                self.fluid.cp0molar() / (self.fluid.cp0molar() - 8.314),
                                P1,
                                self.p_back,
                                0.975,
                                1,
                                1,
                            )
                    else:
                        if self.liquid_level[i - 1] > 0:
                            self.mass_rate[i] = 0
                            # self.fluid.update(CP.PQ_INPUTS, self.P[i], Q)
                            self.fluid.update(
                                CP.DmassUmass_INPUTS,
                                self.rho[i],
                                U_end / self.mass_fluid[i],
                            )
                            self.T_fluid[i] = self.fluid.T()
                            self.P[i] = self.fluid.p()

                        else:
                            self.mass_rate[i] = 0
                            P1, T1, self.U_res[i] = self.UDproblem(
                                U_end / self.mass_fluid[i],
                                self.rho[i],
                                self.P[i - 1],
                                self.T_fluid[i - 1],
                            )

                            self.P[i] = P1
                            self.T_fluid[i] = T1
                            self.fluid.update(CP.PT_INPUTS, self.P[i], self.T_fluid[i])

                else:
                    P1, T1, self.U_res[i] = self.UDproblem(
                        U_end / self.mass_fluid[i],
                        self.rho[i],
                        self.P[i - 1],
                        self.T_fluid[i - 1],
                    )

                    self.P[i] = P1
                    self.T_fluid[i] = T1

                    if len(self.molefracs) == 1 and self.molefracs[0] == 1.0:
                        self.fluid.update(
                            CP.DmassUmass_INPUTS,
                            self.rho[i],
                            self.U_mass[i],
                        )
                    else:
                        try:
                            self.fluid.update(CP.PT_INPUTS, self.P[i], self.T_fluid[i])
                        except:
                            if self.fluid.Q() < 0:
                                self.fluid.update(CP.PQ_INPUTS, self.P[i], 1)
                            else:
                                self.fluid.update(
                                    CP.PQ_INPUTS, self.P[i], self.fluid.Q()
                                )
                    if (
                        self.input["valve"]["flow"] == "discharge"
                        and self.fluid.Q() < 1
                        and self.fluid.Q() >= 0
                    ):
                        self.res_fluid.update(CP.PQ_INPUTS, self.P[i], 1.0)

            else:
                raise NameError("Unknown calculation method: " + self.method)

            Q = self.fluid.Q()
            if Q >= 0 and Q <= 1:
                self.vapour_mole_fraction[i] = Q
            else:
                self.vapour_mole_fraction[i] = 1

            self.H_mass[i] = self.fluid.hmass()
            self.S_mass[i] = self.fluid.smass()
            self.U_mass[i] = self.fluid.umass()

            self.liquid_level[i] = self.calc_liquid_level()

            # Calculating vent temperature (adiabatic) only for discharge problem
            if self.input["valve"]["flow"] == "discharge":
                if "&" in self.species:
                    self.T_vent[i] = self.PHproblem(
                        self.H_mass[i], self.p_back, self.vent_fluid.T()
                    )
                else:
                    try:
                        self.T_vent[i] = PropsSI(
                            "T", "H", self.H_mass[i], "P", self.p_back, self.species
                        )
                    except:
                        self.T_vent[i] = self.vent_fluid.T()

            # For NEM: use gas phase properties for discharge calculations
            # For equilibrium: use main fluid object
            if self.non_equilibrium and self.m_gas[i] > 1e-6:
                # NEM: check if gas phase is two-phase or superheated
                Q_gas = self.fluid_gas.Q()
                if Q_gas >= 0 and Q_gas <= 1:
                    # Two-phase: use saturated vapor properties
                    cpcv = self.fluid_gas.saturated_vapor_keyed_output(CP.iCpmolar) / (
                        self.fluid_gas.saturated_vapor_keyed_output(CP.iCpmolar) - 8.314
                    )
                    Z = self.fluid_gas.saturated_vapor_keyed_output(CP.iZ)
                else:
                    # Superheated gas: use actual gas phase properties
                    cpcv = self.fluid_gas.cp0molar() / (self.fluid_gas.cp0molar() - 8.314)
                    Z = self.fluid_gas.compressibility_factor()
            elif self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                # Equilibrium two-phase: use saturated vapor properties
                cpcv = self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) / (
                    self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) - 8.314
                )
                Z = self.fluid.saturated_vapor_keyed_output(CP.iZ)
            else:
                # Single phase gas: use actual fluid properties
                cpcv = self.fluid.cp0molar() / (self.fluid.cp0molar() - 8.314)
                Z = self.fluid.compressibility_factor()
            # ====================================================================
            # MASS FLOW RATE CALCULATIONS (End of Time Step)
            # ====================================================================
            # Calculate mass flow rate for next time step based on valve type
            # and current thermodynamic state. Uses current pressure, temperature,
            # density to determine compressible flow through valve/orifice.
            #
            # Valve types:
            #   - orifice: Yellow Book compressible flow equation (critical/subcritical)
            #   - controlvalve: Control valve Cv characteristic with time-dependent opening
            #   - psv: Relief valve with API 520/521 sizing equations
            #   - mdot: Constant or time-varying mass flow rate (already set)
            #
            # For two-phase systems: Use saturated vapor properties for mass flow calc
            # Note: "relief" valve type handled separately (not here)
            if input["valve"]["type"] == "orifice":
                if input["valve"]["flow"] == "filling":
                    k = self.res_fluid.cp0molar() / (self.res_fluid.cp0molar() - 8.314)
                    self.mass_rate[i] = -tp.gas_release_rate(
                        self.p_back,
                        self.P[i],
                        self.res_fluid.rhomass(),
                        k,
                        self.CD,
                        self.D_orifice**2 / 4 * math.pi,
                    )
                else:
                    # For NEM: use gas phase density (actual state, not always saturated)
                    # For equilibrium two-phase: use saturated vapor density
                    if self.non_equilibrium and self.m_gas[i] > 1e-6:
                        Q_gas = self.fluid_gas.Q()
                        if Q_gas >= 0 and Q_gas <= 1:
                            rho = self.fluid_gas.saturated_vapor_keyed_output(CP.iDmass)
                        else:
                            rho = self.fluid_gas.rhomass()
                    elif self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                        rho = self.fluid.saturated_vapor_keyed_output(CP.iDmass)
                    else:
                        rho = self.rho[i]
                    self.mass_rate[i] = tp.gas_release_rate(
                        self.P[i],
                        self.p_back,
                        rho,
                        cpcv,
                        self.CD,
                        self.D_orifice**2 / 4 * math.pi,
                    )
            elif input["valve"]["type"] == "hem_release":
                if input["valve"]["flow"] == "filling":
                    raise ValueError(
                        "Filling flow not supported for HEM release valve type"
                    )
                else:
                    self.mass_rate[i] = tp.hem_release_rate(
                        self.P[i],
                        self.p_back,
                        self.CD,
                        self.D_orifice**2 / 4 * math.pi,
                        self.fluid,
                    )
            elif input["valve"]["type"] == "controlvalve":
                Cv = tp.cv_vs_time(
                    self.Cv,
                    self.time_array[i],
                    self.valve_time_constant,
                    self.valve_characteristic,
                )
                if input["valve"]["flow"] == "filling":
                    Z = self.res_fluid.compressibility_factor()
                    MW = self.MW
                    k = self.res_fluid.cp0molar() / (self.res_fluid.cp0molar() - 8.314)
                    self.mass_rate[i] = -tp.control_valve(
                        self.p_back, self.P[i], self.T0, Z, MW, k, Cv
                    )
                else:
                    Z = Z
                    MW = self.MW
                    # For NEM: use gas temperature for discharge
                    T_discharge = self.T_gas[i] if self.non_equilibrium and self.m_gas[i] > 1e-6 else self.T_fluid[i]
                    self.mass_rate[i] = tp.control_valve(
                        self.P[i], self.p_back, T_discharge, Z, MW, cpcv, Cv
                    )
            elif input["valve"]["type"] == "psv":
                # For NEM: use gas temperature for discharge
                T_discharge = self.T_gas[i] if self.non_equilibrium and self.m_gas[i] > 1e-6 else self.T_fluid[i]
                self.mass_rate[i] = tp.relief_valve(
                    self.P[i],
                    self.p_back,
                    self.Pset,
                    self.blowdown,
                    cpcv,
                    self.CD,
                    T_discharge,
                    Z,
                    self.MW,
                    self.D_orifice**2 / 4 * math.pi,
                )
            # Release outflow (thermopack CO2 HEM), additive to any valve rate.
            if self.has_release:
                # Switch a liquid-space release to vapour once the liquid is exhausted.
                if (
                    self.release_type == "liquid"
                    and self.release_phase == "liquid"
                    and self.m_liquid[i] <= 1e-6
                ):
                    self.release_phase = "gas"
                self.mass_rate[i] = self.mass_rate[i] + self.compute_release(self.P[i], i)
            if (
                "end_pressure" in self.input["valve"]
                and self.input['valve']['flow']=='filling' and self.P[i] > self.input["valve"]["end_pressure"]
            ):
                massflow_stop_switch = 1 
            elif (
                "end_pressure" in self.input["valve"]
                and self.input['valve']['flow']=='discharge' and self.P[i] < self.input["valve"]["end_pressure"]
            ):
                massflow_stop_switch = 1
            else:
                massflow_stop_switch = 0 
            if massflow_stop_switch:
                self.mass_rate[i] = 0
        self.isrun = True

        # Cumulative dry-ice mass released to atmosphere (kg): the atmospheric solid mass
        # fraction times the release-hole mass flow, integrated over the run.
        if self.has_release:
            dryice_rate = self.x_solid_atm * self.release_rate  # kg/s of dry ice
            self.m_dryice_cum = np.cumsum(dryice_rate * self.tstep)

        if input["valve"]["type"] == "relief":
            idx_max = self.mass_rate.argmax()
            # Smooth peak value by averaging with neighbors (avoid array index out of bounds)
            if 0 < idx_max < len(self.mass_rate) - 1:
                self.mass_rate[idx_max] = (
                    self.mass_rate[idx_max - 1] + self.mass_rate[idx_max + 1]
                ) / 2
                self.relief_area[idx_max] = (
                    self.relief_area[idx_max - 1] + self.relief_area[idx_max + 1]
                ) / 2
            # print("Relief area:", 2*math.sqrt(max(relief_area[1:])/math.pi), max(self.mass_rate))

    def get_dataframe(self):
        """
        Export simulation results to pandas DataFrame for analysis and archiving.

        Collects all time-series results from the simulation and organizes them
        into a pandas DataFrame with labeled columns. Useful for exporting to
        CSV/Excel, post-processing, or custom plotting.

        Returns
        -------
        pd.DataFrame
            DataFrame with simulation results. Columns include:
            - Time (s)
            - Pressure (bar)
            - Fluid temperature (°C)
            - Wall temperature (°C)
            - Vent temperature (°C)
            - Fluid enthalpy (J/kg)
            - Fluid entropy (J/(kg·K))
            - Fluid internal energy (J/kg)
            - Discharge mass rate (kg/s)
            - Fluid mass (kg)
            - Fluid density (kg/m³)
            - Inner heat transfer coefficient (W/(m²·K))
            - Internal heat flux (W/m²)
            - External heat flux (W/m²)
            - Inner wall temperature (°C)
            - Outer wall temperature (°C)

        Notes
        -----
        Only returns data if simulation has been run (self.isrun == True).
        Temperature values are converted from Kelvin to Celsius.
        Pressure values are converted from Pa to bar.
        """
        if self.isrun == True:
            df = pd.DataFrame(self.time_array, columns=["Time (s)"])

            df.insert(1, "Pressure (bar)", self.P / 1e5, True)
            df.insert(2, "Fluid temperature (oC)", self.T_fluid - 273.15, True)
            df.insert(3, "Wall temperature  (oC)", self.T_vessel - 273.15, True)
            df.insert(4, "Vent temperature  (oC)", self.T_vent - 273.15, True)
            df.insert(5, "Fluid enthalpy (J/kg)", self.H_mass, True)
            df.insert(6, "Fluid entropy (J/kg K)", self.S_mass, True)
            df.insert(7, "Fluid internal energy (J/kg)", self.U_mass, True)
            df.insert(8, "Discharge mass rate (kg/s)", self.mass_rate, True)
            df.insert(9, "Fluid mass (kg)", self.mass_fluid, True)
            df.insert(10, "Fluid density (kg/m3)", self.rho, True)
            df.insert(
                11, "Inner heat transfer coefficient (W/m2 K)", self.h_inside, True
            )
            df.insert(
                12,
                "Internal heat flux (W/m2)",
                self.Q_inner / self.surf_area_inner,
                True,
            )
            df.insert(
                13,
                "External heat flux (W/m2)",
                self.Q_outer / self.surf_area_outer,
                True,
            )
            df.insert(
                14, "Inner wall temperature  (oC)", self.T_inner_wall - 273.15, True
            )
            df.insert(
                13, "Outer wall temperature  (oC)", self.T_outer_wall - 273.15, True
            )
            # CO2 release / dry-ice atmospheric state (appended so column indices are stable)
            if self.has_release:
                # In-vessel inventory breakdown (gas / liquid / solid dry ice)
                df["Vessel gas mass (kg)"] = self.m_gas
                df["Vessel liquid mass (kg)"] = self.m_liquid
                df["Vessel dry-ice mass (kg)"] = self.m_solid
                df["Release mass rate (kg/s)"] = self.release_rate
                df["Atmospheric temperature (oC)"] = self.T_atm - 273.15
                df["Atmospheric vapour mass fraction (-)"] = self.x_vap_atm
                df["Atmospheric dry-ice mass fraction (-)"] = self.x_solid_atm
                df["Throat dry-ice mass fraction (-)"] = self.solid_frac_throat
                df["Cumulative dry-ice mass (kg)"] = self.m_dryice_cum
        return df

    def plot(self, filename=None, verbose=True):
        """
        Creating standard plots for the solved problem

        Parameters
        ----------
        filename : str
            Saving plots to filename if provideed (optional)
        verbose : bool
            Plotting on screen if True (optional)
        """
        import pylab as plt

        if filename != None:
            plt.figure(1, figsize=(12, 7), dpi=300)
        else:
            plt.figure(1, figsize=(8, 6))

        plt.subplot(221)

        # For NEM: plot gas and liquid temperatures separately
        if self.non_equilibrium:
            plt.plot(self.time_array, self.T_gas - 273.15, "r", label="Gas")
            plt.plot(self.time_array, self.T_liquid - 273.15, "b", label="Liquid")
        else:
            plt.plot(self.time_array, self.T_fluid - 273.15, "b", label="Fluid")
        if "thermal_conductivity" not in self.input["vessel"].keys():
            plt.plot(
                self.time_array, self.T_vessel - 273.15, "g", label="Vessel wall dry"
            )
            if self.liquid_level.any() != 0:
                plt.plot(
                    self.time_array,
                    self.T_vessel_wetted - 273.15,
                    # marker="o",
                    color="darkorange",
                    label="Vessel wall wetted",
                )
        if "thermal_conductivity" in self.input["vessel"].keys():
            if "liner_thermal_conductivity" in self.input["vessel"].keys():
                plt.plot(
                    self.time_array,
                    self.T_bonded_wall - 273.15,
                    "g",
                    label="Liner/composite",
                )
                plt.plot(
                    self.time_array,
                    self.T_bonded_wall_wetted - 273.15,
                    color="darkorange",
                    label="Liner/composite wetted",
                )
            plt.plot(
                self.time_array, self.T_inner_wall - 273.15, "g--", label="Inner wall"
            )
            plt.plot(
                self.time_array, self.T_outer_wall - 273.15, "g-.", label="Outer wall"
            )

            plt.plot(
                self.time_array,
                self.T_inner_wall_wetted - 273.15,
                color="darkorange",
                linestyle="--",
                label="Inner wall wetted",
            )
            plt.plot(
                self.time_array,
                self.T_outer_wall_wetted - 273.15,
                color="darkorange",
                linestyle="-.",
                label="Outer wall wetted",
            )

        if self.input["valve"]["flow"] == "discharge":
            plt.plot(self.time_array, self.T_vent - 273.15, "r", label="Vent")
        if "validation" in self.input:
            if "temperature" in self.input["validation"]:
                temp = self.input["validation"]["temperature"]
                if "gas_mean" in temp:
                    plt.plot(
                        np.asarray(temp["gas_mean"]["time"]),
                        np.asarray(temp["gas_mean"]["temp"]) - 273.15,
                        "b.",
                        label="Gas mean",
                    )
                if "gas_high" in temp:
                    plt.plot(
                        np.asarray(temp["gas_high"]["time"]),
                        np.asarray(temp["gas_high"]["temp"]) - 273.15,
                        "b-.",
                        label="Gas high",
                    )
                if "gas_low" in temp:
                    plt.plot(
                        np.asarray(temp["gas_low"]["time"]),
                        np.asarray(temp["gas_low"]["temp"]) - 273.15,
                        "b--",
                        label="Gas low",
                    )
                if "wall_mean" in temp:
                    plt.plot(
                        np.asarray(temp["wall_mean"]["time"]),
                        np.asarray(temp["wall_mean"]["temp"]) - 273.15,
                        "go",
                        label="Wall mean",
                    )
                if "wall_high" in temp:
                    plt.plot(
                        np.asarray(temp["wall_high"]["time"]),
                        np.asarray(temp["wall_high"]["temp"]) - 273.15,
                        "g-.",
                        label="Wall high",
                    )
                if "wall_inner" in temp:
                    plt.plot(
                        np.asarray(temp["wall_inner"]["time"]),
                        np.asarray(temp["wall_inner"]["temp"]) - 273.15,
                        "g+",
                        label="Inner wall",
                    )
                if "wall_low" in temp:
                    plt.plot(
                        np.asarray(temp["wall_low"]["time"]),
                        np.asarray(temp["wall_low"]["temp"]) - 273.15,
                        "g-.",
                        label="Wall high",
                    )
                if "wall_outer" in temp:
                    plt.plot(
                        np.asarray(temp["wall_outer"]["time"]),
                        np.asarray(temp["wall_outer"]["temp"]) - 273.15,
                        "gx",
                        label="Outer wall",
                    )

        plt.legend(loc="best")
        plt.xlabel("Time (seconds)")
        plt.ylabel(r"Temperature ($^\circ$C)")

        plt.subplot(222)
        plt.plot(self.time_array, self.P / 1e5, "b")
        if "validation" in self.input:
            if "pressure" in self.input["validation"]:
                plt.plot(
                    np.asarray(self.input["validation"]["pressure"]["time"]),
                    self.input["validation"]["pressure"]["pres"],
                    "ko",
                    label="Experimental",
                )
            plt.legend(loc="best")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Pressure (bar)")

        plt.subplot(223)
        plt.plot(
            self.time_array, self.q_outer / 1000, "r", label="External heat flux (dry)"
        )
        plt.plot(
            self.time_array,
            self.q_inner / 1000,
            color="darkorange",
            label="Internal heat flux (dry)",
        )

        if self.liquid_level.any() != 0:
            plt.plot(
                self.time_array,
                self.q_outer_wetted / 1000,
                "r--",
                label="External heat flux (wetted)",
            )
            plt.plot(
                self.time_array,
                self.q_inner_wetted / 1000,
                color="darkorange",
                linestyle="--",
                label="Internal heat flux (wetted)",
            )

        plt.legend(loc="best")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Heat flux (kW/m$^2$)")

        plt.subplot(224)
        plt.plot(self.time_array, self.mass_rate, "b", label="Mass flow (kg/s)")
        plt.plot(self.time_array, self.liquid_level, "g", label="Liquid level (m)")
        plt.legend(loc="best")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Vent rate (kg/s) / Liquid level (m)")

        plt.tight_layout()
        if filename != None:
            plt.savefig(filename + "_main.png")

        if verbose:
            plt.show()
        return

    def plot_release(self, filename=None, verbose=True):
        """
        Plot the CO2 release results: release rate + tank pressure, and the atmospheric
        dry-ice / vapour split with cumulative dry-ice mass.

        Only valid for a run with a top-level ``release`` block. Figures are saved as PDF
        (repository convention). Uses the ORS brand palette (navy/red/amber/slate).

        Parameters
        ----------
        filename : str
            If provided, the figure is saved to ``<filename>_release.pdf``.
        verbose : bool
            Show the figure on screen if True.
        """
        if not self.has_release:
            raise ValueError("plot_release requires a run with a 'release' block")

        import pylab as plt

        navy, red, amber, slate = "#002D40", "#D61F39", "#E6A740", "#82979F"
        t = self.time_array

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

        # (1) release rate + tank pressure
        ax1.plot(t, self.release_rate, color=red, label="Release rate (kg/s)")
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Release mass rate (kg/s)", color=red)
        ax1.tick_params(axis="y", labelcolor=red)
        ax1b = ax1.twinx()
        ax1b.plot(t, self.P / 1e5, color=navy, label="Tank pressure (bar)")
        ax1b.axhline(self.release_model.P_TRIPLE / 1e5, color=slate, ls="--",
                     lw=1, label="CO2 triple point")
        ax1b.set_ylabel("Tank pressure (bar)", color=navy)
        ax1b.tick_params(axis="y", labelcolor=navy)
        ax1.set_title("Release rate and tank pressure")

        # (2) atmospheric dry-ice / vapour fractions + cumulative dry-ice mass
        ax2.plot(t, self.x_solid_atm, color=amber, label="Atm. dry-ice fraction")
        ax2.plot(t, self.solid_frac_throat, color=slate, ls=":", label="Throat dry-ice fraction")
        ax2.set_xlabel("Time (s)")
        ax2.set_ylabel("Dry-ice mass fraction (-)")
        ax2.set_ylim(0, 1)
        ax2b = ax2.twinx()
        ax2b.plot(t, self.m_dryice_cum, color=navy, label="Cumulative dry ice (kg)")
        ax2b.set_ylabel("Cumulative dry-ice mass (kg)", color=navy)
        ax2b.tick_params(axis="y", labelcolor=navy)
        ax2.legend(loc="upper left")
        ax2.set_title("Atmospheric dry-ice state")

        plt.tight_layout()
        if filename is not None:
            plt.savefig(filename + "_release.pdf")
        if verbose:
            plt.show()
        return

    def plot_envelope(self, filename=None, verbose=True):
        """
        Creating standard plots for the solved problem

        Parameters
        ----------
        filename : str
            Saving plots to filename if provideed (optional)
        verbose : bool
            Plotting on screen if True (optional)
        """
        import pylab as plt

        if filename != None:
            plt.figure(2, figsize=(12, 7), dpi=300)
        else:
            plt.figure(2, figsize=(8, 6))

        self.fluid.build_phase_envelope("None")
        PE = self.fluid.get_phase_envelope_data()

        plt.plot(PE.T, PE.p, "-", label="HEOS Phase Envelope", color="g")
        plt.plot(self.T_fluid, self.P, "-.", label="P/T fluid trajectory", color="b")
        plt.plot(self.T_fluid[0], self.P[0], "o", label="Start", color="b")
        plt.plot(self.T_fluid[-1], self.P[-1], ".", label="End", color="b")
        plt.xlabel("Temperature [K]")
        plt.ylabel("Pressure [Pa]")
        plt.legend(loc="best")
        plt.tight_layout()

        if filename != None:
            plt.savefig(filename + "_envelope.png")

        if verbose:
            plt.show()

    def plot_tprofile(self, filename=None, verbose=True):
        """
        Creating standard plots for the solved problem

        Parameters
        ----------
        filename : str
            Saving plots to filename if provideed (optional)
        verbose : bool
            Plotting on screen if True (optional)
        """

        # Add some checks if the profile has been constructed
        # return some
        import pylab as plt
        import numpy as np

        if filename != None:
            plt.figure(3, figsize=(8, 6))
        else:
            plt.figure(3, figsize=(8, 6))

        X, Y = np.meshgrid(self.z * 1e3, self.time_array[:-1])
        x0 = self.z[0] * 1e3
        x1 = self.z[-1] * 1e3
        y0 = self.time_array[0]
        y1 = self.time_array[-1]

        # plt.subplot(211)
        plt.contourf(Y, X, np.asarray(self.temp_profile), origin="lower", levels=20)
        # plt.imshow(np.asarray(self.temp_profile).T, aspect = 'auto', extent=(y0,y1,x0,x1), origin='lower')

        plt.colorbar(label="Temperature (K)")
        plt.xlabel("Time (s)")
        plt.ylabel("z (mm)")
        # add title with descriptive text for z axis

        if filename != None:
            plt.savefig(filename + "_tprofile1.png", dpi=300)

        if filename != None:
            plt.figure(4, figsize=(8, 6))
        else:
            plt.figure(4, figsize=(8, 6))
        if verbose:
            plt.show()
        n = math.floor(len(self.time_array) / 15)
        for i in range(len(self.time_array[::n])):
            plt.plot(
                self.temp_profile[::n][i],
                self.z * 1e3,
                label=f"t = {int(self.time_array[::n][i])} s.",
            )
        plt.legend(loc="best")
        plt.ylabel("z (mm)")
        plt.xlabel("Temperature (K)")
        plt.title("Temperature distribution")

        if filename != None:
            plt.savefig(filename + "_tprofile2.png", dpi=300)
        if verbose:
            plt.show()

    def analyze_rupture(self, filename=None):
        """
        Analyze vessel rupture potential under fire exposure conditions.

        Performs a simplified rupture analysis by calculating vessel wall temperatures
        under external fire heat load and comparing von Mises stress (from internal
        pressure) against temperature-dependent allowable tensile stress (ATS).
        Determines if and when vessel rupture may occur.

        The analysis:
        1. Calculates wall temperature evolution under fire exposure
        2. Computes von Mises equivalent stress from internal pressure
        3. Evaluates temperature-dependent material strength (ATS)
        4. Identifies rupture time when stress exceeds strength
        5. Generates diagnostic plots

        Parameters
        ----------
        filename : str, optional
            Base filename for saving plots. If None, plots are displayed on screen.
            Generates two plot files:
            - {filename}_peak_wall_temp.png
            - {filename}_ATS_vonmises.png

        Returns
        -------
        None
            Results are stored in instance variables:
            - self.rupture_time : float or None
                Time when rupture occurs [s], or None if no rupture predicted
            - self.peak_times : ndarray
                Time array for rupture analysis [s]
            - self.von_mises : ndarray
                von Mises equivalent stress history [Pa]
            - self.ATS_wetted : ndarray
                Allowable tensile stress for wetted region [Pa]
            - self.ATS_unwetted : ndarray
                Allowable tensile stress for unwetted region [Pa]
            - self.peak_T_wetted : ndarray
                Wall temperature history for wetted region [K]
            - self.peak_T_unwetted : ndarray
                Wall temperature history for unwetted region [K]

        Notes
        -----
        Requires rupture analysis parameters in input:
        - self.rupture_fire: Fire type (e.g., 'api_pool', 'scandpower_jet')
        - self.rupture_material: Material type for ATS calculation

        The method uses a simplified lumped capacitance model for wall heating.
        Time step for rupture analysis is fixed at 10 seconds.
        """
        from hyddown.materials import steel_Cp, ATS, von_mises
        from hyddown import fire

        pres = lambda x: np.interp(x, self.time_array, self.P)
        q_unwetted = lambda x: np.interp(x, self.time_array, self.q_inner)
        q_wetted = lambda x: np.interp(x, self.time_array, self.q_inner_wetted)

        T0_unwetted = self.T_vessel[0]
        T0_wetted = self.T_vessel_wetted[0]

        thk = self.thickness
        rho = self.vessel_density
        inner_diameter = self.diameter

        dt = 10
        max_time = self.time_array[-1]
        tsteps = int(max_time / dt)

        T_wetted_wall = np.zeros(tsteps + 1)
        T_unwetted_wall = np.zeros(tsteps + 1)
        T_wetted_wall[0] = T0_wetted
        T_unwetted_wall[0] = T0_unwetted
        peak_times = np.zeros(tsteps + 1)
        peak_times[0] = 0

        for i in range(tsteps):
            peak_times[i + 1] = peak_times[i] + dt
            q_fire_wetted = fire.sb_fire(T_wetted_wall[i], self.rupture_fire)
            q_fire_unwetted = fire.sb_fire(T_unwetted_wall[i], self.rupture_fire)
            T_wetted_wall[i + 1] = T_wetted_wall[i] + (
                q_fire_wetted - q_wetted(peak_times[i])
            ) * dt / (thk * rho * steel_Cp(T_wetted_wall[i], self.rupture_material))
            T_unwetted_wall[i + 1] = T_unwetted_wall[i] + (
                q_fire_unwetted - q_unwetted(peak_times[i])
            ) * dt / (thk * rho * steel_Cp(T_unwetted_wall[i], self.rupture_material))

        ATS_wetted = np.array([ATS(T, self.rupture_material, k_s=self.rupture_k_s) for T in T_wetted_wall])
        ATS_unwetted = np.array(
            [ATS(T, self.rupture_material, k_s=self.rupture_k_s) for T in T_unwetted_wall]
        )
        von_mises_wetted = von_mises_unwetted = np.array(
            [von_mises(pres(time), inner_diameter, thk) for time in peak_times]
        )

        self.peak_times = peak_times
        self.von_mises = von_mises_unwetted
        self.ATS_unwetted = ATS_unwetted
        self.ATS_wetted = ATS_wetted
        self.peak_T_wetted = T_wetted_wall
        self.peak_T_unwetted = T_unwetted_wall

        if np.all(ATS_unwetted > von_mises_unwetted) == True:
            self.rupture_time = None
            # print("No rupture")
        elif np.all(ATS_unwetted < von_mises_unwetted) == True:
            self.rupture_time = 0
            # print("Rupture at time=0")
        else:
            rupture_idx = np.where(ATS_unwetted < von_mises_unwetted)[0][0]
            self.rupture_time = (
                peak_times[rupture_idx - 1] + peak_times[rupture_idx]
            ) / 2
            # print("Rupture time +/- 5 s:", self.rupture_time)
            # print("Rupture pressure (bar)", pres(peak_times[rupture_idx - 1]))

        from matplotlib import pyplot as plt

        plt.figure()
        if np.any(self.liquid_level > 0):
            plt.plot(peak_times, T_wetted_wall - 273.15, label="T wetted wall")
        plt.plot(peak_times, T_unwetted_wall - 273.15, label="T unwetted wall")
        plt.xlabel("Time (s)")
        plt.ylabel("Wall temperature (C)")
        plt.legend(loc="best")
        if filename is not None:
            plt.savefig(filename + "_peak_wall_temp.pdf")

        plt.figure()
        plt.plot(
            peak_times,
            np.array([pres(time) for time in peak_times]) / 1e5,
            label="Pressure",
        )
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.legend(loc="best")
        if filename is not None:
            plt.savefig(filename + "_peak_pressure.pdf")

        plt.figure()
        plt.plot(peak_times, von_mises_wetted / 1e6, label="von Mises stress")

        plt.plot(peak_times, ATS_unwetted / 1e6, label="ATS unwetted wall")
        if np.any(self.liquid_level > 0):
            plt.plot(peak_times, ATS_wetted / 1e6, label="ATS wetted wall")
        plt.xlabel("Time (s)")
        plt.ylabel("Allowable Tensile Strength / von Mises Stress (MPa)")
        plt.legend(loc="best")
        if filename is not None:
            plt.savefig(filename + "_ATS_vonmises.pdf")

        if filename is None:
            plt.show()

    def __str__(self):
        return "HydDown vessel filling/depressurization class"

    def generate_report(self):
        """
        Generating a report summarising key features for the problem solved.
        Can be used for e.g. case studies, problem optimisation (external) etc.
        """
        report = {}

        report["start_time"] = self.time_array[0]
        report["end_time"] = self.time_array[-1]

        # Pressure
        report["max_pressure"] = max(self.P)
        report["time_max_pressure"] = self.time_array[np.argmax(self.P)]
        report["min_pressure"] = min(self.P)
        report["time_min_pressure"] = self.time_array[np.argmin(self.P)]

        # Temperatures
        report["max_fluid_temp"] = max(self.T_fluid)
        report["time_max_fluid_temp"] = self.time_array[np.argmax(self.T_fluid)]
        report["min_fluid_temp"] = min(self.T_fluid)
        report["time_min_fluid_temp"] = self.time_array[np.argmin(self.T_fluid)]

        report["max_wall_temp"] = max(self.T_vessel)
        report["time_max_wall_temp"] = self.time_array[np.argmax(self.T_vessel)]
        report["min_wall_temp"] = min(self.T_vessel)
        report["time_min_wall_temp"] = self.time_array[np.argmin(self.T_vessel)]

        report["max_inner_wall_temp"] = max(self.T_inner_wall)
        report["time_max_inner_wall_temp"] = self.time_array[
            np.argmax(self.T_inner_wall)
        ]
        report["min_inner_wall_temp"] = min(self.T_inner_wall)
        report["time_min_inner_wall_temp"] = self.time_array[
            np.argmin(self.T_inner_wall)
        ]

        report["max_outer_wall_temp"] = max(self.T_outer_wall)
        report["time_max_outer_wall_temp"] = self.time_array[
            np.argmax(self.T_outer_wall)
        ]
        report["min_outer_wall_temp"] = min(self.T_outer_wall)
        report["time_min_outer_wall_temp"] = self.time_array[
            np.argmin(self.T_outer_wall)
        ]

        # Mass flows and inventory
        report["max_mass_rate"] = max(self.mass_rate)
        report["time_max_mass_rate"] = self.time_array[self.mass_rate.argmax()]
        report["initial_mass"] = self.mass_fluid[0]
        report["final_mass"] = self.mass_fluid[-1]
        report["volume"] = self.vol

        # CO2 release / dry-ice summary
        if self.has_release:
            report["max_release_rate"] = max(self.release_rate)
            report["total_dryice_mass"] = self.m_dryice_cum[-1]
            report["max_atm_dryice_frac"] = max(self.x_solid_atm)
            report["max_throat_dryice_frac"] = max(self.solid_frac_throat)

        # Heat transfer (Q in W, q in W/m²)
        # Track both max and min to capture extreme values in both directions
        # (discharge: Q_inner > 0, filling: Q_inner < 0)
        report["max_Q_inside"] = max(self.Q_inner)
        report["time_max_Q_inside"] = self.time_array[np.argmax(self.Q_inner)]
        report["min_Q_inside"] = min(self.Q_inner)
        report["time_min_Q_inside"] = self.time_array[np.argmin(self.Q_inner)]

        report["max_Q_outside"] = max(self.Q_outer)
        report["time_max_Q_outside"] = self.time_array[np.argmax(self.Q_outer)]
        report["min_Q_outside"] = min(self.Q_outer)
        report["time_min_Q_outside"] = self.time_array[np.argmin(self.Q_outer)]

        # Heat flux per unit area (use q arrays which store W/m²)
        report["max_heat_flux_inside"] = max(self.q_inner)
        report["time_max_heat_flux_inside"] = self.time_array[np.argmax(self.q_inner)]
        report["min_heat_flux_inside"] = min(self.q_inner)
        report["time_min_heat_flux_inside"] = self.time_array[np.argmin(self.q_inner)]

        report["max_heat_flux_outside"] = max(self.q_outer)
        report["time_max_heat_flux_outside"] = self.time_array[np.argmax(self.q_outer)]
        report["min_heat_flux_outside"] = min(self.q_outer)
        report["time_min_heat_flux_outside"] = self.time_array[np.argmin(self.q_outer)]

        self.report = report
