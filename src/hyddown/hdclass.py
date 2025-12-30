# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

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

        if "thickness" in self.input["vessel"]:
            self.outer_vol = self.inner_vol.add_thickness(
                self.input["vessel"]["thickness"]
            )
        else:
            self.outer_vol = self.inner_vol.add_thickness(0.0)

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

        # Reading valve specific data
        if (
            self.input["valve"]["type"] == "orifice"
            or self.input["valve"]["type"] == "psv"
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

        # valve type
        # - constant_mass
        # - functional mass flow
        self.thickness = 0

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
                if self.input["valve"]["flow"] == "filling":
                    raise ValueError("Filling and Fire heat load not implemented")

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

    def calc_liquid_level(self):
        """
        Calculating liquid level based on current fluid state

        Parameters
        ----------
        fluid : CoolProp AbstractState
            Current fluid state
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

    def PHres(self, T, P, H):
        """
        Residual enthalpy function to be minimised during a PH-problem

        Parameters
        ----------
        H : float
            Enthalpy at initial/final conditions
        P : float
            Pressure at final conditions.
        T : float
            Updated estimate for the final temperature at P,H
        """
        self.vent_fluid.update(CP.PT_INPUTS, P, T)
        return ((H - self.vent_fluid.hmass()) / H) ** 2

    def PHres_relief(self, T, P, H):
        """
        Residual enthalpy function to be minimised during a PH-problem

        Parameters
        ----------
        H : float
            Enthalpy at initial/final conditions
        P : float
            Pressure at final conditions.
        T : float
            Updated estimate for the final temperature at P,H
        """
        self.fluid.update(CP.PT_INPUTS, P, T)
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
        self.T_vessel[0] = self.T0
        self.T_inner_wall[0] = self.T0
        self.T_outer_wall[0] = self.T0
        self.T_vessel_wetted[0] = self.T0
        self.T_inner_wall_wetted[0] = self.T0
        self.T_outer_wall_wetted[0] = self.T0
        self.T_bonded_wall[0] = self.T0
        self.T_bonded_wall_wetted[0] = self.T0
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

        self.time_array[0] = 0

        # setting heat transfer parameters
        T_profile, T_profile2 = 0, 0
        relief_area = []
        for i in tqdm(
            range(1, len(self.time_array)),
            desc="hyddown",
            disable=disable_pbar,
            total=len(self.time_array),
        ):
            self.time_array[i] = self.time_array[i - 1] + self.tstep
            self.mass_fluid[i] = (
                self.mass_fluid[i - 1] - self.mass_rate[i - 1] * self.tstep
            )

            self.rho[i] = self.mass_fluid[i] / self.vol

            # For the simple methods key variable is rho, with either H,S,U or T constant
            # updating T and p is convenient using coolprop.
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

            # This is the **only** complicated part where we have the energy balance included
            elif self.method == "energybalance":
                # Finding the heat supplied/removed by heat transfer from/to outside
                if self.heat_method == "specified_h" or self.heat_method == "detailed":
                    if self.h_in == "calc":
                        if self.vessel_orientation == "horizontal":
                            L = self.diameter
                        else:
                            L = self.length
                        if input["valve"]["flow"] == "filling":
                            # T_film = (self.T_fluid[i - 1] + self.T_vessel[i - 1]) / 2
                            T_film = (
                                self.T_fluid[i - 1] + self.T_inner_wall[i - 1]
                            ) / 2
                            self.transport_fluid.update(
                                CP.PT_INPUTS, self.P[i - 1], T_film
                            )

                            hi = tp.h_inside_mixed(
                                L,
                                # self.T_vessel[i - 1],
                                self.T_inner_wall[i - i],
                                self.T_fluid[i - 1],
                                self.transport_fluid,
                                self.mass_rate[i - 1],
                                self.diameter,
                            )
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
                            T_film = (
                                self.T_fluid[i - 1] + self.T_inner_wall[i - 1]
                            ) / 2
                            try:
                                self.transport_fluid.update(
                                    CP.PT_INPUTS, self.P[i - 1], T_film
                                )
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
                                self.T_fluid[i - 1],
                                self.transport_fluid,
                            )
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

                    self.h_inside[i] = hi
                    self.h_inside_wetted[i] = hiw

                    wetted_area = self.inner_vol.SA_from_h(self.liquid_level[i - 1])
                    if np.isnan(wetted_area):
                        wetted_area = 0

                    self.Q_inner[i] = (
                        (self.surf_area_inner - wetted_area)
                        * hi
                        * (self.T_inner_wall[i - 1] - self.T_fluid[i - 1])
                    )
                    self.q_inner[i] = hi * (
                        self.T_inner_wall[i - 1] - self.T_fluid[i - 1]
                    )

                    self.Q_inner_wetted[i] = (
                        wetted_area
                        * hiw
                        * (self.T_inner_wall_wetted[i - 1] - self.T_fluid[i - 1])
                    )
                    self.q_inner_wetted[i] = self.Q_inner_wetted[i] / wetted_area

                    if np.isnan(self.Q_inner_wetted[i]):
                        self.Q_inner_wetted[i] = 0

                    self.Q_outer[i] = (
                        (self.surf_area_inner - wetted_area)
                        * self.surf_area_outer
                        / self.surf_area_inner
                        * self.h_out
                        * (self.Tamb - self.T_outer_wall[i - 1])
                    )
                    self.q_outer[i] = self.h_out * (
                        self.Tamb - self.T_outer_wall[i - 1]
                    )

                    self.Q_outer_wetted[i] = (
                        wetted_area
                        * self.surf_area_outer
                        / self.surf_area_inner
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

                    if self.Q_outer_wetted[i] == 0 and self.Q_inner_wetted[i] == 0:
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

                    if "thermal_conductivity" in self.input["vessel"].keys():
                        theta = 0.5  # Crank-Nicholson scheme
                        dt = self.tstep / 10
                        k, rho, cp = (
                            self.input["vessel"]["thermal_conductivity"],
                            self.vessel_density,
                            self.vessel_cp,
                        )
                        if (
                            "liner_thermal_conductivity"
                            not in self.input["vessel"].keys()
                        ):
                            nn = 11  # number of nodes
                            z = np.linspace(0, self.thickness, nn)
                            self.z = z
                            mesh = tm.Mesh(
                                z, tm.LinearElement
                            )  # Or `QuadraticElement` to
                            mesh_w = tm.Mesh(
                                z, tm.LinearElement
                            )  # Or `QuadraticElement` to
                            cpeek = tm.isothermal_model(k, rho, cp)
                            cpeek_w = tm.isothermal_model(k, rho, cp)

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
                                t_bonded, T_profile = tm.solve_ht(domain, solver)
                                t_bonded_w, T_profile_w = tm.solve_ht(
                                    domain_w, solver_w
                                )
                            else:
                                bc = [
                                    {
                                        "q": self.Q_outer[i]
                                        / (self.surf_area_inner - wetted_area)
                                    },
                                    {
                                        "q": -self.Q_inner[i]
                                        / (self.surf_area_inner - wetted_area)
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
                                bc_w = [
                                    {"q": self.Q_outer_wetted[i] / wetted_area},
                                    {"q": -self.Q_inner_wetted[i] / wetted_area},
                                ]
                                domain_w = tm.Domain(mesh_w, [cpeek_w], bc_w)
                                domain_w.set_T(T_profile_w[-1, :])
                                solver = {
                                    "dt": dt,
                                    "t_end": self.tstep,
                                    "theta": theta,
                                }
                                t_bonded_w, T_profile_w = tm.solve_ht(
                                    domain_w, solver_w
                                )
                            solver = {"dt": dt, "t_end": self.tstep, "theta": theta}
                            t, T_profile = tm.solve_ht(domain, solver)
                            solver_w = {"dt": dt, "t_end": self.tstep, "theta": theta}
                            t_w, T_profile_w = tm.solve_ht(domain_w, solver_w)

                            self.temp_profile.append(T_profile[-1, :])
                            self.T_outer_wall[i] = T_profile[-1, 0]
                            self.T_inner_wall[i] = T_profile[-1, -1]
                            self.T_outer_wall_wetted[i] = T_profile_w[-1, 0]
                            self.T_inner_wall_wetted[i] = T_profile_w[-1, -1]
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
                                bc = [
                                    {
                                        "q": -self.Q_inner[i]
                                        / (self.surf_area_inner - wetted_area)
                                    },
                                    {
                                        "q": self.Q_outer[i]
                                        / (self.surf_area_outer - wetted_area)
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
                                    {"q": -self.Q_inner_wetted[i] / (wetted_area)},
                                    {"q": self.Q_outer_wetted[i] / wetted_area},
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
                        self.T_inner_wall[i] = self.T_vessel[i]
                        self.T_outer_wall[i] = self.T_vessel[i]
                        self.T_inner_wall_wetted[i] = self.T_vessel_wetted[i]
                        self.T_outer_wall_wetted[i] = self.T_vessel_wetted[i]

                elif self.heat_method == "s-b":
                    if self.vessel_orientation == "horizontal":
                        L = self.diameter
                    else:
                        L = self.length

                    wetted_area = self.inner_vol.SA_from_h(self.liquid_level[i - 1])
                    hi = tp.h_inner(
                        L,
                        self.T_fluid[i - 1],
                        self.T_vessel[i - 1],
                        self.P[i - 1],
                        self.species,
                    )
                    self.h_inside[i] = hi
                    if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                        self.transport_fluid_wet.update(
                            CP.PT_INPUTS, self.P[i - 1], self.T_fluid[i - 1]
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

                    self.Q_inner[i] = (
                        (self.surf_area_inner - wetted_area)
                        * hi
                        * (self.T_vessel[i - 1] - self.T_fluid[i - 1])
                    )

                    self.q_inner[i] = hi * (self.T_vessel[i - 1] - self.T_fluid[i - 1])
                    self.Q_inner_wetted[i] = (
                        wetted_area
                        * hiw
                        * (self.T_inner_wall_wetted[i - 1] - self.T_fluid[i - 1])
                    )
                    self.q_inner_wetted[i] = hiw * (
                        self.T_inner_wall_wetted[i - 1] - self.T_fluid[i - 1]
                    )
                    if np.isnan(self.Q_inner_wetted[i]):
                        self.Q_inner_wetted[i] = 0

                    self.Q_outer[i] = (
                        fire.sb_fire(self.T_vessel[i - 1], self.fire_type)
                        * (self.surf_area_inner - wetted_area)
                        * self.surf_area_outer
                        / self.surf_area_inner
                    )

                    self.q_outer[i] = fire.sb_fire(self.T_vessel[i - 1], self.fire_type)

                    self.Q_outer_wetted[i] = (
                        fire.sb_fire(self.T_vessel_wetted[i - 1], self.fire_type)
                        * wetted_area
                        * self.surf_area_outer
                        / self.surf_area_inner
                    )
                    self.q_outer_wetted[i] = fire.sb_fire(
                        self.T_vessel_wetted[i - 1], self.fire_type
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
                    if self.liquid_level[i - 1] > 0:
                        self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1] + (
                            self.Q_outer_wetted[i] - self.Q_inner_wetted[i]
                        ) * self.tstep / (
                            self.vessel_cp
                            * self.vessel_density
                            * self.vol_solid
                            * wetted_area
                            / self.inner_vol.A
                        )
                    else:
                        # Hack to heat up previous liquid wetted surface
                        self.T_vessel_wetted[i] = self.T_vessel_wetted[i - 1] + (
                            fire.sb_fire(self.T_vessel_wetted[i - 1], self.fire_type)
                            * (self.surf_area_inner - wetted_area)
                            * self.surf_area_outer
                            / self.surf_area_inner
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

                # Not pretty if-statement and a hack for fire relief area estimation. Most cases go directly to the first ...else... clause
                if input["valve"]["type"] == "relief":
                    if self.Pset <= self.P[i - 1]:
                        T1 = self.PHproblem(
                            h_in + self.tstep * self.Q_inner[i] / self.mass_fluid[i],
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

            if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
                cpcv = self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) / (
                    self.fluid.saturated_vapor_keyed_output(CP.iCpmolar) - 8.314
                )
                Z = self.fluid.saturated_vapor_keyed_output(CP.iZ)
            else:
                cpcv = self.fluid.cp0molar() / (self.fluid.cp0molar() - 8.314)
                Z = self.fluid.compressibility_factor()
            # Finally updating the mass rate for the mass balance in the next time step
            # Already done of the valve is "relief" (estimation)
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
                    if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
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
                    self.mass_rate[i] = tp.control_valve(
                        self.P[i], self.p_back, self.T_fluid[i], Z, MW, cpcv, Cv
                    )
            elif input["valve"]["type"] == "psv":
                self.mass_rate[i] = tp.relief_valve(
                    self.P[i],
                    self.p_back,
                    self.Pset,
                    self.blowdown,
                    cpcv,
                    self.CD,
                    self.T_fluid[i],
                    Z,
                    self.MW,
                    self.D_orifice**2 / 4 * math.pi,
                )
            if (
                "end_pressure" in self.input["valve"]
                and self.P[i] > self.input["valve"]["end_pressure"]
            ):
                massflow_stop_switch = 1
            if massflow_stop_switch:
                self.mass_rate[i] = 0
        self.isrun = True

        if input["valve"]["type"] == "relief":
            idx_max = self.mass_rate.argmax()
            self.mass_rate[idx_max] = (
                self.mass_rate[idx_max - 1] + self.mass_rate[idx_max + 1]
            ) / 2
            self.relief_area[idx_max] = (
                self.relief_area[idx_max - 1] + self.relief_area[idx_max + 1]
            ) / 2
            # print("Relief area:", 2*math.sqrt(max(relief_area[1:])/math.pi), max(self.mass_rate))

    def get_dataframe(self):
        """
        Storing relevant results in pandas dataframe for e.g. export
        to csv, excel for archiving or post analysis
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
        plt.ylabel("Temperature ($^\circ$C)")

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
            q_fire_wetted = fire.sb_fire(T_wetted_wall[i], self.sb_peak_fire_type)
            q_fire_unwetted = fire.sb_fire(T_unwetted_wall[i], self.sb_peak_fire_type)
            T_wetted_wall[i + 1] = T_wetted_wall[i] + (
                q_fire_wetted - q_wetted(peak_times[i])
            ) * dt / (thk * rho * steel_Cp(T_wetted_wall[i], self.material))
            T_unwetted_wall[i + 1] = T_unwetted_wall[i] + (
                q_fire_unwetted - q_unwetted(peak_times[i])
            ) * dt / (thk * rho * steel_Cp(T_unwetted_wall[i], self.material))

        ATS_wetted = np.array([ATS(T, self.material) for T in T_wetted_wall])
        ATS_unwetted = np.array([ATS(T, self.material) for T in T_unwetted_wall])
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

        plt.figure(1)
        plt.plot(peak_times, von_mises_wetted / 1e6, label="von Mises stress")

        plt.plot(peak_times, ATS_unwetted / 1e6, label="ATS unwetted wall")
        if sum(self.liquid_dyn_level) > 0:
            plt.plot(peak_times, ATS_wetted / 1e6, label="ATS wetted wall")
        plt.xlabel("Time (s)")
        plt.ylabel("Allowable Tensile Strength / von Mises Stress (MPa)")
        plt.legend(loc="best")
        if filename is not None:
            plt.savefig(
                filename + "_ATS_vonmises.png",
            )
        plt.figure(2)
        if sum(self.liquid_dyn_level) > 0:
            plt.plot(peak_times, T_wetted_wall - 273.15, label="T wetted wall")
        plt.plot(peak_times, T_unwetted_wall - 273.15, label="T unwetted wall")
        plt.xlabel("Time (s)")
        plt.ylabel("Wall temperature (C)")
        plt.legend(loc="best")
        if filename is not None:
            plt.savefig(
                filename + "_peak_wall_temp.png",
            )

        plt.figure(3)
        plt.plot(
            peak_times,
            np.array([pres(time) for time in peak_times]) / 1e5,
            label="Pressure",
        )
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.legend(loc="best")

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

        report["max_Q_inside"] = max(self.Q_inner)
        report["time_max_Q_inside"] = self.time_array[np.argmax(self.Q_inner)]
        report["max_heat_flux_inside"] = max(self.Q_inner / self.surf_area_inner)

        report["max_Q_outside"] = max(self.Q_outer)
        report["time_max_Q_outside"] = self.time_array[np.argmax(self.Q_outer)]
        report["max_heat_flux_outside"] = max(self.Q_outer / self.surf_area_outer)

        self.report = report
