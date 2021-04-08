# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from hyddown import transport as tp

class HydDown:
    def __init__(self, input):
        self.input = input
        self.read_input()
        self.initialize()

    def read_input(self):
        self.length =  self.input['vessel']['length']
        self.diameter =  self.input['vessel']['diameter']

        self.p0 =  self.input['initial']['pressure']
        self.T0 =  self.input['initial']['temperature']
        self.species = 'HEOS::'+ self.input['initial']['fluid'] 

        self.tstep =  self.input['calculation']['time_step']
        self.time_tot =  self.input['calculation']['end_time']
        self.method =  self.input['calculation']['type'] 
        if self.method == "energybalance": self.eta =  self.input['calculation']['eta'] 

        # Reading valve specific data
        if  self.input['valve']['type'] == 'orifice' or  self.input['valve']['type'] == 'psv':
            self.p_back =  self.input['valve']['back_pressure']
            self.D_orifice =  self.input['valve']['diameter']
            self.CD =  self.input['valve']['discharge_coef']
            if  self.input['valve']['type'] == 'psv':
                self.Pset =  self.input['valve']['set_pressure']
                self.blowdown =  self.input['valve']['blowdown']
                self.psv_state = 'closed'
        elif  self.input['valve']['type'] == "controlvalve":
            self.p_back =  self.input['valve']['back_pressure']
            self.Cv =  self.input['valve']['Cv']
            if 'xT' in  self.input['valve']:
                self.xT =  self.input['valve']['xT']
            if 'Fp' in  self.input['valve']:
                self.Fp =  self.input['valve']['Fp']

        # valve type    
        # - constant_mass
        # - functional mass flow
        thickness = 0
        # Reading heat transfer related data/information
        if 'heat_transfer' in  self.input:
            self.heat_method =  self.input['heat_transfer']['type']
            if self.heat_method == "specified_h" or self.heat_method == "specified_U":
                self.Tamb =  self.input['heat_transfer']['temp_ambient']
            if self.heat_method == "specified_U": self.Ufix =  self.input['heat_transfer']['U_fix']
            if self.heat_method == "specified_Q": self.Qfix =  self.input['heat_transfer']['Q_fix']
            if self.heat_method == "specified_h":
                self.vessel_cp =  self.input['vessel']['heat_capacity']
                self.vessel_density =  self.input['vessel']['density']
                self.vessel_orientation =  self.input['vessel']['orientation']
                self.thickness =  self.input['vessel']['thickness']
                self.h_out =  self.input['heat_transfer']['h_outer']
                self.h_in =  self.input['heat_transfer']['h_inner']

    def initialize(self):
        self.vol = self.diameter**2/4 * math.pi * self.length  # m3
        self.vol_tot = (self.diameter + 2 * self.thickness)**2/4 * math.pi * (self.length + 2 * self.thickness)  # m3
        self.vol_solid = self.vol_tot - self.vol
        self.surf_area_outer = (self.diameter + 2 * self.thickness)**2/4 * math.pi * 2 + (self.diameter + 2 * self.thickness) * math.pi * (self.length + 2 * self.thickness)
        self.surf_area_inner = (self.diameter)**2/4 * math.pi * 2 + (self.diameter) * math.pi * self.length

        # data storage
        data_len = int(self.time_tot / self.tstep)
        self.rho = np.zeros(data_len)
        self.T_fluid = np.zeros(data_len)
        self.T_vessel = np.zeros(data_len)
        self.Q_outer = np.zeros(data_len)
        self.Q_inner = np.zeros(data_len)
        self.h_inside = np.zeros(data_len)
        self.T_vent = np.zeros(data_len)
        self.H_mass = np.zeros(data_len)
        self.S_mass = np.zeros(data_len)
        self.U_mass = np.zeros(data_len)
        self.U_tot = np.zeros(data_len)
        self.U_iter = np.zeros(data_len)
        self.m_iter = np.zeros(data_len)
        self.n_iter = np.zeros(data_len)
        self.Qo = np.zeros(data_len)
        self.Qi = np.zeros(data_len)
        self.P = np.zeros(data_len)
        self.mass_fluid = np.zeros(data_len)
        self.mass_rate = np.zeros(data_len)
        self.time_array = np.zeros(data_len)

        self.rho0 = PropsSI('D', 'T', self.T0, 'P', self.p0, self.species)
        self.m0 = self.rho0 * self.vol

    def run(self):
        # Inititialise
        input=self.input
        self.rho[0] = self.rho0
        self.T_fluid[0] = self.T0
        self.T_vessel[0] = self.T0
        self.H_mass[0] = PropsSI('H', 'T', self.T0, 'P', self.p0, self.species)
        self.S_mass[0] = PropsSI('S', 'T', self.T0, 'P', self.p0, self.species)
        self.U_mass[0] = PropsSI('U', 'T', self.T0, 'P', self.p0, self.species)
        self.U_tot[0] = PropsSI('U', 'T', self.T0, 'P', self.p0, self.species) * self.m0
        self.P[0] = self.p0
        self.mass_fluid[0] = self.m0
        cpcv = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p0, self.species) / PropsSI('CVMOLAR', 'T', self.T0, 'P', self.p0, self.species)

        if input['valve']['type'] == 'orifice':
            if input['valve']['flow'] == 'filling':
                k = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p_back, self.species) / PropsSI('CVMOLAR', 'T', self.T0, 'P', self.p_back, self.species)
                self.mass_rate[0] = -tp.gas_release_rate(self.p_back, self.p0, PropsSI('D', 'T', self.T0, 'P', self.p_back, self.species), k, self.CD, self.D_orifice**2/4 * math.pi)
            else:
                self.mass_rate[0] = tp.gas_release_rate(self.p0, self.p_back, self.rho0, cpcv, self.CD, self.D_orifice**2/4 * math.pi)
        elif input['valve']['type'] == 'mdot':
            if input['valve']['flow'] == 'filling':
                self.mass_rate[0] = -input['valve']['mass_flow']
            else:
                self.mass_rate[0] = input['valve']['mass_flow']
        elif input['valve']['type'] == 'controlvalve':
            if input['valve']['flow'] == 'filling':
                Z = PropsSI('Z', 'T', self.T0, 'P', self.p_back, self.species)
                MW = PropsSI('M', 'T', self.T0, 'P', self.p_back, self.species)
                k = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p_back, self.species) / PropsSI('CVMOLAR', 'T', self.T0, 'P', self.p_back, self.species)
                self.mass_rate[0] = -tp.control_valve(self.p_back, self.p0, self.T0, Z, MW, k, self.Cv)
            else:
                Z = PropsSI('Z', 'T', self.T0, 'P', self.p0, self.species)
                MW = PropsSI('M', 'T', self.T0, 'P', self.p0, self.species)
                k = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p0, self.species) / PropsSI('CVMOLAR', 'T',self.T0, 'P', self.p0, self.species)
                self.mass_rate[0] = tp.control_valve(self.p0, self.p_back, self.T0, Z, MW, k, self.Cv)
        elif input['valve']['type'] == 'psv':
            if input['valve']['flow'] == 'filling': raise ValueError("Unsupported valve: ", input['valve']['type'], " for vessel filling.")
            self.mass_rate[0] = tp.relief_valve(self.p0, self.p_back, self.Pset, self.blowdown, self.rho0, cpcv, self.CD, self.D_orifice**2/4 * math.pi)


        self.time_array[0] = 0
        # Run actual integration
        for i in range(1, len(self.time_array)):
            self.time_array[i] = self.time_array[i-1] + self.tstep
            self.mass_fluid[i] = self.mass_fluid[i-1] - self.mass_rate[i-1] * self.tstep
            self.rho[i] = self.mass_fluid[i] / self.vol

            if self.method == "isenthalpic":
                self.T_fluid[i] = PropsSI('T', 'D', self.rho[i], 'H', self.H_mass[i-1], self.species)
                self.P[i] = PropsSI('P', 'D', self.rho[i], 'H', self.H_mass[i-1], self.species)
            elif self.method == "isentropic":
                self.T_fluid[i] = PropsSI('T', 'D', self.rho[i], 'S', self.S_mass[i-1], self.species)
                self.P[i] = PropsSI('P', 'D', self.rho[i], 'S', self.S_mass[i-1], self.species)
            elif self.method == "isothermal":
                self.T_fluid[i] = self.T0
                self.P[i] = PropsSI('P', 'D', self.rho[i], 'T', self.T0, self.species)
            elif self.method == "constantU":
                self.T_fluid[i] = PropsSI('T', 'D', self.rho[i], 'U', self.U_mass[i-1], self.species)
                self.P[i] = PropsSI('P', 'D', self.rho[i], 'U', self.U_mass[i-1], self.species)
            elif self.method == "energybalance":
                P1 = PropsSI('P', 'D', self.rho[i], 'T', self.T_fluid[i-1], self.species)
                T1 = PropsSI('T', 'P', P1, 'H', self.H_mass[i-1], self.species)
        
                if self.heat_method == "specified_h" or self.heat_method == "detailed":
                    if self.h_in == "calc":
                        hi = tp.h_inner(self.length, self.T_fluid[i-1], self.T_vessel[i-1], self.P[i-1], self.species)
                    else:
                        hi = self.h_in
                    self.h_inside[i] = hi
                    self.Q_inner[i] = self.surf_area_inner * hi * (self.T_vessel[i-1] - self.T_fluid[i-1])
                    self.Q_outer[i] = self.surf_area_outer * self.h_out * (self.Tamb - self.T_vessel[i-1])
                    self.T_vessel[i] = self.T_vessel[i-1] + (self.Q_outer[i] - self.Q_inner[i]) * self.tstep / (self.vessel_cp * self.vessel_density * self.vol_solid)
                elif self.heat_method == "specified_U":
                    self.Q_inner[i] = self.surf_area_outer * self.Ufix * (self.Tamb - self.T_fluid[i-1])
                    self.T_vessel[i] = self.T_vessel[0] 
                elif self.heat_method == "specified_Q":
                    self.Q_inner[i] = self.Qfix
                    self.T_vessel[i] = self.T_vessel[0]
                else:
                    self.Q_inner[i] = 0.0
                    self.T_vessel[i] = self.T_vessel[0]

                NMOL = self.mass_fluid[i-1] / PropsSI('M', self.species) 
                NMOL_ADD = (self.mass_fluid[i]-self.mass_fluid[i-1]) / PropsSI('M', self.species) 
                if input['valve']['flow'] == 'filling':
                    T_used = self.T0
                    P_used = self.p_back
                else:
                    T_used = self.T_fluid[i-1]
                    P_used = self.P[i-1]
                U_start = NMOL_ADD * PropsSI('HMOLAR', 'P', P_used, 'T', T_used, self.species) + NMOL * PropsSI('HMOLAR', 'P', self.P[i-1], 'T', self.T_fluid[i-1], self.species) - self.eta *  self.P[i-1] * self.vol + self.Q_inner[i] * self.tstep
                NMOL=NMOL+NMOL_ADD

                U = 0
                nn = 0
                rho1 = 0
                itermax = 1000
                m = 0
                n = 0
                relax = 0.1

                while abs(self.rho[i] - rho1) > 0.01 and m < itermax:
                    m = m + 1
                    rho1 = PropsSI('D', 'T', T1,'P', P1, self.species)
                    dd = self.rho[i] - rho1  #NMOL-nn
                    P1 = P1 + dd * 1e5
                    if m == itermax:
                        raise Exception("Iter max exceeded for rho/P")
                    while abs(U_start - U) / U_start > 0.00001 and n < itermax:
                        n = n + 1
                        U = NMOL * PropsSI('HMOLAR', 'P', P1, 'T', T1, self.species) - self.eta * P1 * self.vol  #Q_inner[i]*tstep
                        d = U_start - U 
                        T1 = T1 + 0.1 * d / U_start * T1
                        if n == itermax:
                            raise Exception("Iter max exceeded for U/T")
            
                self.m_iter[i] = dd
                self.n_iter[i] = d
                self.U_iter[i] = U / NMOL * self.mass_fluid[i]
                self.P[i] = P1
                self.T_fluid[i] = T1
        
            else:
                raise NameError("Unknown calculation method: " + self.method)

            self.H_mass[i] = PropsSI('H', 'T', self.T_fluid[i], 'P', self.P[i], self.species)
            self.S_mass[i] = PropsSI('S', 'T', self.T_fluid[i], 'P', self.P[i], self.species)
            self.U_mass[i] = (self.mass_fluid[i] * PropsSI('H', 'P', self.P[i], 'T', self.T_fluid[i], self.species) - self.P[i] * self.vol) / self.mass_fluid[i]  
            cpcv = PropsSI('CP0MOLAR', 'T', self.T_fluid[i], 'P', self.P[i], self.species) / PropsSI('CVMOLAR', 'T', self.T_fluid[i], 'P', self.P[i], self.species)
    
            if input['valve']['type'] == 'orifice':
                if input['valve']['flow'] == 'filling':
                    k = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p_back, self.species) / PropsSI('CVMOLAR', 'T', self.T0, 'P', self.p_back,self.species)
                    self.mass_rate[i] = -tp.gas_release_rate(self.p_back, self.P[i], PropsSI('D', 'T', self.T0, 'P', self.p_back, self.species), k, self.CD, self.D_orifice**2/4 * math.pi)
                else:
                    self.mass_rate[i] = tp.gas_release_rate(self.P[i], self.p_back, self.rho[i], cpcv, self.CD, self.D_orifice**2/4 * math.pi)
            elif input['valve']['type'] == 'controlvalve':
                if input['valve']['flow'] == 'filling':
                    Z = PropsSI('Z', 'T', self.T0, 'P', self.p_back, self.species)
                    MW = PropsSI('M', 'T', self.T0, 'P', self.p_back, self.species)
                    k = PropsSI('CP0MOLAR', 'T', self.T0, 'P', self.p_back, self.species) / PropsSI('CVMOLAR', 'T', self.T0, 'P', self.p_back, self.species)
                    self.mass_rate[i] = -tp.control_valve(self.p_back, self.P[i], self.T0, Z, MW, k, self.Cv)
                else:
                    Z = PropsSI('Z', 'T', self.T_fluid[i], 'P', self.P[i], self.species)
                    MW = PropsSI('M','T', self.T_fluid[i], 'P', self.P[i], self.species)
                    self.mass_rate[i] = tp.control_valve(self.P[i], self.p_back, self.T_fluid[i], Z, MW, cpcv, self.Cv)
            elif input['valve']['type'] == 'mdot':
                if input['valve']['flow'] == 'filling':
                    self.mass_rate[i] = -input['valve']['mass_flow']
                else:
                    self.mass_rate[i] = input['valve']['mass_flow']
            elif input['valve']['type'] == 'psv':
                self.mass_rate[i] = tp.relief_valve(self.P[i], self.p_back, self.Pset, self.blowdown, self.rho[i], cpcv, self.CD, self.D_orifice**2/4 * math.pi)

    def plot(self):
        import pylab as plt 

        plt.figure()
        plt.subplot(221)
        plt.plot(self.time_array/60, self.T_fluid-273.15, 'b', label="Fluid")
        plt.plot(self.time_array/60, self.T_vessel-273.15, 'g', label="Vessel")
        if 'validation' in self.input:
            if 'temperature' in self.input['validation']:
                temp = self.input['validation']['temperature']
                if 'gas_mean' in temp:
                    plt.plot(np.asarray(temp['gas_mean']['time']) / 60, np.asarray(temp['gas_mean']['temp']) - 273.15, 'b:', label="Gas mean")
                if 'gas_high' in temp:
                    plt.plot(np.asarray(temp['gas_high']['time']) / 60, np.asarray(temp['gas_high']['temp']) - 273.15, 'b-.', label="Gas high")
                if 'gas_low' in temp:
                    plt.plot(np.asarray(temp['gas_low']['time']) / 60, np.asarray(temp['gas_low']['temp']) - 273.15, 'b--', label="Gas low")
                if 'wall_mean' in temp:
                    plt.plot(np.asarray(temp['wall_mean']['time']) / 60, np.asarray(temp['wall_mean']['temp']) - 273.15, 'g:', label="Wall mean")
                if 'wall_high' in temp:
                    plt.plot(np.asarray(temp['wall_high']['time']) / 60, np.asarray(temp['wall_high']['temp']) - 273.15, 'g-.', label="Wall high")
                if 'wall_low' in temp:
                    plt.plot(np.asarray(temp['wall_low']['time']) / 60, np.asarray(temp['wall_low']['temp']) - 273.15, 'g--', label="Wall low")
        plt.legend(loc='best')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Temperature ($^\circ$C)')

        plt.subplot(222)
        plt.plot(self.time_array / 60, self.P / 1e5, 'b', label="Calculated")
        if 'validation' in self.input:
            if 'pressure' in self.input['validation']:
                plt.plot(np.asarray(self.input['validation']['pressure']['time']) / 60, self.input['validation']['pressure']['pres'], 'ko', label="Experimental")
        plt.legend(loc='best')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Pressure (bar)')

        plt.subplot(223)
        plt.plot(self.time_array / 60, self.H_mass, 'b', label='H (J/kg)')
        plt.plot(self.time_array / 60, self.U_mass, 'g', label='U (J/kg)')
        plt.plot(self.time_array / 60, self.S_mass * 100, 'r', label='S*100 (J/kg K)')
        plt.legend(loc='best')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Enthalpy/Internal Energy/Entropy')

        plt.subplot(224)
        plt.plot(self.time_array / 60, self.mass_rate, 'b', label='m_dot')
        plt.xlabel('Time (minutes)')
        plt.ylabel('Vent rate (kg/s)')
        plt.show()

    def report(self):
        pass
