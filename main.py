# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import math
import sys
import yaml
import numpy as np
from CoolProp.CoolProp import PropsSI
from hyddown import transport as tp


class HydDown:
    def __init__(self, input):
        self.input = input
        self.initialize()

    def initialize(self):
        length = input['vessel']['length']
        diameter = input['vessel']['diameter']

        p0 = input['initial']['pressure']
        T0 = input['initial']['temperature']
        species = 'HEOS::'+input['initial']['fluid'] 

        tstep = input['calculation']['time_step']
        time_tot = input['calculation']['end_time']
        method = input['calculation']['type'] 
        if method == "energybalance": eta = input['calculation']['eta'] 

        # Reading valve specific data
        if input['valve']['type'] == 'orifice' or input['valve']['type'] == 'psv':
            p_back = input['valve']['back_pressure']
            D_orifice = input['valve']['diameter']
            CD = input['valve']['discharge_coef']
            if input['valve']['type'] == 'psv':
                Pset = input['valve']['set_pressure']
                blowdown = input['valve']['blowdown']
                psv_state = 'closed'
        elif input['valve']['type'] == "controlvalve":
            p_back = input['valve']['back_pressure']
            Cv = input['valve']['Cv']
            if 'xT' in input['valve']:
                xT = input['valve']['xT']
            if 'Fp' in input['valve']:
                Fp = input['valve']['Fp']

        # valve type    
        # - constant_mass
        # - functional mass flow
        thickness = 0
        # Reading heat transfer related data/information
        if 'heat_transfer' in input:
            heat_method = input['heat_transfer']['type']
            if heat_method == "specified_h" or heat_method == "specified_U":
                Tamb = input['heat_transfer']['temp_ambient']
            if heat_method == "specified_U": Ufix = input['heat_transfer']['U_fix']
            if heat_method == "specified_Q": Qfix = input['heat_transfer']['Q_fix']
            if heat_method == "specified_h":
                vessel_cp = input['vessel']['heat_capacity']
                vessel_density = input['vessel']['density']
                vessel_orientation = input['vessel']['orientation']
                thickness = input['vessel']['thickness']
                h_out = input['heat_transfer']['h_outer']
                h_in = input['heat_transfer']['h_inner']


    def run(self):
        pass

    def plot(self):
        pass

    def report(self):
        pass


if len(sys.argv) > 1:
    input_filename = sys.argv[1]
else:
    input_filename = "input.yml"

with open(input_filename) as infile:
    input = yaml.load(infile, Loader=yaml.FullLoader)

# Intial parameters and setup
length = input['vessel']['length']
diameter = input['vessel']['diameter']

p0 = input['initial']['pressure']
T0 = input['initial']['temperature']
species = 'HEOS::'+input['initial']['fluid'] 

tstep = input['calculation']['time_step']
time_tot = input['calculation']['end_time']
method = input['calculation']['type'] 
if method == "energybalance": eta = input['calculation']['eta'] 

# Reading valve specific data
if input['valve']['type'] == 'orifice' or input['valve']['type'] == 'psv':
    p_back = input['valve']['back_pressure']
    D_orifice = input['valve']['diameter']
    CD = input['valve']['discharge_coef']
    if input['valve']['type'] == 'psv':
        Pset = input['valve']['set_pressure']
        blowdown = input['valve']['blowdown']
        psv_state = 'closed'
elif input['valve']['type'] == "controlvalve":
    p_back = input['valve']['back_pressure']
    Cv = input['valve']['Cv']
    if 'xT' in input['valve']:
        xT = input['valve']['xT']
    if 'Fp' in input['valve']:
        Fp = input['valve']['Fp']

# valve type
# - constant_mass
# - functional mass flow
thickness = 0
# Reading heat transfer related data/information
if 'heat_transfer' in input:
    heat_method = input['heat_transfer']['type']
    if heat_method == "specified_h" or heat_method == "specified_U":
        Tamb = input['heat_transfer']['temp_ambient']
    if heat_method == "specified_U": Ufix = input['heat_transfer']['U_fix']
    if heat_method == "specified_Q": Qfix = input['heat_transfer']['Q_fix']
    if heat_method == "specified_h":
        vessel_cp = input['vessel']['heat_capacity']
        vessel_density = input['vessel']['density']
        vessel_orientation = input['vessel']['orientation']
        thickness = input['vessel']['thickness']
        h_out = input['heat_transfer']['h_outer']
        h_in = input['heat_transfer']['h_inner']


vol = diameter**2/4 * math.pi * length  # m3
vol_tot = (diameter + 2 * thickness)**2/4 * math.pi * (length + 2 * thickness)  # m3
vol_solid = vol_tot-vol
surf_area_outer = (diameter + 2 * thickness)**2/4 * math.pi * 2 + (diameter + 2 * thickness) * math.pi * (length + 2 * thickness)
surf_area_inner = (diameter)**2/4 * math.pi * 2 + (diameter) * math.pi * length

# data storage
data_len = int(time_tot / tstep)
rho = np.zeros(data_len)
T_fluid = np.zeros(data_len)
T_vessel = np.zeros(data_len)
Q_outer = np.zeros(data_len)
Q_inner = np.zeros(data_len)
h_inside = np.zeros(data_len)
T_vent = np.zeros(data_len)
H_mass = np.zeros(data_len)
S_mass = np.zeros(data_len)
U_mass = np.zeros(data_len)
U_tot = np.zeros(data_len)
U_iter = np.zeros(data_len)
m_iter = np.zeros(data_len)
n_iter = np.zeros(data_len)
Qo = np.zeros(data_len)
Qi = np.zeros(data_len)
P = np.zeros(data_len)
mass_fluid = np.zeros(data_len)
mass_rate = np.zeros(data_len)
time_array = np.zeros(data_len)

rho0 = PropsSI('D', 'T', T0, 'P', p0, species)
m0 = rho0 * vol

# Inititialise
rho[0] = rho0
T_fluid[0] = T0
T_vessel[0] = T0
H_mass[0] = PropsSI('H', 'T', T0, 'P', p0, species)
S_mass[0] = PropsSI('S', 'T', T0, 'P', p0, species)
U_mass[0] = PropsSI('U', 'T', T0, 'P', p0, species)
U_tot[0] = PropsSI('U', 'T', T0, 'P', p0, species) * m0
P[0] = p0
mass_fluid[0] = m0
cpcv = PropsSI('CP0MOLAR', 'T', T0, 'P', p0, species) / PropsSI('CVMOLAR', 'T', T0, 'P', p0, species)

if input['valve']['type'] == 'orifice':
    if input['valve']['flow'] == 'filling':
        k = PropsSI('CP0MOLAR', 'T', T0, 'P', p_back, species) / PropsSI('CVMOLAR', 'T', T0, 'P', p_back, species)
        mass_rate[0] = -tp.gas_release_rate(p_back, p0, PropsSI('D', 'T', T0, 'P', p_back, species), k, CD, D_orifice**2/4 * math.pi)
    else:
        mass_rate[0] = tp.gas_release_rate(p0, p_back, rho0, cpcv, CD, D_orifice**2/4 * math.pi)
elif input['valve']['type'] == 'mdot':
    if input['valve']['flow'] == 'filling':
        mass_rate[0] = -input['valve']['mass_flow']
    else:
        mass_rate[0] = input['valve']['mass_flow']
elif input['valve']['type'] == 'controlvalve':
    if input['valve']['flow'] == 'filling':
        Z = PropsSI('Z', 'T', T0, 'P', p_back, species)
        MW = PropsSI('M', 'T', T0, 'P', p_back, species)
        k = PropsSI('CP0MOLAR', 'T', T0, 'P', p_back, species) / PropsSI('CVMOLAR', 'T', T0, 'P', p_back, species)
        mass_rate[0] = -tp.control_valve(p_back, p0, T0, Z, MW, k, Cv)
    else:
        Z = PropsSI('Z', 'T', T0, 'P', p0, species)
        MW = PropsSI('M', 'T', T0, 'P', p0, species)
        k = PropsSI('CP0MOLAR', 'T', T0, 'P', p0, species) / PropsSI('CVMOLAR', 'T', T0, 'P', p0, species)
        mass_rate[0] = tp.control_valve(p0, p_back, T0, Z, MW, k, Cv)
elif input['valve']['type'] == 'psv':
    if input['valve']['flow'] == 'filling': raise ValueError("Unsupported valve: ", input['valve']['type'], " for vessel filling.")
    mass_rate[0] = tp.relief_valve(p0, p_back, Pset, blowdown, rho0, cpcv, CD, D_orifice**2/4 * math.pi)


time_array[0] = 0
# Run actual integration
for i in range(1, len(time_array)):
    time_array[i] = time_array[i-1] + tstep
    mass_fluid[i] = mass_fluid[i-1] - mass_rate[i-1] * tstep
    rho[i] = mass_fluid[i] / vol

    if method == "isenthalpic":
        T_fluid[i] = PropsSI('T', 'D', rho[i], 'H', H_mass[i-1], species)
        P[i] = PropsSI('P', 'D', rho[i], 'H', H_mass[i-1], species)
    elif method == "isentropic":
        T_fluid[i] = PropsSI('T', 'D', rho[i], 'S', S_mass[i-1], species)
        P[i] = PropsSI('P', 'D', rho[i], 'S', S_mass[i-1], species)
    elif method == "isothermal":
        T_fluid[i] = T0
        P[i] = PropsSI('P', 'D', rho[i], 'T', T0, species)
    elif method == "constantU":
        T_fluid[i] = PropsSI('T', 'D', rho[i], 'U', U_mass[i-1], species)
        P[i] = PropsSI('P', 'D', rho[i], 'U', U_mass[i-1], species)
    elif method == "energybalance":
        P1 = PropsSI('P', 'D', rho[i], 'T', T_fluid[i-1], species)
        T1 = PropsSI('T', 'P', P1, 'H', H_mass[i-1], species)
        

        if heat_method == "specified_h" or heat_method == "detailed":
            if h_in == "calc":
                hi = tp.h_inner(length, T_fluid[i-1], T_vessel[i-1], P[i-1], species)
            else:
                hi = h_in
            h_inside[i] = hi
            Q_inner[i] = surf_area_inner * hi * (T_vessel[i-1] - T_fluid[i-1])
            Q_outer[i] = surf_area_outer * h_out * (Tamb - T_vessel[i-1])
            T_vessel[i] = T_vessel[i-1] + (Q_outer[i] - Q_inner[i]) * tstep / (vessel_cp * vessel_density * vol_solid)
        elif heat_method == "specified_U":
            Q_inner[i] = surf_area_outer * Ufix * (Tamb - T_fluid[i-1])
            T_vessel[i] = T_vessel[0] 
        elif heat_method == "specified_Q":
            Q_inner[i] = Qfix
            T_vessel[i] = T_vessel[0]
        else:
            Q_inner[i] = 0.0
            T_vessel[i] = T_vessel[0]

        NMOL = mass_fluid[i-1] / PropsSI('M', species) 
        NMOL_ADD = (mass_fluid[i]-mass_fluid[i-1]) / PropsSI('M', species) 
        if input['valve']['flow'] == 'filling':
            T_used = T0
            P_used = p_back
        else:
            T_used = T_fluid[i-1]
            P_used = P[i-1]
        U_start = NMOL_ADD * PropsSI('HMOLAR', 'P', P_used, 'T', T_used, species) + NMOL * PropsSI('HMOLAR', 'P', P[i-1], 'T', T_fluid[i-1], species) - eta *  P[i-1] * vol + Q_inner[i] * tstep
        NMOL=NMOL+NMOL_ADD

        U = 0
        nn = 0
        rho1 = 0
        itermax = 1000
        m = 0
        n = 0
        relax = 0.1

        while abs(rho[i] - rho1) > 0.01 and m < itermax:
            m = m + 1
            rho1 = PropsSI('D', 'T', T1,'P', P1, species)
            dd = rho[i] - rho1  #NMOL-nn
            P1 = P1 + dd * 1e5
            if m == itermax:
                raise Exception("Iter max exceeded for rho/P")
            while abs(U_start - U) / U_start > 0.00001 and n < itermax:
                n = n + 1
                U = NMOL * PropsSI('HMOLAR', 'P', P1, 'T', T1, species) - eta * P1 * vol  #Q_inner[i]*tstep
                d = U_start - U 
                T1 = T1 + 0.1 * d / U_start * T1
                if n == itermax:
                    raise Exception("Iter max exceeded for U/T")
            
        m_iter[i] = dd
        n_iter[i] = d
        U_iter[i] = U / NMOL * mass_fluid[i]
        P[i] = P1
        T_fluid[i] = T1
        
    else:
        raise NameError("Unknown calculation method: " + method)

    H_mass[i] = PropsSI('H', 'T', T_fluid[i], 'P', P[i], species)
    S_mass[i] = PropsSI('S', 'T', T_fluid[i], 'P', P[i], species)
    U_mass[i] = (mass_fluid[i] * PropsSI('H', 'P', P[i], 'T', T_fluid[i], species) - P[i] * vol) / mass_fluid[i]  #PropsSI('U','T',T_fluid[i],'P',P[i],species)#-(P[i-1]-P[i])*vol/mass_fluid[i]
    cpcv = PropsSI('CP0MOLAR', 'T', T_fluid[i], 'P', P[i], species) / PropsSI('CVMOLAR', 'T', T_fluid[i], 'P', P[i], species)
    
    if input['valve']['type'] == 'orifice':
        if input['valve']['flow'] == 'filling':
            k = PropsSI('CP0MOLAR', 'T', T0, 'P', p_back, species) / PropsSI('CVMOLAR', 'T', T0, 'P', p_back,species)
            mass_rate[i] = -tp.gas_release_rate(p_back, P[i], PropsSI('D', 'T', T0, 'P', p_back, species), k, CD, D_orifice**2/4 * math.pi)
        else:
            mass_rate[i] = tp.gas_release_rate(P[i], p_back, rho[i], cpcv, CD, D_orifice**2/4 * math.pi)
    elif input['valve']['type'] == 'controlvalve':
        if input['valve']['flow'] == 'filling':
            Z = PropsSI('Z', 'T', T0, 'P', p_back, species)
            MW = PropsSI('M', 'T', T0, 'P', p_back, species)
            k = PropsSI('CP0MOLAR', 'T', T0, 'P', p_back,species) / PropsSI('CVMOLAR', 'T', T0, 'P', p_back, species)
            mass_rate[i] = -tp.control_valve(p_back, P[i], T0, Z, MW, k, Cv)
        else:
            Z = PropsSI('Z', 'T', T_fluid[i], 'P', P[i], species)
            MW = PropsSI('M','T', T_fluid[i], 'P', P[i], species)
            mass_rate[i] = tp.control_valve(P[i], p_back, T_fluid[i], Z, MW, cpcv, Cv)
    elif input['valve']['type'] == 'mdot':
        if input['valve']['flow'] == 'filling':
            mass_rate[i] = -input['valve']['mass_flow']
        else:
            mass_rate[i] = input['valve']['mass_flow']
    elif input['valve']['type'] == 'psv':
        mass_rate[i] = tp.relief_valve(P[i], p_back, Pset, blowdown, rho[i], cpcv, CD, D_orifice**2/4 * math.pi)


import pylab as plt 


plt.figure()
plt.subplot(221)
plt.plot(time_array/60, T_fluid-273.15, 'b', label="Fluid")
plt.plot(time_array/60, T_vessel-273.15, 'g', label="Vessel")
if 'validation' in input:
    if 'temperature' in input['validation']:
        temp = input['validation']['temperature']
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
plt.plot(time_array / 60, P / 1e5, 'b', label="Calculated")
if 'validation' in input:
    if 'pressure' in input['validation']:
        plt.plot(np.asarray(input['validation']['pressure']['time']) / 60, input['validation']['pressure']['pres'], 'ko', label="Experimental")
plt.legend(loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Pressure (bar)')

plt.subplot(223)
plt.plot(time_array / 60, H_mass, 'b', label='H (J/kg)')
plt.plot(time_array / 60, U_mass, 'g', label='U (J/kg)')
plt.plot(time_array / 60, S_mass * 100, 'r', label='S*100 (J/kg K)')
plt.legend(loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Enthalpy/Internal Energy/Entropy')

plt.subplot(224)
plt.plot(time_array / 60, mass_rate, 'b', label='m_dot')
plt.xlabel('Time (minutes)')
plt.ylabel('Vent rate (kg/s)')
plt.show()