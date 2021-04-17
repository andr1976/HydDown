---
title: HydDown
subtitle: User guide and technical reference
author: Anders Andreasen
titlepage: true
toc-own-page: true
book: true
reference-section-title: References
bibliography: references.bib
listings: True
---

# Introduction
HydDown is an open source python tool for calculation of Hydrogen (or other pure gas phase species) vessel/container depressurization and filling. The HydDown logo shown in [@Fig:logo] visualizes the key parameters and transport phenomena during gas vessel filling or discharging. The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P). This is caused by change in fluid inventory (density) due to flow of gas either in- or out of the vessel. Further, heat is transfered from or to the surroundings via convective heat transfer on the in- and outside of the vessel - with heat being conducted thorugh the vessel wall. 

![HydDown logo](img/Sketch.png){#fig:logo}

Run the code as simple as: 

    python main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax. 

## Background
This is a small spare time project for calculation of vessel filling and depressurisation behaviour. This is mainly to demonstrate, that although perceived as a very tedious/difficult task to write your own code for such an apparent complex problem, actually a fairly limited amount of code is necessary if you have a good thermodynamic backend. 

A few choices is made to keep things simple to begin with:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Only pure substances are considered
- Gas phase only
- No temperture stratification in the gas phase
- No temperture gradient through vessel wall
- Heat transfer is modelled as simple as possible

These choices makes the problem a lot more simple to solve, First of all the the pure substance Helmholtz energy based equation of state (HEOS) in coolprop offers a lot of convenience in terms of the property pairs/state variables that can be set independently. Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required. That being said the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.


## Requirements

- [Python](http://www.python.org) (3.8 - at least python3)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Coolprop (6.4.1)](http://www.coolprop.org/)
- [cerberus](https://docs.python-cerberus.org/en/latest/)
- [PyYaml](https://pypi.org/project/PyYAML/)
- [pandas](https://pandas.pydata.org/)

The script is running on Windows 10 x64, with stock python installation from python.org and packages installed using pip. Should run om linux (it does on an Ubuntu image on GitHub) or in any conda environment as well, but I haven't checked.

## Units of measure
The SI units are adapted for this project. The following common units are used in the present project and this also applies to the units used in the input files:

Property | Unit | Comment 
----    | ----  | ----
Temperature | K | $^\circ$ C is used in plots
Pressure | Pa    | bar is used in plots
Mass | kg |
Volume | m$^3$ |
Time | s |
Energy | J |
Duty/power | W 
Length | m
Area | m$^2$
Heat flux | W/m$^2$
Heat transfer coefficient | W/(m$^2$ K)
Density | kg/m$^3$
Heat capacity | J/(kg K)

: Unit system {#tbl:units}

As will be noted when presentaing the equations implemented in the code, some of the equations utilise different units than the ones listed in [@tbl:units]. However, it is important to note that the unit conversions are built in to the methods implemented, so the user shall not worry about unit conversion.  

## Credit 
In the making of this document I have sourced a great deal of material (and modified it) from a good collegues M.Sc. thesis [@iskov], co-published papers [@Bjerre2017][@safety4010011] and from on-line material published under permissive licenses (with proper citation). Further, the making of this project would not have possible without the awesome [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999]. I am thankful for enlightning discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark)  and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurisation, nozzle flow and heat transfer considerations.


## License

MIT License

Copyright (c) 2021 Anders Andreasen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Usage 
## Basic usage
Run the code as simple as: 

    python main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax. 

The Yaml input file is edited to reflect the system of interest. 

## Calculation methods 
The following methods are implemented:

- Isothermal i.e. constant temperature of the fluid during depressurisation (for a very slow process with a large heat reservoir)
- Isenthalpic/Adiabatic (no heat transfer with surroundings, no work performed by the expanding fluid)
- Isentropic (no heat transfer with surroundings, PV work performed by the expanding fluid)
- Constant internal energy
- Energy balance. This is the most general case and includes both the ability to transfer heat with surroundings as well as accounting for PV work.

For isothermal/isenthalpic/isentropic/isenergetic calculations the minimal input required are:

- Initial conditions (pressure, temperature)
- vessel dimensions (ID/length)
- valve parameters (Cd, diameter, backpressure)
- Calculation setup (time step, end time)
- Type of gas

If heat transfer is to be considered the calculation type "energybalance" is required. A few options are possible:

- Fixed U (U-value required, and ambient temperature)
- Fixed Q (Q to be applied to the fluid is requried)
- Specified h, the external heat transfer coefficient is provided and either the internal is provided or calculated from assumption of natural convection from a vertical cylinder at high Gr number. Ambient temperature is required.
- Fire (Stefan-Boltzmann equation heat duty)

## Script 

~~~ {.Python}
import yaml
import sys
from hyddown import HydDown

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    hdown.run()
    hdown.verbose=1
    hdown.plot()
~~~

## Module import 

## Input file 

# Theory
In this chapter the basic theory and governing equations for the model implementation in HydDown is presented. The following main topics are covered: 

- thermodynamics
- mass transfer
- heat transfer

## Thermodynamics

### Equation of state
The equation of state used by HydDown is the Helmholtz energy formulation as implemented in [CoolProp](http://www.coolprop.org/) [@doi:10.1021/ie4033999]. Most of the text in the present section has been sourced from the CoolProp documentation to be as accurate and true to the source as possible. The Helmholtz energy formulation is a convenient construction of the equation of state because all the thermodynamic properties of interest can be obtained directly from partial derivatives of the Helmholtz energy.

It should be noted that the EOS are typically valid over the entire range of the fluid, from subcooled liquid to superheated vapor, to supercritical fluid.
In general, the EOS are based on non-dimensional terms $\delta$ and $\tau$, where these terms are defined by

$$ \delta=\rho/\rho_c $$    
$$ \tau=T_c/T $$
    
where $\rho_c$ and $T_c$ are the critical density of the fluid if it is a pure fluid. For pseudo-pure mixtures, the critical point is typically not used as the reducing state point, and often the maximum condensing temperature on the saturation curve is used instead.

The non-dimensional Helmholtz energy of the fluid is given by

$$ \alpha=\alpha^0+\alpha^r $$
    
where $\alpha^0$ is the ideal-gas contribution to the Helmholtz energy, and $\alpha^r$ is the residual Helmholtz energy contribution which accounts for non-ideal behavior.  For a given set of $\delta$ and $\tau$, each of the terms $\alpha^0$ and $\alpha^r$ are known.  The exact form of the Helmholtz energy terms is fluid dependent, but a relatively simple example is that of Nitrogen, which has the ideal-gas Helmholtz energy of

$$  \alpha^0=\ln\delta+a_1\ln\tau+a_2+a_3\tau+a_4\tau^{-1}+a_5\tau^{-2}+a_6\tau^{-3}+a_7\ln[1-\exp(-a_8\tau)] $$
    
and the non-dimensional residual Helmholtz energy of

$$    \alpha^r=\sum_{k=1}^{6}{N_k\delta^{i_k}\tau^{j_k}}+\sum_{k=7}^{32}{N_k\delta^{i_k}\tau^{j_k}\exp(-\delta^{l_k})}+\sum_{k=33}^{36}{N_k\delta^{i_k}\tau^{j_k}\exp(-\phi_k(\delta-1)^2-\beta_k(\tau-\gamma_k)^2)} $$
    
and all the terms other than $\delta$ and $\tau$ are fluid-dependent correlation parameters.

The other thermodynamic parameters can then be obtained through analytic derivatives of the Helmholtz energy terms.  For instance, the pressure is given by

$$    p=\rho RT\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau} \right] $$
    
and the specific internal energy by

$$    \frac{u}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] $$

and the specific enthalpy by


$$    \frac{h}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] +\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+1 $$

which can also be written as

$$    \frac{h}{RT}=\frac{u}{RT}+\frac{p}{\rho RT} $$
    
The specific entropy is given by

$$    \frac{s}{R}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]-\alpha^0-\alpha^r $$
    
and the specific heats at constant volume and constant pressure respectively are given by

$$    \frac{c_v}{R}=-\tau^2 \left[\left(\frac{\partial^2 \alpha^0}{\partial \tau^2}\right)_{\delta}+ \left(\frac{\partial^2 \alpha^r}{\partial \tau^2}\right)_{\delta} \right] $$
    
$$ \frac{c_p}{R}=\frac{c_v}{R}+\dfrac{\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}-\delta\tau\left(\frac{\partial^2 \alpha^r}{\partial \delta\partial\tau}\right)\right]^2}{\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]} $$
    
The EOS is set up with temperature and density as the two independent properties, but often other inputs are known, most often temperature and pressure because they can be directly measured.  As a result, if the density is desired for a known temperature and pressure, it can be obtained iteratively.

### First law for flow process
The control volume sketched in [@Fig:firstlaw], seprated from the surrounding by a control surface, is used as a basis for the analysis of an open thermodynamic system with flowing streams (fs) in and out, according to [@sva]

A general mass balance or continuity equation can be written: 

![Control volume with one entrance and one exit. The image has been sourced from [@firstlaw].](img/First_law_open_system.png){#fig:firstlaw}

$$ \frac{m_{cv}}{dt} + \Delta \left( \dot{m} \right)_{fs}= 0 $$ {#eq:continuity}

The first term is the accumulation term i.e. the rate of change of the mass inside the control volume , $m_cv$ , and the $\Delta$ in the second term represents the diference between the outflow and the inflow

$$ \Delta \left( \dot{m} \right) _{fs} = \dot{m}_2 - \dot{m}_1 $$

An energy balance for the control volume, with the first law of thermodynamics applied, needs to account for all the energy modes with can cross the control surface. Each stream has a total energy 

$$ U + \frac{1}{2}u^2 + zg $$

where the first term is the specfic internal energy, the second term is the kinetic energy and the last term is the potential energy. The rate at which each stream transports energy in or out of the control volume is given by 

$$ \dot{m} (U + \frac{1}{2}u^2 + zg) $$

and in total 

$$  \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs}$$

Further work (not to be confused with shaft work) is also associated with each stream in order to move the stream into or out from the control volume (one can think of a hypothetical piston pushing the fluid at constant pressure), and the work is $PV$ on the basis of the specific fluid volume. The work rate for each stream is 

$$ \dot{m}(PV) $$

and in total 

$$ \Delta\left[ \dot{m}(PV) \right]_{fs} $$

Further, heat may be transfered to (or from) the control volume at a rate $\dot{Q}$ and shaft work may be applied, $\dot{W}_{shaft}$. Combining all this with the accumulation term given by the change in total internal energy the following general energy balance can be written 

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs} + \Delta \left[ \dot{m}(PV) \right]_{fs} = \dot{Q} +\dot{W}_{shaft}   $$

Applying the relation $H = U + PV$, setting $\dot{W}_{shaft} = 0$ since no shaft work is applied to the vessel, and assmuning that kinetic and potential energy changes can be omitted the energy balance simplifies to 

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} H \right]_{fs} = \dot{Q}  $$

The equation can be further simplified if only a single port acting as either inlet or outlet is present 

$$ \frac{d(mU)_{cv}}{dt} + \dot{m} H  = \dot{Q}  $$ {#eq:energybalanmce}

where the sign of $\dot{m}$ determines if the control volume is either emptied of filled. The continuity equation [@eq:continuity] and the energy balance [@eq:energybalanmce] combined with the equation of state are the key equations that shall be solved/intergrated in order to calculate the change in temperature and pressure as a function of time. 

## Flow devices
### Restriction Orifice
When a fluid flows through a constriction or opening such as an orifice, the velocity will be affected by conditions upstream and downstream.
If the upstream pressure is high enough, relative to the downstream pressure, the velocity will reach the speed of sound (Ma = 1) and the flow rate obtained will be the critical flow rate. The maximum downstream pressure for the flow to still be sonic (Ma = 1), is when $P_d = P_c$. The ratio of the critical and upstream pressure is defined by equation [@eq:P_critical].

$$ \frac{P_{c}}{P_u}=\left (\frac{2}{k+1} \right)^\frac{k}{k-1} $$ {#eq:P_critical}

- $P_c$ is the critical pressure. [kPa]
- $P_d$ is the downstream pressure. [kPa]
- $P_u$ is the upstream pressure. [kPa]
- k is the isentropic expansion factor, approximated by the ideal gas heat capacity ratio $C_p/C_v$. [-]

In order to calculated the mass flow rate through an orifice equation 5.49 is used
based on literature from the Committee for the Prevention of Disasters [@yellowbook]. It is
assumed that only gas exits the orifice, as the PSV is positioned on the top of the
vessel

To account for the difference in choked and non-choked flow a set limit pressure is introduced as in equation [@eq:plimit]. If the downstream pressure, $P_{down}$, is below the pressure limit, $ P_{limit}$, then the flow is choked, and the pressure used, $P_{used}$, in equation [@eq:massfloworifice] should be the pressure limit, $P_{limit}$. Otherwise if the downstream pressure, $P_{down}$, is greater than or equal to the pressure limit, $P_{limit}$, the flow is no longer choked and the pressure used should be the downstream pressure, $P_{down}$ [@yellowbook].

$$ P_{limit}=P_{up} \cdot \left ( \frac{2}{k+1} \right ) ^{\frac{k}{k-1}} $$ {#eq:plimit}

$$ \dot{m}_{flow}= C_d  \cdot A \cdot\sqrt{\left ( \frac{2 k}{k-1}\right )  \cdot  P_{up} \cdot \rho \cdot \left (  \frac{P_{used}}{P_{up}} \right )^{\frac{2}{k}} \left (1-\left ( \frac{P_{used}}{P_{up}} \right ) ^{\frac{k-1}{k}} \right )} $$ {#eq:massfloworifice}

- $\rho$ is the density of the gas upstream. $[kg/m^3]$
- $P_{limit}$ is the pressure limit of the upstream absolute pressure. $[bara]$
- $P_{up}$ is the absolute pressure upstream of the orifice. $[bara]$
- $k$ is the ratio of the heat capacities at constant pressure, $C_p$, and at constant volume, $C_v$.
- $P_{down}$ is the absolute pressure downstream of the orifice. $[bara]$
- $P_{used}$ is the pressure used in the mass flow equation based on choked or non-coked conditions. $[bara]$
- $\dot{m}_{flow}$ is the mass flow through the orifice. $[kg/s]$
- $C_d$ is the discharge coefficient of the orifice opening. $[-]$
- $A$ is the cross sectional area of the orifice. $[m^2]$

### Pressure safety valve / Relief valve
A PSV / relief valve is a mechanical device actuated by the static pressure in the vessel and
a conventional PSV is often used for gas/vapor systems. A conventional PSV is
a spring-loaded device which will activate at a predetermined opening pressure,
and relieve the vessel pressure until a given reseat pressure has been reached. Both
the opening pressure and the reset pressure is above the vessel operating pressure,
and the PSV will remain closed until the pressure inside the vessel increases to the
opening pressure.
The operation of a conventional spring-loaded PSV is based on a force balance.
A conventional PSV can be seen in [@Fig:psv]. A spring exerts a force on a disc
blocking the inlet of the PSV. When the pressure inside the vessels reaches the
opening pressure, the force exerted on the disc by the gas, will be larger than the
force exerted by the spring and the PSV will open and allow the gas to flow out of
the vessels. The flow of gas out of the vessels will lower the pressure and thereby
also the force exerted on the disc. When the pressure in the vessels is reduced to
the reset pressure, the PSV will close and the disc will again hinder the gas flow.

![Conventional/pop action PSV adapted from [@iskov] and [@API520]](img/PSV.pdf){#fig:psv}

The relief valve model implemented in HydDown is the API 520 equations [@API520] for gas relief for both sonic/critical as well as subcritical flow. No corrections factors are implemeted in HydDown. 

For sonic flow (critical flow), as indicated in equation, the  mass flow thorugh the PSV can be determined by equation [@eq:Sationary_sizing_sonic].

$$ W = \frac{A {C \cdot K_d \cdot  K_b \cdot  K_c \cdot  P_1}}{\sqrt{\frac{T\cdot Z}{M}}} $$ {#eq:Sationary_sizing_sonic}

- A is the effective discharge area. [mm$^2$]
- W is the mass flow through the device. [kg/h]
- C is a coefficient, as a function of k, as defined in equation [@eq:Sationary_sizing_C_function].
- $K_d$, $K_b$, and $K_c$ are correction factors.
- $P_1$ is the allowable upstream absolute pressure. [kPa]
- T is the temperature of the inlet gas at relieving conditions. [K]
- M is the molecular mass of the gas at relieving conditions. [kg/kmol]
- Z is the compressibility factor for the gas.

$K_d$ is the effective coefficient of discharge, with a typical value of 0.975, for an installed PSV. $K_b$ is a back-pressure correction factor between 0 and 1, assumed to be 1. $K_c$ is a correction factor used when a rupture disk is installed upstream, otherwise it is 1. In the present implementation a value of 1 is assumed. 

$$ C=0.03948 \sqrt{k \left (\frac{2}{k+1}\right)^{\left (\frac{k+1}{k-1}\right)}} $$ {#eq:Sationary_sizing_C_function}

For subsonic flow (subcritical flow), the effective discharge area of the PSV is determined by equation [@eq:Sationary_sizing_subsonic].

$$ W = \frac{A \cdot {F_2 \cdot K_d  \cdot  K_c }}{17.9 \sqrt{\frac{T\cdot Z}{M\cdot P_1\cdot (P_1-P_2)}}} $$ {#eq:Sationary_sizing_subsonic}

F$_2$ is the coefficient of subcritical flow which can be determined from [@eq:Sationary_sizing_F2].

$$ F_2=\sqrt{\left ( \frac{k}{k-1} \right ) r^{\left (\frac{2}{k} \right )} \left (  \frac{1-r^{\left ( \frac{k-1}{k} \right )}}{1-r}\right )} $$ {#eq:Sationary_sizing_F2}

where r is the ratio of backpressure to upstream relieving pressure, $P_2 / P_1$.

When modelling a pop action PSV/relief valve under dynamic conditions, the valve will go from closed to fully open in a short
period of time when the set pressure, $P_{set}$, is reached. The pop action is illustrated in [@Fig:psv_hyst] which shows the opening
and closing hysteresis of the PSV as a function of pressure. In order to close the shall be reduced below the reseat pressure. 

![Relief valve hysteresis adapted from [@iskov]](img/hysteresis.pdf){#fig:psv_hyst}

When specifying PSV's it is comon to use standrd API sizes as shown in [@tbl:psv_sizes]

Size    | Area [in$^2$]     |   Area [m$^2$] 
--------|-------------------|----------------------------
D       |   0.110           |   7.09676 $\cdot$ 10$^{-5}$
E       |   0.196           |   1.26451 $\cdot$ 10$^{-4}$
F       |   0.307           |   1.98064 $\cdot$ 10$^{-4}$
G       |   0.503           |   3.24515 $\cdot$ 10$^{-4}$
H       |   0.785           |   5.06450 $\cdot$ 10$^{-4}$
J       |   1.287           |   8.30320 $\cdot$ 10$^{-4}$
K       |   1.838           |   1.18580 $\cdot$ 10$^{-3}$
L       |   2.853           |   1.84064 $\cdot$ 10$^{-3}$
M       |   3.600           |   2.32257 $\cdot$ 10$^{-3}$
N       |   4.340           |   2.79999 $\cdot$ 10$^{-3}$
P       |   6.380           |   4.11612 $\cdot$ 10$^{-3}$
Q       |   11.050          |   7.12901 $\cdot$ 10$^{-3}$
R       |   16.000          |   1.03225 $\cdot$ 10$^{-2}$
T       |   26.000          |   1.67741 $\cdot$ 10$^{-2}$

: Standard PSV orifice sizes according to API {#tbl:psv_sizes}

### Control Valve
For calculating the the mass flow through a control valve the ANSI/ISA [@borden][@ISA] methodology also described in IEC 60534 [@IEC60534]. 

The flow model for a compressible fluid in the turbulent regime is 

$$ W = C N_6 F_P Y \sqrt{x_{sizing} p_1 \rho_1} $$ 

or equivalent 

$$ W = C N_8 F_P Y \sqrt{\frac{x_{sizing}M}{T_1 Z_1}} $$ 

- C is the flow coefficient ($C_v$ or $K_v$)
- $N_8$ is a unit specific constant, 94.8 for $C_v$ and bar as pressure unit
- $F_P$ is a piping geometry factor [-]
- Y is the expansion factor [-]
- $x_{sizing}$ is is the pressure drop used for sizing[-]
- $p_1$ is the upstream pressure [bar]
- $\rho_1$ is the upstream density [kg/m$^3$]
- M is the molecular weight [kg/kmol]
- $T_1$ is the upstream temperature [K]
- $Z_1$ is the upstream compressibility [-]

In HydDown the piping geometry factor is not yet implemented and asuumed to be 1. The pressure drop ratio $x_{sizing}$ used for sizing is determined as the lesser of the actual pressure drop ratio, $x$, and the choked pressure drop ratio $x_{choked}$. The actual pressure drop ratio is given by:

$$ x \frac{\Delta p}{p_1}$$

The pressure drop ratio at which flow no longer increases with increased value in pressure
drop ratio, is the choked pressure drop ratio, given by the following equation

$$ x_{choked} = F_\gamma x_{TP} $$

The factor $x_T$ is based on air near atmospheric pressure as the flowing fluid with a specific
heat ratio of 1.40. If the specific heat ratio for the flowing fluid is not 1.40, the factor $F_\gamma$ is used to adjust $x_T$. Use the following equation to calculate the specific heat ratio factor:

$$ F_\gamma = \frac{\gamma}{1.4} $$

where $\gamma$ is the ideal gas $C_p/C_v$. It should be noted that the above equation has been derived from perfect gas behaviour and externsion of an orifice model with $\gamma$ in the range of 1.08 to 1.65. If used outside the assumptions flow calculations may become inaccurate. 

The expansion factor Y accounts for the change in density as the fluid passes from the valve
inlet to the vena contracta. It also accounts for the change in the vena contracta area as the
pressure differential is varied.

$$ Y = 1 -  \frac{x_{sizing}}{3x_{choked}}$$

## Heat transfer

### Natural convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to depressurisation can be well approximated by that of natural convection as found from measured Nusselt numbers being well correlated with Rayleigh number, with no apparent improvement in model performance by included the Reynold number [@woodfield].

To determine the heat transfer for the gas-wall interface, Newton's law of cooling is applied, as given in equation [@eq:newton].

$$  \frac{dQ}{dt} = -h A ( T_{s} - T_{gas} ) $$  {#eq:newton}

- $d Q$ is the change in thermal energy due to convective heat transfer. [J]
- $d t$ is the change in time during the heat transfer. [s]
- $h$ is the convective heat transfer. [W/m$^2$ $\cdot$ K ]
- $A$ is the area normal to the direction of the heat transfer. [m$^2$]
- $T_{s}$ is the surface temperature of the geometry. [K]
- $T_{gas}$ is the temperature of the surrounding gas. [K]

Equation [@eq:NewtonsLawOfCooling} indicates the heat transfer to the surface of any geometry such as a plate or cylindrical wall, by means of heat convection from its surroundings [@cengel][@geankoplis].

The convective heat transfer will need to be estimated for the the gas-wall interface, by the use of empirical relations for the Nusselt number.
The Nusselt number describes the ratio of convective heat transfer to conductive heat transfer, normal to a surface area, as given in equation [@eq:Nu].

$$ Nu=\frac{hL}{k} $$ {#eq:Nu}

- $Nu$ is the Nusselt number. [-]
- $h$ is the convective heat transfer. [W/m$^2$$\cdot$K]
- $L$ is a characteristic length of the geometry. [m]
- $k$ is the thermal conductivity of the gas. [W/m$\cdot$K]

The characteristic length $L$ used is the height of the gas volume.  

The empirical correlations used to calculate the Nusselt number of the gas-wall interface is a function of the Rayleigh number, which can be defined by the Grashof number and Prandtl number, as in equation [@eq:rayleigh_gas].

$$ Ra=Gr \cdot Pr $$ {#eq:rayleigh_gas}

- $Ra$ is the Rayleigh number. [-]
- $Gr$ is the Grashof  number. [-]
- $Pr$ is the Prandtl number. [-]

The Grashof number is a dimensionless number which approximates the ratio of the buoyancy forces to viscous forces, as given in equation [@eq:grashof_gas}. The Prandtl number is a dimensionless number defined as the ratio of the momentum diffusivity to thermal diffusivity, as given in equation [@eq:prandtl_gas].

$$ Gr=\frac{\beta g\rho^2 L^3 \Delta T }{\mu^2} $$ {#eq:grashof_gas}

$$ Pr=\frac{c_p \mu}{k} $$ {#eq:prandtl_gas}

- $\beta$ is the coefficient of volume expansion. [1/K]
- $g$ is the standard acceleration of gravity. [m/s$^2$]
- $\rho$ is the gas density. [kg/m$^3$]
- $L$ is the characteristic length. [m]
- $\Delta T$ is the temperature difference of the surface and gas. [K] 
- $\mu$ is the dynamic viscosity. [kg/m$\cdot$s]
- $c_p$ is the heat capacity of gas. [J/kg$\cdot$K]
- $k$ is the thermal conductivity of gas. [J/m$\cdot$K]

It is important to note that the properties in the above equations shall be evaluated at the fluid film temperature which can be approximated by the average of the the fluid bulk temperature and the vessel wall temperature [@geankoplis].

### Mixed convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to filling can be well approximated by that of combined natural convection and forced convenction as found from measured Nusselt numbers being well correlated with Rayleigh and Reynolds number [@woodfield].

For mixed convection the effective Nusselt number, $Nu$, can be approximated by 

$$ Nu = (Nu_{forced}^n + Nu_{natural}^n)^{\frac{1}{n}} $$

During charging with different gases (H$_2$, N$_2$ and Argon), Woodfield *et al.* demonstrated that in roder to provide a good fit to the experimentally determined Nusselt number a correlation based on both Reynolds and Rayleigh number was necessary. They found a good fit with the following formula (n=1)

$$ Nu = Nu_{forced} +  Nu_{natural} = 0.56Re_d^0.67 + 0.104Ra_H^0.352 $$

### Conduction


### Fire heat loads
The heat transfer from the flame to the shell is modelled using the recommended approach from Scandpower [@scandpower]. The heat transfer from the flame to the vessel shell is divided into radiation, convection and reradiation as seen in equation [@eq:flame].

$$ q_f=\underbrace{{\alpha}_s \cdot {\varepsilon}_f \cdot \sigma \cdot T_f^4}_\text{Radiation}+\underbrace{h_f \cdot (T_f-T_s(t))}_\text{Convection}-\underbrace{{\varepsilon}_s \cdot \sigma \cdot T_s(t)^4 }_\text{Reradiation} $$ {#eq:flame}

- $q_f$ is the flame heat flux. [W/m$^2$]
- ${\alpha}_s$ is the vessel surface absorptivity. [-]
- ${\varepsilon}_f$ is the flame emissivity. [-]
- $\sigma$ is the Stefan-Boltzmann constant, $\sigma$ = $5.67 \cdot 10 ^ {-8}$  [W/m$^2 \cdot$ K$^4$]
- $T_f$ is the flame temperature. [K]
- $h_f$ is the convection heat transfer coefficient between the flame and the surface. [W/m$^2 \cdot$ K]
- $T_s(t)$ is the time dependent surface temperature. [K]
- ${\varepsilon}_s$ is the surface emissivity. [-]

This model assumes that the pressure vessel is fully engulfed by the flame. This means that the view factor for the radiation is unity and is therefore not taken into consideration.
The convective heat transfer coefficients for a jet fire and a pool fire, and recommended values for the emissivity and absorptivity, are given by Scandpower as [@scandpower]

- $h_{jet~fire}$ = 100 [W/m$^2 \cdot$K]
- $h_{pool~fire}$ = 30 [W/m$^2 \cdot$K]
- ${\alpha}_s$ = 0.85
- ${\varepsilon}_s$ = 0.85
- ${\varepsilon}_f$ = 1.0 (optical thick flames, thickness > 1 m)

The flame temperature is found by solving equation [@eq:flame2] for the incident heat flux in relation to the ambient conditions. The flame temperature is kept constant throughout the simulation.

$$ q_{total}=\sigma \cdot T_f^4 + h_f \cdot (T_f-T_{amb})$$ {#eq:flame2}

- $q_{total}$ is the incident flame heat flux as given in table [@tbl:heatfluxes1]. [W/m$^2$]
- $T_{amb}$ is the ambient temperature $\approx$ 293 K (20$^\circ$ C)

The heat flux used to calculate the flame temperature is given in table [@tbl:heatfluxes1].

|                        | Small jet fire  [kW/m$^2$]  |  Large jet fire  [kW/m$^2$] |  Pool fire  [kW/m$^2$]
| ----                   |  ----           |  ----           | ----
| Peak heat load         |  250            |  350            | 150         
| Background heat load   |   0             |   100           |  100         

: Incident heat fluxes for various fire scenarios given by Scandpower [@scandpower] {#tbl:heatfluxes1}


## Model implementation

A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation. For each step, the mass relief/ left in the vessel is known. Since the volume is fixed the mass density is directly given. For the calculation methods (isentropic,isenthalpic,isenergetic etc), Coolprop allows specifying density and either H,S or U directly - this is very handy and normally only TP, PH, TS property pairs are implemented, and you would need to code a second loop to make it into am UV, VH or SV calculation. Coolprop is very convenient for this, however for a cubic EOS and for multicomponent Helmholtz energy EOS coolprop only supports a subset of state variables to be specified directly (T,P,quality). For this reason single component HEOS is the main target of this project.  



### Isothermal process

### Isentropic process

### Isenthalpic process


# Validation
The code is provided as-is. However, comparisons have been made to a few experiments from the literature.

The following gases and modes are considered:

- High pressure nitrogen discharge
- High pressure hydrogen filling
- High pressure hydrogen discharge
- Low pressure air discharge 
- Low pressure air filling

## Nitrogen discharge
Calculations with HydDown is compared to  experiment I1 from ref. [@Haque1992b]. The experiment is a blowdown of a vertically oriented cylindrical vessel with flat ends. The vessel length is 1.524 m, the inside diameter is 0.273 m and the wall thickness is 25 mm. The vessel is filled with N$_2$ at 150 bar, at 15$^\circ$C. Ambient temperature is 15$^\circ$C. The blowdown orifice diameter is 6.35 mm. The results are shown in [@Fig:N2val]. The didirscharge coefficient of the orifice has been set to 0.8 in order to match the vessel pressure profile. The back pressure is set to atmospheric conditions.

![Calculations of nitrogen discharge emulating experiment I1 from [@Haque1992b]. The figure shows calculated gas an wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](img/N2_filling.png){#fig:N2val}

As seen from [@Fig:N2val], the calculations compare well with the experimental results. The calculated temperature of the bulk vapor is within the experimental range of measured temperature at all times during the simulation. It is also noted that the minimum temperature is reached at approx. the same time as in the experiments. The calculated vessel inner wall temperature does not decline as rapidly as the experiments—but from around a calculation time of 60 s, the temperature is within the experimentally observed inner wall temperature. The main reason for the inability to match the vessel wall temperature is that the model ignores the temperature gradient from the outer to the inner wall surface and uses an average material temperture. Especially at the beginning of the discharge it is considered likely that a significant temperature gradient will exist. 

## Hydrogen filling 

## Air discharge/filling
