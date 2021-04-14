---
title: HydDown
subtitle: User guide
author: Anders Andreasen
titlepage: true
toc-own-page: true
book: true
logo: img/Sketch.png
logo-width: 100mm
reference-section-title: References
bibliography: references.bib
---

# HydDown

![HydDown logo](img/Sketch.png){#fig:logo}

The HydDown logo shown in [@Fig:logo] visualizes the key parameters and transport phenomena during gas vessel filling or discharging  [@wrigstad2017mastery]. 

Hydrogen (or other pure gas phase species) depressurization calculations

This code is published under an MIT license.

Run the code as simple as: 

    python main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax. 

## Background
This is a small spare time project for calculation of vessel depressurisation behaviour. This is mainly to demonstrate, that although perceived as a very tedious/difficult task to write your own code for such an apparent complex problem, actually a fairly limited amount of code is necessary if you have a good thermodynamic backend. 

A few choices is made to keep things simple to begin with:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Only pure substances are considered
- Gas phase only
- No temperture stratification in the gas phase
- No temperture gradient through vessel wall
- Heat transfer is modelled as simple as possible

Code will be as simple as possible - a single file to start with - no fancy module structure, no GUI, nothing - just a simple script with an input file.
These choices makes the problem a lot more simple to solve, First of all the the pure substance Helmholtz energy based equation of state (HEOS) in coolprop offers a lot of convenience in terms of the property pairs/state variables that can be set independently. Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required. That being said the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.

Over time the code could evolve to implement the Stefan-Bolzmann fire duty equation as described in e.g. API521 / Scandpower guideline, pressure relief valve instead of blow down valve, pressurisation behaviour etc.

## Description
The following methods are implemented:

- Isothermal i.e. constant temperature of the fluid during depressurisation (for a very slow process with a large heat reservoir)
- Isenthalpic/Adiabatic (no heat transfer with surroundings, no work performed by the expanding fluid)
- Isentropic (no heat transfer with surroundings, PV work performed by the expanding fluid)
- Constant internal energy
- Energy balance. This is the most general case and includes both the ability to transfer heat with surroundings as well as accounting for PV work (PV work efficiency can be specied, but 100% or close to is recommended)

A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation. For each step, the mass relief/ left in the vessel is known. Since the volume is fixed the mass density is directly given. For the simple methods (isentropic,isenthalpic,isenergetic etc), Coolprop allows specifying density and either H,S or U directly - this is very handy and normally only TP, PH, TS property pairs are implemented, and you would need to code a second loop to make it into am UV, VH or SV calculation. Coolprop is very convenient for this, however for a cubic EOS and for multicomponent Helmholtz energy EOS coolprop only supports a subset of state variables to be specified directly (T,P,quality). For this reason single component HEOS is the main target of this small project.  

In case the "Energy balance" method is applied, the heat added from convection and work is accounted for. In this case some of the Coolprop convenience is lost and P,T for the new time step is iteratively solved in a nested loop, with pressure being iterated to match the density, and an inner loop with an energy balance solved to estimate the final temperture. The two iteration loops can most likely be improved massively. 

## Basic usage
The Yaml input file is edited to reflect the system of interest. For isothermal/isenthalpic/isentropic/isenergetic calculations the minimal input required are:

- Initial conditions (pressure, temperature)
- vessel dimensions (ID/length)
- valve parameters (Cd, diameter, backpressure)
- Calculation setup (time step, end time)
- Type of gas

If heat transfer is to be considered the calculation type "energybalance" is required. A few options are possible:

- Fixed U (U-value required, and ambient temperature)
- Fixed Q (Q to be applied to the fluid is requried)
- Specified h, the external heat transfer coefficient is provided and either the internal is provided or calculated from assumption of natural convection from a vertical cylinder at high Gr number. Ambient temperature is required.
- Detailed (not yet implememted)
- Fire (Not yet implemented)

## Validation
The code is mainly for demonstration, and is provided as-is. However, a comparison has been made to experiment I1 from Haque et al. (see full ref below). The results are shown below.

## Requirements

- [Python](http://www.python.org) (3.8 - at least python3)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Coolprop (6.4.1)](http://www.coolprop.org/)
- Yaml 

The script is running on Windows 10 x64, with stock python installation from python.org and packages installed using pip. Should run om linux or in any conda environment as well, but I haven't checked.

## TO-DO
- Addition of heat transfer correlations for different geometries and also considering mixed convection
- Fire heat load using Stefan-Boltzmann equation as per API521/Scandpower guideline
- Vessel loading/filling
- Relief valve (spring loaded) as an alternative to orifice depressurisation
- Additional plots
- More description
- Manual
- Example input files
- More validation cases with different gases
- Multi-component / multi-phase with cubic EOS (long term goal)

## References

- [C. J. Geankoplis; Transport Processes and Unit Operations, International Edition, Prentice-Hall, 1993](https://www.amazon.co.uk/Transport-Processes-Unit-Operations-International/dp/013045253X)
- Haque, M.A.; Richardson, S.M.; Saville, G.; Chamberlain, G.; Shirvill, L. Blowdown of pressure vessels II. Experimental Validation of Computer Model and Case Studies. Trans. IChemE 1992, 70, 10–17.
- [Woodfield, Monde, Mitsutake; Measurement of Averaged Heat Transfer Coefficients in High-Pressure Vessel during Charging with Hydrogen, Nitrogen or Argon Gas; Journal of Thermal Science and Technology; 2007 Volume 2 Issue 2 Pages 180-191 ](https://doi.org/10.1299/jtst.2.180)
- [Woodfield, Monde, Takano; Heat Transfer Characteristics for Practical Hydrogen Pressure Vessels Being Filled at High Pressure, Journal of Thermal Science and Technology 3(2):241-253](http://dx.doi.org/10.1299/jtst.3.241)
- [On the Adequacy of API 521 Relief-Valve Sizing Method for Gas-Filled Pressure Vessels Exposed to Fire](https://doi.org/10.3390/safety4010011)
- American Petroluem Institute. Pressure-Relieving and Depressuring Systems, API Standard 521; American Petroleum Institute: Washington, DC, USA, 2014.
- Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.
