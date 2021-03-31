# HydDown
Hydrogen (or other pure gas phase species) depressurization calculations

## Background
This is a small spare time project for calculation of vessel depressurisation behaviour. This is mainly to demonstrate, that although perceived as a very tedious/difficult task to write your own code for such an apparent complex problem, actually a fairly limited amount of code is necessary if you have a good thermodynamic backend. 

A few choices is made to keep things simple to begin with:

- Coolprop is used as thermodynamic backend
- Only pure substances are considered
- Gas phase only
- Heat transfer is modelled as simple as possible

Code will be as simple as possible - likely a single file to start with - no fancy module structure, manual input for starters - no input file parser, no GUI, nothing - just a simple script.
These choices make the problem a lot more simple to solve, First of all the the pure substance Helmholtz energy based equation of state in coolprop offers a lot of convenience in terms of the property pairs/state variables that can be set independtly. Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required. That being said the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.

Over time the code could evolve to account for convective heat transfer inside/outside vessel, also with a heat balance of the vessel walls, acting to attenuate the temperature drop during depressurisation. Implenentation of the Stefan-Bolzmann fire duty equation as described in e.g. API521 / Scandpower guideline, pressure relief valve instead of blow down valve, pressurisation behaviour etc.

## Description


## Requirements

- [Python](http://www.python.org) (3.8 - at least python3)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Coolprop (6.4.1)](http://www.coolprop.org/)

The script is running on Windows 10 x64, with stock python installation from python.org and packages installed using pip. Should run om linux or in any conda environment as well, but I haven't checked.

## TO-DO

- Internal/external heat transfer
- Fire heat load using Stefan-Boltzmann equation as per API521/Scandpower guideline
- Vessel wall temperature
- Vessel loading/filling
- Relief valve (spring loaded) as an alternative to orifice depressurisation
- Config file
- Additional plots
- More description

## References

- [C. J. Geankoplis; Transport Processes and Unit Operations, International Edition, Prentice-Hall, 1993](https://www.amazon.co.uk/Transport-Processes-Unit-Operations-International/dp/013045253X)
- HYSYS/Unisim
- [Woodfield, Monde, Takano; Heat Transfer Characteristics for Practical Hydrogen Pressure Vessels Being Filled at High Pressure, Journal of Thermal Science and Technology 3(2):241-253](http://dx.doi.org/10.1299/jtst.3.241)
- [On the Adequacy of API 521 Relief-Valve Sizing Method for Gas-Filled Pressure Vessels Exposed to Fire](https://doi.org/10.3390/safety4010011)
- American Petroluem Institute. Pressure-Relieving and Depressuring Systems, API Standard 521; American Petroleum Institute: Washington, DC, USA, 2014.
- Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.
