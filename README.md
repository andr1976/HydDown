# HydDown
Hydrogen (or other pure gas phase species) depressurization calculations

## Background

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
