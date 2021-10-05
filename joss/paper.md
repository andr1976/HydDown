---
title: 'HydDown: A Python package for calculation of hydrogen (or other gas) pressure vessel filling and discharge'
tags:
  - Python
  - Gas storage
  - Depressurisation
  - Blow-down
  - Pressure cylinder
  - Energy storage
  - Hydrogen
authors:
  - name: Anders Andreasen
    orcid: 0000-0003-0475-323X
    affiliation: 1
affiliations:
 - name: Ramboll Energy, Field Development, Bavnehøjvej 5, DK-6700 Esbjerg, Denmark
   index: 1
date: 04 August 2021
bibliography: paper.bib

---

# Summary
HydDown [@anders_andreasen_2021_5154096] is a Python package for calculation of pressure vessel behaviour during filling (pressurisation) or discharge (depressurisation/blow-down). More specifically, the software allows calculation of vessel pressure, fluid inventory temperature as well as vessel wall temperature as a function of time during either filling or discharge operations. The applications are manifold and some examples are: 

* Rapid filling of vehicle high-pressure hydrogen cylinders
* Rapid emergency discharge from high-pressure gas cylinders
* Dynamic pressure safety valve behaviour for gas-filled pressure vessels subject to fire heat load

A typical system modelled is shown in \autoref{fig:sketch} and it visualizes the key parameters and transport phenomena during gas vessel filling or discharging. The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P). This is caused by a change in fluid inventory (density) due to flow of gas either in or out of the vessel. Furthermore, heat is transfered from or to the surroundings via convective heat transfer on the in- and outside of the vessel - with heat being conducted through the vessel wall. The easiest way to explore the basic capabilities fo the code is via the [HydDown Streamlit app](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py).  

![Gas-filled pressure vessel subject to gas discharge and heat transfer between vessel and gas inventory. \label{fig:sketch}](../docs/img/Sketch.png){ width=70% }

In its essence, the code solves the mass and energy balances with gas thermodynamics calculated using the [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999] including both real gas equation of state and gas transport properties. The energy balance is the first law of thermodynamics for an open system exchanging both heat and mass with the surroundings [@sva]. Heat transfer between gas inventory and vessel wall is accounted for using either natural convection or mixed forced convection/natural convection [@woodfield][@geankoplis]. The mass balance is closed using an applicable flow equation e.g. orifice  [@yellowbook], relief valve [@API520],  control valve [@borden][@ISA][@IEC60534], a fixed mass rate or a predefined mass rate in or out of the pressure vessel. 
The code also allows an external heat load to be applied using the Stefan-Boltzmann equation for fire heat load mimicking background heat load from both pool and jet fire [@scandpower][@API521].

A few choices have been made to keep things simple:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Gas phase only
- No temperature stratification in vessel inventory
- No temperature gradient through vessel wall (applicable for high heat conductivity / thin-walled vessels)

Still, the code allows a variety of mass flow devices and heat transfer both with ambient air and also considering fire heat loads, with the vessel fluid inventory being rigorously described by a Helmholtz energy formulation of a real gas equation of state. Typical calculation output is shown in \autoref{fig:N2discharge} and \autoref{fig:H2filling} with experimental data included for comparison. 

![Calculations of nitrogen discharge emulating experiment I1 from [@Haque1992b]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right). \label{fig:N2discharge}](../docs/img/N2_filling.png)

![Simulation of hydrogen cylinder pressurisation using a pressurisation rate of 10 MPa/min. Comparison between calculated (full line) and measured gas temperature [@STRIEDNIG] (stipulated line) is shown in the upper left graph. \label{fig:H2filling}](../docs/img/Striednig_fillingH2_10MPa_min.png)

# Statement of need
With an increasing demand of clean(er) energy and associated storage such as e.g. on-board hydrogen storage for hydrogen powered vehicles, compressed air energy storage (CAES), compressed biogas, compressed natural gas (CNG) etc., the need for tools, which can simulate the filling and emptying of the pressure containers holding the fluid is indeed present both for operational reasons, but even more importantly for safety reasons. 

Apparently, no free (as in speach, as in beer) tool exists today which accomplishes the same tasks as HydDown. It has been one of the main goals in producing this software to provide a free tool, which can calculate the pressure vessel response during filling and discharge of gas. There are, however, many tools available which can do a subset of HydDown capabilities, the same as HydDown, or even more. Most of these tools are commercial/proprietary (and closed source) and come with a significant license cost. A single tool has been identified, which is free (as in beer), but still closed source [@h2fills]. Below is a partial list of tools, which have comparable capabilities as HydDown - for review of more codes, please refer to [@SHAFIQ2020104].  

| Software                      | Cost                  |  Filling/discharge    | Fire heat laod    |
|-------------------------------|-----------------------|-----------------------|-------------------|
| [Aspen HYSYS](https://www.aspentech.com/en/products/engineering/aspen-hysys)                  |  Fee           |  Both                 | Yes               |
| [Unisim Design](https://www.honeywellprocess.com/en-US/online_campaigns/connected_plant/Pages/process-simulation.html)                 |  Fee           |  Both                 | Yes               |
| VBsim [@DALESSANDRO2015719]   |  N/A             |  Discharge            | No                |
| BLOWSIM [@blowsim]            |  N/A             |  Discharge            | Yes               |
| H2FillS [@h2fills]            |  Free            |  Filling              | No                |

Compared to other software codes, HydDown has a number of advantages: It is open source, computations are fast - especially for single component gas, the code is simple to install, the code is easy to run using a simple input file format, both filling and discharge can be handled and fire heat load can be handled for calculations of emergency response. In this way, the code has advantages over existing codes on one or more parameters. However, the code is still limited to gas phase behaviour and heat transfer through the vessel wall is not rigorously modelled. Hence composite cylinders with different layers of material may give inaccurate results, especially for low heat conductivity materials. Some of the above codes can handle this.   

It is the intent that the present software can be of use for a number of tasks including (but not limited to):

* Design aid for experimental facilities for e.g. mass flow device sizing e.g. orifice, control valve etc. 
* Benchmarking of other software codes. 
* Safety device sizing and review (blow-down valve, pressure safety valve), for instance how fast shall blow-down be? How large should an installed pressure safety valve be?
* Thermal response of vessel during filling/discharge. Is the design (lower/upper) temperatures exceeded? Should the rate be lower? Is any pre-cooling/heating required? 

# Acknowledgements
The making of this project would not have been possible without the great [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999]. The author is also thankful for enlightning discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark) and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurisation, nozzle flow and heat transfer considerations. This work has also benefitted significantly from the understanding of thermo-mechanical pressure vessel behaviour obtained in previous works [@Bjerre2017],[@safety4010011],[@iskov]. Furthermore, this project relies on high quality open source Python packages: NumPy [@Walt:2011:NAS:1957373.1957466][@harris2020array], SciPy [@virtanen2020scipy], matplotlib [@Hunter:2007], pandas [@mckinney-proc-scipy-2010], [PyYaml](https://pyyaml.org/wiki/PyYAMLDocumentation), [cerberus](https://docs.python-cerberus.org/en/stable/), [Streamlit](https://streamlit.io/), tqdm [@tqdmref] and pytest [@pytest]. Sussanne Tolstrup, language secretary at Ramboll Energy, Field Development, is acknowledged for proofreading this manuscript. 

# References
