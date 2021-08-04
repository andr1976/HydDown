---
title: 'HydDown: A Python package for calculation of hydrogen (or other gas) pressure vessel filling and discharge'
tags:
  - Python
  - Gas storage
  - Depressurisation
  - Blow-down
  - Pressure cylinder
  - Energy storage
authors:
  - name: Anders Andreasen
    orcid: 0000-0003-0475-323X
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Ramboll Energy, Field Development, Bavnehøjvej 5, DK-6700 Esbjerg, Denmark
   index: 1
date: 04 August 2021
bibliography: paper.bib

---

# Summary
HydDown is a Python package for calculation of pressure vessel behaviour during filling (pressurisation) or discharge (depressurisation/blow-down). More specifically the software allows calculation of vessel pressure, fluid inventory temperature as well as vessel wall temperature as a function of time during either filling or discharge operations. The applications are manifold and some examples are: 

* Rapid filling of vehicle high-pressure hydrogen cylinders
* Rapid emergency discharge from high pressure gas cylinders
* Dynamic pressure safety valve behaviour for gas filled pressure vessels subject to fire heat load

A typical system modelled is shown in \autoref{fig:sketch} and it visualizes the key parameters and transport phenomena during gas vessel filling or discharging. The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P). This is caused by change in fluid inventory (density) due to flow of gas either in- or out of the vessel. Further, heat is transfered from or to the surroundings via convective heat transfer on the in- and outside of the vessel - with heat being conducted thorugh the vessel wall.   

![Gas filled pressure vessel subject to gas discharge and heat transfer between vessel and gas inventory. \label{fig:sketch}](../docs/img/Sketch.png)

In its essence the code solves the mass and energy balances, with gas thermodynamics calculated using the [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999] including both real gas equation of state and gas transport properties. The energy balance is the first law of thermodynamics for an open system exchanging both heat and mass with the surroundings [@sva]. Heat transfer between gas inventory and vessel wall is accounted for using either natural convection or mixed forced convenction/natural convection [@woodfield][@geankoplis].T he mass balance is closed using an applicable flow equation e.g. orifice  [@yellowbook], relief valve [@API520],  control valve [@borden][@ISA][@IEC60534], a fixed mass rate or a predefined mass rate in or out of the pressure vessel. 
The code also allows an external heat load to be applied using the Stefan-Boltzmann equation for fire heat load minimcking both background heat load from pool and jet fire [@scandpower][@API521].


# Statement of need



# Acknowledgements
The making of this project would not have possible without the great [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999]. I am also thankful for enlightning discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark) and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurisation, nozzle flow and heat transfer considerations. This work has also befitted significantly from the understanding of thermo mechanical pressure vessel behaviour obtained in previous works [@Bjerre2017],[@safety4010011],[@iskov]. 

# References