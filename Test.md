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
HydDown is an open source python3 tool for calculation of hydrogen (or other pure gas phase species) vessel/container depressurization and filling.
The HydDown logo shown in [@Fig:logo] visualizes the key parameters and transport phenomena during gas vessel filling or discharging.
The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P).
This is caused by change in fluid inventory (density) due to flow of gas either in or out of the vessel.
Further, heat is transferred from or to the surroundings via convective heat transfer on the in- and outside of the vessel with heat being conducted through the vessel wall.

![HydDown logo](img/Sketch.png){#fig:logo}

Running the code is as simple as: 

    python main.py input.yml

where `main.py` is the main script and `input.yml` is the input file in Yaml syntax. 
