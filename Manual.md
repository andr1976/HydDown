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

![HydDown logo](docs/img/Sketch.png){#fig:logo}

Running the code is as simple as:

    python main.py input.yml

where `main.py` is the main script and `input.yml` is the input file in Yaml syntax.

## Citing HydDown
If you use HydDown please cite the following reference [@Andreasen2021]:

Andreasen, A., (2021). HydDown: A Python package for calculation of hydrogen (or other gas) pressure vessel filling and discharge. Journal of Open Source Software, 6(66), 3695, https://doi.org/10.21105/joss.03695

    @article{Andreasen2021,
      doi = {10.21105/joss.03695},
      url = {https://doi.org/10.21105/joss.03695},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {66},
      pages = {3695},
      author = {Anders Andreasen},
      title = {HydDown: A Python package for calculation of hydrogen (or other gas) pressure vessel filling and discharge},
      journal = {Journal of Open Source Software}
    }

This manual can be cited as [@Andreasen2024]:

Anders Andreasen. HydDown - User guide and technical reference. 2024. ⟨hal-04858235⟩

    @report{Andreasen2024,
      url = {https://hal.science/hal-04858235},
      year = {2024},
      publisher = {HAL open science},
      author = {Anders Andreasen},
      title = {HydDown -- user guide and technical reference},
    }

## Background
HydDown started as a small spare-time project for calculation of vessel filling and depressurization behaviour.
This was mainly to demonstrate that although perceived as a very tedious and difficult task to write your own code for such an apparently complex problem, actually only a fairly limited amount of code is necessary if you have a good thermodynamic backend.
The code has evolved to a point where additional features will increase the complexity to the point where it is no longer a simple tool.

A few choices has been made to keep things simple:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Only pure substances are considered (limited multi-component capabilities are included)
- Gas phase only
- No temperature stratification in the gas phase
- A default of of no temperature gradient through vessel wall. 
- Heat transfer is modelled as constant or simplified using empirical
  correlations

These choices make the problem much simpler to solve:
First of all, the the pure substance Helmholtz energy based equation of state (HEOS) in `coolprop` offers a lot of convenience in terms of the property pairs/state variables that can be set independently.
Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required.
That being said, the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.

In the latest revision 1-D transient heat conduction through the vessel wall is now an option if required for low thermal conductivity materials and e.g. type III/IV vessels. 

## Getting the software
The source code can be obtained either from GitHub (via `git` or via the latest tar-ball release) or via **pip** . No packaged releases have currently been planned for **conda**.  

The main branch is located here:

[`https://github.com/andr1976/HydDown`](https://github.com/andr1976/HydDown)

Clone the repo by:

    git clone https://github.com/andr1976/HydDown.git

Running from source via `git` the dependencies must be installed manually from the repo root dir:

    pip install -r requirements.txt

Installation of latest release via **pip** also installe dependecies automatically:

    pip install hyddown

In case `pip` links to a v2.7 of python you will get an error. If so try the following:

    python3 -m pip install hyddown

where python3 is the symlink or full path to the python3 executable installed on your system.

## Requirements

- [Python](http://www.python.org) (Check CoolProp for latest suppeorted python version)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Coolprop (6.4.1)](http://www.coolprop.org/)
- [cerberus](https://docs.python-cerberus.org/en/latest/)
- [PyYaml](https://pypi.org/project/PyYAML/)
- [pandas](https://pandas.pydata.org/)
- [Scipy](https://www.scipy.org/)
- [tqdm](https://tqdm.github.io/)

The script is running on Windows 10 x64, with stock python installation from python.org and packages installed using pip.
It should also run on Linux (it does on an Ubuntu image on GitHub) or in any conda environment as well, but this hasn't been checked.

## Testing
Although testing is mainly intended for automated testing (CI) during development using github actions, testing of the installed package can be done for source install by:

    python -m pytest

run from the root folder. In writing 33+ test should pass.

For the package installed with pip navigate to the install directory `../site-packages/hyddown` and run:

    python -m pytest

## Units of measure
The SI unit system is adapted for this project.
The following common units are used in the present project and this also applies to the units used in the input files:

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
Thermal conductivity | W/(m K)
Density | kg/m$^3$
Heat capacity | J/(kg K)

: Unit system {#tbl:units}

As will be noted when presenting the equations implemented in the code, some of the equations utilise different units than the ones listed in [@tbl:units].
However, it is important to note that unit conversions are built in to the methods implemented, so the user shall not worry about unit conversion.  

## Credit
In the making of this document a great deal of material has been sourced (and modified) from a good colleague's M.Sc. thesis [@iskov], from co-published papers [@Bjerre2017][@safety4010011] and from on-line material published under permissive licenses (with proper citation).
Further, the making of this project would not have possible without the awesome [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999].
I am also thankful for enlightening discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark)  and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurization, nozzle flow and heat transfer considerations.

The present document is typeset using Markdown + [pandoc](https://pandoc.org/) with the [Eisvogel](https://github.com/Wandmalfarbe/pandoc-latex-template) template.

## License

MIT License

Copyright (c) 2021-2025 Anders Andreasen

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
From the source repository, running the code is as simple as:

    python src/main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax.

The Yaml input file is edited to reflect the system of interest.

Further, a usable copy of a main script is installed in the python installation Scripts/ folder:

    /python3x/Scripts/hyddown_main.py

or from GitHub/in the source folder:

    HydDown/Scripts/hyddown_main.py

Various example files are included for inspiration under site-packages:

    /site-packages/hyddown/examples/

or in the git/source directory tree

    HydDown/src/hyddown/examples/

## Demos
A few demonstrations of the codes capability are available as [`streamlit`](https://streamlit.io/) apps:

- Vessel gas pressurisation/depressurization for various gases at [`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py)
- Response to external fire of vessel equipped with a PSV using the Stefan-Boltzmann fire equation at[`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_sbapp.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_sbapp.py)
- Vessel depressurization including slow opening of valve at [`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_cvapp.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_cvapp.py)

The various streamlit apps are also included in the scripts folder for running the scripts on a local machine.

## Calculation methods
The following methods are implemented:

- Isothermal: constant temperature of the fluid during depressurization (for a very slow process with a large heat reservoir)
- Isenthalpic: constant enthalpy of the fluid, no heat transfer with surroundings, no work performed by the expanding fluid for depressurization, no work added from inlet stream
- Isentropic: constant entropy of the fluid, no heat transfer with surroundings, PV work performed by the expanding fluid
- Isenergetic: constant internal energy of the fluid
- Energy balance: this is the most general case and is based on the first law of thermodynamics applied to a flow process.

For `isothermal`/`isenthalpic`/`isentropic`/`isenergetic` calculations the minimal input required are:

- Initial conditions (pressure, temperature)
- vessel dimensions (ID, length)
- valve parameters (Cd, diameter, backpressure)
- Calculation setup (time step, end time)
- Type of gas

If heat transfer is to be considered the calculation type `energybalance` is required. A few options are possible:

- Fixed U: U-value, i.e. overall heat transfer coefficient, is required and ambient temperature required.
Based on the ambient (external) temperature, the fluid temperature, and the U-value, the heat flux Q is calculated for each time step.
In this calculation method the temperature of the vessel wall is not calculated:

$$ Q = UA(T_{ambient} - T_{fluid})  $$

- Fixed Q: The external heat flux, Q, applied to the fluid is specified and constant.
The ambient temperature is not required.
Further, the vessel wall temperature is redundant and not calculated.
- Specified h: The external heat transfer coefficient must be specified and either the internal heat transfer coefficient is provided or calculated from the assumption of natural convection from a vertical cylinder at high Gr number.
Ambient temperature is required. Using this method the wall temperature is calculated from an energy balance over the vessel wall taking in and out flux to/from the external ambient plenum as well as heat flux to/from the fluid inventory into account.
- Fire: The Stefan-Boltzmann equation is applied for estimating the external heat duty.
The fire heat flux depends on the vessel wall temperature, and the wall temperature is continuously updated as in the method with specified h.

More elaborate description of the required input for the different calculation types are provided in [@Sec:input].

## Script
HydDown comes with a script which can be used as a command-line tool to start calculations.
If an input filename (path) is given as the first argument, this input file will be used.
If no arguments are passed, the script will look for an input file with the name `input.yml`.
The content of the main script is shown below:


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
To use HydDown simple import the main calculation class `HydDown`.

~~~ {.Python}
from hyddown import HydDown
~~~

## Input file examples
When using HydDown a dictionary holding all relevant input in order for HydDown to do vessel calculations shall be provided when the class is initialized. One way is to read an input file. For HydDown a Yaml format is chosen, but if JSON is a preference this should also work with the Yaml parser being substituted with a JSON parser.

An example of a minimal file for an isentropic vessel depressurization (no heat transfer) is shown below

~~~ {.Yaml}
vessel:
  length: 1.524
  diameter: 0.273
initial:
  temperature: 388.0
  pressure: 15000000.
  fluid: "N2"
calculation:
  type: "isentropic"
  time_step: 0.05
  end_time: 100.
valve:
  flow: "discharge"
  type: "orifice"
  diameter: 0.00635
  discharge_coef: 0.8
  back_pressure: 101300.
~~~

A more elaborate example which includes heat transfer and with validation data (some data points dropped for simplicity) included:

~~~ {.Yaml}
vessel:
  length: 1.524
  diameter: 0.273
  thickness: 0.025
  heat_capacity: 500
  density: 7800.
  orientation: "vertical"
initial:
  temperature: 288.0
  pressure: 15000000.
  fluid: "N2"
calculation:
  type: "energybalance"
  time_step: 0.05
  end_time: 100.
valve:
  flow: "discharge"
  type: "orifice"
  diameter: 0.00635
  discharge_coef: 0.8
  back_pressure: 101300.
heat_transfer:
  type: "specified_h"
  temp_ambient: 288.
  h_outer: 5
  h_inner: 'calc'
validation:
  temperature:
    gas_high:
      time: [0.050285, ... , 99.994]
      temp: [288.93, ... ,241.29]
    gas_low:
      time: [0.32393, ..., 100.11]
      temp: [288.67, ... ,215.28]
    wall_low:
      time: [0.32276, ... , 100.08]
      temp: [288.93, ... ,281.72]
    wall_high:
      time: [0.049115, ... ,100.06]
      temp: [289.18, ... ,286.09]
  pressure:
    time: [0.28869, ... , 98.367]
    pres: [150.02, ... ,1.7204]
~~~

## Input fields and hierarchy {#sec:input}
In the following the full hierarchy of input for the different calculation types is summarised.

At the top level the following fields are accepted, with the last being optional and the second last dependant on calculation type:

~~~ {.Yaml}
initial: mandatory
vessel: mandatory
calculation: mandatory
valve: mandatory
heat_transfer: depends on calculation type
validation: optional
~~~

### Calculation

The subfields under `calculation`, with value formats and options are:

~~~ {.Yaml}
calculation:
  type: "isothermal", "isentropic", "isenthalpic", "constantU", "energybalance"
  time_step: number
  end_time: number
~~~

The simulation end time is specified as well as the fixed time step used in the integration of the differential equations to be solved.
The four main calculation types are shown as well.

### Vessel

~~~ {.Yaml}
vessel:
  length: number, mandatory
  diameter: number, mandatory
  thickness: number, required when heat transfer is calculated
  heat_capacity: number, required when heat transfer is calculated
  density: number, required when heat transfer is calculated
  thermal_conductivity: number, required for 1-D transient heat transfer
  liner_thickness: number, required only for bi-material 1-D transient heat transfer
  liner_heat_capacity: number, required only for bi-material 1-D transient heat transfer 
  liner_density: number, required only for bi-material 1-D transient heat transfer
  liner_thermal_conductivity: number, required only for bi-material 1-D transient heat transfer  
  orientation: string, required when heat transfer is calculated
~~~

### Initial
~~~ {.Yaml}
initial:
  temperature: number, mandatory
  pressure: number, mandatory
  fluid: string, mandatory , e.g. "N2", "H2"
~~~

### Valve
The `valve` field determines the mode of the process is it for depressurization/discharge from the vessel using the value `discharge` or if mass flow is entering the vessel, filling it up/increasing pressure with the value `filling`.

Different types of mass flow devices can be specified:

- Restriction `orifice`
- Relief valve/pressure safety valve (only for `discharge` not `filling`)
- Control valve (`controlvalve`)
- A specified mass flow (`mdot`)

For the physical devices a `back_pressure` is required for the flow calculations.
The value of the `back_pressure` **is also used to specify the reservoir pressure when the vessel is filled**.
See also [@Sec:flow] for details about the calculation of flow rates.

~~~ {.Yaml}
valve:
  flow: string, mandatory "discharge" or "filling"
  type: string, mandatory "orifice", "controlvalve", "psv", "mdot"
  back_pressure: number, required for type "orifice", "controlvalve" and "psv"
  diameter: number, required for "orifice" and "psv"
  discharge_coef: number, required for "orifice" and "psv"
  Cv: number, required for "control_valve"
  characteristic:  string, optional for "control_valve"
  time_constant: number, optional for "control_valve"
~~~

For the calculations using a control valve as mass transfer device a linear rate for the actuator can be specified using the `time_constant` field. An associated valve characteristic is associated in order to translate the linear actuator rate to an opening Cv (of max. Cv). Three caracteristics can be specified:

- Linear
- Equal percentage (exponential)
- Quick opening/fast (square root)

### Heat transfer
For more information about the actual estimation of heat transfer see also [@Sec:heat]. For the case of `fire` heat input predetermined parameters for the Stefan-Boltzmann fire equation are used to calculate different background heat loads cf. [@Tbl:fire]

| Source          |  Fire type      |       Heat load ($kW/m^2$) |
|-----------------|-----------------|----------------------------|
| API521          | Pool            |                   60       |
| API521          | Jet             |                  100       |
| Scandpower      | Pool            |                  100       |
| Scandpower      | Jet             |                  100       |

: Fire heat loads {#tbl:fire}

~~~ {.Yaml}
heat_transfer:
  type: string, mandatory, "specified_h", "specified_Q", "specified_U", "s-b"
  temp_ambient: number, required for type "specified_h", "specified_U"
  h_outer: number, required for type "specified_h"
  h_inner: number or 'calc', required for type "specified_h"
  U_fix: number, required for type "specified_U"
  Q_fix: number, required for type "specified_Q"
  fire: string, required for type "s-b", 'api_pool','api_jet','scandpower_pool','scandpower_jet'
  D_thoat: number, required for flow type "filling", set to vessel ID as a starting point
~~~

### Validation
In order to plot measured data against simulated data the field `validation` is included.

The following arrays (one or more) are supported in the sub-field `temperature`:

- `gas_high`: highest measured values of the bulk gas
- `gas_low`: lowest measured values of the bulk gas
- `gas_mean`: average measured values of the bulk gas
- `wall_mean`: average measured values of the vessel (inner) wall
- `wall_high`  i.e. highest measured values of the vessel (inner) wall
- `wall_low`  i.e. lowest measured values of the vessel (inner) wall
- `wall_outer`  i.e. external measured values of the vessel (outer) wall
- `wall_inner`  i.e. internal measured values of the vessel (inner) wall

For each of the above fields arrays for `time`and `temp`shall be supplied with matching length.

A field for measured vessel pressure is also possible, where `time`and `pres` shall be identical length arrays. See also example below.

~~~ {.Yaml}
validation:
  temperature:
    gas_high:
      time: [0.050285, ... , 99.994]
      temp: [288.93, ... ,241.29]
    gas_low:
      time: [0.32393, ..., 100.11]
      temp: [288.67, ... ,215.28]
    wall_low:
      time: [0.32276, ... , 100.08]
      temp: [288.93, ... ,281.72]
    wall_high:
      time: [0.049115, ... ,100.06]
      temp: [289.18, ... ,286.09]
  pressure:
    time: [0.28869, ... , 98.367]
    pres: [150.02, ... ,1.7204]
~~~


# Theory
In this chapter the basic theory and governing equations for the model implementation in HydDown is presented.
The following main topics are covered:

- thermodynamics
- mass transfer
- heat transfer

## Thermodynamics

### Equation of state
The equation of state used by HydDown is the Helmholtz energy formulation as implemented in [CoolProp](http://www.coolprop.org/) [@doi:10.1021/ie4033999].
Most of the text in the present section has been sourced from the CoolProp documentation to be as accurate and true to the source as possible.
The Helmholtz energy formulation is a convenient construction of the equation of state because all the thermodynamic properties of interest can be obtained directly from partial derivatives of the Helmholtz energy.

It should be noted that the EOS are typically valid over the entire range of the fluid, from subcooled liquid to superheated vapor, to supercritical fluid.
In general, the EOS are based on non-dimensional terms $\delta$ and $\tau$, where these terms are defined by

$$ \delta=\rho/\rho_c $$
$$ \tau=T_c/T $$

where $\rho_c$ and $T_c$ are the critical density of the fluid if it is a pure fluid.
For pseudo-pure mixtures, the critical point is typically not used as the reducing state point, and often the maximum condensing temperature on the saturation curve is used instead.

The non-dimensional Helmholtz energy of the fluid is given by

$$ \alpha=\alpha^0+\alpha^r $$

where $\alpha^0$ is the ideal-gas contribution to the Helmholtz energy, and $\alpha^r$ is the residual Helmholtz energy contribution which accounts for non-ideal behavior.
For a given set of $\delta$ and $\tau$, each of the terms $\alpha^0$ and $\alpha^r$ are known.
The exact form of the Helmholtz energy terms is fluid dependent, but a relatively simple example is that of Nitrogen, which has the ideal-gas Helmholtz energy of

$$  \alpha^0=\ln\delta+a_1\ln\tau+a_2+a_3\tau+a_4\tau^{-1}+a_5\tau^{-2}+a_6\tau^{-3}+a_7\ln[1-\exp(-a_8\tau)] $$

and the non-dimensional residual Helmholtz energy of

$$    \alpha^r=\sum_{k=1}^{6}{N_k\delta^{i_k}\tau^{j_k}}+\sum_{k=7}^{32}{N_k\delta^{i_k}\tau^{j_k}\exp(-\delta^{l_k})}+\sum_{k=33}^{36}{N_k\delta^{i_k}\tau^{j_k}\exp(-\phi_k(\delta-1)^2-\beta_k(\tau-\gamma_k)^2)} $$

and all the terms other than $\delta$ and $\tau$ are fluid-dependent correlation parameters.

The other thermodynamic parameters can then be obtained through analytic derivatives of the Helmholtz energy terms.
For instance, the pressure is given by

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

The EOS is set up with temperature and density as the two independent properties, but often other inputs are known, most often temperature and pressure because they can be directly measured.
As a result, if the density is desired for a known temperature and pressure, it can be obtained iteratively.

### First law for flow process {#sec:firstlaw}
The control volume sketched in [@Fig:firstlaw], separated from the surrounding by a control surface, is used as a basis for the analysis of an open thermodynamic system with flowing streams (fs) in and out, according to [@sva]

A general mass balance or continuity equation can be written:

![Control volume with one entrance and one exit. The image has been sourced from [@firstlaw].](docs/img/First_law_open_system.png){#fig:firstlaw}

$$ \frac{m_{cv}}{dt} + \Delta \left( \dot{m} \right)_{fs}= 0 $$ {#eq:continuity}

The first term is the accumulation term i.e. the rate of change of the mass inside the control volume, $m_cv$, and the $\Delta$ in the second term represents the difference between the outflow and the inflow

$$ \Delta \left( \dot{m} \right) _{fs} = \dot{m}_2 - \dot{m}_1 $$

An energy balance for the control volume, with the first law of thermodynamics applied, needs to account for all the energy modes with can cross the control surface. Each stream has a total energy

$$ U + \frac{1}{2}u^2 + zg $$

where the first term is the specific internal energy, the second term is the kinetic energy and the last term is the potential energy.
The rate at which each stream transports energy in or out of the control volume is given by

$$ \dot{m} (U + \frac{1}{2}u^2 + zg) $$

and in total

$$  \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs}$$

Furthermore, work (not to be confused with shaft work) is also associated with each stream in order to move the stream into or out from the control volume (one can think of a hypothetical piston pushing the fluid at constant pressure), and the work is equal to $PV$ on the basis of the specific fluid volume.
The work rate for each stream is

$$ \dot{m}(PV) $$

and in total

$$ \Delta\left[ \dot{m}(PV) \right]_{fs} $$

Further, heat may be transferred to (or from) the control volume at a rate $\dot{Q}$ and shaft work may be applied, $\dot{W}_{shaft}$.
Combining all this with the accumulation term given by the change in total internal energy the following general energy balance can be written:

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs} + \Delta \left[ \dot{m}(PV) \right]_{fs} = \dot{Q} +\dot{W}_{shaft}   $$

Applying the relation $H = U + PV$, setting $\dot{W}_{shaft} = 0$ since no shaft work is applied to or extracted from the vessel, and assuming that kinetic and potential energy changes can be omitted the energy balance simplifies to:

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} H \right]_{fs} = \dot{Q}  $$

The equation can be further simplified if only a single port acting as either inlet or outlet is present:

$$ \frac{d(mU)_{cv}}{dt} + \dot{m} H  = \dot{Q}  $$ {#eq:energybalance}

where the sign of $\dot{m}$ determines if the control volume is either emptied or filled.
The continuity equation [@Eq:continuity] and the energy balance [@Eq:energybalance] combined with the equation of state are the key equations that shall be solved/integrated in order to calculate the change in temperature and pressure as a function of time.

## Flow devices {#sec:flow}
### Restriction Orifice
When a fluid flows through a restriction or opening such as an orifice, the velocity will be affected by conditions upstream and downstream.
If the upstream pressure is high enough, relative to the downstream pressure, the velocity will reach the speed of sound (Ma = 1) and the flow rate obtained will be the critical flow rate.
This condition is referred to as choked flow.
The maximum downstream pressure for the flow to still be sonic (Ma = 1), is when $P_d = P_c$.
The ratio of the critical and upstream pressure is defined by equation [@Eq:P_critical].

$$ \frac{P_{c}}{P_u}=\left (\frac{2}{k+1} \right)^\frac{k}{k-1} $$ {#eq:P_critical}

- $P_c$ is the critical pressure. [kPa]
- $P_d$ is the downstream pressure. [kPa]
- $P_u$ is the upstream pressure. [kPa]
- k is the isentropic coefficient, approximated by the ideal gas heat capacity ratio $C_p/C_v$. [-]

In order to calculate the mass flow rate through an orifice equations are used based on literature from the Committee for the Prevention of Disasters [@yellowbook].

To account for the difference in choked and non-choked flow a set limit pressure is introduced as in equation [@Eq:plimit].
If the downstream pressure, $P_{down}$, is below the pressure limit, $ P_{limit}$, then the flow is choked, and the pressure used, $P_{used}$, in equation [@Eq:massfloworifice] should be the pressure limit, $P_{limit}$.
Otherwise if the downstream pressure, $P_{down}$, is greater than or equal to the pressure limit, $P_{limit}$, the flow is no longer choked and the pressure used should be the downstream pressure, $P_{down}$ [@yellowbook].

$$ P_{limit}=P_{up} \cdot \left ( \frac{2}{k+1} \right ) ^{\frac{k}{k-1}} $$ {#eq:plimit}

$$ \dot{m}_{flow}= C_d  \cdot A \cdot\sqrt{\left ( \frac{2 k}{k-1}\right )  \cdot  P_{up} \cdot \rho \cdot \left (  \frac{P_{used}}{P_{up}} \right )^{\frac{2}{k}} \left (1-\left ( \frac{P_{used}}{P_{up}} \right ) ^{\frac{k-1}{k}} \right )} $$ {#eq:massfloworifice}

- $\rho$ is the density of the gas upstream. $[kg/m^3]$
- $P_{limit}$ is the pressure limit of the upstream absolute pressure. $[bara]$
- $P_{up}$ is the absolute pressure upstream of the orifice. $[bara]$
- $k$ is the ratio of the heat capacities at constant pressure, $C_p$, and at constant volume, $C_v$.
- $P_{down}$ is the absolute pressure downstream of the orifice. $[bara]$
- $P_{used}$ is the pressure used in the mass flow equation based on choked or non-choked conditions. $[bara]$
- $\dot{m}_{flow}$ is the mass flow through the orifice. $[kg/s]$
- $C_d$ is the discharge coefficient of the orifice opening. $[-]$
- $A$ is the cross sectional area of the orifice. $[m^2]$

### Pressure safety valve / Relief valve
A PSV / relief valve is a mechanical device actuated by the static pressure in the vessel and a conventional PSV is often used for gas/vapor systems.
A conventional PSV is a spring-loaded device which will activate at a predetermined opening pressure, and relieve the vessel pressure until a given reseat pressure has been reached.
Both the opening pressure and the reset pressure is above the vessel operating pressure and the PSV will remain closed until the pressure inside the vessel increases to the opening pressure.
The operation of a conventional spring-loaded PSV is based on a force balance.
A conventional PSV can be seen in [@Fig:psv].
A spring exerts a force on a disc blocking the inlet of the PSV.
When the pressure inside the vessels reaches the opening pressure, the force exerted on the disc by the gas will be larger than the force exerted by the spring and the PSV will open and allow the gas to flow out of the vessel.
The flow of gas out of the vessels will lower the pressure and thereby also the force exerted on the disc.
When the pressure in the vessels is reduced to the reset pressure, the PSV will close and the disc will again hinder the gas flow.

![Conventional/pop action PSV adapted from [@iskov] and [@API520]](docs/img/PSV.pdf){#fig:psv width=30%}

The relief valve model implemented in HydDown is the API 520 equations [@API520] for gas relief for both sonic/critical as well as subcritical flow. No corrections factors are implemented in HydDown.

For sonic flow (critical flow), as indicated in equation, the  mass flow through the PSV can be determined by equation [@eq:Stationary_sizing_sonic].

$$ W = \frac{A {C \cdot K_d \cdot  K_b \cdot  K_c \cdot  P_1}}{\sqrt{\frac{T\cdot Z}{M}}} $$ {#eq:Stationary_sizing_sonic}

- A is the effective discharge area. [mm$^2$]
- W is the mass flow through the device. [kg/h]
- C is a coefficient, as a function of k, as defined in equation [@eq:Stationary_sizing_C_function].
- $K_d$, $K_b$, and $K_c$ are correction factors.
- $P_1$ is the allowable upstream absolute pressure. [kPa]
- T is the temperature of the inlet gas at relieving conditions. [K]
- M is the molecular mass of the gas at relieving conditions. [kg/kmol]
- Z is the compressibility factor for the gas.

$K_d$ is the effective coefficient of discharge, with a typical value of 0.975, for an installed PSV.
$K_b$ is a back-pressure correction factor between 0 and 1, assumed to be 1.
$K_c$ is a correction factor used when a rupture disk is installed upstream, otherwise it is 1.
In the present implementation a value of 1 is assumed.

$$ C=0.03948 \sqrt{k \left (\frac{2}{k+1}\right)^{\left (\frac{k+1}{k-1}\right)}} $$ {#eq:Stationary_sizing_C_function}

For subsonic flow (subcritical flow), the effective discharge area of the PSV is determined by equation [@eq:Stationary_sizing_subsonic].

$$ W = \frac{A \cdot {F_2 \cdot K_d  \cdot  K_c }}{17.9 \sqrt{\frac{T\cdot Z}{M\cdot P_1\cdot (P_1-P_2)}}} $$ {#eq:Stationary_sizing_subsonic}

F$_2$ is the coefficient of subcritical flow which can be determined from [@eq:Stationary_sizing_F2].

$$ F_2=\sqrt{\left ( \frac{k}{k-1} \right ) r^{\left (\frac{2}{k} \right )} \left (  \frac{1-r^{\left ( \frac{k-1}{k} \right )}}{1-r}\right )} $$ {#eq:Stationary_sizing_F2}

where $r$ is the ratio of backpressure to upstream relieving pressure, $P_2 / P_1$.

When modelling a pop action PSV/relief valve under dynamic conditions, the valve will go from closed to fully open in a short period of time when the set pressure, $P_{set}$, is reached.
The pop action is illustrated in [@Fig:psv_hyst] which shows the opening and closing hysteresis of the PSV as a function of pressure.
In order to close the pressure shall be reduced below the reseat pressure.

![Relief valve hysteresis adapted from [@iskov]](docs/img/hysteresis.pdf){#fig:psv_hyst width=80%}

When specifying PSVs it common to use standard API sizes as shown in [@tbl:psv_sizes]

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
For calculating the mass flow through a control valve, the ANSI/ISA [@borden][@ISA] methodology also described in IEC 60534 [@IEC60534] is applied.

The flow model for a compressible fluid in the turbulent regime is:

$$ W = C N_6 F_P Y \sqrt{x_{sizing} p_1 \rho_1} $$

or equivalently:

$$ W = C N_8 F_P p_1 Y \sqrt{\frac{x_{sizing}M}{T_1 Z_1}} $$

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

In HydDown, the piping geometry factor is not yet implemented and assumed to be 1.
The pressure drop ratio $x_{sizing}$ used for sizing is determined as the lesser of the actual pressure drop ratio, $x$, and the choked pressure drop ratio $x_{choked}$.
The actual pressure drop ratio is given by:

$$ x \frac{\Delta p}{p_1}$$

The pressure drop ratio at which flow no longer increases with increased value in pressure drop ratio is the choked pressure drop ratio, given by the following equation:

$$ x_{choked} = F_\gamma x_{TP} $$

The factor $x_T$ is based on air near atmospheric pressure as the flowing fluid with a specific heat ratio of 1.40.
If the specific heat ratio for the flowing fluid is not 1.40, the factor $F_\gamma$ is used to adjust $x_T$.
Use the following equation to calculate the specific heat ratio factor:

$$ F_\gamma = \frac{\gamma}{1.4} $$

where $\gamma$ is the ideal gas $C_p/C_v$.
It should be noted that the above equation has been derived from perfect gas behaviour and extension of an orifice model with $\gamma$ in the range of 1.08 to 1.65.
If used outside the assumptions flow calculations may become inaccurate.

The expansion factor Y accounts for the change in density as the fluid passes from the valve inlet to the vena contracta.
It also accounts for the change in the vena contracta area as the pressure differential is varied.

$$ Y = 1 -  \frac{x_{sizing}}{3x_{choked}}$$


## Heat transfer {#sec:heat}

### Natural convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to depressurization can be well approximated by that of natural convection as found from measured Nusselt numbers being well correlated with Rayleigh number, with no apparent improvement in model performance by inclusion of the Reynolds number in the model [@woodfield].

To determine the heat transfer for the gas-wall interface, the following is applied cf. equation [@Eq:newton]:

$$  \frac{dQ}{dt} = \dot{Q} = h A ( T_{s} - T_{gas} ) $$  {#eq:newton}

- $d Q$ is the change in thermal energy due to convective heat transfer. [J]
- $d t$ is the change in time during the heat transfer. [s]
- $h$ is the convective heat transfer. [W/m$^2$ $\cdot$ K ]
- $A$ is the area normal to the direction of the heat transfer. [m$^2$]
- $T_{s}$ is the inner surface temperature of the geometry. [K]
- $T_{gas}$ is the temperature of the bulk gas inside the vessel. [K]

The convective heat transfer will need to be estimated for the the gas-wall interface, by the use of empirical relations for the Nusselt number.
The Nusselt number describes the ratio of convective heat transfer to conductive heat transfer, normal to a surface area, as given in equation [@eq:Nu].

$$ Nu=\frac{hL}{k} $$ {#eq:Nu}

- $Nu$ is the Nusselt number. [-]
- $h$ is the convective heat transfer. [W/m$^2$$\cdot$K]
- $L$ is a characteristic length of the geometry. [m]
- $k$ is the thermal conductivity of the gas. [W/m$\cdot$K]

The characteristic length $L$ used is the height of the gas volume i.e. the length of a vertical vessel or the diameter of a horizontal vessel.

The empirical correlations used to calculate the Nusselt number of the gas-wall interface is a function of the Rayleigh number, which can be defined by the Grashof number and Prandtl number, as in equation [@Eq:rayleigh_gas]:

$$ Ra=Gr \cdot Pr $$ {#eq:rayleigh_gas}

- $Ra$ is the Rayleigh number. [-]
- $Gr$ is the Grashof  number. [-]
- $Pr$ is the Prandtl number. [-]

The Grashof number is a dimensionless number which approximates the ratio of the buoyancy forces to viscous forces, as given in equation [@Eq:grashof_gas}:

$$ Gr=\frac{\beta g\rho^2 L^3 \Delta T }{\mu^2} $$ {#eq:grashof_gas}

The Prandtl number is a dimensionless number defined as the ratio of the momentum diffusivity to thermal diffusivity, as given in equation [@Eq:prandtl_gas]:

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

The mechanism of heat transfer on the outside of the vessel at ambient conditions is similar to the above.
In HydDown the external heat transfer is not modelled currently.
A heat transfer coefficient shall be provided.

### Mixed convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to filling can be well approximated by that of combined natural convection and forced convection as found from measured Nusselt numbers being well correlated with Rayleigh and Reynolds number [@woodfield].

For mixed convection the effective Nusselt number, $Nu$, can be approximated by

$$ Nu = (Nu_{forced}^n + Nu_{natural}^n)^{\frac{1}{n}} $$

During charging with different gases (H$_2$, N$_2$ and Argon), Woodfield *et al.* demonstrated that in order to provide a good fit to the experimentally determined Nusselt number a correlation based on both Reynolds and Rayleigh number was necessary.
The following formula was found to provide a good fit (n=1):

$$ Nu = Nu_{forced} +  Nu_{natural} = 0.56Re_d^{0.67} + 0.104Ra_H^{0.352} $$

### Conduction
For accurate prediction of the outer and especially the inner wall temperature for correct estimation of internal convective heat transfer and the average material temperature, the general equation of 1-D unsteady heat transfer shall be solved:

$$ \frac{\delta T}{\delta t} = \frac{k}{\rho c_p} \frac{\delta^2 T}{\delta x^2} $$

- T is temperature
- x is the spatial (1-D) coordinate
- k is the thermal conductivity
- $\rho$ is the material density
- $C_p$ is the specific heat capacity

To be solved, the initial values and boundary values must be specified.
In its default state (if thermal cobductivity is not applied for the vessel), HydDown does not include the unsteady heat transfer model, i.e. the assumption is that the temperature from outer to inner surface is uniform and equal to the average temperature.
This is obviously a crude approximation, but might be justified depending in the Biot number:

$$ Bi = \frac{hL}{k} $$

The Biot number is a simple measure of the ratio of the heat transfer resistances at the surface of a body to the inside of a body.
The ratio gives an indication to which extent the temperature will vary in space (gradient) when the body is subject to a displacement in temperature at the surface boundary layer.
Striednig *et al.* [@STRIEDNIG] concluded that for a type I (steel) cylinder the Biot number was approx. 0.03 and hence the error in assuming a uniform temperature in the vessel wall was low.

With a typical thermal conductivity of 45 $W/m K$ for steel and a heat transfer coefficient up to 600 $W/m^2 K$ [@woodfield] the Biot number for a vessel with a wall thickness of 2 cm is 0.27.
This is significantly higher than that approximated by [@STRIEDNIG].
However, the Biot number is significantly lower than 1, and the assumption of a uniform temperature is reasonable.
However, for increased wall thickness, and/or for different materials with lower thermal conductivity, the error may grow to an unacceptable level.

Especially for vessels with low conductivity materials (or very thick walls) accurate estimation of the vessel wall temperatures requires the 1-D transient heat transfer problem to be solved. HydDown incorporates the [*thermesh*](https://github.com/wjbg/thermesh) code provided under an MIT license by Wouter Grouve [@thermesh]. The implemented model applies Cartesian coodinates as also applied in the h2fills program by NREL [@KUROKI], i.e. the curved vessel wall is assumed a flat plate. This may to some extent be justied by the relatively large radius of a cylindrical storage container compared to the vessel wall thichness (see also justification references in [@KUROKI]).

In HydDown, the 1-D transient heat conductivity problem is solved with *thermesh* assuming temperature independent vessel material density, heat capacity and thermal conductivity applying the Crank-Nicholson scheme using Neumann boundary conditions:

$$  -k_{\text{z}}\frac{\partial T}{\partial z}\Biggr|_{z=0} = q_\text{left}(t), \qquad
    -k_{\text{z}}\frac{\partial T}{\partial z}\Biggr|_{z=L} = q_{\text{right}}(t) $$

The heat flux at the outer (left) and inner (right) is calculated/updated for each major time step for the explicit solver for mass and energy balances. For each major time step *thermesh* is used to update the outer and inner wall temperature by solving the 1-D transient heat conductivity using a *minor* time step equal to 1/10th of the major time step and 11 nodes. For each time *thermesh* is called during major time steps the mesh is initialised with the temperature profile at the end of the former major time step.

### Fire heat loads
The heat transfer from the flame to the shell is modelled using the recommended approach from Scandpower [@scandpower] and API [@API521].
The heat transfer from the flame to the vessel shell is divided into radiation, convection, and reradiation as seen in equation [@Eq:flame]:

$$ q_f=\underbrace{{\alpha}_s \cdot {\varepsilon}_f \cdot \sigma \cdot T_f^4}_\text{Radiation}+\underbrace{h_f \cdot (T_f-T_s(t))}_\text{Convection}-\underbrace{{\varepsilon}_s \cdot \sigma \cdot T_s(t)^4 }_\text{Reradiation} $$ {#eq:flame}

- $q_f$ is the flame heat flux. [W/m$^2$]
- ${\alpha}_s$ is the vessel surface absorptivity. [-]
- ${\varepsilon}_f$ is the flame emissivity. [-]
- $\sigma$ is the Stefan-Boltzmann constant, $\sigma$ = $5.67 \cdot 10 ^ {-8}$  [W/m$^2 \cdot$ K$^4$]
- $T_f$ is the flame temperature. [K]
- $h_f$ is the convection heat transfer coefficient between the flame and the surface. [W/m$^2 \cdot$ K]
- $T_s(t)$ is the time dependent surface temperature. [K]
- ${\varepsilon}_s$ is the surface emissivity. [-]

This model assumes that the pressure vessel is fully engulfed by the flame.
This means that the view factor for the radiation is unity and is therefore not taken into consideration.
The convective heat transfer coefficients for a jet fire and a pool fire, and recommended values for the emissivity and absorptivity, are given by Scandpower as [@scandpower]:

- $h_{jet~fire}$ = 100 [W/m$^2 \cdot$K]
- $h_{pool~fire}$ = 30 [W/m$^2 \cdot$K]
- ${\alpha}_s$ = 0.85
- ${\varepsilon}_s$ = 0.85
- ${\varepsilon}_f$ = 1.0 (optical thick flames, thickness > 1 m)

The flame temperature is found by solving equation [@Eq:flame2] for the incident heat flux in relation to the ambient conditions.
The flame temperature is kept constant throughout the simulation:

$$ q_{total}=\sigma \cdot T_f^4 + h_f \cdot (T_f-T_{amb})$$ {#eq:flame2}

- $q_{total}$ is the incident flame heat flux as given in table [@Tbl:heatfluxes1]. [W/m$^2$]
- $T_{amb}$ is the ambient temperature $\approx$ 293 K (20$^\circ$ C)

The heat flux used to calculate the flame temperature is given in table [@tbl:heatfluxes1].

|                        | Small jet fire  [kW/m$^2$]  |  Large jet fire  [kW/m$^2$] |  Pool fire  [kW/m$^2$]
| ----                   |  ----           |  ----           | ----
| Peak heat load         |  250            |  350            | 150         
| Background heat load   |   0             |   100           |  100         

: Incident heat fluxes for various fire scenarios given by Scandpower [@scandpower] {#tbl:heatfluxes1}


## Model implementation
A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation as described in [@Sec:flow]:

$$ m_{cv}(i+1) =  m_{cv}(i) + \dot{m}(i) \Delta t  $$ {#eq:euler_mass}

$$\dot{m}(i) = f(P,T,) $$

For each step, the mass relief/ left in the vessel is known.
Since the volume is fixed the mass density is directly given.

For the calculation methods (isentropic, isenthalpic, isenergetic, etc.), CoolProp allows specifying density and either H, S, or U directly - this is very handy and normally only TP, PH, or TS property pairs are implemented and you would need to code a second loop to make it into am UV, VH, or SV calculation.

In the following, the high-level steps in the solution procedure are outlined for the different calculation types.
To illustrate calls to the equation of state in CoolProp, the notation $EOS(D,T)$ corresponds to calculation of the state at known density and temperature, for example.

### Isothermal process
For an isothermal process, the solution procedure for each calculation step is the following with the mass in the control volume being calculated from [@Eq:euler_mass]:

$$ T(i+1) = T(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ P(i+1) = EOS(D(i+1),T(i+1)) $$

The notation $EOS(T,P)$ refers to a call to CoolProp either via PropsSI or the AbstractState solving the equation of state for specified temperature and pressure. Different state pairs are applied includeing combinations of density (D), specific enthalpy (H), specific entropy (S), specific internal energy (U).  

### Isenthalpic process
For an isenthalpic process:

$$ H(i+1) = H(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ P(i+1) = EOS(D(i+1),H(i+1)) $$
$$ T(i+1) = EOS(D(i+1),H(i+1)) $$

### Isentropic process
For an isentropic process:

$$ S(i+1) = S(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V}$$
$$ P(i+1) = EOS(D(i+1),S(i+1)) $$
$$ T(i+1) = EOS(D(i+1),S(i+1)) $$

### Isenergetic process
For an isenergetic process:

$$ U(i+1) = U(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V}$$
$$ P(i+1) = EOS(D(i+1),U(i+1)) $$
$$ T(i+1) = EOS(D(i+1),U(i+1)) $$

### Energy balance
The general first law applied to a flow process as outlined in [@Sec:firstlaw] subject to an explicit Euler scheme is:

$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ U_{cv}(i+1) = \frac{m_{cv}(i)U_{cv}(i) - \left( \dot{m}(i) H (i) +  \dot{Q}(i) \right) \Delta t}{m_{cv}(i+1)}  $$ {#eq:firstlaw_euler}

The above assumes that mass flow is positive when leaving the control volume and heat rate is positive when transferred to the control volume.
$H(i)$ is the specific enthalpy of the fluid in the control volume for a discharging process and it is equal to the energy of the entering stream for a filling process.
The heat rate is calculated as outlined in [#Sec:heat].

For the vessel wall a simple heat balance is also made:

$$ T_{wall}(i+1) = T_{wall}(i)  + \frac{\dot{Q}_{outer} - \dot{Q}_{inner} } {c_p m_{vessel}} \Delta t $$

where $\dot{Q}_{outer}$ is the convective heat transfer to or from the ambient surroundings from the outer surface of the vessel, with positive values indicating that heat is transferred from the surroundings to the vessel wall.
This is either a fixed heat transfer coefficient with a specified ambient temperature (outer surface) or a calculated fire heat load.
$\dot{Q}_{inner}$ is the internal heat transfer, either a fixed number or calculated as natural convection (for discharge) or mixed natural and forced convection (filling).

When the temperature profile of the vessel wall is calculated using the 1-D transient conductivity equation the vessel temperature as calculated above represent the mean vessel temperature and addition to this the outer and inner wall temperatures are calculated.

If a fixed $\dot{Q}$ is provided or a fixed overall heat transfer coefficient is provided the vessel wall temperature heat balance is not solved.

The remaining steps are update of temperature and pressure:

$$ P(i+1) = EOS(D(i+1),U(i+1)) $$
$$ T(i+1) = EOS(D(i+1),U(i+1)) $$

## Multicomponent mixtures
Although not an initial requirement, the code can handle multi-component gas mixtures. When calculations are done for multicomponent mixtures, the code runs significantly slower.
Please note, that in case that liquid condensate is formed during discharge calculations and even if the calculations does not stop, the results cannot be considered reliable.
This is because component balances are not made and it is always assumed that the discharge composition is the same as the global composition inside the vessel.
When liquid condensate is formed the composition of the vapour phase will differ from the global composition.

There are a few examples of multicomponent mixtures included with HydDown.
In order to specify multicomponent mixtures the below example can be used as guidance:

    "Methane[9.01290129e-01]&Ethane[6.35063506e-02]&N2[7.80078008e-03]&CO2[2.34023402e-02]&Propane[3.50035004e-03]&Butane[5.00050005e-04]"

For component names please refer to the [Coolprop manual](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids).

An example calculation for the above mixture is shown in [@Fig:ng_all].

The pressure and temperature trajectory is visualised along with the fluid mixture phase envelope in [@Fig:ng_pe].
As seen from the figure, this case is borderline and the pressure/temperature trajectory just coincides with the dew line on the phase envelope.
This plot is included in the example main script that comes with HydDown and serves as an important quality control.

For multi-component mixtures CoolProp does not provide solver for PH and UV-problems and these solvers have been included in HydDown by wrapping the CoolProp PT solver with a Nelder-Mead algorithm for solving for internal energy and density or pressure and enthalpy.

![Gas discharge of multicomponent mixture](docs/img/ng_all.png){#fig:ng_all width=90%}

![Multicomponent mixture phase envelope and pressure and temperature trajectory of the vessel inventory during discharge](docs/img/ng_pe.png){#fig:ng_pe width=70%}

# Validation
The code is provided as-is.
However, comparisons have been made to some experiments from the literature.

The following gases and modes are considered:

- High pressure hydrogen discharge
- High pressure hydrogen filling
- High pressure nitrogen discharge
- High pressure multi-component gas discharge

These cases are all for Type I (steel) cylinders and the simple heat transfer / heat balance model has been applied, ignoring the thermal gradient in the cylinder wall. It has been checked by applying the 1D transient heat conduction model that the results are basically unaffected by this assumption. In addition to these, a few cases for Type IV cylinders have also been included in order to validate the 1D transient heat conduction model applied. This includes

- High hydrogen discharge - comparison with commercial software.
- High pressure helium discharge experiment (Karlshuhe Institute of Technology)
- High pressure hydrogen discharge (GasTeF experiment)

Details on Type III experiments are sparse but some general observations are made for the difference between type III and type IV cylinders and a single type III filling experiment is modelled.

Furthermore, HydDown has also been benchmarked against external code by Ruiz and Moscardelli [@maraggi], who compared their GeoH2 app against HydDown and found excellent agreement for both pressure, mass flow rate and gas temperature for three different cases of discharge and filling.

## Hydrogen discharge (Type I)
Calculations with HydDown for high pressure hydrogen vessel discharge have been compared to three experiments from [@byrnes].
Byrnes *et al.* [@byrnes] have reported a number of rapid depressurization experiments using both nitrogen and hydrogen, with blowdown times being varied between 14 seconds to 33 minutes.
The vessel used was a type "K" gas cylinder with a volume of 0.0499 m$^3$, an ID of 8.56 inches and shell length of 55 inches.
The vessel is placed in a vertical position.
A control valve is mounted on the neck of the cylinder with a fixed opening percentage during each experiment, set to provide the required depressurization rate.
It is mentioned that the gas cylinder is placed inside a weather proof cabinet with perlite insulation in order to emulate adiabatic conditions.
However the simulations run here are with a finite heat transfer coefficient on the outside of the cylinder.
In their paper, detailed pressure and temperature traces are provided for three experimental runs: 7, 8, and 9, with depressurization from 2000 psia in 30 seconds, 480 seconds and 14 seconds, respectively.
These experiments are simulated in [@Fig:byrnes_run7], [@Fig:byrnes_run8] and [@Fig:byrnes_run9].

![Calculations of hydrogen discharge experiment run 7 from [@byrnes]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](docs/img/Byrnes_run7.png){#fig:byrnes_run7}

The figure layout is general for all the validation studies conducted and will be described in brief for eased readability.
Each experiment figure is divided in 4 subplots: the upper left shows the simulated bulk gas temperature and vessel wall (average) temperature with measured values included with stipulated lines (if available), the upper right shows the simulated and measured pressure (if available, shown with points), the lower left shows the specific enthalpy, specific internal energy, and the specific entropy of the bulk gas content as a function of time, the lower right shows the simulated mass rate leaving or entering the vessel.

As seen from [@Fig:byrnes_run7], [@Fig:byrnes_run8], and [@Fig:byrnes_run9], the calculations generally compare well with the experimental results.
It appears as though the faster depressurization times are predicted the best, as seen from run 7 and 9.
The much slower run 8 does not have as good a fit for the bulk gas temperature.
However, on an absolute scale the results are still within 5$^\circ$C of the experimentally determined gas temperature.
The trend in *goodness of fit* could indicate that the heat transfer from the outside of the vessel to the gas lacks some details, which become detectable at the very slow depressurization.
For the fast depressurization the heat transfer from the outside is of less importance and the inside heat transfer coefficient may be predicted better.

![Simulation of run 8 from [@byrnes].](docs/img/Byrnes_run8.png){#fig:byrnes_run8}

![Simulation of run 9 from [@byrnes].](docs/img/Byrnes_run9.png){#fig:byrnes_run9}

## Hydrogen filling (Type I)
In order to compare HydDown to experimental values of hydrogen filling operation the experiments reported by Striednig *et al.* [@STRIEDNIG] are used.
Experimental results were obtained using a type I tank (steel) with a volume of 0.0235 m$^3$.
The filling of hydrogen of gas into the vessel was carried out with an upstream reservoir kept at 350 bar and the flow was controlled with an electronically controlled dispenser.
Vessel details are provided below:

| Parameter       | Value       |
|-----------------|-------------|
| Mass            | 69..5 kg    |
| Length          | 0.61 m      |
| OD              | 0.280 m     |
| Wall thickness  | 0.0129 m    |
| Density         | 7740 kg/m$^3$ |
| Heat capacity   | 470 J/(kg K) |

: Key hydrogen vessel data [@STRIEDNIG] {#tbl:h2filling}

[@STRIEDNIG] reports many different experiments.
In this comparison we will use experiments applying different pressurisation rates / durations using identical initial conditions.
The experiments are reported in Fig. 6 and Fig. 7 in [@STRIEDNIG].

In its current stage HydDown does not directly support specifying a fixed pressurisation rate expressed as a constant rate of change in pressure / pressure ramp rate (PRR).
In order to emulate the experiments, an orifice model is used with an artificially high reservoir pressure in order to ensure that the flow is critical during the entire pressurisation.
This gives a constant mass rate and translates to a fairly constant pressure ramp rate, although with real gas effects being clearly observable.

The results of the simulations with HydDown and the experiments from [@STRIEDNIG] are shown in [@Fig:striednig5], [@Fig:striednig10], [@Fig:striednig30].
As seen from the figures the experimental gas temperature appears to be simulated very accurately.
The simulation with the fastest pressurisation under-predicts the experimental temperature slightly.
The main features of the pressurisation temperature profile is a rapid increase in temperature, at the beginning of the pressurisation, where the effects of filling dominate.
At approximately 20-30 s, the temperature increase levels off to a plateau, where the length of the plateau depends on the pressurisation rate.
The slower the depressurization, the longer the plateau.
During this temperature plateau, the heat transfer from the gas to the vessel and surroundings balance the temperature increase that would occur if performed at adiabatic (isolated) conditions.
When the gas flow has ceased upon reaching the final pressure, the temperature starts to drop due to cooling of the gas by the colder vessel wall driven by convective heat transfer.
The temperature stabilises at a temperature which is higher than the initial ambient temperature due to the heating of the steel vessel.
The cooling of the vessel by ambient air is slower and would require a longer run time to be clearly visible.

![Simulation of $H_2$ pressurization using 5 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_5MPa_min.png){#fig:striednig5}

![Simulation of $H_2$ pressurization using 10 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_10MPa_min.png){#fig:striednig10}

![Simulation of $H_2$ pressurization using 30 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_30MPa_min.png){#fig:striednig30}

## Nitrogen discharge (Type I)
Calculations with HydDown are compared to experiment I1 from ref. [@Haque1992b].
The experiment is a blowdown of a vertically oriented cylindrical vessel with flat ends.
The vessel length is 1.524 m, the inside diameter is 0.273 m, and the wall thickness is 25 mm.
The vessel is filled with N$_2$ at 150 bar and  15$^\circ$C.
Ambient temperature is 15$^\circ$C.
The blowdown orifice diameter is 6.35 mm.
The results are shown in [@Fig:N2val].
The discharge coefficient of the orifice has been set to 0.8 in order to match the vessel pressure profile.
The back pressure is set to atmospheric conditions.


As seen from [@Fig:N2val], the calculations compare well with the experimental results.
The calculated temperature of the bulk vapor is within the experimental range of measured temperature at all times during the simulation.
It is also noted that the minimum temperature is reached at approximately the same time as in the experiments.
The calculated vessel inner wall temperature does not decline as rapidly as the experiments, but from around a calculation time of 60 s, the temperature is within the experimentally observed inner wall temperature.

![Calculations of nitrogen discharge emulating experiment I1 from [@Haque1992b]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](docs/img/N2_filling.png){#fig:N2val}


## Multi-component discharge (Type I)
The experiments were conducted by Haque *et al.* [@Haque1992b] and the experimental data has been extracted from the Ph.d. thesis of Wong [@WONG].

The experimental vessel is a full-size suction scrubber for a gas compressor which has a total volume of 2.78 m3, with length and inside diameter of 3.240 m (2.75 m tan-to-tan)
and 1.130 m respectively. The wall thickness is 59 mm. For modelling the length of a flat-ended cylinder has been adjusted to give a total volume of 2.78 m3. An orifice diameter of 6.35 mm was used for discharge and for HydDown modelling the discharge coefficient was adjusted to 0.97 to match the experimental pressure. The results of HydDown and comparison with the experimental results are shown in [@fig:Multi-val].

![Calculations of discharge of a non-condensable mixture (0.91/0.09 molefraction of methane/ethane) emulating experiment I1 from [@Haque1992b]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](docs/img/NG_validation.png){#fig:Multi-val}

The minimum gas temperature simulated by HydDown is -6$^\circ$C, compared to the minimum measured average gas temperature of
-10.3$^\circ$C. The difference in measured vs calculated wall temperature is 2.3$^\circ$C, with a lower measured wall temperature.  


## 1-D transient heat transfer
A few notes about vessels with poor thermal conductivity and composite materials. When modelling systems with high Biot number and in particular composite materials the complexity increases significantly. Not just because of the more difficult numerical problem, but even more so because of uncertainty in key parameters such as thermal conductivity, density and heat capacity. Composite materials such as carbon fibre or glass fibre reinforced epoxy systems can be manufactured in many different ways (fibre orientation etc.) which effects the previously mentioned properties. These properties, in particular the thermal conductivity, will influence the results significantly. If these properties are not accurately informed for the system to be analysed, sourcing data from literature shall be done with caution.

### Validation against commercial simulation tool
The first validation case is made against the commercial tool Honeywell Unisim Design in dynamics mode. The initial conditions and vessel details are summarised below cf. [@tbl:1D-valid-1].


| Parameter       | Value       |
|-----------------|-------------|
| Length          | 6.4 m      |
| ID              | 0.619 m     |
| Wall thickness  | 0.024 m    |
| Density         | 2000 kg/m$^3$ |
| Heat capacity   | 962 J/(kg K) |
| Thermal conductivity | 0.5 W/(m K) |
| Outside $h$ | 2 W/(m$^2$ K) |
| Initial pressure| 182 bar|
| Initial temperature| 6 $^\circ$C |
|Discharge mass flow | 0.02 kg/s|
| Gas | Hydrogen |

: Key vessel data and intial conditions {#tbl:1D-valid-1}

The comparison between HydDown results and the corresponding simulations using Unisim Design Dynamics are shown in Figure [@fig:unisim_val]: As seen from the results the results obtained using HydDown closely resembles the Unisim Design results. The outer wall temperature is matched perfectly and the final gas temperature and inner vessel temperature deviates marginally by 1.5 and 1.9 $^\circ$C, respectively. Unisim predicts a sligtly lower gas temeprature and HydDown predicts a slightly lower inner wall temperature.  

![Calculations of vessel wall temperature (inner/outer) with 1D transient heat conduction during hydrogen discharge for comparison with simulation performed with Unisim Design Dynamics. The figure shows calculated gas and wall temperature (lines) compared to Unisim simulations (left) with Unisim results shown with points, calculated and Unisim simulation pressure (denoted "Experimental") (right).](docs/img/unisim_validation.png){#fig:unisim_val}


### Validation against KIT experiment (Type IV)
Using data from Molkov *et al.* [@MOLKOV1][@MOLKOV2] experiments performed at HYKA-HyJet research facility at Karlsruhe Institute of Technology (KIT) for a 19 liter type IV pressure vessel are simulated. The storage vessel was initially charged to 700 bar with helium gas and then cooled down to a normal room temperature (293 K ) before start of the depressurisation. Tank characteristics were not available and the required parameters were extracted from a similar tank by Molkov *et al.*, see references within refs. [@MOLKOV1][@MOLKOV2]. Discharge was thorugh 1 mm nozzle and a discharge coefficient of 0.9 was applied in HydDown as in the study by Molkov *et al.*. The initial conditions and vessel details are summarised below cf. [@tbl:1D-KIT].

| Type IV tank        |                |
|---------------------|----------------|
| External length     | 904 mm         |
| ID                  | 180 mm         |
| OD                  | 228 mm         |
| **HDPE liner**      |                |
| Thickness           | 7 mm           |
| Density             | 945 kg/m$^3$   |
| Heat capacity       |1584 J/(kg K)   |
| Thermal conductivity|0.385 W/(m K)   |
| **CFRP shell**      |                |
| Thickness           |17 mm           |
| Density             |1360 kg/m$^3$   |
| Heat capacity       | 1020 J/(kg K)  |
| Thermal conductivity| 0.5 W/(m K)    |
| **Initial conditions** |             |
| Outside $h$         | 8 W/(m$^2$ K)  |
| Initial pressure    |  700 bar       |
| Initial temperature | 20 $^\circ$C   |
| Gas                 | Helium         |

: Key vessel data and intial conditions {#tbl:1D-KIT}


In HydDown the type IV pressure vessel geometry as assumed that of a flat ended cylinder, with the length adjusted to give a total inventory volume of 19 liter.  In order to simulate the system a lumped heat capacity and density for the type IV cylinder is made. The thermal conductivty is set to that of the liner material, and the outer heat transfer coefficient is set to 8 W/(m$^2$K) as applied in the H2fills software [@KUROKI]. For a type IV vessel, the application of an average value for density, heat capacity and thermal conductivity can be an acceptable approximation since the liner and composite shell material are not too dissimilar.

![Calculations of vessel wall temperature (inner/outer) with 1D transient heat conduction during helium discharge for comparison with KIT experiments. The figure shows calculated gas and wall temperature (lines) compared to experimental gas temperature (left) and calculated and measured pressure (right).](docs/img/KIT_1.png){#fig:KIT_val}

HydDown simulations are compared with the KIT experiement in Figure [@fig:KIT_val]. As seen from the results the depressurisation pressure is matched very well. The experimental gas temperature is also matched fairly well. The lowest simulated gas temperature is 182.3 K and the lowest experimental temperature is 177.5 K a difference of 4.8 K. The minimum gas temperature occurs at an earlier time compared to the experimental results. The former at approx. 75 sec. and the latter at approx. 100 sec. This is also in agreement with the simulation model presented by Molkov *et al.* [@MOLKOV1]. The final gas temperature of the simulations is 237 K compared to the experimental value of 216 K.

Molkov *et al.* [@MOLKOV1] also made a thermal analysis of the thermocouple arrangement used in KIT experiment in order to estimate thermal lag in the temperature measurement. Incorporating the thermal model of the thermocouple arrangement displayed an improved prediction of the time of the minimum measured gas temperature as well as the final measured temperature.


## Validation against GasTeF experiments (Type IV)
The experimental setup is described in more detail in refs. [@ACOSTA][@DEMIGUEL]. Different tank types where tested both for filling and discharge including a 29 liter type IV tank and a 40 liter type III tank. In the following the type IV tank is modelled. Key details such as the initial conditions and vessel details are summarised below cf. [@tbl:1D-GasTeF]. The geometry and thermal properties of the liner and composite shell material was not informed in the cited references and it is assumed to be similar to the KIT vessel descibed in the previous section, with the liner and composite shell thickness scaled to match the total thickness of the GasTeF type IV vessel.  

| Type IV tank        |                |
|---------------------|----------------|
| External length     | 827 mm         |
| ID                  | 230 mm         |
| OD                  | 279 mm         |
| Liner               | HDPE           |
| Composite shell     | G&CRPE         |
| **Initial conditions** |             |
| Initial pressure    |  700 bar       |
| Initial temperature | 25$\pm 2^\circ$C   |
| Gas                 | Hydrogen       |
| Discharge rate      | 1.8 g/s        |

: Key vessel data and initial conditions. The mass flow during discharge is constant at 1.8 g/s until the pressure drops below 5 MPa after which it drops {#tbl:1D-GasTeF}

The simulations with HydDown and comparison against the GasTeF experient is shown in Figure [@fig:GasTeF_val]. The experimental setup included several thermocouples mounted in different positions inside the test vessel, both at the center line and nearer the top and bottom. The experimental setup did not include a direct measurement of the internal liner temperature interface towards the gas nor the composite shell. The experimental points for the gas temperature as shown in the figure includes the lowest measured temperatures (near the bottom), the gas temperature measured in the middle and towards the top (highest temperatures). The experiments display significant temperature stratification during discharge experiments.


![Calculations of vessel wall temperature (inner/outer) with 1D transient heat conduction during hydrogen discharge for comparison with GasTeF experiments. The figure shows calculated gas and wall temperature (lines) compared to experimental gas temperature.](docs/img/demiguel.png){#fig:GasTeF_val width=80%}

As seen the prediction of the external surface temperature of the composite shell is matched very well. The average gas temperature as calculated with HydDown matches the temperatures recorded in the lower half of the test vessel. Further, the minimum in gas temperature between 450 and 500 seconds are also matched well. The measured pressure was not reported in ref. [@DEMIGUEL]

## Validation against Type III cylinder discharge experiments
As shown previously type IV cylinder thermal behaviour is matched fairly well when the liner/composite thermal material properties are lumped. For a Type III vessel which has very dissimilar thermal properties of the liner and the composite shell, this approach is not advisable. Luckily, the implemented 1D transient heat conduction code *thermesh* allows modelling bi-materials and this is also implemented in HydDown.

Unfortunately, it has not been possible to find reported discharge experiments which include detailed information about the exact liner and shell composite thickness.
In lack of a good and complete validation cases, an example is made using the KIT cylinder in the previous section and replacing  the HPDE liner with an aluminum liner.


| Type III tank        |                |
|---------------------|----------------|
| External length     | 904 mm         |
| ID                  | 180 mm         |
| OD                  | 228 mm         |
| **Aluminum liner**      |                |
| Thickness           | 7 mm           |
| Density             | 2700 kg/m$^3$   |
| Heat capacity       |900 J/(kg K)   |
| Thermal conductivity|237 W/(m K)   |
| **CFRP shell**      |                |
| Thickness           |17 mm           |
| Density             |1360 kg/m$^3$   |
| Heat capacity       | 1020 J/(kg K)  |
| Thermal conductivity| 0.5 W/(m K)    |
| **Initial conditions** |             |
| Outside $h$         | 8 W/(m$^2$ K)  |
| Initial pressure    |  700 bar       |
| Initial temperature | 20 $^\circ$C   |
| Gas                 | Helium         |

: Key vessel data and initial conditions for Type III cylinder simulations. The mass flow during discharge is constant at 1.8 g/s until the pressure drops below 5 MPa after which it drops. {#tbl:1D-TypeIII}

The results of HydDown simulation with the *artificial* type III cylinder is shown in Figure [@fig:KIT_typeIII]. As seen the immediate effect of replacing the HDPE liner with aluminum (compared to Figure [@fig:KIT_val]) is a lower temperature drop in the gas phase and a liner temperature which approached the gas temperature very closely. This is in general agreement with the observations in ref. [@DEMIGUEL].

![Calculations of vessel wall temperature (inner/outer) with 1D transient heat conduction during helium discharge for a hypothetical type III cylinder with dimensions similar to the KIT type IV cylinder.](docs/img/KIT_type_III.png){#fig:KIT_typeIII}

In Figures [@fig:figureIV] and [@fig:figureIII] the simulated temperature across the vessel liner/composite wall is shown for the real type IV cylinder from the KIT experiment and for the *artificial* type III cylinder, respectively. As seen there is a significant difference between the two bi-materials. The type III cylinder has almost zero temperature gradient across the aluminum liner, which can be rationalised by the very high thermal conductivity. Thus, the entire thermal gradient is over the composite shell material. These observations are in-line with simulation results from ref. [@MELIDEO20177304].

<div id="fig:figureIV">
![Temperature thorugh wall at different times](docs/img/KIT_II_tprofile2.png){#fig:KIT_typeIV_temp1 width=70%}
![Spatioal and temporal temprature contour](docs/img/KIT_II_tprofile1.png){#fig:KIT_typeIV_temp2 width=80%}

Calculations of type IV vessel wall temperature profile with 1D transient heat conduction. z=0 is the bonding interfance between liner and composite, z<0 is the liner and z>0 is the composite shell
</div>

<div id="fig:figureIII">
![Temperature thorugh wall at different times](docs/img/KIT_III_tprofile2.png){#fig:KIT_typeIII_temp1 width=70%}
![Spatioal and temporal temprature contour](docs/img/KIT_III_tprofile1.png){#fig:KIT_typeIII_temp2 width=80%}

Calculations of type III vessel wall temperature profile with 1D transient heat conduction. z=0 is the bonding interfance between liner and composite, z<0 is the liner and z>0 is the composite shell
</div>

## Validation against Type III cylinder filling experiments
For HydDown type III cylinder filling experimental validation, the work of Dicken and Mérida [@Dicken] is used. They conducted experiments of filling a 74 liter type III cylinder and also compared measurements with CFD calculations. Details of type III cyinder tested is shown in Table [@tbl:Dicken-TypeIII]. For HydDown simulation the length of a flat-ended cylinder has been adjusted to give a total cylinder volume of 74 liter. The time-dependent filling mass flow  has been sourced from ref. [@Dicken] and scaled in order to match the final experimental pressure.


| Type III tank        |                |
|---------------------|----------------|
| External length     | 893 mm         |
| ID                  | 358  mm         |
| OD                  | 396 mm         |
| **Aluminum liner**      |                |
| Thickness           | 4 mm           |
| Density             | 2700 kg/m$^3$   |
| Heat capacity       |900 J/(kg K)   |
| Thermal conductivity|167 W/(m K)   |
| **CFRP shell**      |                |
| Thickness           |15 mm           |
| Density             |938 kg/m$^3$   |
| Heat capacity       | 1494 J/(kg K)  |
| Thermal conductivity| 1.0 W/(m K)    |
| **Initial conditions** |             |
| Outside $h$         | 8 W/(m$^2$ K)  |
| Initial pressure    |  93 bar       |
| Initial temperature | 20.25 $^\circ$C   |
| Gas                 | H2         |

: Key vessel data and initial conditions for Type III cylinder simulations.  {#tbl:Dicken-TypeIII}

Simulation results from HydDown is compared to the experimental data of Dicken and Mérida in Figure [@fig:Dicken_typeIII]. As seen from the results the gas temperature calculated with HydDown is generally higher than the experimental results. The final calculated temperature is 74.2$^\circ$C, compared to the final measured temperature of 68.9$^\circ$C. The CFD calculations performed by Dicken and Mérida [@Dicken] was also higher than the experimental value (71.3$^\circ$C). Another CFD analysis of the investigated system by Hall and Ramasamy [@Hall] revealed a final mean gas temperature of 71.4$^\circ$C, and a zero-dimensional model by the same authors gave an end temperature of 72.4$^\circ$C.   

![Calculations of vessel pressure, gas temperature and wall temperature (inner/outer) with 1D transient heat conduction model during hydrogen filling. Comparison is made against experimental results from Dicken and Mérida [@Dicken] for pressure and gas temperature. Inner wall (liner) temperature is compared to CFD simulations [@Dicken] ](docs/img/dicken_typeIII.png){#fig:Dicken_typeIII}

In lack of measured vessel inner wall temperature (liner) the values calculated at the end of filling by Dicken and Mérida using their CFD model is used for comparison. The CFD simulations revealed very large variations in the heat transfer coefficient and resulting wall temperatures. For comparison with HydDown the average CFD calculated temperature is used. This is taken as the arithmetric average of all seven positions from front to back of the cylinder. As seen from Figure [@fig:Dicken_typeIII] the average value from HydDown (52.5$^\circ$C)compares very well with the  average CFD results at the end of filling (52.1$^\circ$C). However, estimated inner wall temperature spanned temperatures from 26 to 67 $^\circ$C, with the highest temperatures recorded opposite to the entry of hydrogen.   


# Similar software
To be written

# Future developmet of HydDown
To be written 

# Example use cases 
To be written
