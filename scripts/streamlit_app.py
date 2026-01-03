# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license


import streamlit as st
from PIL import Image
import os

st.set_page_config(
    page_title="Hyddown webapp",
    page_icon="👋",
)

with st.sidebar:
    try:
        image_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            "..",
            "docs",
            "img",
            "Sketch.png",
        )
        icon = Image.open(image_path)
        st.image(icon, use_container_width=True, caption="HydDown")
    except:
        pass

st.write("# Welcome to HydDown!")

st.markdown(
    """
    Hyddown is an open-source tool for hydrogen (or other pure species) pressure vessel filling and discharge calculations.
    The tool is based on rigorous thermodynamics and heat transfer models to provide accurate predictions of pressure, temperature, and other relevant parameters during the depressurisation/pressurisation process.

    **👈 Select a demo from the sidebar** to see some examples
    of what Hyddown can do!
    ### Want to learn more?
    - Check out [HydDown](https://github.com/andr1976/HydDown)
    - Jump into our [documentation](https://github.com/andr1976/HydDown/blob/main/docs/MANUAL.pdf)
    - Read the paper in JOSS: [![DOI](https://joss.theoj.org/papers/10.21105/joss.03695/status.svg)](https://doi.org/10.21105/joss.03695)

    ### About Hyddown
    Thermodynamics are implemented using the [**CoolProp** library](http://www.coolprop.org/) which provides high accuracy real gas properties for a wide range of fluids including hydrogen, nitrogen, helium, methane, carbon dioxide, and many others.
    See also the original [CoolProp paper](https://pubs.acs.org/doi/full/10.1021/ie4033999).

    Hyddown provides a number of different calculation models, including:
    - Isentropic models for quick estimates
    - Energy balance models for more accurate predictions
    - Heat transfer models to account for thermal effects such as vessel wall temperature    
    The tool is designed to be user-friendly and accessible, with a web-based interface that allows users to input their vessel and fluid parameters and obtain results quickly and easily.
    
    A number of advanced features are also available, including:
    - 1-D transient heat transfer models for vessel wall including dual layer walls for modelling e.g. composite vessels or insulated vessels using [**thermesh** library](https://github.com/wjbg/thermesh)
    - Customizable valve models for different flow conditions
    - External heat input for fire heat loads via Stefan-Boltzmann equation
    - Single component two-phase conditions can be managed for e.g. LPG tanks, CO2 storage vessels etc. Both for vessel response to fire or ambient boil-off calculations
    - Vessel rupture calculations based on von Misses stress criteria according to [Scandpower guidelines](https://www.iomosaic.com/diersweb/docs/Scandpower%20Fire%20Guidelines%20Version%202.pdf). 
    
   """
)
