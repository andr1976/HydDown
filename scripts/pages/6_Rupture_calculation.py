# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import streamlit as st
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import base64

try:
    from hyddown import HydDown
except:
    import sys
    import os

    hyddown_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "..", "..", "src"
    )
    sys.path.append(os.path.abspath(hyddown_path))
    from hyddown import HydDown


def get_table_download_link(df, filename):
    """
    Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(
        csv.encode()
    ).decode()  # some strings <-> bytes conversions necessary here
    filename = filename + ".csv"
    return f'<a href="data:application/octet-stream;base64,{b64}" download={filename}>Download csv file</a>'


def read_input():
    sideb = st.sidebar

    with sideb:
        try:
            icon = Image.open("../docs/img/Sketch.png")
            st.image(icon, use_column_width=True, caption="HydDown")
        except:
            pass

        with st.form(key="my_form"):
            submit_button = st.form_submit_button(label="Run calculation")
            heattran = True
            c1, c2 = st.columns(2)

            with c2:
                length = st.text_input("Vessel length (m):", 10.0)
                diam = st.text_input("Vessel diam (m):", 3.0)
                thk = st.text_input("Vessel thichness (m):", 0.050)
                orientation = st.selectbox(
                    "Vessel orientation", ("horizontal", "vertical")
                )
                orifice_diam = st.text_input("Orifice diam (mm):", 13)
                orifice_diam = float(orifice_diam) / 1000
                tstep = st.text_input("Time step (s):", 1.0)
                end_time = st.text_input("End time (s):", 1500)

            with c1:
                pres = st.text_input("Initial pressure (bar):", 100.0)
                pres = float(pres) * 1e5

                back_pressure = st.text_input("Back pres. (bar):", 1.013)
                back_pressure = float(back_pressure) * 1e5

                set_pressure = st.text_input("PSV S/P (bar):", 121.0)
                set_pressure = float(set_pressure) * 1e5

                blowdown = st.text_input("Blowdown (%):", 10.0)
                blowdown = float(blowdown) / 100.0

                fluid = st.selectbox("Select fluid", ("H2", "CH4", "N2", "air"))

                mode = "discharge"

                temp = st.text_input("Initial temp. (C):", 25)
                temp = float(temp) + 273.15

                fire_type = st.selectbox(
                    "Fire type",
                    ("api_jet", "api_pool", "scandpower_pool", "scandpower_jet"),
                )

            material = st.selectbox(
                "Select vessel material",
                ("CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"),
            )
            peak_fire_type = st.selectbox(
                "Select peak fire type",
                (
                    "scandpower_jet_peak_large",
                    "scandpower_jet_peak_small",
                    "scandpower_pool_peak",
                ),
            )

            density = st.text_input("Vessel material density (kg/m3):", 7740)
            density = float(density)

            cp = st.text_input("Vessel material heat capacity (J/kg K):", 470)
            cp = float(cp)

    input = {}
    input["calculation"] = {}
    input["vessel"] = {}
    input["initial"] = {}
    input["valve"] = {}
    input["heat_transfer"] = {}

    input["calculation"]["type"] = "energybalance"
    input["calculation"]["time_step"] = float(tstep)
    input["calculation"]["end_time"] = float(end_time)

    input["vessel"]["length"] = float(length)
    input["vessel"]["diameter"] = float(diam)
    input["vessel"]["heat_capacity"] = cp
    input["vessel"]["density"] = density
    input["vessel"]["orientation"] = orientation
    input["vessel"]["thickness"] = float(thk)

    input["initial"]["pressure"] = pres
    input["initial"]["temperature"] = temp
    input["initial"]["fluid"] = fluid
    input["valve"]["flow"] = mode
    input["valve"]["type"] = "psv"
    input["valve"]["diameter"] = float(orifice_diam)
    input["valve"]["discharge_coef"] = 0.975
    input["valve"]["set_pressure"] = set_pressure
    input["valve"]["back_pressure"] = back_pressure
    input["valve"]["blowdown"] = blowdown
    # input['valve']['end_pressure']=end_pressure

    input["heat_transfer"]["type"] = "s-b"
    input["heat_transfer"]["fire"] = fire_type
    input["rupture"] = {}
    input["rupture"]["fire"] = peak_fire_type
    input["rupture"]["material"] = material

    return input


if __name__ == "__main__":
    st.set_page_config(layout="wide")
    input = read_input()
    hdown = HydDown(input)
    with st.spinner("Calculating, please wait...."):
        hdown.run(disable_pbar=True)
    st.title("HydDown rigorous gas vessel fire PSV discharge calculation")
    st.subheader(r"https://github.com/andr1976/HydDown")
    my_expander = st.expander("Description")
    my_expander.write(
        "Real gas vessel fire PSV relief with fire modelled via Stefan-Boltzmann equation. Orifice size (Cd = 0.975) is specified."
        "Rupture evaluation is performed based on the selected peak fire type and vessel material using von Mises stress criteria."
    )
    my_expander.write(
        "For more information about the calculations and validation of the code please refer to the [manual](https://github.com/andr1976/HydDown/raw/main/docs/MANUAL.pdf)"
    )

    df = hdown.get_dataframe()
    file_name = st.text_input("Filename for saving data:", "saved_data")

    st.markdown(get_table_download_link(df, file_name), unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        hdown.plot(verbose=False)
        st.pyplot(plt.gcf())
        plt.tight_layout()

    with col2:
        hdown.analyze_rupture()
        st.pyplot(plt.gcf())
        plt.tight_layout()
