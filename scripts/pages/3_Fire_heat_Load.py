# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import streamlit as st
import pandas as pd
from PIL import Image
import base64
import matplotlib.pyplot as plt

try:
    from hyddown import HydDown
except:
    import sys
    import os

    hyddown_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "src")
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
            image_path = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                "..",
                "docs",
                "img",
                "Sketch.png",
            )
            icon = Image.open(image_path)
            st.image(icon, use_column_width=True, caption="HydDown")
        except:
            pass

        with st.form(key="my_form"):
            submit_button = st.form_submit_button(label="Run calculation")
            heattran = True
            c1, c2 = st.columns(2)

            with c2:
                length = st.text_input("Vessel length (m):", 0.463)
                diam = st.text_input("Vessel diam (m):", 0.254)
                thk = st.text_input("Vessel thickness (m):", 0.016)
                orientation = st.selectbox(
                    "Vessel orientation", ("horizontal", "vertical")
                )
                orifice_diam = st.text_input("Orifice diam (mm):", 0.40)
                orifice_diam = float(orifice_diam) / 1000
                tstep = st.text_input("Time step (s):", 1.0)

            with c1:
                pres = st.text_input("Initial pressure (bar):", 100.0)
                pres = float(pres) * 1e5

                back_pressure = st.text_input("Back pres. (bar):", 1)
                back_pressure = float(back_pressure) * 1e5

                # fluid = st.selectbox('Select fluid', ('H2', 'He', 'N2', 'air', 'CH4', 'O2'))
                fluid = st.selectbox(
                    "Select fluid", ("H2", "NG", "NG1", "He", "N2", "air", "CH4", "O2")
                )
                if fluid == "NG":
                    fluid = "Methane[0.89571]&Ethane[5.6739e-02]&Propane[2.30395e-02]&Butane[1.03E-02]&Pentane[2.67E-03]&CO2[0.84e-02]&N2[0.3080e-2]"
                if fluid == "NG1":
                    fluid = "Methane[0.860231]&Ethane[0.078217]&Propane[0.033786]&Butane[9.210E-03]&Pentane[2.573E-03]&Hexane[3.560E-04]&CO2[1.206E-02]&N2[3.701E-03]"

                mode = "discharge"
                temp = st.text_input("Initial temp. (C):", 25)
                temp = float(temp) + 273.15
                end_time = st.text_input("End time (s):", 240)

            fire_type = st.selectbox(
                "Select fire heat load type",
                ("scandpower_pool", "scandpower_jet", "api_pool", "api_jet"),
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
    input["valve"]["type"] = "orifice"
    input["valve"]["diameter"] = float(orifice_diam)
    input["valve"]["discharge_coef"] = 0.84
    input["valve"]["back_pressure"] = back_pressure
    # input['valve']['end_pressure']=end_pressure

    input["heat_transfer"]["type"] = "s-b"
    input["heat_transfer"]["fire"] = fire_type
    return input


if __name__ == "__main__":
    st.set_page_config(layout="wide")

    input = read_input()
    hdown = HydDown(input)

    with st.spinner("Calculating, please wait...."):
        hdown.run(disable_pbar=True)

    st.title("HydDown - Fire heat load on gas filled vessel")
    st.subheader(r"https://github.com/andr1976/HydDown")
    my_expander = st.expander("Description")

    my_expander.write(
        "Gas filled vessel subject to external fire heat load during blowdown.  Orifice size (Cd = 0.84) is specified for desired depressurisation rate."
    )
    my_expander.write(
        "For more information about the calculations and validation of the code please refer to the [manual](https://github.com/andr1976/HydDown/raw/main/docs/MANUAL.pdf)"
    )

    df = hdown.get_dataframe()
    file_name = st.text_input("Filename for saving data:", "saved_data")

    st.markdown(get_table_download_link(df, file_name), unsafe_allow_html=True)
    hdown.plot(verbose=False)
    st.pyplot(plt.gcf())
    plt.tight_layout()
