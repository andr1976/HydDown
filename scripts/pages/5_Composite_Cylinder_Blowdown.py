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
            heattran = True  # st.checkbox("Include heat transfer", value=True)
            c1, c2 = st.columns(2)

            with c2:
                length = st.text_input("Vessel length (m):", 0.463)
                diam = st.text_input("Vessel diam (m):", 0.254)
                orientation = st.selectbox(
                    "Vessel orientation", ("horizontal", "vertical")
                )
                orifice_diam = st.text_input("Orifice diam (mm):", 1)
                orifice_diam = float(orifice_diam) / 1000
                tstep = st.text_input("Time step (s):", 1.0)

            with c1:
                pres = st.text_input("Initial pressure (bar):", 100.0)
                pres = float(pres) * 1e5

                back_pressure = 1.013e5

                # fluid = st.selectbox('Select fluid', ('H2', 'He', 'N2', 'air', 'CH4', 'O2'))
                fluid = st.selectbox(
                    "Select fluid", ("CH4", "NG", "NG1", "He", "N2", "air", "H2", "O2")
                )
                if fluid == "NG":
                    fluid = "Methane[0.89571]&Ethane[5.6739e-02]&Propane[2.30395e-02]&Butane[1.03E-02]&Pentane[2.67E-03]&CO2[0.84e-02]&N2[0.3080e-2]"
                if fluid == "NG1":
                    fluid = "Methane[0.860231]&Ethane[0.078217]&Propane[0.033786]&Butane[9.210E-03]&Pentane[2.573E-03]&Hexane[3.560E-04]&CO2[1.206E-02]&N2[3.701E-03]"

                mode = (
                    "discharge"  # st.selectbox("Select mode", ("filling", "discharge"))
                )
                temp = st.text_input("Initial temp. (C):", 25)
                temp = float(temp) + 273.15
                end_time = st.text_input("End time (s):", 240)

            liner_density = st.text_input("Liner/shell material density (kg/m3):", 945)
            liner_density = float(liner_density)

            liner_cp = st.text_input(
                "Liner/shell material heat capacity (J/kg K):", 1584
            )
            liner_cp = float(liner_cp)

            liner_k = float(
                float(
                    st.text_input(
                        "Liner/shell material thermal conductivity (W/m K):", 0.385
                    )
                )
            )
            liner_thk = float(st.text_input("Liner/shell thickness (mm):", 10)) / 1000

            density = st.text_input(
                "Insulation/outer shell material density (kg/m3):", 1360.0
            )
            density = float(density)

            cp = st.text_input(
                "Insulation/outer shell material heat capacity (J/kg K):", 1020
            )
            cp = float(cp)

            k = float(
                st.text_input(
                    "Insulation/outer shell material thermal conductivity (W/m K):",
                    0.5,
                )
            )
            thk = float(st.text_input("Vessel thickness (mm):", 20)) / 1000

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
    input["vessel"]["thermal_conductivity"] = k
    input["vessel"]["liner_heat_capacity"] = liner_cp
    input["vessel"]["liner_density"] = liner_density
    input["vessel"]["liner_thermal_conductivity"] = liner_k

    input["vessel"]["orientation"] = orientation
    input["vessel"]["thickness"] = float(thk)
    input["vessel"]["liner_thickness"] = float(liner_thk)

    input["initial"]["pressure"] = pres
    input["initial"]["temperature"] = temp
    input["initial"]["fluid"] = fluid
    input["valve"]["flow"] = mode
    input["valve"]["type"] = "orifice"
    input["valve"]["diameter"] = float(orifice_diam)
    input["valve"]["discharge_coef"] = 0.84
    input["valve"]["back_pressure"] = back_pressure
    # input['valve']['end_pressure']=end_pressure

    input["heat_transfer"]["type"] = "specified_h"
    input["heat_transfer"]["temp_ambient"] = 298
    input["heat_transfer"]["h_outer"] = 8
    input["heat_transfer"]["h_inner"] = "calc"
    return input


if __name__ == "__main__":
    st.set_page_config(layout="wide")

    input = read_input()
    hdown = HydDown(input)

    with st.spinner("Calculating, please wait...."):
        hdown.run(disable_pbar=True)

    st.title("HydDown - Composite cylinder blowdown")
    st.subheader(r"https://github.com/andr1976/HydDown")
    my_expander = st.expander("Description")

    my_expander.write(
        "Real gas vessel depressurisation of composite cylider e.g. Type III/IV hydrogen storage cylinder or e.g. an insulated tankwith heat transfer from gas to vessel and ambient and vice versa. Orifice size (Cd = 0.84) is specified for desired pressurisation/depressurisation rate."
    )
    my_expander.write(
        "For more information about the calculations and validation of the code please refer to the [manual](https://github.com/andr1976/HydDown/raw/main/docs/MANUAL.pdf)"
    )

    df = hdown.get_dataframe()
    file_name = st.text_input("Filename for saving data:", "saved_data")

    st.markdown(get_table_download_link(df, file_name), unsafe_allow_html=True)
    hdown.plot(verbose=False)
    st.pyplot(plt.gcf())

    hdown.plot_tprofile(verbose=False)
    st.pyplot(plt.gcf())
    plt.tight_layout()
