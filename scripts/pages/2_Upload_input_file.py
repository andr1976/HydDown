# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import streamlit as st
import pandas as pd
from PIL import Image
import base64
import matplotlib.pyplot as plt
import yaml

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


def run_calculation(hdown):
    hdown.run()


def read_input():

    sideb = st.sidebar

    with sideb:

        with st.form(key="my_form"):
            # Upload input yaml file
            infile = st.file_uploader(
                "Upload your input YAML file", type=["yml", "yaml"]
            )
            input = yaml.load(infile, Loader=yaml.FullLoader)

            submit_button = st.form_submit_button(label="Run calculation")
        hdown = HydDown(input)

        if submit_button:
            run_calculation(hdown)
            # submit_button = st.form_submit_button(label="Run calculation")

            st.success("Calculation completed!")
            st.markdown("### Download results:")
            df = pd.DataFrame(hdown.get_dataframe())
            st.markdown(
                get_table_download_link(df, "hyddown_results"), unsafe_allow_html=True
            )
    if submit_button:
        hdown.plot(verbose=False)
        st.pyplot(plt.gcf())

        plt.tight_layout()
        if "thermal_conductivity" in input["vessel"].keys():
            hdown.plot_tprofile(verbose=False)
            st.pyplot(plt.gcf())


read_input()
