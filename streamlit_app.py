# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import streamlit as st
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as  pd
from hyddown import HydDown

if __name__ == "__main__":
    #matplotlib.use('TkAgg')

    sideb = st.sidebar
    length = sideb.text_input('Vessel length (m):',10)
    diam = sideb.text_input('Vessel diam (m):',3) 
    orifice_diam = sideb.text_input('Orifice diam (mm):',100) 
    orifice_diam = float(orifice_diam)/1000
    pres = sideb.text_input('Initial pressure (bar):',100)
    pres = float(pres)*1e5
    temp = sideb.text_input('Initial temperature (C):',25)
    temp = float(temp)+273.15
    fluid = sideb.selectbox(
        'Select fluid',
        ('N2', 'He', 'H2', 'air', 'CH4'))
    option = sideb.selectbox(
        'Select calculation type',
        ('isothermal', 'isenthalpic', 'isentropic'))
    
    tstep = sideb.text_input('Calculation time step (s):',0.5) 
    end_time = sideb.text_input('Calculation end time (s):',100) 
    
    input={}
    input['calculation'] = {}
    input['vessel'] = {}
    input['initial'] = {}
    input['valve'] = {}

    input['calculation']['type'] = option
    input['calculation']['time_step'] = float(tstep)
    input['calculation']['end_time'] = float(end_time)
    input['vessel']['length'] = float(length)
    input['vessel']['diameter'] = float(diam)
    input['initial']['pressure'] = pres
    input['initial']['temperature'] = temp
    input['initial']['fluid'] = fluid
    input['valve']['flow'] = 'discharge'
    input['valve']['type'] = 'orifice'
    input['valve']['diameter'] = float(orifice_diam)
    input['valve']['discharge_coef'] = 0.84
    input['valve']['back_pressure'] = 1e5

    
    col = st.beta_columns(1)
    st.title('HydDown adiabatic demo')
    st.subheader(r'https://github.com/andr1976/HydDown')
    my_expander = st.beta_expander("Description")
    my_expander.write('Real gas vessel depressurisation for pure and pseudo-pure components. No heat transfer is enabled in this demo version.')

    col1, col2= st.beta_columns(2)
    hdown=HydDown(input)
    hdown.run()
    
    temp_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Temperature (C)': hdown.T_fluid-273.15})
    pres_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Pressure (bar)': hdown.P/1e5})

    col1.line_chart(pres_data.rename(columns={'Time (s)':'index'}).set_index('index'))
    col1.text('Time (s)')
    col2.line_chart(temp_data.rename(columns={'Time (s)':'index'}).set_index('index'))
    col2.text('Time (s)')
    
    
    mdot_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Mass rate (kg/s)': hdown.mass_rate})
    mass_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Fluid inventory (kg)': hdown.mass_fluid})
    col1.line_chart(mdot_data.rename(columns={'Time (s)':'index'}).set_index('index'))
    col1.text('Time (s)')
    col2.line_chart(mass_data.rename(columns={'Time (s)':'index'}).set_index('index'))
    col2.text('Time (s)')
    
