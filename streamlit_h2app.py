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
    length = sideb.text_input('Vessel length (m):',0.463)
    diam = sideb.text_input('Vessel diam (m):',0.254) 
    thk = sideb.text_input('Vessel thichness (m):',0.016)

    orifice_diam = sideb.text_input('Orifice diam (mm):',0.40) 
    orifice_diam = float(orifice_diam)/1000
    pres = sideb.text_input('Initial pressure (bar):', 50.)
    pres = float(pres)*1e5

    end_pressure = sideb.text_input('Target pressure (bar):',240) 
    end_pressure= float(end_pressure)*1e5

    temp = sideb.text_input('Initial temperature (C):',25)
    temp = float(temp)+273.15
    fluid = 'H2'

    tstep = sideb.text_input('Calculation time step (s):',1.0) 
    end_time = sideb.text_input('Calculation end time (s):',240) 
   


    input={}
    input['calculation'] = {}
    input['vessel'] = {}
    input['initial'] = {}
    input['valve'] = {}
    input['heat_transfer'] = {}

    input['calculation']['type'] = 'energybalance'
    input['calculation']['time_step'] = float(tstep)
    input['calculation']['end_time'] = float(end_time)
    
    input['vessel']['length'] = float(length)
    input['vessel']['diameter'] = float(diam)
    input['vessel']['heat_capacity']=470
    input['vessel']['density']=7740. 
    input['vessel']['orientation']="horizontal"
    input['vessel']['thickness']=float(thk)

    
    input['initial']['pressure'] = pres
    input['initial']['temperature'] = temp
    input['initial']['fluid'] = fluid
    input['valve']['flow'] = 'filling'
    input['valve']['type'] = 'orifice'
    input['valve']['diameter'] = float(orifice_diam)
    input['valve']['discharge_coef'] = 0.84
    input['valve']['back_pressure'] = 2*end_pressure
    input['valve']['end_pressure']=end_pressure

    input['heat_transfer']['type']='specified_h'
    input['heat_transfer']['temp_ambient']=298
    input['heat_transfer']['h_outer']=5
    input['heat_transfer']['h_inner']='calc'
    input['heat_transfer']['D_throat']=float(diam)

    col = st.beta_columns(1)
    st.title('HydDown Type I H2 cylinder filling')
    st.subheader(r'https://github.com/andr1976/HydDown')
    my_expander = st.beta_expander("Description")
    my_expander.write('Real gas vessel pressurisation at constant mass flow with heat transfer from gas to vessel and ambient. Orifice size is arbitrary and used for pressurisation rate.')

    col1, col2= st.beta_columns(2)
    hdown=HydDown(input)
    hdown.run()
    
    temp_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Fluid temperature (C)': hdown.T_fluid-273.15, 'Wall temperature (C)': hdown.T_vessel-273.15})
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
    