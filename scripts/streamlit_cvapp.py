# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import streamlit as st
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as  pd
from PIL import Image
import base64

try:
    from hyddown import HydDown
except:
    import sys
    import os
    hyddown_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),"..","src")
    sys.path.append(os.path.abspath(hyddown_path))
    from hyddown import HydDown

def get_table_download_link(df,filename):
    """
    Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    filename=filename+'.csv'
    return f'<a href="data:application/octet-stream;base64,{b64}" download={filename}>Download csv file</a>'

def read_input():
    sideb = st.sidebar
    
    with sideb:
        try:
            image_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),"..","docs","img","Sketch.png")
            icon = Image.open(image_path)
            st.image(icon, use_column_width=True, caption="HydDown")
        except:
            pass
        
        with st.form(key='my_form'):
            submit_button = st.form_submit_button(label='Run calculation')
            heattran = st.checkbox("Include heat transfer",value=True)
            c1,c2 = st.columns(2)
            
            with c2:
                length = st.text_input('Vessel length (m):',1015.4)
            
                diam = st.text_input('Vessel diam (m):',0.45563) 
                
                orientation = st.selectbox('Vessel orientation', ('horizontal', 'vertical'))
                temp_amb = float(st.text_input('Ambient T(C):',25) ) + 273.15
                cv = float(st.text_input('Valve Cv:',110.4) )
                
                tstep = st.text_input('Time step (s):',1.0) 

            with c1:
                pres = st.text_input('Initial pressure (bar):', 77.)
                pres = float(pres)*1e5

                back_pressure = st.text_input('Back pressure (bar):',1) 
                back_pressure= float(back_pressure)*1e5

                fluid = st.selectbox('Select fluid', ('CH4', 'NG', 'He', 'N2', 'air', 'H2','O2'))
                if fluid == 'NG':
                    fluid = "Methane[0.89571]&Ethane[5.6739e-02]&Propane[2.30395e-02]&Butane[1.03E-02]&Pentane[2.67E-03]&CO2[0.84e-02]&N2[0.3080e-2]"
                mode = 'discharge' #st.selectbox('Select mode', ('filling', 'discharge'))
    
                temp = st.text_input('Initial temp. (C):',1)
                temp = float(temp)+273.15
                end_time = st.text_input('End time (s):',900) 

                valve_tconst = float(st.text_input('Valve time const (s):',300))
               
            density = st.text_input('Vessel material density (kg/m3):',7740) 
            density= float(density)
            thk = st.text_input('Vessel thichness (m):',0.02619)
            thk = float (thk)
            cp = st.text_input('Vessel material heat capacity (J/kg K):',470) 
            cp = float(cp)


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
    input['vessel']['heat_capacity']=cp
    input['vessel']['density']=density
    input['vessel']['orientation']=orientation
    input['vessel']['thickness']=float(thk)

    
    input['initial']['pressure'] = pres
    input['initial']['temperature'] = temp
    input['initial']['fluid'] = fluid
    input['valve']['flow'] = mode
    input['valve']['type'] = 'controlvalve'
    input['valve']['characteristic'] = 'linear'
    input['valve']['time_constant'] = valve_tconst
    input['valve']['Cv'] = cv
    input['valve']['back_pressure'] = back_pressure
    #input['valve']['end_pressure']=end_pressure


    input['heat_transfer']['type']='specified_h'
    input['heat_transfer']['temp_ambient']=temp_amb
    input['heat_transfer']['h_outer']=5
    if heattran == True:
        input['heat_transfer']['h_inner']='calc'
    else:
        input['heat_transfer']['h_inner']=0.0
    input['heat_transfer']['D_throat']=float(diam)
    return input

    

if __name__ == "__main__":
    #matplotlib.use('TkAgg')
    st.set_page_config(layout='wide')

    input = read_input()
    hdown=HydDown(input)

    with st.spinner('Calculating, please wait....'):
        hdown.run(disable_pbar=True) 
    
    st.title('HydDown rigorous demo')
    st.subheader(r'https://github.com/andr1976/HydDown')
    my_expander = st.expander("Description")

    my_expander.write('Real gas vessel depressurisation with heat transfer from gas to vessel and ambient and vice versa. Flow is controlled via an ANSI/ISA control valve.')
    my_expander.write('For more information about the calculations and validation of the code please refer to the [manual](https://github.com/andr1976/HydDown/raw/main/docs/MANUAL.pdf)')

    df=hdown.get_dataframe()
    file_name=st.text_input('Filename for saving data:','saved_data') 
    
    st.markdown(get_table_download_link(df,file_name), unsafe_allow_html=True)

    col1, col2= st.columns(2)

    if input['valve']['flow']=='discharge':
        temp_data = pd.DataFrame({'Time (s)': hdown.time_array, 'Fluid temperature (C)': hdown.T_fluid-273.15, 'Wall temperature (C)': hdown.T_vessel-273.15, 'Vent temperature (C)': hdown.T_vent-273.15})
    else:
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
    
