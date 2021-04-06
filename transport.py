# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import math 
from CoolProp.CoolProp import PropsSI

def Gr(L,Tfluid,Tvessel,P,species):
    """
    Calculation of Grasshof number. See eq. 4.7-4 in 
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition, 
    Prentice-Hall, 1993
    """
    T=(Tfluid+Tvessel)/2
    beta=PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T,'P',P,species)
    nu=PropsSI('V','T',T,'P',P,species)/PropsSI('D','T',T,'P',P,species)
    Gr = 9.81 * beta * abs(Tvessel-Tfluid)*L**3/nu**2
    return Gr

def Pr(T,P,species):
    """
    Calculation of Prandtl number, eq. 4.5-6 in 
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition, 
    Prentice-Hall, 1993
    """
    Pr = PropsSI('C','T',T,'P',P,species)*PropsSI('V','T',T,'P',P,species)/PropsSI('L','T',T,'P',P,species)
    return Pr 

def Nu(Ra,Pr):
    """
    Calculation of Nusselt number for natural convection. See eq. 4.7-4  and Table 4.7-1 in 
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition, 
    Prentice-Hall, 1993
    """
    if Ra>=1e9:
        NNu = 0.13*Ra**0.333
    elif Ra<1e9 and Ra>1e4:
        NNu = 0.59*Ra**0.25
    else:
        NNu = 1.36*Ra**0.20
    return NNu

def h_inner(L,Tfluid,Tvessel,P,species):
    """
    Calculation of heat transfer coefficient from Nusselt number
    """
    NPr=Pr((Tfluid+Tvessel)/2,P,species)
    NGr=Gr(L,Tfluid,Tvessel,P,species)
    NRa=NPr*NGr
    NNu=Nu(NRa,NPr)
    return NNu*PropsSI('L','T',(Tfluid+Tvessel)/2,'P',P,species)/L

def gas_release_rate(P1,P2,rho,k,CD,area):
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996
    """
    if P1>P2:
        if P1 / P2 > ((k + 1) / 2)**((k) / (k - 1)):
            flow_coef = 1
        else:
            flow_coef = 2 / (k - 1) * (((k + 1) / 2)**((k + 1) / (k - 1))) * ((P2 / P1)**(2 / k)) * (1 - (P2 / P1)**((k - 1) / k))

        return math.sqrt(flow_coef) * CD * area * math.sqrt(rho * P1 * k * (2 / (k + 1))**((k + 1) / (k - 1)))
    else:
        return 0 

def relief_valve(P1,Pback,Pset,blowdown,rho,k,CD,area):
    """
    Pop action relief valve model including hysteresis. 
    The pressure shall rise above P_set to open and 
    decrease below P_reseat (P_set*(1-blowdown)) to close
    """
    global psv_state
    if P1>Pset:
        eff_area=area
        psv_state="open"
    elif P1<Pset*(1-blowdown):
        eff_area=0
        psv_state="closed"
    else:
        if psv_state=="open":
            eff_area=area
        elif psv_state=="closed":
            eff_area=0
        else:
            raise ValueError("Unknown PSV open/close state.")

    if eff_area > 0:
        return gas_release_rate(P1,Pback,rho,k,CD,area)
    else:
        return 0.0

def control_valve(P1,P2,T,Z,MW,gamma,Cv,xT=0.75,FP=1):
    """
    Flow calculated from ANSI/ISA control valve equations for single phase gas flow.
    Equation 19 pp. 132 in
    Control Valves / Guy Borden, editor; Paul Friedmann, style editor
    """
    P1=P1/1e5
    P2=P2/1e5
    MW=MW*1000
    N8=94.8
    Fk=gamma/1.4
    x = (P1-P2)/P1
    Y= 1. - min(x,Fk*xT) / (3. * Fk * xT)
    mass_flow = N8 * FP * Cv * P1 * Y * (MW  * min(x,xT*Fk) / T / Z)**0.5
    return mass_flow/3600 # kg/s

    
if __name__=='__main__':
    gamma=1.2
    P1=10.e5          # bar
    P2=5.5e5           # bar
    xT=0.75         # default
    MW=20.    
    T1=20.+273.15   # Kelvin
    Z1=0.9          
    print(Gr(0.305,311,505.4,1e5,'HEOS::air'))
    print(PropsSI('D','T',(311+504)/2,'P',1e5,'HEOS::air'))
    print(PropsSI('V','T',(311+504)/2,'P',1e5,'HEOS::air')/PropsSI('D','T',(311+504)/2,'P',1e5,'HEOS::air'))
    print(control_valve(P1,P2,T1,Z1,MW,gamma,500))
    print(gas_release_rate(P1,P2,9.12,1.1,0.85,0.01))
    print(Pr((311+504)/2,1e5,'HEOS::air'))