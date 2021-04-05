import math 

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

def control_valve(P1,P2,T,Z,MW,gamma,Cv,xT=0.75,FP=1):
    """
    Flow calculated from ANSI/ISA control valve equations for single phase gas flow.
    Equation 19 pp. 132 in
    Control Valves / Guy Borden, editor; Paul Friedmann, style editor
    """
    P1=P1/1e5
    P2=P2/1e5
    N8=94.8
    Fk=gamma/1.4
    x = (P1-P2)/P1
    Y= 1. - min(x,Fk*xT) / (3. * Fk * xT)
    mass_flow = N8 * FP * Cv * P1 * Y * (MW * min(x,xT*Fk) / T / Z)**0.5
    return mass_flow/3600 # kg/s

    
if __name__=='__main__':
    gamma=1.2
    P1=10.e5          # bar
    P2=5.5e5           # bar
    xT=0.75         # default
    MolW=20.    
    T1=20.+273.15   # Kelvin
    Z1=0.9          

    print(control_valve(P1,P2,T1,Z1,MolW,gamma,500))
    print(gas_release_rate(P1,P2,9.12,1.1,0.85,0.01))
    