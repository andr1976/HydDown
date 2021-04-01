from CoolProp.CoolProp import PropsSI
import math
import numpy as np

def gas_release_rate(P1,P2,T,rho,MW,k,CD,area):
    p_limit = (P1 * (2 / (1 + k)) ** (k / (k - 1)))
    if P2 < p_limit:
        p_used = p_limit
    else:
        p_used=P2
    if P1>P2:
        retval = CD * area * math.sqrt((2 * k / (k - 1)) * P1 * rho * (p_used / P1) **  (2 / k) * (1 - (p_used / P1) ** ((k - 1) / k)))
    else: 
        retval = 0 
    return retval # kg/s

def Gr(L,Tfluid,Tvessel,P,species):
    T=(Tfluid+Tvessel)/2
    beta=PropsSI('ISOBARIC_EXPANSION_COEFFICIENT','T',T,'P',P,species)
    nu=PropsSI('V','T',T,'P',P,species)/PropsSI('D','T',T,'P',P,species)
    Gr = 9.81 * beta * (Tvessel-Tfluid)*L**3/nu**2
    return Gr

def Pr(T,P,species):
    Pr = PropsSI('C','T',T,'P',P,species)*PropsSI('V','T',T,'P',P,species)/PropsSI('L','T',T,'P',P,species)
    return Pr 

def Nu(Ra):
    return 0.104*Ra**0.352

def h_inner(L,Tfluid,Tvessel,P,species):
    NPr=Pr(Tfluid,P,species)
    print("Pr:", NPr)
    NGr=Gr(L,Tfluid,Tvessel,P,species)
    print("Gr:", NGr)
    NRa=NPr*NGr
    print("Ra:", NRa)
    NNu=Nu(NRa)
    print("Nu:", NNu)
    return NNu*PropsSI('L','T',Tfluid,'P',P,species)/L


# Intial parameters and setup
length=.212#10 #internal
diameter=0.075#3 #internal
thickness=0.0030 # m

p0=1.0e7 #Pa
T0=273+36 #K
tstep=0.05 # sec
D_orifice=0.0008 #m
CD=0.65
p_back=1e5 # Pa
time_tot = 50 #s
species='HEOS::H2'
method="energybalance"
eta=1 

heat_method="specified_h"
Tamb=298.
Ufix=20.
ho=5.
hi=30.
Qfix = 0.0 #1000e3#1000. #"W"
vessel_cp=500 # J/kg K
vessel_density=7800 # kg/m3

vol=diameter**2/4*3.1415*length #m3
vol_tot=(diameter+2*thickness)**2/4*3.1415*(length+2*thickness) #m3
vol_solid=vol_tot-vol
surf_area_outer=(diameter+2*thickness)**2/4*3.1415*2+(diameter+2*thickness)*3.1415*(length+2*thickness)
surf_area_inner=(diameter)**2/4*3.1415*2+(diameter)*3.1415*length

# data storage
data_len = int(time_tot / tstep)
rho = np.zeros(data_len)
T_fluid = np.zeros(data_len)
T_vessel = np.zeros(data_len)
Q_outer = np.zeros(data_len)
Q_inner = np.zeros(data_len)
h_inside = np.zeros(data_len)
T_vent = np.zeros(data_len)
H_mass = np.zeros(data_len)
S_mass = np.zeros(data_len)
U_mass = np.zeros(data_len)
U_tot = np.zeros(data_len)
U_iter=np.zeros(data_len)
m_iter=np.zeros(data_len)
n_iter=np.zeros(data_len)
Qo=np.zeros(data_len)
Qi=np.zeros(data_len)
P = np.zeros(data_len)
mass_vessel = np.zeros(data_len)
mass_rate = np.zeros(data_len)
time_array = np.zeros(data_len)

rho0 = PropsSI('D','T',T0,'P',p0,species)
m0 = rho0*vol


# Inititialise 
rho[0] = rho0
T_fluid[0] = T0
T_vessel[0] = T0
H_mass[0] = PropsSI('H','T',T0,'P',p0,species)
S_mass[0] = PropsSI('S','T',T0,'P',p0,species)
U_mass[0] = PropsSI('U','T',T0,'P',p0,species)
U_tot[0] = PropsSI('U','T',T0,'P',p0,species)*m0
P[0] = p0
mass_vessel[0] = m0
cpcv=PropsSI('CP0MOLAR','T',T0,'P',p0,species)/PropsSI('CVMOLAR','T',T0,'P',p0,species)

mass_rate[0] = gas_release_rate(p0,p_back,T0,rho0,PropsSI('M',species),cpcv,CD,D_orifice**2/4*3.1415)
time_array[0] = 0

for i in range(1,len(time_array)):
    time_array[i]=time_array[i-1]+tstep
    mass_vessel[i]=mass_vessel[i-1]-mass_rate[i-1]*tstep
    rho[i]=mass_vessel[i]/vol

    if method == "isenthalpic":
        T_fluid[i]=PropsSI('T','D',rho[i],'H',H_mass[i-1],species)
        P[i]=PropsSI('P','D',rho[i],'H',H_mass[i-1],species)
    elif method=="isentropic":
        T_fluid[i]=PropsSI('T','D',rho[i],'S',S_mass[i-1],species)
        P[i]=PropsSI('P','D',rho[i],'S',S_mass[i-1],species)
    elif method=="isothermal":
        T_fluid[i]=T0
        P[i]=PropsSI('P','D',rho[i],'T',T0,species)
    elif method=="constantU":
        T_fluid[i]=PropsSI('T','D',rho[i],'U',U_mass[i-1],species)
        P[i]=PropsSI('P','D',rho[i],'U',U_mass[i-1],species)
    elif method=="energybalance":
        P1 = PropsSI('P','D',rho[i],'T',T_fluid[i-1],species)
        #T1 = PropsSI('T','P',P1,'H',H_mass[i-1]+Q_inner[i-1]*tstep/mass_vessel[i],species)
        T1 = PropsSI('T','P',P1,'H',H_mass[i-1],species)
        NMOL=mass_vessel[i]/PropsSI('M',species) #vol*PropsSI('D','T',T_fluid[i-1],'P',P[i-1],species)/PropsSI('M',species)
        #Q=Uheat*surf_area*(298-T_fluid[i-1])
        hi=h_inner(length,T_fluid[i-1],T_vessel[i-1],P[i-1],species)
        h_inside[i]=hi
        if heat_method=="specified_h" or heat_method=="detailed":
            Q_inner[i]=surf_area_inner*hi*(T_vessel[i-1]-T_fluid[i-1])
            Q_outer[i]=surf_area_outer*ho*(Tamb-T_vessel[i-1])
            T_vessel[i]=T_vessel[i-1]+(Q_outer[i]-Q_inner[i])*tstep/(vessel_cp*vessel_density*vol_solid)
        elif heat_method=="specified_U":
            Q_inner[i]=surf_area_outer*Ufix*(Tamb-T_fluid[i-1])
            T_vessel[i]=T_vessel[0] 
        elif heat_method=="specified_Q":
            Q_inner[i]=Qfix
            T_vessel[i]=T_vessel[0]
        else:
            Q_inner[i]=0.0
            T_vessel[i]=T_vessel[0]
        print("i: ",i," Time: ",time_array[i]," Qinner: ",Q_inner[i]," Qouter: ", Q_outer[i], " h_inner: ",h_inner(length,T_fluid[i-1],T_vessel[i-1],P[i-1],species))
        U_start=NMOL*PropsSI('HMOLAR','P',P[i-1],'T',T_fluid[i-1],species)-eta*P[i-1]*vol+Q_inner[i]*tstep
        #U_start=NMOL*PropsSI('UMOLAR','P',P[i-1],'T',T_fluid[i-1],species)
        
        U=0
        nn=0
        rho1=0
        itermax=1000
        m=0
        n=0
        relax=0.1

        while abs(rho[i]-rho1)>0.01 and m<itermax:
            m=m+1
            rho1=PropsSI('D','T',T1,'P',P1,species)
            #nn=vol*PropsSI('D','T',T1,'P',P1,species)/PropsSI('M',species)
            dd=rho[i]-rho1#NMOL-nn
            P1 = P1 + dd*1e4
            if m==itermax:
                raise Exception("Iter max exceeded for rho/P")
            while abs(U_start-U)/U_start>0.00001 and n<itermax:
                n=n+1
                U=NMOL*PropsSI('HMOLAR','P',P1,'T',T1,species)-eta*P1*vol#Q_inner[i]*tstep
                d=U_start - U 
                T1 = T1 + 0.1* d / U_start* T1
                if n==itermax:
                    raise Exception("Iter max exceeded for U/T")
            
                #print(n,d,T1)
        m_iter[i]=dd
        n_iter[i]=d
        U_iter[i]=U/NMOL*mass_vessel[i]
        P[i]=P1
        T_fluid[i]=T1
        
    else:
        raise NameError("Unknown calculation method: "+method)


    H_mass[i]=PropsSI('H','T',T_fluid[i],'P',P[i],species)
    S_mass[i]=PropsSI('S','T',T_fluid[i],'P',P[i],species)
    U_mass[i]=(mass_vessel[i]*PropsSI('H','P',P[i],'T',T_fluid[i],species)-P[i]*vol)/mass_vessel[i]#PropsSI('U','T',T_fluid[i],'P',P[i],species)#-(P[i-1]-P[i])*vol/mass_vessel[i]
    cpcv=PropsSI('CP0MOLAR','T',T_fluid[i],'P',P[i],species)/PropsSI('CVMOLAR','T',T_fluid[i],'P',P[i],species)
    mass_rate[i] = gas_release_rate(P[i],p_back,T_fluid[i],rho[i],PropsSI('M',species),cpcv,CD,D_orifice**2/4*3.1415)


import pylab as plt 

plt.figure()
plt.subplot(221)
plt.plot(time_array/60, T_fluid-273.15,label="Fluid")
plt.plot(time_array/60, T_vessel-273.15,label="Vessel")
plt.legend(loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Temperature ($^\circ$C)')

plt.subplot(222)
plt.plot(time_array/60,P/1e5)
plt.xlabel('Time (minutes)')
plt.ylabel('Pressure (bar)')

plt.subplot(223)
plt.plot(time_array/60,H_mass,label='H (J/kg)')
plt.plot(time_array/60,U_mass,label='U (J/kg)')
plt.plot(time_array/60, S_mass*100,label='S*100 (J/kg K)')
plt.plot(time_array/60,U_iter)
plt.legend(loc='best')
plt.xlabel('Time (minutes)')
plt.ylabel('Enthalpy/Internal Energy/Entropy/')

plt.subplot(224)
plt.plot(time_array/60,mass_rate,label='m_dot')
plt.xlabel('Time (minutes)')
plt.ylabel('Vent rate (kg/s)')
plt.show()



