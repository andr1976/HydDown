from CoolProp.CoolProp import PropsSI
import time
import numpy as np

import CoolProp.CoolProp as CP

HEOS=CP.AbstractState("PR","CH4&Ethane&N2&CO2&Propane&Butane")
x=np.asarray([0.9012,0.0635,0.0078,0.0234,0.0035,0.0005])
x=x/sum(x)
x=[9.01290129e-01,6.35063506e-02,7.80078008e-03,2.34023402e-02,3.50035004e-03,5.00050005e-04]

#"CH4[9.01290129e-01]&Ethane[6.35063506e-02]&N2[7.80078008e-03]&CO2[2.34023402e-02]&Propane[3.50035004e-03]&Butane[5.00050005e-04]"

HEOS.set_mole_fractions(x)
print(x)

def U(species):
    return PropsSI("U","P",150e5,"T",298,species)


species = "HEOS::CH4"

t1=time.time()

for i in range(100):
    u = U(species) 

t2 = time.time()

print("Pure time: ", t2-t1)

species = "PR::CH4[0.9]&Ethane[0.1]"

t1=time.time()

for i in range(100):
    u = U(species) 

t2 = time.time()

print("Mixture time: ", t2-t1)


t1=time.time()

HEOS.update(CP.PT_INPUTS,150e5,298)
for i in range(100):
    HEOS.update(CP.PT_INPUTS,150e5,298)
    u = HEOS.umass()

t2 = time.time()

print("Low level HEOS mixture time: ", t2-t1)

HEOS.build_phase_envelope("None")
PE=HEOS.get_phase_envelope_data()

import pylab as plt

plt.plot(PE.T, PE.p, '-', label = 'SRK with transformations in multi-fluid', color = 'g')

plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.legend(loc='best')
plt.tight_layout()
plt.show()