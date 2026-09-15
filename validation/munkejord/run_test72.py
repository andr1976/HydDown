import os, yaml, numpy as np
from hyddown.hdclass import HydDown
H=os.path.dirname(os.path.abspath(__file__))
hd=HydDown(yaml.safe_load(open(os.path.join(H,"test72_dense.yml")))); hd.run(disable_pbar=True)
t=hd.time_array
np.savez(os.path.join(H,"test72_out.npz"),t=t,P=hd.P,m=hd.mass_fluid,ms=hd.m_solid,Tg=hd.T_gas,Tl=hd.T_liquid,mdot=hd.mass_rate)
print("m0=%.1f kg (meas W0 51.4)"%hd.mass_fluid[0],flush=True)
print("retained dry ice=%.1f kg (measured ~11.7)"%hd.m_solid[-1],flush=True)
print("end P=%.2f bar end m=%.1f kg"%(hd.P[-1]/1e5,hd.mass_fluid[-1]),flush=True)
print("min Tl=%.1f min Tg=%.1f, max Tg>150s=%.1f C"%(np.min(hd.T_liquid)-273.15,np.min(hd.T_gas)-273.15,np.max(hd.T_gas[t>150])-273.15),flush=True)
print("DONE",flush=True)
