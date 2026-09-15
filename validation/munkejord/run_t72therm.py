import os, yaml, numpy as np
from hyddown.hdclass import HydDown
H=os.path.dirname(os.path.abspath(__file__))
hd=HydDown(yaml.safe_load(open(os.path.join(H,"test72_thermesh.yml")))); hd.run(disable_pbar=True)
t=hd.time_array
np.savez(os.path.join(H,"t72therm_out.npz"),t=t,P=hd.P,m=hd.mass_fluid,ms=hd.m_solid,Tg=hd.T_gas,Tl=hd.T_liquid)
print("THERMESH Test72: retained dry ice=%.1f kg (lumped 7.8, measured 11.7)"%hd.m_solid[-1],flush=True)
print("end P=%.2f bar m=%.1f kg minTl=%.1f"%(hd.P[-1]/1e5,hd.mass_fluid[-1],np.min(hd.T_liquid)-273.15),flush=True)
print("DONE",flush=True)
