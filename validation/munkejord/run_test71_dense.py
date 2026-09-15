import os, yaml, numpy as np
from hyddown.hdclass import HydDown
H=os.path.dirname(os.path.abspath(__file__))
inp=yaml.safe_load(open(os.path.join(H,"test71_dense.yml")))
hd=HydDown(inp); hd.run(disable_pbar=True)
t=hd.time_array
np.savez(os.path.join(H,"test71_dense_out.npz"), t=t,P=hd.P,m=hd.mass_fluid,mg=hd.m_gas,ml=hd.m_liquid,ms=hd.m_solid,Tg=hd.T_gas,Tl=hd.T_liquid,mdot=hd.mass_rate,Tvw=hd.T_vessel_wetted,Tv=hd.T_vessel)
print("m0=%.1f kg (measured 49.4)"%hd.mass_fluid[0], flush=True)
print("retained dry ice=%.1f kg (measured 8.4)"%hd.m_solid[-1], flush=True)
print("end P=%.2f bar  end m=%.1f kg"%(hd.P[-1]/1e5, hd.mass_fluid[-1]), flush=True)
print("min Tl=%.1f  min Tg=%.1f C"%(np.min(hd.T_liquid)-273.15, np.min(hd.T_gas)-273.15), flush=True)
print("max Tg after 120s=%.1f C (measured top ~-25)"%(np.max(hd.T_gas[t>120])-273.15), flush=True)
print("DONE", flush=True)
