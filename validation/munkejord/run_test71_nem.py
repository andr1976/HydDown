"""Run the Munkejord Test 71 NEM case through HydDown (full blowdown) and save results."""
import os, yaml, numpy as np
from hyddown.hdclass import HydDown
HERE = os.path.dirname(os.path.abspath(__file__))
inp = yaml.safe_load(open(os.path.join(HERE, "test71_nem.yml")))
inp["calculation"]["time_step"] = 0.02
inp["calculation"]["end_time"] = 350.0
hd = HydDown(inp)
hd.run(disable_pbar=True)

t = hd.time_array
np.savez(os.path.join(HERE, "test71_nem_out.npz"),
         t=t, P=hd.P, m=hd.mass_fluid, mg=hd.m_gas, ml=hd.m_liquid, ms=hd.m_solid,
         Tg=hd.T_gas, Tl=hd.T_liquid, mdot=hd.mass_rate,
         Tvw=hd.T_vessel_wetted, Tv=hd.T_vessel)
print(f"backend {type(hd.release_model).__name__}", flush=True)
print(f"m0 = {hd.mass_fluid[0]:.1f} kg", flush=True)
print(f"retained dry ice = {hd.m_solid[-1]:.1f} kg   (measured 8.4 kg / 17%)", flush=True)
print(f"end P = {hd.P[-1]/1e5:.2f} bar   end m = {hd.mass_fluid[-1]:.1f} kg", flush=True)
print(f"min T_gas = {np.min(hd.T_gas)-273.15:.1f} C   min T_liquid = {np.min(hd.T_liquid)-273.15:.1f} C", flush=True)
print(f"max T_gas after plateau = {np.max(hd.T_gas[t>120])-273.15:.1f} C  (measured top warms to ~-25 C)", flush=True)
print("DONE", flush=True)
