"""Compare 1D vs 0D heat transfer models for thin-wall rupture case."""
import matplotlib.pyplot as plt
import numpy as np
import yaml
from hyddown import HydDown

# Run 0D model (lumped capacitance) - thin wall
print("Running 0D model (rupture_thin.yml)...")
with open("src/hyddown/examples/rupture_thin.yml") as f:
    input_0d = yaml.safe_load(f)
hdown_0d = HydDown(input_0d)
hdown_0d.run()

# Run 1D model (detailed heat transfer) - thin wall
print("Running 1D model (rupture_1d_thin.yml)...")
with open("src/hyddown/examples/rupture_1d_thin.yml") as f:
    input_1d = yaml.safe_load(f)
hdown_1d = HydDown(input_1d)
hdown_1d.run()

# Create comparison plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Pressure comparison
ax = axes[0, 0]
ax.plot(hdown_0d.time_array, hdown_0d.P/1e6, 'b-', label='0D model', linewidth=2)
ax.plot(hdown_1d.time_array, hdown_1d.P/1e6, 'r--', label='1D model', linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Pressure (MPa)')
ax.set_title('Pressure Comparison (Thin Wall: 30mm)')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Gas temperature comparison
ax = axes[0, 1]
ax.plot(hdown_0d.time_array, hdown_0d.T_fluid-273.15, 'b-', label='0D model', linewidth=2)
ax.plot(hdown_1d.time_array, hdown_1d.T_fluid-273.15, 'r--', label='1D model', linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Gas Temperature (°C)')
ax.set_title('Gas Temperature Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Wall temperature comparison
ax = axes[1, 0]
ax.plot(hdown_0d.time_array, hdown_0d.T_vessel-273.15, 'b-', label='0D: Bulk wall', linewidth=2)
ax.plot(hdown_1d.time_array, hdown_1d.T_vessel-273.15, 'r-', label='1D: Mean wall', linewidth=2)
ax.plot(hdown_1d.time_array, hdown_1d.T_inner_wall-273.15, 'g--', label='1D: Inner wall', linewidth=1.5)
ax.plot(hdown_1d.time_array, hdown_1d.T_outer_wall-273.15, 'm--', label='1D: Outer wall', linewidth=1.5)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Wall Temperature (°C)')
ax.set_title('Wall Temperature Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 4: Temperature difference
ax = axes[1, 1]
temp_diff = hdown_1d.T_vessel - hdown_0d.T_vessel
ax.plot(hdown_1d.time_array, temp_diff, 'k-', linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Temperature Difference (K)')
ax.set_title('1D Mean Wall - 0D Bulk Wall')
ax.axhline(y=0, color='r', linestyle='--', alpha=0.5)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('comparison_thin_0d_vs_1d.png', dpi=150)
print("\nPlot saved as 'comparison_thin_0d_vs_1d.png'")

# Print statistics
print("\n" + "="*60)
print("COMPARISON STATISTICS - THIN WALL (30 mm)")
print("="*60)
print(f"\nFinal values at t={hdown_0d.time_array[-1]:.0f}s:")
print(f"  Pressure:")
print(f"    0D model: {hdown_0d.P[-1]/1e6:.3f} MPa")
print(f"    1D model: {hdown_1d.P[-1]/1e6:.3f} MPa")
print(f"    Difference: {(hdown_1d.P[-1]-hdown_0d.P[-1])/1e6:.3f} MPa")
print(f"\n  Gas Temperature:")
print(f"    0D model: {hdown_0d.T_fluid[-1]-273.15:.2f} °C")
print(f"    1D model: {hdown_1d.T_fluid[-1]-273.15:.2f} °C")
print(f"    Difference: {hdown_1d.T_fluid[-1]-hdown_0d.T_fluid[-1]:.2f} K")
print(f"\n  Wall Temperature:")
print(f"    0D model (bulk):   {hdown_0d.T_vessel[-1]-273.15:.2f} °C")
print(f"    1D model (mean):   {hdown_1d.T_vessel[-1]-273.15:.2f} °C")
print(f"    1D model (inner):  {hdown_1d.T_inner_wall[-1]-273.15:.2f} °C")
print(f"    1D model (outer):  {hdown_1d.T_outer_wall[-1]-273.15:.2f} °C")
print(f"    Temperature gradient: {hdown_1d.T_outer_wall[-1]-hdown_1d.T_inner_wall[-1]:.2f} K")
print(f"    Difference (mean): {hdown_1d.T_vessel[-1]-hdown_0d.T_vessel[-1]:.2f} K")

print(f"\n  Thermal resistance effect:")
print(f"    Thick wall (136mm) gradient: 151.46 K")
print(f"    Thin wall (30mm) gradient:   {hdown_1d.T_outer_wall[-1]-hdown_1d.T_inner_wall[-1]:.2f} K")
print(f"    Reduction factor: {151.46/(hdown_1d.T_outer_wall[-1]-hdown_1d.T_inner_wall[-1]):.2f}x")

plt.show()
