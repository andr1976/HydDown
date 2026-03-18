"""Test that wetted wall equals unwetted wall when no liquid present."""
import yaml
from hyddown import HydDown

# Run 1D model
print("Running 1D model (rupture_1d.yml)...")
with open("src/hyddown/examples/rupture_1d.yml") as f:
    input_1d = yaml.safe_load(f)
hdown = HydDown(input_1d)
hdown.run()

# Check final time step
print(f"\nFinal time step (t={hdown.time_array[-1]:.0f}s):")
print(f"  Liquid level: {hdown.liquid_level[-1]:.6f} m")
print(f"  Wetted area: {hdown.inner_vol.SA_from_h(hdown.liquid_level[-1]):.6f} m²")
print(f"\n  Unwetted wall temperatures:")
print(f"    Inner:  {hdown.T_inner_wall[-1]-273.15:.2f} °C")
print(f"    Outer:  {hdown.T_outer_wall[-1]-273.15:.2f} °C")
print(f"    Mean:   {hdown.T_vessel[-1]-273.15:.2f} °C")
print(f"\n  Wetted wall temperatures (should equal unwetted when no liquid):")
print(f"    Inner:  {hdown.T_inner_wall_wetted[-1]-273.15:.2f} °C")
print(f"    Outer:  {hdown.T_outer_wall_wetted[-1]-273.15:.2f} °C")
print(f"    Mean:   {hdown.T_vessel_wetted[-1]-273.15:.2f} °C")

# Check if they're equal
if hdown.liquid_level[-1] == 0:
    inner_match = abs(hdown.T_inner_wall[-1] - hdown.T_inner_wall_wetted[-1]) < 0.01
    outer_match = abs(hdown.T_outer_wall[-1] - hdown.T_outer_wall_wetted[-1]) < 0.01
    mean_match = abs(hdown.T_vessel[-1] - hdown.T_vessel_wetted[-1]) < 0.01

    print(f"\n  Temperature matching check:")
    print(f"    Inner walls match: {inner_match}")
    print(f"    Outer walls match: {outer_match}")
    print(f"    Mean temps match:  {mean_match}")

    if inner_match and outer_match and mean_match:
        print("\n✓ PASS: Wetted walls correctly equal unwetted walls when no liquid")
    else:
        print("\n✗ FAIL: Wetted walls differ from unwetted walls despite no liquid!")
        print(f"    Difference (inner): {hdown.T_inner_wall[-1] - hdown.T_inner_wall_wetted[-1]:.2f} K")
        print(f"    Difference (outer): {hdown.T_outer_wall[-1] - hdown.T_outer_wall_wetted[-1]:.2f} K")
        print(f"    Difference (mean):  {hdown.T_vessel[-1] - hdown.T_vessel_wetted[-1]:.2f} K")
