#!/usr/bin/env python3
"""
Run all validation cases with visual plots for inspection.

This script runs all validation YAML files and generates:
- Detailed report of results
- Time-dependent plots for each successful case
- Summary statistics
"""

import os
import sys
import time
import yaml
from pathlib import Path
import traceback
import matplotlib.pyplot as plt

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from hyddown import HydDown


def run_validation_case_with_plot(yaml_file, output_dir):
    """Run a single validation case and generate plots."""

    result = {
        "file": os.path.basename(yaml_file),
        "status": "unknown",
        "error": None,
        "runtime": 0,
        "final_pressure": None,
        "final_temperature": None,
        "final_mass": None,
        "initial_pressure": None,
        "initial_temperature": None,
        "time_steps": 0,
        "fluid": None,
        "calc_type": None,
        "valve_type": None,
        "plot_file": None,
    }

    try:
        start_time = time.time()

        # Load input file
        with open(yaml_file, 'r') as f:
            input_data = yaml.load(f, Loader=yaml.FullLoader)

        # Extract metadata
        result["fluid"] = input_data.get("initial", {}).get("fluid", "Unknown")
        result["calc_type"] = input_data.get("calculation", {}).get("type", "Unknown")
        result["valve_type"] = input_data.get("valve", {}).get("type", "Unknown")
        result["initial_pressure"] = input_data.get("initial", {}).get("pressure", None)
        result["initial_temperature"] = input_data.get("initial", {}).get("temperature", None)

        # Run simulation
        hdown = HydDown(input_data)
        hdown.run()

        # Extract results
        result["runtime"] = time.time() - start_time
        result["time_steps"] = len(hdown.time_array)
        result["final_pressure"] = hdown.P[-1] if len(hdown.P) > 0 else None
        result["final_temperature"] = hdown.T_fluid[-1] if len(hdown.T_fluid) > 0 else None
        result["final_mass"] = hdown.mass_fluid[-1] if len(hdown.mass_fluid) > 0 else None
        result["status"] = "SUCCESS"

        # Generate plots
        plot_filename = output_dir / f"{Path(yaml_file).stem}_plot.png"

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f"Validation Case: {result['file']}", fontsize=14, fontweight='bold')

        # Plot 1: Pressure vs Time
        ax1 = axes[0, 0]
        ax1.plot(hdown.time_array, hdown.P / 1e5, 'b-', linewidth=2)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Pressure (bar)')
        ax1.set_title('Pressure vs Time')
        ax1.grid(True, alpha=0.3)

        # Plot 2: Temperature vs Time
        ax2 = axes[0, 1]
        ax2.plot(hdown.time_array, hdown.T_fluid - 273.15, 'r-', linewidth=2, label='Fluid')
        if hasattr(hdown, 'T_wall') and len(hdown.T_wall) > 0:
            ax2.plot(hdown.time_array, hdown.T_wall - 273.15, 'g--', linewidth=2, label='Wall')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel(r'Temperature ($^\circ$C)')
        ax2.set_title('Temperature vs Time')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        # Plot 3: Mass vs Time
        ax3 = axes[1, 0]
        ax3.plot(hdown.time_array, hdown.mass_fluid, 'g-', linewidth=2)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Mass (kg)')
        ax3.set_title('Mass vs Time')
        ax3.grid(True, alpha=0.3)

        # Plot 4: Mass Flow Rate vs Time
        ax4 = axes[1, 1]
        ax4.plot(hdown.time_array, hdown.mass_rate, 'm-', linewidth=2)
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Mass Flow Rate (kg/s)')
        ax4.set_title('Mass Flow Rate vs Time')
        ax4.grid(True, alpha=0.3)

        # Add metadata text
        metadata_text = (
            f"Fluid: {result['fluid']}\n"
            f"Type: {result['calc_type']}\n"
            f"Valve: {result['valve_type']}\n"
            f"Runtime: {result['runtime']:.3f}s\n"
            f"Steps: {result['time_steps']}"
        )
        fig.text(0.02, 0.02, metadata_text, fontsize=9, verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
        plt.close()

        result["plot_file"] = str(plot_filename)

    except Exception as e:
        result["status"] = "FAILED"
        result["error"] = str(e)
        result["traceback"] = traceback.format_exc()

    return result


def generate_summary(results):
    """Generate summary report."""

    total = len(results)
    success = sum(1 for r in results if r["status"] == "SUCCESS")
    failed = sum(1 for r in results if r["status"] == "FAILED")

    print()
    print("=" * 100)
    print("VALIDATION SUMMARY")
    print("=" * 100)
    print(f"Total Cases:     {total}")
    print(f"Successful:      {success} ({100*success/total:.1f}%)")
    print(f"Failed:          {failed} ({100*failed/total:.1f}%)")
    print()

    if failed > 0:
        print("FAILED CASES:")
        print("-" * 100)
        for r in results:
            if r["status"] == "FAILED":
                print(f"  ✗ {r['file']}: {r['error']}")
        print()

    # Fluid type statistics
    successful_results = [r for r in results if r["status"] == "SUCCESS"]
    if successful_results:
        print("SUCCESS BY FLUID TYPE:")
        print("-" * 100)
        fluids = set(r['fluid'] for r in successful_results)
        for fluid in sorted(fluids):
            count = sum(1 for r in successful_results if r['fluid'] == fluid)
            print(f"  ✓ {fluid}: {count} cases")
        print()

    # Runtime statistics
    if successful_results:
        runtimes = [r['runtime'] for r in successful_results]
        print("PERFORMANCE:")
        print("-" * 100)
        print(f"  Average Runtime: {sum(runtimes)/len(runtimes):.3f}s")
        print(f"  Total Runtime:   {sum(runtimes):.3f}s")
        print()


def main():
    """Main execution."""

    print("=" * 100)
    print("HydDown Validation Suite with Visual Inspection")
    print("=" * 100)
    print()

    # Find validation files
    script_dir = Path(__file__).parent
    validation_dir = script_dir / "validation"
    output_dir = script_dir / "validation_plots"

    # Create output directory
    output_dir.mkdir(exist_ok=True)

    if not validation_dir.exists():
        print(f"Error: Validation directory not found: {validation_dir}")
        return

    yaml_files = sorted(validation_dir.glob("*.yml"))

    print(f"Found {len(yaml_files)} validation cases")
    print(f"Plots will be saved to: {output_dir}")
    print()

    # Run all cases
    results = []

    for i, yaml_file in enumerate(yaml_files, 1):
        print(f"[{i}/{len(yaml_files)}] Running: {yaml_file.name}...", end=" ", flush=True)
        result = run_validation_case_with_plot(yaml_file, output_dir)
        results.append(result)

        if result["status"] == "SUCCESS":
            print(f"✓ SUCCESS ({result['runtime']:.2f}s) - Plot saved")
        else:
            print(f"✗ FAILED - {result['error']}")

    # Generate summary
    generate_summary(results)

    print(f"All plots saved to: {output_dir}")
    print()


if __name__ == "__main__":
    main()
