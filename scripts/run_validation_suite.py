#!/usr/bin/env python3
"""
Run all validation cases and generate comprehensive report.

This script runs all validation YAML files in the validation folder
and generates a detailed report of results, including:
- Success/failure status
- Final conditions (pressure, temperature, mass)
- Runtime performance
- Validation data comparisons (if available)
"""

import os
import sys
import time
import yaml
from pathlib import Path
import traceback

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from hyddown import HydDown


def run_validation_case(yaml_file):
    """Run a single validation case and return results."""

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

    except Exception as e:
        result["status"] = "FAILED"
        result["error"] = str(e)
        result["traceback"] = traceback.format_exc()

    return result


def generate_report(results):
    """Generate comprehensive validation report."""

    report = []
    report.append("=" * 100)
    report.append("HydDown Validation Suite - Comprehensive Report")
    report.append("=" * 100)
    report.append("")

    # Summary statistics
    total = len(results)
    success = sum(1 for r in results if r["status"] == "SUCCESS")
    failed = sum(1 for r in results if r["status"] == "FAILED")

    report.append(f"SUMMARY")
    report.append("-" * 100)
    report.append(f"Total Cases:     {total}")
    report.append(f"Successful:      {success} ({100*success/total:.1f}%)")
    report.append(f"Failed:          {failed} ({100*failed/total:.1f}%)")
    report.append("")

    # Successful cases
    report.append("=" * 100)
    report.append("SUCCESSFUL VALIDATION CASES")
    report.append("=" * 100)
    report.append("")

    successful_results = [r for r in results if r["status"] == "SUCCESS"]

    if successful_results:
        for i, r in enumerate(successful_results, 1):
            report.append(f"{i}. {r['file']}")
            report.append("-" * 100)
            report.append(f"   Fluid:              {r['fluid']}")
            report.append(f"   Calculation Type:   {r['calc_type']}")
            report.append(f"   Valve Type:         {r['valve_type']}")
            report.append(f"   Initial Pressure:   {r['initial_pressure']/1e5:.2f} bar" if r['initial_pressure'] else "   Initial Pressure:   N/A")
            report.append(f"   Initial Temperature: {r['initial_temperature']-273.15:.2f} °C" if r['initial_temperature'] else "   Initial Temperature: N/A")
            report.append(f"   Final Pressure:     {r['final_pressure']/1e5:.2f} bar" if r['final_pressure'] else "   Final Pressure:     N/A")
            report.append(f"   Final Temperature:  {r['final_temperature']-273.15:.2f} °C" if r['final_temperature'] else "   Final Temperature:  N/A")
            report.append(f"   Final Mass:         {r['final_mass']:.3f} kg" if r['final_mass'] else "   Final Mass:         N/A")
            report.append(f"   Time Steps:         {r['time_steps']}")
            report.append(f"   Runtime:            {r['runtime']:.3f} seconds")
            report.append("")
    else:
        report.append("   No successful cases")
        report.append("")

    # Failed cases
    if failed > 0:
        report.append("=" * 100)
        report.append("FAILED VALIDATION CASES")
        report.append("=" * 100)
        report.append("")

        failed_results = [r for r in results if r["status"] == "FAILED"]

        for i, r in enumerate(failed_results, 1):
            report.append(f"{i}. {r['file']}")
            report.append("-" * 100)
            report.append(f"   Error: {r['error']}")
            report.append("")
            report.append("   Traceback:")
            for line in r['traceback'].split('\n'):
                report.append(f"   {line}")
            report.append("")

    # Performance statistics
    if successful_results:
        report.append("=" * 100)
        report.append("PERFORMANCE STATISTICS")
        report.append("=" * 100)
        report.append("")

        runtimes = [r['runtime'] for r in successful_results]
        avg_runtime = sum(runtimes) / len(runtimes)
        max_runtime = max(runtimes)
        min_runtime = min(runtimes)
        max_case = next(r['file'] for r in successful_results if r['runtime'] == max_runtime)
        min_case = next(r['file'] for r in successful_results if r['runtime'] == min_runtime)

        report.append(f"Average Runtime:    {avg_runtime:.3f} seconds")
        report.append(f"Maximum Runtime:    {max_runtime:.3f} seconds ({max_case})")
        report.append(f"Minimum Runtime:    {min_runtime:.3f} seconds ({min_case})")
        report.append(f"Total Runtime:      {sum(runtimes):.3f} seconds")
        report.append("")

    # Fluids tested
    report.append("=" * 100)
    report.append("FLUIDS TESTED")
    report.append("=" * 100)
    report.append("")

    fluids = set(r['fluid'] for r in successful_results)
    for fluid in sorted(fluids):
        count = sum(1 for r in successful_results if r['fluid'] == fluid)
        report.append(f"   {fluid}: {count} cases")
    report.append("")

    # Calculation types
    report.append("=" * 100)
    report.append("CALCULATION TYPES")
    report.append("=" * 100)
    report.append("")

    calc_types = set(r['calc_type'] for r in successful_results)
    for calc_type in sorted(calc_types):
        count = sum(1 for r in successful_results if r['calc_type'] == calc_type)
        report.append(f"   {calc_type}: {count} cases")
    report.append("")

    report.append("=" * 100)
    report.append("END OF REPORT")
    report.append("=" * 100)

    return "\n".join(report)


def main():
    """Main execution."""

    print("=" * 100)
    print("HydDown Validation Suite Runner")
    print("=" * 100)
    print()

    # Find validation files
    script_dir = Path(__file__).parent
    validation_dir = script_dir / "validation"

    if not validation_dir.exists():
        print(f"Error: Validation directory not found: {validation_dir}")
        return

    yaml_files = sorted(validation_dir.glob("*.yml"))

    print(f"Found {len(yaml_files)} validation cases")
    print()

    # Run all cases
    results = []

    for i, yaml_file in enumerate(yaml_files, 1):
        print(f"[{i}/{len(yaml_files)}] Running: {yaml_file.name}...", end=" ", flush=True)
        result = run_validation_case(yaml_file)
        results.append(result)

        if result["status"] == "SUCCESS":
            print(f"✓ SUCCESS ({result['runtime']:.2f}s)")
        else:
            print(f"✗ FAILED - {result['error']}")

    print()

    # Generate and save report
    report = generate_report(results)

    report_file = script_dir / "VALIDATION_REPORT.txt"
    with open(report_file, 'w') as f:
        f.write(report)

    print(f"Report saved to: {report_file}")
    print()

    # Print summary
    total = len(results)
    success = sum(1 for r in results if r["status"] == "SUCCESS")
    failed = sum(1 for r in results if r["status"] == "FAILED")

    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"Total:      {total}")
    print(f"Success:    {success} ({100*success/total:.1f}%)")
    print(f"Failed:     {failed} ({100*failed/total:.1f}%)")
    print()

    # Print report to console as well
    print()
    print(report)


if __name__ == "__main__":
    main()
