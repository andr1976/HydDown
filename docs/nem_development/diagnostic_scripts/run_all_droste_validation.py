"""
Complete validation script for Droste & Schoen (1988) experiments 1-3
Runs NEM simulations, performs rupture analysis, and generates validation plots
"""

import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_NEM_DIR = os.path.abspath(os.path.join(_SCRIPT_DIR, '..'))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'src'))

import yaml
from hyddown import HydDown
import matplotlib.pyplot as plt
import numpy as np

# Experimental data from Droste & Schoen (1988)
EXPERIMENTAL_DATA = {
    'exp1': {
        'name': 'Experiment 1 (Oct 1982)',
        'initial_T': 10,  # °C
        'initial_P': 5.5,  # bar
        'psv_opening_time': 340,  # s
        'rupture_time': 720,  # s
        'pressure_data': {
            'time': [0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720],
            'pres': [5.5, 7.0, 9.8, 12.5, 14.3, 15.6, 17.0, 18.5, 20.0, 21.0, 22.0, 23.0, 24.5]
        },
        'liquid_temp_data': {
            'time': [0, 120, 240, 360, 480, 600, 720],
            'temp': [10, 25, 45, 55, 62, 68, 72]
        }
    },
    'exp2': {
        'name': 'Experiment 2 (Nov 1983) - Corrected blowdown',
        'initial_T': 37,  # °C
        'initial_P': 13.5,  # bar
        'psv_opening_time': 100,  # s
        'rupture_time': 440,  # s
        'pressure_data': {
            'time': [0, 60, 120, 180, 240, 300, 360, 420, 440],
            'pres': [13.5, 15.6, 17.5, 19.5, 21.5, 23.5, 25.5, 27.0, 27.5]
        },
        'liquid_temp_data': {
            'time': [0, 120, 240, 360, 440],
            'temp': [37, 60, 72, 80, 85]
        }
    },
    'exp3': {
        'name': 'Experiment 3 (Dec 1983)',
        'initial_T': 26,  # °C
        'initial_P': 9.8,  # bar
        'psv_opening_time': 150,  # s
        'rupture_time': 540,  # s
        'pressure_data': {
            'time': [0, 60, 120, 180, 240, 300, 360, 420, 480, 540],
            'pres': [9.8, 12.5, 15.6, 17.8, 20.0, 22.0, 24.0, 25.5, 26.5, 27.5]
        },
        'liquid_temp_data': {
            'time': [0, 120, 240, 360, 480, 540],
            'temp': [26, 48, 65, 75, 82, 86]
        }
    }
}

def run_experiment(yaml_file, exp_key):
    """Run a single experiment and return results"""
    print(f"\n{'='*80}")
    print(f"Running {EXPERIMENTAL_DATA[exp_key]['name']}")
    print(f"YAML file: {yaml_file}")
    print(f"{'='*80}")

    # Load YAML file
    with open(yaml_file) as infile:
        input_dict = yaml.load(infile, Loader=yaml.FullLoader)

    # Initialize and run simulation
    hd = HydDown(input_dict)
    hd.run()

    # Perform rupture analysis
    print("\nPerforming rupture analysis...")
    hd.analyze_rupture()

    # Extract key metrics
    results = {
        'hd': hd,
        'psv_opening_time': None,
        'final_pressure': hd.P[-1] / 1e5,  # bar
        'final_temp_liquid': hd.T_liquid[-1] - 273.15,  # °C
        'final_temp_gas': hd.T_gas[-1] - 273.15,  # °C
        'rupture_time': hd.rupture_time if hasattr(hd, 'rupture_time') else None
    }

    # Find PSV opening time
    for i, mdot in enumerate(hd.mass_rate):
        if abs(mdot) > 1e-6:  # Non-zero mass flow (discharge is negative)
            results['psv_opening_time'] = hd.time_array[i]
            break

    # Print summary
    print(f"\nResults Summary:")
    print(f"  PSV opening time: {results['psv_opening_time']:.1f} s (exp: {EXPERIMENTAL_DATA[exp_key]['psv_opening_time']} s)")
    if results['rupture_time']:
        print(f"  Rupture time: {results['rupture_time']:.1f} s (exp: {EXPERIMENTAL_DATA[exp_key]['rupture_time']} s)")
    else:
        print(f"  No rupture detected (exp: {EXPERIMENTAL_DATA[exp_key]['rupture_time']} s)")
    print(f"  Final pressure: {results['final_pressure']:.1f} bar")
    print(f"  Final liquid temp: {results['final_temp_liquid']:.1f} °C")
    print(f"  Final gas temp: {results['final_temp_gas']:.1f} °C")

    return results

def create_validation_plot(exp_key, results, exp_data):
    """Create validation plot for a single experiment"""
    hd = results['hd']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    fig.suptitle(f"Droste & Schoen (1988) - {exp_data['name']}\nNEM Validation",
                 fontsize=14, fontweight='bold')

    # Pressure plot
    ax1.plot(hd.time_array, hd.P / 1e5, 'b-', linewidth=2, label='NEM Simulation')
    ax1.plot(exp_data['pressure_data']['time'], exp_data['pressure_data']['pres'],
             'ro', markersize=8, label='Experimental Data')
    ax1.axhline(y=15.6, color='g', linestyle='--', linewidth=1.5, label='PSV Set Pressure (15.6 bar)')

    # Mark PSV opening
    if results['psv_opening_time']:
        ax1.axvline(x=results['psv_opening_time'], color='orange', linestyle=':',
                   linewidth=2, label=f"PSV Opening ({results['psv_opening_time']:.0f}s)")
    ax1.axvline(x=exp_data['psv_opening_time'], color='orange', linestyle='--',
               linewidth=1.5, alpha=0.5, label=f"Exp PSV Opening ({exp_data['psv_opening_time']}s)")

    # Mark rupture
    if results['rupture_time']:
        ax1.axvline(x=results['rupture_time'], color='red', linestyle=':',
                   linewidth=2, label=f"Sim Rupture ({results['rupture_time']:.0f}s)")
    ax1.axvline(x=exp_data['rupture_time'], color='red', linestyle='--',
               linewidth=1.5, alpha=0.5, label=f"Exp Rupture ({exp_data['rupture_time']}s)")

    ax1.set_xlabel('Time [s]', fontsize=12)
    ax1.set_ylabel('Pressure [bar]', fontsize=12)
    ax1.set_title('Pressure Evolution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='best', fontsize=9)

    # Temperature plot
    ax2.plot(hd.time_array, hd.T_liquid - 273.15, 'b-', linewidth=2, label='NEM Liquid')
    ax2.plot(hd.time_array, hd.T_gas - 273.15, 'c--', linewidth=1.5, label='NEM Gas')
    ax2.plot(exp_data['liquid_temp_data']['time'], exp_data['liquid_temp_data']['temp'],
             'ro', markersize=8, label='Exp Liquid')

    # Mark PSV opening
    if results['psv_opening_time']:
        ax2.axvline(x=results['psv_opening_time'], color='orange', linestyle=':',
                   linewidth=2, label=f"PSV Opening ({results['psv_opening_time']:.0f}s)")

    # Mark rupture
    if results['rupture_time']:
        ax2.axvline(x=results['rupture_time'], color='red', linestyle=':',
                   linewidth=2, label=f"Sim Rupture ({results['rupture_time']:.0f}s)")
    ax2.axvline(x=exp_data['rupture_time'], color='red', linestyle='--',
               linewidth=1.5, alpha=0.5, label=f"Exp Rupture ({exp_data['rupture_time']}s)")

    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Temperature [°C]', fontsize=12)
    ax2.set_title('Temperature Evolution', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='best', fontsize=9)

    plt.tight_layout()

    return fig

def calculate_errors(exp_key, results, exp_data):
    """Calculate validation errors"""
    errors = {}

    # PSV opening time error
    if results['psv_opening_time']:
        errors['psv_opening'] = ((results['psv_opening_time'] - exp_data['psv_opening_time']) /
                                 exp_data['psv_opening_time'] * 100)
    else:
        errors['psv_opening'] = None

    # Rupture time error
    if results['rupture_time']:
        errors['rupture_time'] = ((results['rupture_time'] - exp_data['rupture_time']) /
                                  exp_data['rupture_time'] * 100)
    else:
        errors['rupture_time'] = None

    # Final pressure error (at experimental rupture time)
    idx_rupture = np.argmin(np.abs(results['hd'].time_array - exp_data['rupture_time']))
    sim_pressure_at_rupture = results['hd'].P[idx_rupture] / 1e5
    exp_pressure_at_rupture = exp_data['pressure_data']['pres'][-1]
    errors['pressure_at_rupture'] = ((sim_pressure_at_rupture - exp_pressure_at_rupture) /
                                      exp_pressure_at_rupture * 100)

    # Final liquid temperature error
    sim_temp_at_rupture = results['hd'].T_liquid[idx_rupture] - 273.15
    exp_temp_at_rupture = exp_data['liquid_temp_data']['temp'][-1]
    errors['temp_at_rupture'] = ((sim_temp_at_rupture - exp_temp_at_rupture) /
                                  (exp_temp_at_rupture + 273.15) * 100)

    return errors

def main():
    """Main validation script"""

    experiments = {
        'exp1': os.path.join(_NEM_DIR, 'nem_droste_exp1.yml'),
        'exp2': os.path.join(_NEM_DIR, 'nem_droste_exp2.yml'),
        'exp3': os.path.join(_NEM_DIR, 'nem_droste_exp3.yml'),
    }

    all_results = {}
    all_errors = {}

    # Run all experiments
    for exp_key, yaml_file in experiments.items():
        results = run_experiment(yaml_file, exp_key)
        all_results[exp_key] = results

        # Calculate errors
        errors = calculate_errors(exp_key, results, EXPERIMENTAL_DATA[exp_key])
        all_errors[exp_key] = errors

        # Create validation plot
        fig = create_validation_plot(exp_key, results, EXPERIMENTAL_DATA[exp_key])
        filename = f"droste_{exp_key}_validation.png"
        fig.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"\nValidation plot saved: {filename}")
        plt.close(fig)

    # Print summary table
    print(f"\n{'='*100}")
    print("VALIDATION SUMMARY TABLE")
    print(f"{'='*100}")
    print(f"{'Experiment':<15} {'PSV Open [s]':<15} {'Rupture [s]':<15} {'P Error [%]':<12} {'T Error [%]':<12}")
    print(f"{'':<15} {'Sim / Exp / Err%':<15} {'Sim / Exp / Err%':<15} {'(at rupture)':<12} {'(at rupture)':<12}")
    print(f"{'-'*100}")

    for exp_key in ['exp1', 'exp2', 'exp3']:
        exp_data = EXPERIMENTAL_DATA[exp_key]
        results = all_results[exp_key]
        errors = all_errors[exp_key]

        psv_str = f"{results['psv_opening_time']:.0f} / {exp_data['psv_opening_time']} / {errors['psv_opening']:.1f}%" if errors['psv_opening'] else "N/A"

        if results['rupture_time']:
            rupt_str = f"{results['rupture_time']:.0f} / {exp_data['rupture_time']} / {errors['rupture_time']:.1f}%"
        else:
            rupt_str = f"No / {exp_data['rupture_time']} / N/A"

        print(f"{exp_data['name']:<15} {psv_str:<15} {rupt_str:<15} {errors['pressure_at_rupture']:>10.1f}% {errors['temp_at_rupture']:>10.1f}%")

    print(f"{'='*100}")

    # Tuned parameters summary
    print(f"\nTUNED PARAMETERS:")
    print(f"  PSV diameter: 0.016 m (DN16, ~5/8 inch)")
    print(f"  Fire model: scandpower_pool (100 kW/m² base)")
    print(f"  Fire scaling: 0.5 (50% engulfment, effective 50 kW/m²)")
    print(f"  Time step: 0.1 s")
    print(f"  Non-equilibrium model: Enabled (h_gas_liquid: calc)")

    print(f"\n{'='*100}")
    print("Validation complete! All plots saved.")
    print(f"{'='*100}\n")

if __name__ == "__main__":
    main()
