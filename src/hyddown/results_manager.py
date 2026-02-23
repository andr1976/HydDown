# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Results storage and management for HydDown simulations.

This module provides the ResultsManager class which handles all time-series
data storage, retrieval, and export for HydDown simulations. It encapsulates
the numpy array management and provides clean interfaces for accessing results.

Key Features:
- Automatic array initialization based on simulation duration
- Dictionary-like access to time-series data
- DataFrame export for analysis
- Timestep tracking and validation
- Memory-efficient numpy array storage

Typical usage:
    results = ResultsManager(time_total=100, time_step=0.1)
    results['pressure'][i] = 1e6
    results['temperature'][i] = 300
    df = results.to_dataframe()
"""

import numpy as np
import pandas as pd


class ResultsManager:
    """
    Manage time-series results storage for HydDown simulations.

    This class encapsulates all result arrays and provides clean interfaces
    for storing, accessing, and exporting simulation results. It handles:
    - Numpy array initialization and management
    - Dictionary-like access to results
    - DataFrame export
    - Timestep tracking

    Parameters
    ----------
    time_total : float
        Total simulation time [s]
    time_step : float
        Time step size [s]

    Attributes
    ----------
    data_len : int
        Number of time steps in simulation
    current_step : int
        Current time step index (for validation)

    Examples
    --------
    >>> results = ResultsManager(time_total=100, time_step=0.1)
    >>> results['pressure'][0] = 1e6
    >>> results['temperature'][0] = 300
    >>> df = results.to_dataframe()
    """

    def __init__(self, time_total, time_step):
        """
        Initialize results storage arrays.

        Parameters
        ----------
        time_total : float
            Total simulation time [s]
        time_step : float
            Time step size [s]
        """
        self.time_total = time_total
        self.time_step = time_step
        self.data_len = int(time_total / time_step)
        self.current_step = 0

        # ====================================================================
        # THERMODYNAMIC STATE ARRAYS
        # ====================================================================
        self.rho = np.zeros(self.data_len)  # Fluid density [kg/m³]
        self.P = np.zeros(self.data_len)  # Pressure [Pa]
        self.T_fluid = np.zeros(self.data_len)  # Fluid temperature [K]
        self.H_mass = np.zeros(self.data_len)  # Specific enthalpy [J/kg]
        self.S_mass = np.zeros(self.data_len)  # Specific entropy [J/(kg·K)]
        self.U_mass = np.zeros(self.data_len)  # Specific internal energy [J/kg]
        self.U_tot = np.zeros(self.data_len)  # Total internal energy [J]
        self.U_res = np.zeros(self.data_len)  # Reservoir internal energy [J]

        # ====================================================================
        # MASS FLOW AND INVENTORY ARRAYS
        # ====================================================================
        self.mass_fluid = np.zeros(self.data_len)  # Fluid mass [kg]
        self.mass_rate = np.zeros(self.data_len)  # Mass flow rate [kg/s]
        self.relief_area = np.zeros(self.data_len)  # Relief valve area [m²]

        # ====================================================================
        # TWO-PHASE PROPERTIES ARRAYS
        # ====================================================================
        self.vapour_mole_fraction = np.zeros(self.data_len)  # Vapor mole fraction [-]
        self.vapour_mass_fraction = np.zeros(self.data_len)  # Vapor mass fraction [-]
        self.vapour_volume_fraction = np.zeros(self.data_len)  # Vapor volume fraction [-]
        self.liquid_level = np.zeros(self.data_len)  # Liquid level height [m]

        # ====================================================================
        # WALL TEMPERATURE ARRAYS
        # ====================================================================
        self.T_vessel = np.zeros(self.data_len)  # Vessel wall temp (unwetted) [K]
        self.T_vessel_wetted = np.zeros(self.data_len)  # Vessel wall temp (wetted) [K]
        self.T_inner_wall = np.zeros(self.data_len)  # Inner wall temp (unwetted) [K]
        self.T_inner_wall_wetted = np.zeros(self.data_len)  # Inner wall temp (wetted) [K]
        self.T_outer_wall = np.zeros(self.data_len)  # Outer wall temp (unwetted) [K]
        self.T_outer_wall_wetted = np.zeros(self.data_len)  # Outer wall temp (wetted) [K]
        self.T_bonded_wall = np.zeros(self.data_len)  # Bonded interface temp (unwetted) [K]
        self.T_bonded_wall_wetted = np.zeros(self.data_len)  # Bonded interface temp (wetted) [K]

        # ====================================================================
        # HEAT TRANSFER ARRAYS
        # ====================================================================
        self.Q_inner = np.zeros(self.data_len)  # Internal heat flux (unwetted) [W]
        self.Q_inner_wetted = np.zeros(self.data_len)  # Internal heat flux (wetted) [W]
        self.Q_outer = np.zeros(self.data_len)  # External heat flux (unwetted) [W]
        self.Q_outer_wetted = np.zeros(self.data_len)  # External heat flux (wetted) [W]
        self.q_inner = np.zeros(self.data_len)  # Internal heat flux density (unwetted) [W/m²]
        self.q_inner_wetted = np.zeros(self.data_len)  # Internal heat flux density (wetted) [W/m²]
        self.q_outer = np.zeros(self.data_len)  # External heat flux density (unwetted) [W/m²]
        self.q_outer_wetted = np.zeros(self.data_len)  # External heat flux density (wetted) [W/m²]
        self.h_inside = np.zeros(self.data_len)  # Internal heat transfer coef (unwetted) [W/(m²·K)]
        self.h_inside_wetted = np.zeros(self.data_len)  # Internal heat transfer coef (wetted) [W/(m²·K)]

        # ====================================================================
        # OTHER ARRAYS
        # ====================================================================
        self.T_vent = np.zeros(self.data_len)  # Vent temperature [K]
        self.time_array = np.zeros(self.data_len)  # Time array [s]
        self.temp_profile = []  # Temperature profiles through wall thickness

    def __getitem__(self, key):
        """
        Dictionary-like access to result arrays.

        Parameters
        ----------
        key : str
            Name of the result array

        Returns
        -------
        np.ndarray
            The requested result array

        Examples
        --------
        >>> results['pressure'][0] = 1e6
        >>> P = results['pressure']
        """
        return getattr(self, key)

    def __setitem__(self, key, value):
        """
        Dictionary-like setting of result arrays.

        Parameters
        ----------
        key : str
            Name of the result array
        value : array-like
            Values to set

        Examples
        --------
        >>> results['pressure'] = np.array([1e6, 1.1e6, 1.2e6])
        """
        setattr(self, key, value)

    def get_all_arrays(self):
        """
        Get dictionary of all result arrays.

        Returns
        -------
        dict
            Dictionary mapping array names to numpy arrays

        Examples
        --------
        >>> arrays = results.get_all_arrays()
        >>> print(arrays.keys())
        """
        return {
            # Thermodynamic state
            'rho': self.rho,
            'P': self.P,
            'T_fluid': self.T_fluid,
            'H_mass': self.H_mass,
            'S_mass': self.S_mass,
            'U_mass': self.U_mass,
            'U_tot': self.U_tot,
            'U_res': self.U_res,
            # Mass flow and inventory
            'mass_fluid': self.mass_fluid,
            'mass_rate': self.mass_rate,
            'relief_area': self.relief_area,
            # Two-phase properties
            'vapour_mole_fraction': self.vapour_mole_fraction,
            'vapour_mass_fraction': self.vapour_mass_fraction,
            'vapour_volume_fraction': self.vapour_volume_fraction,
            'liquid_level': self.liquid_level,
            # Wall temperatures
            'T_vessel': self.T_vessel,
            'T_vessel_wetted': self.T_vessel_wetted,
            'T_inner_wall': self.T_inner_wall,
            'T_inner_wall_wetted': self.T_inner_wall_wetted,
            'T_outer_wall': self.T_outer_wall,
            'T_outer_wall_wetted': self.T_outer_wall_wetted,
            'T_bonded_wall': self.T_bonded_wall,
            'T_bonded_wall_wetted': self.T_bonded_wall_wetted,
            # Heat transfer
            'Q_inner': self.Q_inner,
            'Q_inner_wetted': self.Q_inner_wetted,
            'Q_outer': self.Q_outer,
            'Q_outer_wetted': self.Q_outer_wetted,
            'q_inner': self.q_inner,
            'q_inner_wetted': self.q_inner_wetted,
            'q_outer': self.q_outer,
            'q_outer_wetted': self.q_outer_wetted,
            'h_inside': self.h_inside,
            'h_inside_wetted': self.h_inside_wetted,
            # Other
            'T_vent': self.T_vent,
            'time_array': self.time_array,
        }

    def to_dataframe(self, surf_area_inner=None, surf_area_outer=None):
        """
        Export results to pandas DataFrame.

        Converts all time-series results into a pandas DataFrame with
        labeled columns. Temperature values converted from K to °C.
        Pressure values converted from Pa to bar.

        Parameters
        ----------
        surf_area_inner : float, optional
            Inner surface area [m²]. Required for heat flux density calculations.
        surf_area_outer : float, optional
            Outer surface area [m²]. Required for heat flux density calculations.

        Returns
        -------
        pd.DataFrame
            DataFrame with simulation results

        Notes
        -----
        Only includes data up to current_step if simulation is incomplete.
        """
        # Determine how much data to include
        end_idx = self.current_step if self.current_step > 0 else self.data_len

        # Create DataFrame with time as first column
        df = pd.DataFrame(self.time_array[:end_idx], columns=["Time (s)"])

        # Add thermodynamic state
        df.insert(1, "Pressure (bar)", self.P[:end_idx] / 1e5, True)
        df.insert(2, "Fluid temperature (oC)", self.T_fluid[:end_idx] - 273.15, True)
        df.insert(3, "Wall temperature  (oC)", self.T_vessel[:end_idx] - 273.15, True)
        df.insert(4, "Vent temperature  (oC)", self.T_vent[:end_idx] - 273.15, True)
        df.insert(5, "Fluid enthalpy (J/kg)", self.H_mass[:end_idx], True)
        df.insert(6, "Fluid entropy (J/kg K)", self.S_mass[:end_idx], True)
        df.insert(7, "Fluid internal energy (J/kg)", self.U_mass[:end_idx], True)

        # Add mass flow and inventory
        df.insert(8, "Discharge mass rate (kg/s)", self.mass_rate[:end_idx], True)
        df.insert(9, "Fluid mass (kg)", self.mass_fluid[:end_idx], True)
        df.insert(10, "Fluid density (kg/m3)", self.rho[:end_idx], True)

        # Add heat transfer coefficients
        df.insert(
            11, "Inner heat transfer coefficient (W/m2 K)", self.h_inside[:end_idx], True
        )

        # Add heat flux densities (if surface areas provided)
        if surf_area_inner is not None:
            df.insert(
                12,
                "Internal heat flux (W/m2)",
                self.Q_inner[:end_idx] / surf_area_inner,
                True,
            )

        if surf_area_outer is not None:
            df.insert(
                13,
                "External heat flux (W/m2)",
                self.Q_outer[:end_idx] / surf_area_outer,
                True,
            )

        # Add wall temperatures
        df.insert(
            14, "Inner wall temperature  (oC)", self.T_inner_wall[:end_idx] - 273.15, True
        )
        df.insert(
            15, "Outer wall temperature  (oC)", self.T_outer_wall[:end_idx] - 273.15, True
        )

        return df

    def get_array_names(self):
        """
        Get list of all array names.

        Returns
        -------
        list of str
            Names of all result arrays
        """
        return list(self.get_all_arrays().keys())

    def trim_to_current_step(self):
        """
        Trim all arrays to current step (for incomplete simulations).

        This is useful when simulation ends early or is interrupted.
        """
        for name, array in self.get_all_arrays().items():
            if isinstance(array, np.ndarray):
                setattr(self, name, array[:self.current_step])
