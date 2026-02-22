# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Runtime safety checks for HydDown simulations.

This module provides defensive checks for issues that can ONLY be detected
during simulation execution, not from static input validation. Cerberus
validator handles all input file validation - this module handles runtime
numerical and physical errors.

Runtime-only checks:
- Negative values from numerical errors (mass, pressure, temperature)
- NaN/Inf detection during calculation
- Array bounds checking
- CFL stability condition (depends on computed flow rates)
- Triple point violations during discharge
- Thermodynamic solver convergence
- Division by zero protection
"""

import warnings
import numpy as np
from CoolProp.CoolProp import PropsSI
from hyddown.exceptions import (
    NegativeMassError,
    TriplePointViolation,
    InvalidStateError,
    NumericalInstabilityError,
    NumericalStabilityWarning,
    ThermodynamicConvergenceError,
    ConvergenceWarning,
)


def check_positive_runtime(value, name, time=None, step=None):
    """
    Check that a computed value is positive during simulation.

    This checks for negative values that arise from numerical errors during
    time integration, NOT from invalid input (Cerberus handles input validation).

    Parameters
    ----------
    value : float
        Computed value to check
    name : str
        Name of the variable (for error messages)
    time : float, optional
        Simulation time when check occurs [s]
    step : int, optional
        Time step number when check occurs

    Raises
    ------
    NegativeMassError
        If computed value is negative or zero

    Examples
    --------
    >>> check_positive_runtime(100000, "pressure", time=5.0, step=100)
    >>> check_positive_runtime(-50, "mass", time=10.0)  # Raises error
    """
    if value <= 0:
        msg = f"{name} became non-positive during simulation: {value:.3e}"
        if time is not None:
            msg += f" at t={time:.3f}s"
        if step is not None:
            msg += f" (step {step})"
        msg += ". This indicates numerical instability - try reducing time step."

        raise NegativeMassError(
            msg,
            variable=name,
            value=value,
            time=time
        )


def check_array_bounds(index, array_length, context=""):
    """
    Check array index is within bounds during simulation.

    This prevents IndexError crashes from array access bugs.
    Particularly important for relief valve smoothing and similar operations.

    Parameters
    ----------
    index : int
        Index to check
    array_length : int
        Length of the array
    context : str, optional
        Description of operation (for error messages)

    Raises
    ------
    IndexError
        If index is out of bounds

    Examples
    --------
    >>> check_array_bounds(5, 10, "mass_rate smoothing")
    >>> check_array_bounds(10, 10)  # Raises IndexError
    """
    if index < 0 or index >= array_length:
        msg = f"Index {index} out of bounds for array of length {array_length}"
        if context:
            msg += f" in {context}"
        msg += f". Valid range is [0, {array_length-1}]"
        raise IndexError(msg)


def check_nan_inf(value, name, time=None, step=None):
    """
    Check for NaN or Inf in computed values.

    NaN/Inf indicate serious numerical problems (division by zero,
    overflow, etc.) and must be caught immediately.

    Parameters
    ----------
    value : float or np.ndarray
        Computed value(s) to check
    name : str
        Name of the variable
    time : float, optional
        Simulation time [s]
    step : int, optional
        Time step number

    Raises
    ------
    ValueError
        If value contains NaN or Inf

    Examples
    --------
    >>> check_nan_inf(5.0, "temperature")
    >>> check_nan_inf(np.nan, "pressure", time=10.0)  # Raises ValueError
    """
    has_nan = np.any(np.isnan(value))
    has_inf = np.any(np.isinf(value))

    if has_nan or has_inf:
        problem = "NaN" if has_nan else "Inf"
        msg = f"{name} contains {problem}: {value}"
        if time is not None:
            msg += f" at t={time:.3f}s"
        if step is not None:
            msg += f" (step {step})"
        msg += ". This indicates numerical instability."
        raise ValueError(msg)


def check_triple_point_runtime(temperature, pressure, fluid_name, time=None):
    """
    Check if state is above triple point during simulation.

    Triple point violations can occur during discharge even if initial
    conditions are valid. This is a physical constraint violation.

    Parameters
    ----------
    temperature : float
        Current temperature [K]
    pressure : float
        Current pressure [Pa]
    fluid_name : str
        CoolProp fluid name (e.g., "CO2", "N2")
    time : float, optional
        Simulation time [s]

    Raises
    ------
    TriplePointViolation
        If temperature or pressure is below triple point

    Warnings
    --------
    Issues warning if within 5% of triple point

    Examples
    --------
    >>> check_triple_point_runtime(300, 101325, "N2")  # OK
    >>> check_triple_point_runtime(200, 101325, "CO2")  # Raises error
    """
    try:
        # Get triple point properties (not all fluids have this data)
        T_triple = PropsSI("Ttriple", fluid_name)
        P_triple = PropsSI("ptriple", fluid_name)

        # Check temperature
        if temperature < T_triple:
            msg = (f"{fluid_name} temperature {temperature:.2f}K dropped below "
                   f"triple point {T_triple:.2f}K during simulation")
            if time is not None:
                msg += f" at t={time:.3f}s"
            msg += ". Solid phase would form - EOS invalid."

            raise TriplePointViolation(
                msg,
                fluid=fluid_name,
                temperature=temperature,
                pressure=pressure,
                T_triple=T_triple,
                P_triple=P_triple
            )

        # Check pressure
        if pressure < P_triple:
            msg = (f"{fluid_name} pressure {pressure:.2e}Pa dropped below "
                   f"triple point {P_triple:.2e}Pa during simulation")
            if time is not None:
                msg += f" at t={time:.3f}s"
            msg += ". Solid phase would form - EOS invalid."

            raise TriplePointViolation(
                msg,
                fluid=fluid_name,
                temperature=temperature,
                pressure=pressure,
                T_triple=T_triple,
                P_triple=P_triple
            )

        # Warning if close to triple point (within 5%)
        if temperature < T_triple * 1.05:
            msg = (f"{fluid_name} temperature {temperature:.2f}K is within 5% "
                   f"of triple point {T_triple:.2f}K")
            if time is not None:
                msg += f" at t={time:.3f}s"
            msg += ". Results may be inaccurate near phase boundary."
            warnings.warn(msg, NumericalStabilityWarning)

    except ValueError:
        # Fluid doesn't have triple point data in CoolProp, skip check
        pass


def check_cfl_stability(mass_flow_rate, vessel_mass, time_step,
                        characteristic_fraction=0.1):
    """
    Check CFL (Courant-Friedrichs-Lewy) stability condition.

    For explicit Euler integration, the time step must be small compared
    to the characteristic time scale (mass inventory / mass flow rate).
    This check can only be done at runtime with actual computed flow rates.

    Parameters
    ----------
    mass_flow_rate : float
        Current mass flow rate [kg/s] (absolute value)
    vessel_mass : float
        Current mass in vessel [kg]
    time_step : float
        Integration time step [s]
    characteristic_fraction : float, default 0.1
        Fraction of characteristic time for stability (typically 0.05 to 0.2)

    Warnings
    --------
    Issues NumericalStabilityWarning if time step is too large

    Returns
    -------
    float
        Recommended maximum time step [s], or None if no limit

    Examples
    --------
    >>> check_cfl_stability(0.1, 50.0, 1.0)  # May warn
    >>> check_cfl_stability(0.001, 50.0, 0.01)  # OK
    """
    # Avoid division by zero for very small flow rates
    if abs(mass_flow_rate) < 1e-10:
        return None  # No significant flow, any time step is OK

    # Characteristic time: how long to change vessel mass significantly
    # tau = mass / dmdt
    characteristic_time = vessel_mass / abs(mass_flow_rate)

    # Recommended maximum time step
    dt_max_recommended = characteristic_fraction * characteristic_time

    # Warn if current time step is too large
    if time_step > dt_max_recommended:
        warnings.warn(
            f"Time step dt={time_step:.3f}s may be too large for stability. "
            f"Recommended: dt < {dt_max_recommended:.3f}s "
            f"(based on characteristic time {characteristic_time:.3f}s). "
            f"Current mass flow rate = {abs(mass_flow_rate):.3e} kg/s, "
            f"vessel mass = {vessel_mass:.3e} kg. "
            f"Consider reducing time step if results show oscillations.",
            NumericalStabilityWarning,
            stacklevel=2
        )

    return dt_max_recommended


def check_optimization_convergence(result, solver_name, state_vars=None,
                                   tolerance=1e-4):
    """
    Check if scipy optimization converged successfully.

    This checks convergence of thermodynamic solvers (PHproblem, UDproblem)
    which use numerical optimization for multicomponent fluids. Convergence
    can only be checked at runtime.

    Parameters
    ----------
    result : scipy.optimize.OptimizeResult
        Result object from scipy.optimize.minimize or root_scalar
    solver_name : str
        Name of the solver (for error messages)
    state_vars : dict, optional
        Input state variables being solved for
    tolerance : float, default 1e-4
        Acceptable residual tolerance

    Raises
    ------
    ThermodynamicConvergenceError
        If optimization did not converge

    Warnings
    --------
    Issues ConvergenceWarning if convergence was marginal

    Examples
    --------
    >>> from scipy.optimize import minimize
    >>> result = minimize(lambda x: x**2, x0=1.0)
    >>> check_optimization_convergence(result, "test_solver")
    """
    # Check if optimization succeeded
    if hasattr(result, 'success'):
        if not result.success:
            msg = f"{solver_name} failed to converge: {result.message}"
            if state_vars:
                msg += f"\nInput state: {state_vars}"
            raise ThermodynamicConvergenceError(
                msg,
                solver=solver_name,
                state_vars=state_vars,
                iterations=getattr(result, 'nit', None)
            )

    # Check converged flag (for root_scalar)
    if hasattr(result, 'converged'):
        if not result.converged:
            msg = f"{solver_name} did not converge"
            if state_vars:
                msg += f" for state: {state_vars}"
            raise ThermodynamicConvergenceError(
                msg,
                solver=solver_name,
                state_vars=state_vars,
                iterations=getattr(result, 'iterations', None)
            )

    # Warn if residual is large (marginal convergence)
    if hasattr(result, 'fun'):
        if isinstance(result.fun, (list, np.ndarray)):
            residual = np.max(np.abs(result.fun))
        else:
            residual = abs(result.fun)

        if residual > tolerance:
            warnings.warn(
                f"{solver_name} converged but residual {residual:.2e} "
                f"exceeds tolerance {tolerance:.2e}. Results may be inaccurate.",
                ConvergenceWarning,
                stacklevel=2
            )

    # Warn if required many iterations
    nit = getattr(result, 'nit', getattr(result, 'iterations', None))
    if nit is not None and nit > 100:
        warnings.warn(
            f"{solver_name} required {nit} iterations. "
            "Consider improving initial guess for better performance.",
            ConvergenceWarning,
            stacklevel=2
        )


def safe_divide(numerator, denominator, name="division", default=0.0):
    """
    Safely divide two numbers, returning default if denominator is zero.

    This prevents division by zero errors that occur at runtime (e.g.,
    when wetted_area = 0 in two-phase calculations).

    Parameters
    ----------
    numerator : float
        Numerator value
    denominator : float
        Denominator value
    name : str, optional
        Description of the division operation (for warnings)
    default : float, default 0.0
        Value to return if denominator is zero or very small

    Returns
    -------
    float
        numerator / denominator, or default if division would fail

    Examples
    --------
    >>> safe_divide(10.0, 2.0)
    5.0
    >>> safe_divide(10.0, 0.0, default=0.0)
    0.0
    >>> safe_divide(10.0, 1e-20, name="heat flux")
    0.0
    """
    # Check for zero or near-zero denominator
    if abs(denominator) < 1e-15:
        # Silently return default - this is expected behavior
        # (e.g., wetted_area = 0 is normal for fully vapor state)
        return default

    return numerator / denominator


def check_bounds_for_smoothing(idx, array_length, context=""):
    """
    Check that index and neighbors are within bounds for smoothing operations.

    Relief valve smoothing and similar operations need idx-1, idx, idx+1.
    This checks all three indices are valid.

    Parameters
    ----------
    idx : int
        Center index to smooth
    array_length : int
        Array length
    context : str, optional
        Description of operation

    Returns
    -------
    bool
        True if smoothing is safe (idx-1, idx, idx+1 all valid)

    Examples
    --------
    >>> check_bounds_for_smoothing(5, 10)
    True
    >>> check_bounds_for_smoothing(0, 10)  # Can't access idx-1
    False
    >>> check_bounds_for_smoothing(9, 10)  # Can't access idx+1
    False
    """
    # Need idx-1, idx, and idx+1 to be valid
    if idx < 1 or idx >= array_length - 1:
        return False
    return True
