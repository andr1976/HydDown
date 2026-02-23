# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Custom exceptions for HydDown package.

This module defines a hierarchy of exceptions for different error conditions
that can occur during pressure vessel simulations. Using specific exception
types allows calling code to handle different error cases appropriately.

Exception Hierarchy:
    HydDownError (base)
    ├── ThermodynamicError
    │   ├── ThermodynamicConvergenceError
    │   ├── InvalidStateError
    │   └── PhaseEquilibriumError
    ├── NumericalError
    │   ├── NumericalInstabilityError
    │   ├── IntegrationError
    │   └── ConvergenceError
    ├── ConfigurationError
    │   ├── ValveConfigurationError
    │   ├── VesselConfigurationError
    │   └── HeatTransferConfigurationError
    └── PhysicalConstraintError
        ├── TriplePointViolation
        ├── CriticalPointViolation
        └── NegativeMassError

Usage:
    from hyddown.exceptions import ThermodynamicConvergenceError

    if not result.success:
        raise ThermodynamicConvergenceError(
            f"Failed to converge: {result.message}"
        )
"""


class HydDownError(Exception):
    """
    Base exception for all HydDown errors.

    All custom exceptions in HydDown inherit from this base class,
    allowing calling code to catch all HydDown-specific errors with
    a single except clause if desired.
    """
    pass


# ============================================================================
# Thermodynamic Errors
# ============================================================================

class ThermodynamicError(HydDownError):
    """
    Base class for thermodynamic calculation errors.

    Raised when CoolProp property calculations fail or produce
    non-physical results.
    """
    pass


class ThermodynamicConvergenceError(ThermodynamicError):
    """
    Failed to converge in thermodynamic property solver.

    Raised when numerical optimization fails to find a valid thermodynamic
    state (e.g., during PHproblem or UDproblem calculations). This typically
    occurs with multicomponent mixtures when scipy.optimize cannot find
    a solution.

    Attributes
    ----------
    solver : str
        Name of the solver that failed ('PHproblem', 'UDproblem', etc.)
    state_vars : dict
        Input state variables that caused failure
    iterations : int, optional
        Number of iterations attempted
    """

    def __init__(self, message, solver=None, state_vars=None, iterations=None):
        super().__init__(message)
        self.solver = solver
        self.state_vars = state_vars
        self.iterations = iterations


class InvalidStateError(ThermodynamicError):
    """
    Thermodynamic state is outside valid range for equation of state.

    Raised when requested state conditions are outside the valid range
    of the CoolProp equation of state (e.g., below triple point, above
    maximum temperature, or at unphysical densities).

    Attributes
    ----------
    variable : str
        The state variable that is out of range ('temperature', 'pressure', etc.)
    value : float
        The invalid value
    valid_range : tuple
        (min, max) valid range for the variable
    """

    def __init__(self, message, variable=None, value=None, valid_range=None):
        super().__init__(message)
        self.variable = variable
        self.value = value
        self.valid_range = valid_range


class PhaseEquilibriumError(ThermodynamicError):
    """
    Error in phase equilibrium calculations for two-phase systems.

    Raised when calculations involving vapor-liquid equilibrium fail,
    such as when trying to determine liquid level or saturation properties.
    """
    pass


# ============================================================================
# Numerical Errors
# ============================================================================

class NumericalError(HydDownError):
    """
    Base class for numerical integration and solver errors.
    """
    pass


class NumericalInstabilityError(NumericalError):
    """
    Time step too large for numerical stability.

    Raised when the explicit Euler integration scheme becomes unstable,
    typically because the time step is too large relative to the
    characteristic time scales of the problem. This can manifest as
    negative masses, temperatures, or pressures.

    Attributes
    ----------
    time_step : float
        The time step that caused instability [s]
    recommended_dt : float, optional
        Recommended smaller time step [s]
    characteristic_time : float, optional
        Characteristic time scale of the problem [s]
    """

    def __init__(self, message, time_step=None, recommended_dt=None,
                 characteristic_time=None):
        super().__init__(message)
        self.time_step = time_step
        self.recommended_dt = recommended_dt
        self.characteristic_time = characteristic_time


class IntegrationError(NumericalError):
    """
    Error during time integration of mass/energy balances.

    Raised when the numerical integration produces invalid results,
    such as NaN or Inf values in state variables.
    """
    pass


class ConvergenceError(NumericalError):
    """
    General convergence failure in iterative solvers.

    Raised when iterative algorithms (not thermodynamic solvers)
    fail to converge within specified tolerance or iteration limit.
    """
    pass


# ============================================================================
# Configuration Errors
# ============================================================================

class ConfigurationError(HydDownError):
    """
    Base class for invalid input configuration errors.

    These errors indicate problems with the input YAML file beyond
    what the Cerberus schema validator can catch (e.g., physically
    impossible values or inconsistent parameters).
    """
    pass


class ValveConfigurationError(ConfigurationError):
    """
    Invalid valve parameters.

    Raised when valve configuration is physically impossible or
    inconsistent (e.g., negative diameter, discharge coefficient > 1,
    relief valve set pressure below operating pressure).

    Attributes
    ----------
    valve_type : str
        Type of valve ('orifice', 'relief_valve', 'control_valve')
    parameter : str
        The problematic parameter name
    value : float
        The invalid value
    """

    def __init__(self, message, valve_type=None, parameter=None, value=None):
        super().__init__(message)
        self.valve_type = valve_type
        self.parameter = parameter
        self.value = value


class VesselConfigurationError(ConfigurationError):
    """
    Invalid vessel geometry or material parameters.

    Raised when vessel configuration is physically impossible
    (e.g., negative dimensions, wall thickness greater than diameter,
    invalid material properties).
    """
    pass


class HeatTransferConfigurationError(ConfigurationError):
    """
    Invalid heat transfer parameters.

    Raised when heat transfer configuration is inconsistent or
    physically impossible (e.g., negative heat transfer coefficient,
    invalid fire scenario).
    """
    pass


# ============================================================================
# Physical Constraint Errors
# ============================================================================

class PhysicalConstraintError(HydDownError):
    """
    Base class for violations of physical constraints.

    These errors indicate that the simulation has reached a state
    that violates fundamental physical constraints.
    """
    pass


class TriplePointViolation(PhysicalConstraintError):
    """
    Fluid state below triple point temperature or pressure.

    Raised when calculations result in conditions below the fluid's
    triple point, where the equation of state is not valid and
    solid phase would form.

    Attributes
    ----------
    fluid : str
        Name of the fluid
    temperature : float
        Current temperature [K]
    pressure : float
        Current pressure [Pa]
    T_triple : float
        Triple point temperature [K]
    P_triple : float
        Triple point pressure [Pa]
    """

    def __init__(self, message, fluid=None, temperature=None, pressure=None,
                 T_triple=None, P_triple=None):
        super().__init__(message)
        self.fluid = fluid
        self.temperature = temperature
        self.pressure = pressure
        self.T_triple = T_triple
        self.P_triple = P_triple


class CriticalPointViolation(PhysicalConstraintError):
    """
    Fluid state exceeds equation of state validity limits.

    Raised when calculations produce conditions far beyond the
    critical point or other EOS validity limits.
    """
    pass


class NegativeMassError(PhysicalConstraintError):
    """
    Mass, pressure, or temperature became negative.

    Raised when numerical integration produces unphysical negative
    values for extensive properties. This usually indicates numerical
    instability or time step too large.

    Attributes
    ----------
    variable : str
        The variable that became negative ('mass', 'pressure', 'temperature')
    value : float
        The negative value
    time : float
        Simulation time when error occurred [s]
    """

    def __init__(self, message, variable=None, value=None, time=None):
        super().__init__(message)
        self.variable = variable
        self.value = value
        self.time = time


# ============================================================================
# Warning Classes (for non-fatal issues)
# ============================================================================

class HydDownWarning(UserWarning):
    """
    Base warning class for HydDown.

    Warnings are issued for conditions that may affect accuracy but
    don't prevent calculation from completing.
    """
    pass


class NumericalStabilityWarning(HydDownWarning):
    """
    Warning that time step may be too large for optimal stability.

    Issued when CFL condition or other stability criteria suggest
    the time step should be reduced, but calculation can continue.
    """
    pass


class AccuracyWarning(HydDownWarning):
    """
    Warning that results may have reduced accuracy.

    Issued when approximations or numerical issues may affect
    accuracy but not prevent calculation completion.
    """
    pass


class ConvergenceWarning(HydDownWarning):
    """
    Warning about slow or marginal convergence.

    Issued when iterative solvers converge but required many
    iterations or achieved marginal tolerance.
    """
    pass
