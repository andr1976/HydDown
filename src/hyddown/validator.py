# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
Input validation for HydDown YAML configuration files.

This module uses the Cerberus validation library to enforce schema-based validation
of user-provided input dictionaries parsed from YAML files. It ensures that all
required parameters are present, values are of correct types, and fall within
acceptable ranges before calculations begin.

The validation is hierarchical and mode-dependent:
- Mandatory ruleset: vessel geometry, initial conditions, calculation setup
- Heat transfer validation: depends on calculation type (energybalance, etc.)
- Valve validation: depends on valve type (orifice, control valve, relief valve)

Validation errors are reported with descriptive messages to help users correct
their input files. The module prevents invalid calculations that would lead to
runtime errors or physically meaningless results.

Main functions:
- validation(): Top-level validation function called by HydDown class
- validate_mandatory_ruleset(): Validates core required parameters
- heat_transfer_validation(): Validates heat transfer settings
- valve_validation(): Validates valve parameters based on valve type
"""

from cerberus import Validator
from cerberus.errors import ValidationError
from hyddown.exceptions import (
    ConfigurationError,
    VesselConfigurationError,
    ValveConfigurationError,
    HeatTransferConfigurationError,
)


# ============================================================================
# REUSABLE SCHEMA COMPONENTS (Phase 1.5 Refactoring)
# ============================================================================

# Material property constraints (consistent across all uses)
MATERIAL_PROPERTIES = {
    "heat_capacity": {
        "required": False,
        "type": "number",
        "min": 1,
        "max": 10000,  # Sanity check [J/(kg·K)]
    },
    "thermal_conductivity": {
        "required": False,
        "type": "number",
        "min": 0.0001,
        "max": 500,  # Max for metals ~400 W/(m·K)
    },
    "density": {
        "required": False,
        "type": "number",
        "min": 1,
        "max": 25000,  # Max: osmium ~22,000 kg/m³
    },
    "thickness": {
        "required": False,
        "type": "number",
        "min": 0.0001,  # Wall thickness must be > 0 when specified
    },
}


def create_liner_properties():
    """Create liner property schema (same constraints as base materials)."""
    return {
        "liner_thickness": {"required": False, "type": "number", "min": 0.0},
        "liner_heat_capacity": {"required": False, "type": "number", "min": 1},
        "liner_thermal_conductivity": {
            "required": False,
            "type": "number",
            "min": 0.0001,
        },
        "liner_density": {"required": False, "type": "number", "min": 1},
    }


# Time series validation schema (for validation data)
TIME_SERIES_SCHEMA = {
    "required": False,
    "type": "dict",
    "contains": ["time", "data"],
    "schema": {
        "data": {
            "required": False,
            "type": "list",
            "schema": {"type": "number"},
        },
        "time": {
            "required": False,
            "type": "list",
            "schema": {"type": "number"},
        },
    },
}


def create_validation_data_schema():
    """Create schema for validation data section (eliminates massive duplication)."""
    # Pressure validation data
    pressure_schema = {
        "type": "dict",
        "required": False,
        "contains": ["time", "pres"],
        "schema": {
            "pres": {
                "required": False,
                "type": "list",
                "schema": {"type": "number"},
            },
            "time": {
                "required": False,
                "type": "list",
                "schema": {"type": "number"},
            },
        },
    }

    # Temperature data schema (same for all temperature types)
    temp_schema = {
        "required": False,
        "type": "dict",
        "contains": ["time", "temp"],
        "schema": {
            "temp": {
                "required": False,
                "type": "list",
                "schema": {"type": "number"},
            },
            "time": {
                "required": False,
                "type": "list",
                "schema": {"type": "number"},
            },
        },
    }

    # All temperature types use the same schema
    temp_types = [
        "gas_high",
        "gas_low",
        "gas_mean",
        "wall_high",
        "wall_mean",
        "wall_low",
        "wall_inner",
        "wall_outer",
    ]

    temperature_schema = {
        "type": "dict",
        "required": False,
        "allowed": temp_types,
        "schema": {temp_type: temp_schema.copy() for temp_type in temp_types},
    }

    return {
        "pressure": pressure_schema,
        "temperature": temperature_schema,
    }


def create_vessel_schema(required_fields=None):
    """
    Create vessel schema with specified required fields.

    Parameters
    ----------
    required_fields : list, optional
        Fields that should be required. If None, all are optional.

    Returns
    -------
    dict
        Vessel schema with specified requirements
    """
    required_fields = required_fields or []

    schema = {
        "length": {
            "type": "number",
            "min": 0.001,  # Length must be > 0 (m)
        },
        "diameter": {
            "type": "number",
            "min": 0.001,  # Diameter must be > 0 (m)
        },
        "orientation": {
            "type": "string",
            "allowed": ["vertical", "horizontal"],
        },
        "type": {
            "type": "string",
            "allowed": ["Flat-end", "ASME F&D", "DIN", "Hemispherical"],
        },
        "liquid_level": {
            "type": "number",
            "min": 0,
        },
    }

    # Add material properties
    schema.update(MATERIAL_PROPERTIES.copy())
    schema.update(create_liner_properties())

    # Mark specified fields as required
    for field in required_fields:
        if field in schema:
            schema[field]["required"] = True

    # Set all non-required fields
    for field in schema:
        if "required" not in schema[field]:
            schema[field]["required"] = False

    return schema


def create_top_level_schema():
    """Create the common top-level section schema used across all validations."""
    return {
        "initial": {"required": True},
        "calculation": {"required": True},
        "validation": {"required": False},
        "rupture": {"required": False},
    }


def create_valve_schema(valve_type):
    """
    Create valve schema for specific valve type.

    Parameters
    ----------
    valve_type : str
        Valve type: 'relief', 'psv', 'orifice', 'controlvalve', 'mdot'

    Returns
    -------
    dict
        Valve schema with type-specific requirements
    """
    # Base valve schema (all possible fields)
    base_schema = {
        "type": {
            "required": True,
            "type": "string",
            "allowed": ["orifice", "psv", "controlvalve", "mdot", "relief"],
        },
        "flow": {
            "required": True,
            "type": "string",
            "allowed": ["discharge", "filling"],
        },
        "diameter": {"type": "number", "min": 0},
        "discharge_coef": {
            "type": "number",
            "min": 0.001,
            "max": 1.0,  # Physical limit
        },
        "set_pressure": {"type": "number", "min": 0},
        "end_pressure": {"type": "number", "min": 0},
        "blowdown": {"type": "number", "min": 0, "max": 1},
        "back_pressure": {"type": "number", "min": 0},
        "Cv": {"type": "number", "min": 0},
        "mdot": {"type": ["number", "list"]},
        "time": {"type": ["number", "list"]},
        "time_constant": {"type": "number", "min": 0},
        "characteristic": {
            "type": "string",
            "allowed": ["linear", "eq", "fast"],
        },
    }

    # Type-specific requirements
    if valve_type == "relief":
        base_schema["set_pressure"]["required"] = True
        base_schema["back_pressure"]["required"] = True

    elif valve_type == "psv":
        base_schema["diameter"]["required"] = True
        base_schema["discharge_coef"]["required"] = True
        base_schema["set_pressure"]["required"] = True
        base_schema["blowdown"]["required"] = True
        base_schema["back_pressure"]["required"] = True

    elif valve_type == "orifice":
        base_schema["diameter"]["required"] = True
        base_schema["discharge_coef"]["required"] = True
        base_schema["back_pressure"]["required"] = True

    elif valve_type == "controlvalve":
        base_schema["Cv"]["required"] = True
        base_schema["back_pressure"]["required"] = True

    elif valve_type == "mdot":
        base_schema["mdot"]["required"] = True
        base_schema["back_pressure"]["required"] = True

    return base_schema


def format_cerberus_errors(errors, section="input"):
    """
    Convert Cerberus error dict to user-friendly ConfigurationError.

    Parameters
    ----------
    errors : dict
        Cerberus validation errors
    section : str
        Section name (vessel, valve, heat_transfer, etc.)

    Raises
    ------
    ConfigurationError
        With formatted, actionable error message
    """
    error_messages = []

    def flatten_errors(err_dict, prefix=""):
        """Recursively flatten nested error dictionaries."""
        for field, error_list in err_dict.items():
            if isinstance(error_list, list):
                for error in error_list:
                    if isinstance(error, dict):
                        # Nested errors - recurse
                        flatten_errors(error, f"{prefix}{field}.")
                    else:
                        # Leaf error - format it
                        error_messages.append(f"  {prefix}{field}: {error}")
            else:
                error_messages.append(f"  {prefix}{field}: {error_list}")

    flatten_errors(errors)

    if not error_messages:
        error_messages = [f"  Unknown validation error: {errors}"]

    message = f"Invalid {section} configuration:\n" + "\n".join(error_messages)

    # Raise specific exception based on section
    if section == "vessel":
        raise VesselConfigurationError(message)
    elif section == "valve":
        raise ValveConfigurationError(message)
    elif section == "heat_transfer":
        raise HeatTransferConfigurationError(message)
    else:
        raise ConfigurationError(message)


# ============================================================================
# CROSS-FIELD VALIDATION (Physical Constraints)
# ============================================================================


def validate_relief_valve_physics(input):
    """
    Validate physical constraints for relief valves.

    Parameters
    ----------
    input : dict
        Structure holding input

    Raises
    ------
    ValveConfigurationError
        If physical constraints violated
    """
    if input["valve"]["type"] not in ["relief", "psv"]:
        return  # Not a relief valve

    valve = input["valve"]
    set_p = valve.get("set_pressure")
    back_p = valve.get("back_pressure")

    # Check set pressure > back pressure
    if set_p is not None and back_p is not None:
        if set_p <= back_p:
            raise ValveConfigurationError(
                f"Relief valve set_pressure ({set_p} Pa) must be greater than "
                f"back_pressure ({back_p} Pa).\n"
                f"The valve cannot open if set pressure is at or below back pressure."
            )

    # Check blowdown is reasonable
    blowdown = valve.get("blowdown")
    if blowdown is not None and blowdown >= 1.0:
        raise ValveConfigurationError(
            f"Relief valve blowdown ({blowdown}) must be < 1.0.\n"
            f"Typical values are 0.1-0.2 (10%-20%)."
        )


def validate_composite_vessel(input):
    """
    Validate composite vessel liner constraints.

    Parameters
    ----------
    input : dict
        Structure holding input

    Raises
    ------
    VesselConfigurationError
        If liner thickness >= wall thickness
    """
    vessel = input.get("vessel", {})
    thickness = vessel.get("thickness")
    liner_thickness = vessel.get("liner_thickness")

    if thickness and liner_thickness:
        if liner_thickness >= thickness:
            raise VesselConfigurationError(
                f"liner_thickness ({liner_thickness} m) must be less than "
                f"thickness ({thickness} m).\n"
                f"The liner cannot be thicker than the entire wall."
            )


def validate_time_step_sanity(input):
    """
    Validate time step is reasonable for simulation duration.

    Parameters
    ----------
    input : dict
        Structure holding input

    Raises
    ------
    ConfigurationError
        If time step is too large relative to end time
    """
    calc = input.get("calculation", {})
    time_step = calc.get("time_step")
    end_time = calc.get("end_time")

    if time_step and end_time:
        if time_step > end_time:
            raise ConfigurationError(
                f"time_step ({time_step} s) is greater than end_time ({end_time} s).\n"
                f"The simulation would complete in a single step."
            )

        # Warn if very few steps
        num_steps = end_time / time_step
        if num_steps < 10:
            import warnings

            warnings.warn(
                f"Only {num_steps:.1f} time steps will be simulated. "
                f"Consider reducing time_step for better resolution.",
                UserWarning,
            )


# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================


def validate_mandatory_ruleset(input):
    """
    Validate input

    Parameters
    ----------
    input : dict
        Structure holding input

    Return
    ----------
    retval : bool
        True for success, False for failure
    """
    schema_general = {
        "initial": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "temperature": {
                    "required": True,
                    "type": "number",
                    "min": 0.001,  # Kelvin scale, effectively > 0
                },
                "pressure": {
                    "required": True,
                    "type": "number",
                    "min": 0.001,  # Pressure must be > 0 (Pa)
                },
                "fluid": {"required": True, "type": "string"},
            },
        },
        "calculation": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": [
                        "energybalance",
                        "isenthalpic",
                        "isentropic",
                        "isothermal",
                        "specified_U",
                    ],
                },
                "time_step": {
                    "required": True,
                    "type": "number",
                    "min": 0.000001,
                    "max": 3600,  # Sanity check: max 1 hour per step
                },
                "end_time": {
                    "required": True,
                    "type": "number",
                    "min": 0.001,  # Must have positive duration
                },
            },
        },
        "vessel": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": create_vessel_schema(required_fields=["length", "diameter"]),
        },
        "rupture": {
            "required": False,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "material": {
                    "required": True,
                    "type": "string",
                    "allowed": ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"],
                },
                "fire": {"required": False, "type": "string"},
            },
        },
        "valve": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": ["orifice", "psv", "controlvalve", "mdot", "relief"],
                },
                "flow": {
                    "required": True,
                    "type": "string",
                    "allowed": ["discharge", "filling"],
                },
                "diameter": {"type": "number", "min": 0},  # Allow 0 for no-flow
                "discharge_coef": {
                    "type": "number",
                    "min": 0.001,  # > 0
                    "max": 1.0,  # Physical limit for ideal orifice
                },
                "set_pressure": {"type": "number", "min": 0},
                "end_pressure": {"type": "number", "min": 0},
                "blowdown": {"type": "number", "min": 0, "max": 1},
                "back_pressure": {"type": "number", "min": 0},
                "Cv": {"type": "number", "min": 0},
                "mdot": {"type": ["number", "list"]},
                "time": {"type": "list"},
                "time_constant": {"type": "number", "min": 0},
                "characteristic": {
                    "type": "string",
                    "allowed": ["linear", "eq", "fast"],
                },
            },
        },
        "heat_transfer": {
            "required": False,
            "type": "dict",
            "allow_unknown": False,
            "allowed": [
                "Q_fix",
                "U_fix",
                "h_inner",
                "h_outer",
                "temp_ambient",
                "type",
                "fire",
                "s-b",
                "D_throat",
            ],
            "schema": {
                "type": {
                    "type": "string",
                    "allowed": ["specified_Q", "specified_h", "specified_U", "s-b"],
                },
                "Q_fix": {"required": False, "type": "number"},
                "U_fix": {
                    "required": False,
                    "type": "number",
                    "min": 0.001,  # > 0
                    "max": 1000,  # Sanity check [W/(m²·K)]
                },
                "temp_ambient": {
                    "required": False,
                    "type": "number",
                    "min": 0.001,  # Kelvin scale > 0
                    "max": 2000,  # Max for fire scenarios
                },
                "h_outer": {
                    "required": False,
                    "type": "number",
                    "min": 0,
                    "max": 10000,  # Sanity check [W/(m²·K)]
                },
                "h_inner": {"required": False, "type": ["number", "string"]},
                "fire": {
                    "required": False,
                    "type": "string",
                    "allowed": [
                        "api_pool",
                        "api_jet",
                        "scandpower_pool",
                        "scandpower_jet",
                    ],
                },
                "D_throat": {"required": False, "type": "number", "min": 0},
            },
        },
        "validation": {
            "required": False,
            "type": "dict",
            "allow_unknown": False,
            "allowed": ["pressure", "temperature"],
            "schema": create_validation_data_schema(),
        },
    }

    v = Validator(schema_general)
    retval = v.validate(input)

    if not retval:
        # Raise specific exceptions instead of printing
        if "vessel" in v.errors:
            format_cerberus_errors({"vessel": v.errors["vessel"]}, "vessel")
        elif "valve" in v.errors:
            format_cerberus_errors({"valve": v.errors["valve"]}, "valve")
        elif "heat_transfer" in v.errors:
            format_cerberus_errors({"heat_transfer": v.errors["heat_transfer"]}, "heat_transfer")
        elif "initial" in v.errors:
            format_cerberus_errors({"initial": v.errors["initial"]}, "initial conditions")
        elif "calculation" in v.errors:
            format_cerberus_errors({"calculation": v.errors["calculation"]}, "calculation")
        else:
            # Generic configuration error for other cases
            format_cerberus_errors(v.errors, "input")

    return retval


def heat_transfer_validation(input):
    """
    Validate input['heat_transfer'] deeper than cerberus.

    Parameters
    ----------
    input : dict
        Structure holding input

    Return
    ----------
    bool
        True for success, False for failure
    """
    retval = True

    if input["calculation"]["type"] == "energybalance":
        ht_type = input["heat_transfer"]["type"]

        # Build schema using factory functions instead of massive duplication
        base_schema = create_top_level_schema()
        base_schema["vessel"] = {"required": True}
        base_schema["valve"] = {"required": True}
        base_schema["heat_transfer"] = {"required": True}

        if ht_type == "specified_h":
            # Vessel requirements for specified_h
            required_vessel_fields = [
                "length",
                "diameter",
                "thickness",
                "heat_capacity",
                "density",
                "orientation",
            ]

            schema_heattransfer = base_schema.copy()
            schema_heattransfer["vessel"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "schema": create_vessel_schema(required_fields=required_vessel_fields),
            }
            schema_heattransfer["heat_transfer"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "allowed": ["h_inner", "h_outer", "temp_ambient", "type", "D_throat"],
                "schema": {
                    "type": {"type": "string", "allowed": ["specified_h"]},
                    "temp_ambient": {
                        "required": True,
                        "type": "number",
                        "min": 0.001,
                        "max": 2000,
                    },
                    "h_outer": {
                        "required": True,
                        "type": "number",
                        "min": 0,
                        "max": 10000,
                    },
                    "h_inner": {"required": True, "type": ["number", "string"]},
                    "D_throat": {"required": False, "type": "number", "min": 0},
                },
            }

        elif ht_type == "specified_Q":
            # Vessel requirements for specified_Q (minimal)
            required_vessel_fields = ["length", "diameter"]

            schema_heattransfer = base_schema.copy()
            schema_heattransfer["vessel"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "schema": create_vessel_schema(required_fields=required_vessel_fields),
            }
            schema_heattransfer["heat_transfer"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "allowed": ["Q_fix", "type"],
                "schema": {
                    "type": {"type": "string", "allowed": ["specified_Q"]},
                    "Q_fix": {"required": False, "type": "number"},
                },
            }

        elif ht_type == "specified_U":
            # Vessel requirements for specified_U (minimal)
            required_vessel_fields = ["length", "diameter"]

            schema_heattransfer = base_schema.copy()
            schema_heattransfer["vessel"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "schema": create_vessel_schema(required_fields=required_vessel_fields),
            }
            schema_heattransfer["heat_transfer"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "allowed": ["U_fix", "type", "temp_ambient"],
                "schema": {
                    "type": {"type": "string", "allowed": ["specified_U"]},
                    "U_fix": {
                        "required": False,
                        "type": "number",
                        "min": 0.001,
                        "max": 1000,
                    },
                    "temp_ambient": {
                        "required": True,
                        "type": "number",
                        "min": 0.001,
                        "max": 2000,
                    },
                },
            }

        elif ht_type == "s-b":
            # Vessel requirements for s-b (Stefan-Boltzmann fire)
            required_vessel_fields = [
                "length",
                "diameter",
                "thickness",
                "heat_capacity",
                "density",
                "orientation",
            ]

            schema_heattransfer = base_schema.copy()
            schema_heattransfer["vessel"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "schema": create_vessel_schema(required_fields=required_vessel_fields),
            }
            schema_heattransfer["heat_transfer"] = {
                "required": True,
                "type": "dict",
                "allow_unknown": False,
                "allowed": ["fire", "type"],
                "schema": {
                    "type": {
                        "required": True,
                        "type": "string",
                        "allowed": ["s-b"],
                    },
                    "fire": {
                        "required": True,
                        "type": "string",
                        "allowed": [
                            "api_pool",
                            "api_jet",
                            "scandpower_pool",
                            "scandpower_jet",
                        ],
                    },
                },
            }
        # Validate with the schema
        v = Validator(schema_heattransfer)
        retval = v.validate(input)

        if not retval:
            # Raise specific exceptions for heat transfer errors
            if "heat_transfer" in v.errors:
                format_cerberus_errors(
                    {"heat_transfer": v.errors["heat_transfer"]}, "heat_transfer"
                )
            elif "vessel" in v.errors:
                format_cerberus_errors({"vessel": v.errors["vessel"]}, "vessel")
            else:
                format_cerberus_errors(v.errors, f"heat_transfer ({ht_type})")

    return retval


def valve_validation(input):
    """
    Validate input['valve'] deeper than cerberus.

    Parameters
    ----------
    input : dict
        Structure holding input

    Return
    ----------
    bool
        True for success, False for failure
    """
    valve_type = input["valve"]["type"]

    # Build schema using factory function
    base_schema = create_top_level_schema()
    base_schema["vessel"] = {"required": True}

    # Heat transfer requirements vary by valve type
    if valve_type == "relief":
        base_schema["heat_transfer"] = {"required": True}
    else:
        base_schema["heat_transfer"] = {"required": False}

    # Create valve schema with type-specific requirements
    base_schema["valve"] = {
        "required": True,
        "type": "dict",
        "allow_unknown": False,
        "schema": create_valve_schema(valve_type),
    }

    # Validate
    v = Validator(base_schema)
    retval = v.validate(input)

    if not retval:
        # Raise specific exceptions for valve errors
        if "valve" in v.errors:
            format_cerberus_errors({"valve": v.errors["valve"]}, "valve")
        else:
            format_cerberus_errors(v.errors, f"valve ({valve_type})")

    return retval


def validation(input):
    """
    Aggregate validation using cerberus and cross-field checks.

    Parameters
    ----------
    input : dict
        Structure holding input

    Return
    ----------
    bool
        True for success, False for failure
    """
    # Schema-based validation
    schema_valid = (
        validate_mandatory_ruleset(input)
        and valve_validation(input)
        and heat_transfer_validation(input)
    )

    if not schema_valid:
        return False

    # Cross-field validation (physical constraints)
    try:
        validate_relief_valve_physics(input)
        validate_composite_vessel(input)
        validate_time_step_sanity(input)
    except (
        ConfigurationError,
        VesselConfigurationError,
        ValveConfigurationError,
    ):
        # Re-raise configuration errors
        raise

    return True
