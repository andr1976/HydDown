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
                "temperature": {"required": True, "type": "number"},
                "pressure": {"required": True, "type": "number"},
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
                "time_step": {"required": True, "type": "number", "min": 0.000001},
                "end_time": {"required": True, "type": "number", "min": 0},
            },
        },
        "vessel": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "length": {"required": True, "type": "number"},
                "diameter": {"required": True, "type": "number"},
                "thickness": {"required": False, "type": "number", "min": 0.0},
                "heat_capacity": {"required": False, "type": "number", "min": 1},
                "thermal_conductivity": {
                    "required": False,
                    "type": "number",
                    "min": 0.0001,
                },
                "density": {"required": False, "type": "number", "min": 1},
                "liner_thickness": {"required": False, "type": "number", "min": 0.0},
                "liner_heat_capacity": {"required": False, "type": "number", "min": 1},
                "liner_thermal_conductivity": {
                    "required": False,
                    "type": "number",
                    "min": 0.0001,
                },
                "liner_density": {"required": False, "type": "number", "min": 1},
                "orientation": {
                    "required": False,
                    "type": "string",
                    "allowed": ["vertical", "horizontal"],
                },
                "type": {
                    "required": False,
                    "type": "string",
                    "allowed": {"Flat-end", "ASME F&D", "DIN", "Hemispherical"},
                },
                "liquid_level": {
                    "required": False,
                    "type": "number",
                    "min": 0,
                },
            },
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
                "diameter": {"type": "number", "min": 0},
                "discharge_coef": {"type": "number", "min": 0},
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
                "U_fix": {"required": False, "type": "number", "min": 0},
                "temp_ambient": {"required": False, "type": "number", "min": 0},
                "h_outer": {"required": False, "type": "number", "min": 0},
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
            "schema": {
                "pressure": {
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
                },
                "temperature": {
                    "type": "dict",
                    "required": False,
                    "allowed": [
                        "wall_high",
                        "wall_low",
                        "wall_mean",
                        "wall_inner",
                        "wall_outer",
                        "gas_high",
                        "gas_low",
                        "gas_mean",
                    ],
                    "schema": {
                        "gas_high": {
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
                        },
                        "gas_low": {
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
                        },
                        "gas_mean": {
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
                        },
                        "wall_high": {
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
                        },
                        "wall_mean": {
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
                        },
                        "wall_low": {
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
                        },
                        "wall_inner": {
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
                        },
                        "wall_outer": {
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
                        },
                    },
                },
            },
        },
    }

    v = Validator(schema_general)
    retval = v.validate(input)
    if v.errors:
        print(v.errors)

    return retval


def heat_transfer_validation(input):
    """
    Validate input['heat_transfer'] deeper than cerberus

    Parameters
    ----------
    input : dict
        Structure holding input


    Return
    ----------
        : bool
        True for success, False for failure
    """
    retval = True

    if input["calculation"]["type"] == "energybalance":
        if input["heat_transfer"]["type"] == "specified_h":
            schema_heattransfer = {
                "initial": {"required": True},
                "calculation": {"required": True},
                "validation": {"required": False},
                "rupture": {"required": False},
                "valve": {"required": True},
                "vessel": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "schema": {
                        "length": {"required": True, "type": "number"},
                        "diameter": {"required": True, "type": "number"},
                        "thickness": {"required": True, "type": "number", "min": 0.0},
                        "heat_capacity": {"required": True, "type": "number", "min": 1},
                        "thermal_conductivity": {
                            "required": False,
                            "type": "number",
                            "min": 0.0001,
                        },
                        "density": {"required": True, "type": "number", "min": 1},
                        "liner_thickness": {
                            "required": False,
                            "type": "number",
                            "min": 0.0,
                        },
                        "liner_heat_capacity": {
                            "required": False,
                            "type": "number",
                            "min": 1,
                        },
                        "liner_thermal_conductivity": {
                            "required": False,
                            "type": "number",
                            "min": 0.0001,
                        },
                        "liner_density": {
                            "required": False,
                            "type": "number",
                            "min": 1,
                        },
                        "orientation": {
                            "required": True,
                            "type": "string",
                            "allowed": ["vertical", "horizontal"],
                        },
                        "type": {
                            "required": False,
                            "type": "string",
                            "allowed": ["Flat-end", "DIN", "ASME F&D", "Hemispherical"],
                        },
                        "liquid_level": {
                            "required": False,
                            "type": "number",
                            "min": 0,
                        },
                    },
                },
                "heat_transfer": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "allowed": [
                        "h_inner",
                        "h_outer",
                        "temp_ambient",
                        "type",
                        "D_throat",
                    ],
                    "schema": {
                        "type": {"type": "string", "allowed": ["specified_h"]},
                        "temp_ambient": {"required": True, "type": "number", "min": 0},
                        "h_outer": {"required": True, "type": "number", "min": 0},
                        "h_inner": {"required": True, "type": ["number", "string"]},
                        "D_throat": {"required": False, "type": "number", "min": 0},
                    },
                },
            }

        elif input["heat_transfer"]["type"] == "specified_Q":
            schema_heattransfer = {
                "initial": {"required": True},
                "calculation": {"required": True},
                "validation": {"required": False},
                "rupture": {"required": False},
                "valve": {"required": True},
                "vessel": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "schema": {
                        "length": {"required": True, "type": "number"},
                        "diameter": {"required": True, "type": "number"},
                        "thickness": {"required": False, "type": "number", "min": 0.0},
                        "heat_capacity": {
                            "required": False,
                            "type": "number",
                            "min": 1,
                        },
                        "density": {"required": False, "type": "number", "min": 1},
                        "orientation": {
                            "required": False,
                            "type": "string",
                            "allowed": ["vertical", "horizontal"],
                        },
                        "type": {
                            "required": False,
                            "type": "string",
                            "allowed": ["Flat-end", "DIN", "ASME F&D", "Hemispherical"],
                        },
                        "liquid_level": {
                            "required": False,
                            "type": "number",
                            "min": 0,
                        },
                    },
                },
                "heat_transfer": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "allowed": ["Q_fix", "type"],
                    "schema": {
                        "type": {"type": "string", "allowed": ["specified_Q"]},
                        "Q_fix": {"required": False, "type": "number"},
                    },
                },
            }

        elif input["heat_transfer"]["type"] == "specified_U":
            schema_heattransfer = {
                "initial": {"required": True},
                "calculation": {"required": True},
                "validation": {"required": False},
                "valve": {"required": True},
                "rupture": {"required": False},
                "vessel": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "schema": {
                        "length": {"required": True, "type": "number"},
                        "diameter": {"required": True, "type": "number"},
                        "thickness": {"required": False, "type": "number", "min": 0.0},
                        "heat_capacity": {
                            "required": False,
                            "type": "number",
                            "min": 1,
                        },
                        "density": {"required": False, "type": "number", "min": 1},
                        "orientation": {
                            "required": False,
                            "type": "string",
                            "allowed": ["vertical", "horizontal"],
                        },
                        "type": {
                            "required": False,
                            "type": "string",
                            "allowed": ["Flat-end", "DIN", "ASME F&D", "Hemispherical"],
                        },
                        "liquid_level": {
                            "required": False,
                            "type": "number",
                            "min": 0,
                        },
                    },
                },
                "heat_transfer": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "allowed": ["U_fix", "type", "temp_ambient"],
                    "schema": {
                        "type": {"type": "string", "allowed": ["specified_U"]},
                        "U_fix": {"required": False, "type": "number", "min": 0.0},
                        "temp_ambient": {"required": True, "type": "number", "min": 0},
                    },
                },
            }

        elif input["heat_transfer"]["type"] == "s-b":
            schema_heattransfer = {
                "initial": {"required": True},
                "calculation": {"required": True},
                "validation": {"required": False},
                "valve": {"required": True},
                "rupture": {"required": False},
                "vessel": {
                    "required": True,
                    "type": "dict",
                    "allow_unknown": False,
                    "schema": {
                        "length": {"required": True, "type": "number"},
                        "diameter": {"required": True, "type": "number"},
                        "thickness": {"required": True, "type": "number", "min": 0.0},
                        "heat_capacity": {"required": True, "type": "number", "min": 1},
                        "density": {"required": True, "type": "number", "min": 1},
                        "orientation": {
                            "required": True,
                            "type": "string",
                            "allowed": ["vertical", "horizontal"],
                        },
                        "type": {
                            "required": False,
                            "type": "string",
                            "allowed": ["Flat-end", "DIN", "ASME F&D", "Hemispherical"],
                        },
                        "liquid_level": {
                            "required": False,
                            "type": "number",
                            "min": 0,
                        },
                    },
                },
                "heat_transfer": {
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
                },
            }
        v = Validator(schema_heattransfer)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)

    return retval


def valve_validation(input):
    """
    Validate input['valve'] deeper than cerberus

    Parameters
    ----------
    input : dict
        Structure holding input


    Return
    ----------
        : bool
        True for success, False for failure
    """
    schema_relief = {
        "initial": {"required": True},
        "calculation": {"required": True},
        "validation": {"required": False},
        "vessel": {"required": True},
        "rupture": {"required": False},
        "heat_transfer": {"required": True},
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
                "diameter": {"required": False, "type": "number", "min": 0},
                "discharge_coef": {"required": False, "type": "number", "min": 0},
                "set_pressure": {"required": True, "type": "number", "min": 0},
                "end_pressure": {"type": "number", "min": 0},
                "blowdown": {"required": False, "type": "number", "min": 0, "max": 1},
                "back_pressure": {"required": True, "type": "number", "min": 0},
                "Cv": {"type": "number", "min": 0},
                "mdot": {"type": ["number", "list"]},
                "time": {"type": "list"},
            },
        },
    }

    schema_psv = {
        "initial": {"required": True},
        "calculation": {"required": True},
        "validation": {"required": False},
        "vessel": {"required": True},
        "rupture": {"required": False},
        "heat_transfer": {"required": False},
        "valve": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": ["orifice", "psv", "controlvalve", "mdot"],
                },
                "flow": {
                    "required": True,
                    "type": "string",
                    "allowed": ["discharge", "filling"],
                },
                "diameter": {"required": True, "type": "number", "min": 0},
                "discharge_coef": {"required": True, "type": "number", "min": 0},
                "set_pressure": {"required": True, "type": "number", "min": 0},
                "end_pressure": {"type": "number", "min": 0},
                "blowdown": {"required": True, "type": "number", "min": 0, "max": 1},
                "back_pressure": {"required": True, "type": "number", "min": 0},
                "Cv": {"type": "number", "min": 0},
                "mdot": {"type": ["number", "list"]},
                "time": {"type": "list"},
            },
        },
    }

    schema_orifice = {
        "initial": {"required": True},
        "calculation": {"required": True},
        "validation": {"required": False},
        "vessel": {"required": True},
        "rupture": {"required": False},
        "heat_transfer": {"required": False},
        "valve": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": ["orifice", "psv", "controlvalve", "mdot"],
                },
                "flow": {
                    "required": True,
                    "type": "string",
                    "allowed": ["discharge", "filling"],
                },
                "diameter": {"required": True, "type": "number", "min": 0},
                "discharge_coef": {"required": True, "type": "number", "min": 0},
                "set_pressure": {"type": "number", "min": 0},
                "end_pressure": {"type": "number", "min": 0},
                "blowdown": {"type": "number", "min": 0, "max": 1},
                "back_pressure": {"required": True, "type": "number", "min": 0},
                "Cv": {"type": "number", "min": 0},
                "mdot": {"type": ["number", "list"]},
                "time": {"type": "list"},
            },
        },
    }

    schema_control_valve = {
        "initial": {"required": True},
        "rupture": {"required": False},
        "calculation": {"required": True},
        "validation": {"required": False},
        "vessel": {"required": True},
        "heat_transfer": {"required": False},
        "valve": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": ["orifice", "psv", "controlvalve", "mdot"],
                },
                "flow": {
                    "required": True,
                    "type": "string",
                    "allowed": ["discharge", "filling"],
                },
                "back_pressure": {"required": True, "type": "number", "min": 0},
                "Cv": {"required": True, "type": "number", "min": 0},
                "time_constant": {"type": "number", "min": 0},
                "characteristic": {
                    "type": "string",
                    "allowed": ["linear", "eq", "fast"],
                },
            },
        },
    }

    schema_mdot = {
        "initial": {"required": True},
        "calculation": {"required": True},
        "validation": {"required": False},
        "vessel": {"required": True},
        "rupture": {"required": False},
        "heat_transfer": {"required": False},
        "valve": {
            "required": True,
            "type": "dict",
            "allow_unknown": False,
            "schema": {
                "type": {
                    "required": True,
                    "type": "string",
                    "allowed": ["orifice", "psv", "controlvalve", "mdot"],
                },
                "flow": {
                    "required": True,
                    "type": "string",
                    "allowed": ["discharge", "filling"],
                },
                "mdot": {"required": True, "type": ["number", "list"]},
                "time": {"type": ["number", "list"]},
                "back_pressure": {"required": True, "type": "number", "min": 0},
            },
        },
    }

    if input["valve"]["type"] == "relief":
        v = Validator(schema_relief)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)

    if input["valve"]["type"] == "psv":
        v = Validator(schema_psv)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)
    elif input["valve"]["type"] == "orifice":
        v = Validator(schema_orifice)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)
    elif input["valve"]["type"] == "controlvalve":
        v = Validator(schema_control_valve)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)
    elif input["valve"]["type"] == "mdot":
        v = Validator(schema_mdot)
        retval = v.validate(input)
        if v.errors:
            print(v.errors)

    return retval


def validation(input):
    """
    Aggregate validation using cerberus

    Parameters
    ----------
    input : dict
        Structure holding input


    Return
    ----------
        : bool
        True for success, False for failure
    """
    return (
        validate_mandatory_ruleset(input)
        and valve_validation(input)
        and heat_transfer_validation(input)
    )
