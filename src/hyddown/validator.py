# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

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
        'initial': {
            'required' : True,
            'type': 'dict',
            'allow_unknown': False,
            'schema': {
                'temperature':{'required' : True,'type':'number'},
                'pressure':{'required' : True,'type': 'number'},
                'fluid': {'required' : True,'type': 'string'},
                },
            },
        'calculation': {
            'required' : True,
            'type': 'dict',
            'allow_unknown': False,
            'schema': {
                'type': {
                    'required' : True,
                    'type': 'string', 
                    'allowed': ['energybalance',
                                'isenthalpic',
                                'isentropic',
                                'isothermal',
                                'constant_U'],
                    },
                'time_step': {'required' : True,'type': 'number', 'min': 0.000001},
                'end_time': {'required' : True,'type': 'number', 'min': 0},
            }
        },
        'vessel': {
            'required' : True,
            'type': 'dict',
            'allow_unknown': False,  
            'schema': {             
                'length': {'required' : True,'type': 'number'},
                'diameter': {'required' : True,'type': 'number'},
                'thickness': {'required': False,'type': 'number', 'min': 0.0},
                'heat_capacity': {'required': False, 'type': 'number', 'min': 1},
                'density': {'required': False, 'type': 'number', 'min': 1},
                'orientation': {
                    'required': False, 
                    'type': 'string', 
                    'allowed': ['vertical', 'horizontal']
                }
            }
        },
        'valve':{
            'required' : True,
            'type': 'dict',
            'allow_unknown': False,
            'schema':{
                'type': {'required' : True,'type': 'string', 'allowed': ['orifice', 'psv', 'controlvalve', 'mdot']},
                'flow': {'required' : True,'type': 'string', 'allowed': ['discharge', 'filling']},
                'diameter': {'type': 'number', 'min': 0},
                'discharge_coef': {'type': 'number', 'min': 0},
                'set_pressure': {'type': 'number', 'min': 0},
                'end_pressure': {'type': 'number', 'min': 0},
                'blowdown': {'type': 'number', 'min': 0, 'max': 1},
                'back_pressure': {'type': 'number', 'min': 0},
                'Cv': {'type': 'number', 'min': 0},
                'mdot': {'type': ['number','list']},
                'time' : {'type': 'list'},
            }
        },
        'heat_transfer':{
            'required': False,
            'type': 'dict',
            'allow_unknown': False,
            'allowed': ['Q_fix','h_inner','h_outer','temp_ambient','type','fire','s-b','D_throat'],
            'schema':{
                'type': {'type': 'string','allowed': ['specified_Q','specified_h','specified_U','s-b']}, 
                'Q_fix': {'required': False, 'type': 'number'},
                'U_fix': {'required': False, 'type': 'number', 'min': 0},
                'temp_ambient': {'required': False, 'type': 'number', 'min': 0},
                'h_outer': {'required': False, 'type': 'number', 'min': 0},
                'h_inner': {'required': False, 'type': ['number','string']},
                'fire': {'required': False, 'type': 'string','allowed': ['api_pool','api_jet','scandpower_pool','scandpower_jet']},
                'D_throat' : {'required': False, 'type': 'number', 'min': 0},
            }
        },
        'validation':{
            'required': False,
            'type': 'dict',
            'allow_unknown': False,
            'allowed': ['pressure','temperature'],
            'schema':{
                'pressure': {'type': 'dict',
                            'required': False, 
                            'contains': ['time','pres'], 
                            'schema': {
                                'pres': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                }
                
                    },
                'temperature': {'type': 'dict',
                        'required': False,
                        'allowed': ['wall_high','wall_low','wall_mean','gas_high','gas_low','gas_mean'],
                        'schema': {'gas_high':{'required': False, 'type': 'dict', 'contains': ['time','temp'], 
                                    'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    }, 
                                    'gas_low': {'required': False, 'type':'dict', 'contains': ['time','temp'],
                                    'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    }, 
                                    'gas_mean':{'required': False, 'type': 'dict', 'contains': ['time','temp'],
                                     'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    },
                                    'wall_high':{'required': False, 'type': 'dict', 'contains': ['time','temp'],
                                     'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    },
                                    'wall_mean':{'required': False, 'type': 'dict', 'contains': ['time','temp'],
                                     'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    },
                                    'wall_low':{'required': False, 'type': 'dict', 'contains': ['time','temp'],
                                     'schema': {
                                        'temp': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}, 
                                        'time': {'required': False, 'type': 'list', 'schema': {'type': 'number'}}
                                       }
                                    }
                                }
                } 
            }
        }
    }

    v = Validator(schema_general)
    retval = v.validate(input)
    if v.errors:
        print(v.errors)
    
    if input['calculation']['type'] == 'energybalance':
        if input['heat_transfer']['type'] == 'specified_h':
            schema_heattransfer = {
                'initial' : {'required' : True},
                'calculation' : {'required' : True},
                'validation' : {'required' : False},
                'valve' : {'required' : True},
                'vessel': {
                    'required' : True,
                    'type': 'dict',
                    'allow_unknown': False,  
                    'schema': {             
                        'length': {'required' : True,'type': 'number'},
                        'diameter': {'required' : True,'type': 'number'},
                        'thickness': {'required': True,'type': 'number', 'min': 0.0},
                        'heat_capacity': {'required': True, 'type': 'number', 'min': 1},
                        'density': {'required': True, 'type': 'number', 'min': 1},
                        'orientation': {
                            'required': True, 
                            'type': 'string', 
                            'allowed': ['vertical', 'horizontal']
                        }
                    }
                },
                'heat_transfer':{
                    'required': True,
                    'type': 'dict',
                    'allow_unknown': False,
                    'allowed': ['h_inner','h_outer','temp_ambient','type','D_throat'],
                    'schema':{
                        'type': {'type': 'string','allowed': ['specified_h']}, 
                        'temp_ambient': {'required': True, 'type': 'number', 'min': 0},
                        'h_outer': {'required': True, 'type': 'number', 'min': 0},
                        'h_inner': {'required': True, 'type': ['number','string']},
                        'D_throat' : {'required': False, 'type': 'number', 'min': 0},
                    }
                },
            }  
            v = Validator(schema_heattransfer)
            retval = v.validate(input)
            if v.errors:
                print(v.errors)

        elif input['heat_transfer']['type']=='specified_Q':
            schema_heattransfer = {
                'initial' : {'required' : True},
                'calculation' : {'required' : True},
                'validation' : {'required' : False},
                'valve' : {'required' : True},
                'vessel': {
                    'required' : True,
                    'type': 'dict',
                    'allow_unknown': False,  
                    'schema': {             
                        'length': {'required' : True,'type': 'number'},
                        'diameter': {'required' : True,'type': 'number'},
                        'thickness': {'required': False,'type': 'number', 'min': 0.0},
                        'heat_capacity': {'required': False, 'type': 'number', 'min': 1},
                        'density': {'required': False, 'type': 'number', 'min': 1},
                        'orientation': {
                            'required': False, 
                            'type': 'string', 
                            'allowed': ['vertical', 'horizontal']
                        }
                    }
                },
                'heat_transfer':{
                    'required': True,
                    'type': 'dict',
                    'allow_unknown': False,
                    'allowed': ['Q_fix','type'],
                    'schema':{
                        'type': {'type': 'string','allowed': ['specified_Q']}, 
                        'Q_fix' : {'required': False, 'type': 'number'}, 
                    }
                },
            }
            v = Validator(schema_heattransfer)
            retval = v.validate(input)
            if v.errors:
                print(v.errors)

        elif input['heat_transfer']['type']=='specified_U':
            schema_heattransfer = {
                'initial' : {'required' : True},
                'calculation' : {'required' : True},
                'validation' : {'required' : False},
                'valve' : {'required' : True},
                'vessel': {
                    'required' : True,
                    'type': 'dict',
                    'allow_unknown': False,  
                    'schema': {             
                        'length': {'required' : True,'type': 'number'},
                        'diameter': {'required' : True,'type': 'number'},
                        'thickness': {'required': False,'type': 'number', 'min': 0.0},
                        'heat_capacity': {'required': False, 'type': 'number', 'min': 1},
                        'density': {'required': False, 'type': 'number', 'min': 1},
                        'orientation': {
                            'required': False, 
                            'type': 'string', 
                            'allowed': ['vertical', 'horizontal']
                        }
                    }
                },
                'heat_transfer':{
                    'required': True,
                    'type': 'dict',
                    'allow_unknown': False,
                    'allowed': ['U_fix','type','temp_ambient'],
                    'schema':{
                        'type': {'type': 'string','allowed': ['U_fix']}, 
                        'U_fix': {'required': False, 'type': 'number', 'min' : 0.0},
                        'temp_ambient': {'required': True, 'type': 'number', 'min': 0},
                    }
                },
            }
            v = Validator(schema_heattransfer)
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
    
    if input['calculation']['type'] == 'energybalance':
        if 'heat_transfer' in input:
            if 'specified_h' in input['heat_transfer']['type']:
                if 'h_inner' in input['heat_transfer'] and 'h_outer' in input['heat_transfer'] and 'temp_ambient' in input['heat_transfer']:
                    if 'orientation' in input['vessel'] and 'thickness' in input['vessel'] and 'heat_capacity' in input['vessel'] and 'density' in input['vessel']:
                        pass
                    else:
                        return False
                else:
                    return False
            if 'specified_Q' in input['heat_transfer']['type']:
                if 'Q_fix' in input['heat_transfer']:
                    pass
                else:
                    return False
            if 'specified_U' in input['heat_transfer'] and 'temp_ambient' in input['heat_transfer']:
                if 'Q_fix' in input['heat_transfer']['type']:
                    pass
                else:
                    return False
            if 's-b' in input['heat_transfer']['type']:
                if 'fire' in input['heat_transfer'] and 'orientation' in input['vessel'] and 'thickness' in input['vessel'] and 'heat_capacity' in input['vessel'] and 'density' in input['vessel']:
                    pass
                else:
                    return False
        else:
            return False
    else:
        return True

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
    if input['valve']['type'] == 'psv':
        if ('diameter' in input['valve'] 
            and 'discharge_coef' in input['valve'] 
            and 'set_pressure' in input['valve'] 
            and 'blowdown' in input['valve'] 
            and 'back_pressure' in input['valve']):
            pass 
        else:
            return False
    if input['valve']['type'] == 'orifice':
        if ('diameter' in input['valve'] and 'back_pressure' in input['valve'] and 'discharge_coef' in input['valve']):
            pass
        else: 
            return False
    if input['valve']['type'] == 'mdot':
        pass
    if input['valve']['type'] == 'controlvalve':
        if ('Cv' in input['valve'] and 'back_pressure' in input['valve']):
            pass
        else: 
            return False
    
    return True

def validation(input):
    """
    Aggregate validation using cerberus and homebrew validation 

    Parameters
    ----------
    input : dict 
        Structure holding input 

    
    Return
    ----------
        : bool 
        True for success, False for failure
    """
    return validate_mandatory_ruleset(input) and valve_validation(input) and heat_transfer_validation(input)
