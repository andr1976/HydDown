from cerberus import Validator
from cerberus.errors import ValidationError 


def validate_mandatory_ruleset(input):
    schema = {
        'initial': {
            'type': 'dict',
            'allow_unknown': False,
            'schema': {
                'temperature':{'type':'number'},
                'pressure':{'type': 'number'},
                'fluid': {'type': 'string'},
                },
            },
        'calculation': {
            'type': 'dict',
            'allow_unknown': False,
            'schema': {
                'type': {
                    'type': 'string', 
                    'allowed': ['energybalance',
                                'isenthalpic',
                                'isentropic',
                                'isothermal',
                                'constant_U'],
                    },
                'time_step': {'type': 'number', 'min': 0.000001},
                'end_time': {'type': 'number', 'min': 0},
                'eta': {'type': 'number', 'min': 0, 'max': 1}
            }
        },
        'vessel': {
            'type': 'dict',
            'allow_unknown': False,  
            'schema': {             
                'length': {'type': 'number'},
                'diameter': {'type': 'number'},
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
            'type': 'dict',
            'allow_unknown': False,
            'schema':{
                'type': {'type': 'string', 'allowed': ['orifice', 'psv', 'controlvalve', 'mdot']},
                'flow': {'type': 'string', 'allowed': ['discharge', 'filling']},
                'diameter': {'type': 'number', 'min': 0},
                'discharge_coef': {'type': 'number', 'min': 0},
                'set_pressure': {'type': 'number', 'min': 0},
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
                        'allowed': ['wall_high','wall_low','gas_high','gas_low','gas_average'],
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
                                    'gas_average':{'required': False, 'type': 'dict', 'contains': ['time','temp'],
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

    v = Validator(schema)
    retval = v.validate(input)
    return retval


def heat_transfer_validation(input):
    if input['calculation']['type'] == 'energybalance':
        if 'heat_transfer' in input:
            if 'specified_h' in input['heat_transfer']['type']:
                if 'h_inner' in input['heat_transfer'] and 'h_outer' in input['heat_transfer'] and 'temp_ambient' in input['heat_transfer']:
                    if 'orientation' in input['vessel'] and 'thickness' in input['vessel'] and 'heat_capacity' in input['vessel'] and 'density' in input['vessel']:
                        return True
                    else:
                        return False
                else:
                    return False
            if 'specified_Q' in input['heat_transfer']['type']:
                if 'Q_fix' in input['heat_transfer']:
                    return True
                else:
                    return False
            if 'specified_U' in input['heat_transfer'] and 'temp_ambient' in input['heat_transfer']:
                if 'Q_fix' in input['heat_transfer']['type']:
                    return True
                else:
                    return False
            if 's-b' in input['heat_transfer']['type']:
                if 'fire' in input['heat_transfer'] and 'orientation' in input['vessel'] and 'thickness' in input['vessel'] and 'heat_capacity' in input['vessel'] and 'density' in input['vessel']:
                    return True
                else:
                    return False
        else:
            return False
    else:
        return True

def valve_validation(input):
    if input['valve']['type'] == 'psv':
        if ('diameter' in input['valve'] 
            and 'discharge_coef' in input['valve'] 
            and 'set_pressure' in input['valve'] 
            and 'blowdown' in input['valve'] 
            and 'back_pressure' in input['valve']):
            return True 
        else:
            return False
    if input['valve']['type'] == 'orifice':
        if ('diameter' in input['valve'] and 'back_pressure' in input['valve'] and 'discharge_coef' in input['valve']):
            return True
        else: 
            return False
    if input['valve']['type'] == 'mdot':
        return True
    if input['valve']['type'] == 'controlvalve':
        if ('Cv' in input['valve'] and 'back_pressure' in input['valve']):
            return True
        else: 
            return False

def validation(input):
    return validate_mandatory_ruleset(input) and valve_validation(input) and heat_transfer_validation(input)