from cerberus import Validator 

def get_example_input(fname):
    import os
    import yaml

    if "C:\\Users\\ANRA" in os.getcwd():
        fname = r"C:\\Users\\ANRA\\Documents\\GitHub\\HydDown\\examples\\" + fname
    else:
        fname = r"//home//runner//work//HydDown//HydDown//examples//" + fname

    with open(fname) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)

    return input


def define_mandatory_ruleset():
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
            'allow_unknown': True,
            'schema': {
                'type': {
                    'type': 'string', 
                    'allowed': ['energybalance',
                                'isenthalpic',
                                'isentropic',
                                'isothermal',
                                'constant_U']
                    },
                'time_step': {'type': 'number', 'min': 0.000001},
                'end_time': {'type': 'number', 'min': 0},
            }
        },
        'vessel': {
            'type': 'dict',
            'allow_unknown': False,  # this overrides the behaviour for
            'schema': {             # the validation of this definition
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
            'allow_unknown': True,
            'schema':{
                'type': {'type': 'string', 'allowed': ['orifice', 'psv', 'controlvalve', 'mdot']},
                'flow': {'type': 'string', 'allowed': ['discharge', 'filling']},
            }
        },
        'heat_transfer':{
            'required': False,
            'type': 'dict',
            'allow_unknown': True,
            'schema':{
                'type': {'type': 'string','allowed': ['specified_Q','specified_h','specified_U']}, 
            }
        }
    }
    return schema

if __name__=="__main__":
    import os 
    schema=define_mandatory_ruleset()
    v = Validator(schema,allow_unknown=True)
    
    for fname in os.listdir("..//examples/"):
        input = get_example_input("..//examples//"+fname)
        print(v.validate(input))
        print(v.errors)