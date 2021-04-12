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
            'allow_unknown': True,
            'schema': {
                'temperature':{'type':'number'},
                'pressure':{'type': 'number'},
                'fluid': {'type': 'string'},
                'orientation': {'type': 'string', 'required': 'false'},
            }
        },
        'calculation': {
            'type': 'dict',
            'allow_unknown': True,
            'schema': {
                'type': {'type': 'string'},
                'time_step': {'type': 'number'},
                'end_time': {'type': 'number'},
            }
        },
        'vessel': {
            'type': 'dict',
            'allow_unknown': True,  # this overrides the behaviour for
            'schema': {             # the validation of this definition
                'length': {'type': 'number'},
                'diameter': {'type': 'number'},
            }
        },
        'valve':{
            'type': 'dict',
            'allow_unknown': True,
            'schema':{
                'type': {'type': 'string'},
                'flow': {'type': 'string'},
            }
        }
    }
    return schema

if __name__=="__main__":
    scema=define_mandatory_ruleset()
    v = Validator(scema,allow_unknown=True)
        
    input = get_example_input("..//examples//psv.yml")
    print(v.validate(input))