# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import yaml
import sys

try:
    from hyddown import HydDown
except:
    import sys
    import os
    hyddown_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),"..","src")
    sys.path.append(os.path.abspath(hyddown_path))

    from hyddown import HydDown


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    
    hdown.run(disable_pbar=False)
    hdown.plot()