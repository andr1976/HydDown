# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import yaml
import sys
from hyddown import HydDown

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    hdown.run()
    hdown.plot('docs/img/'+ input_filename.split('\\')[1].split('.')[0] +'.png')