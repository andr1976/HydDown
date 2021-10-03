# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import yaml
import sys
from hyddown import HydDown
import time
from tqdm import tqdm


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    
    for i in tqdm(hdown.run(generator=True),desc='hyddown',total=len(hdown.time_array)):
        pass
        
    hdown.plot()