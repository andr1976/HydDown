# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import yaml
import sys
import time

try:
    from hyddown import HydDown
except:
    import sys
    import os
    sys.path.append(os.path.abspath("../src"))
    from hyddown import HydDown


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    start = time.time()
    hdown.run()
    end = time.time()
    print('Elapsed time: ',end-start,' sec.')
    hdown.plot()

    hdown.fluid.build_phase_envelope("None")
    PE=hdown.fluid.get_phase_envelope_data()

    import pylab as plt

    plt.figure(2)
    plt.plot(PE.T, PE.p, '-', label = 'HEOS Phase Envelope', color = 'g')
    plt.plot(hdown.T_fluid,hdown.P,'-.',label = 'P/T fluid trajectory', color = 'b' )
    plt.plot(hdown.T_fluid[0],hdown.P[0],'o',label = 'Start', color = 'b' )
    plt.plot(hdown.T_fluid[-1],hdown.P[-1],'.',label = 'End', color = 'b' )
    plt.xlabel('Temperature [K]')
    plt.ylabel('Pressure [Pa]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()