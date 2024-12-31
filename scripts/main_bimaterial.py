# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

import yaml
import sys
from hyddown import HydDown
import numpy as np

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "mdot_filling.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)

    # mdot = np.loadtxt("mass_flow.csv", delimiter=";")

    hdown = HydDown(input)
    # hdown.p_back = 50e5
    hdown.run(disable_pbar=False)

    hdown.plot()

    import matplotlib.pyplot as plt
    import math

    fig = plt.figure()
    ax = plt.subplot(111)
    times = hdown.time_array[0::100]
    profile = hdown.temp_profile[0::100]
    z = hdown.z
    for i in range(len(times)):
        ax.plot(profile[i], z * 1e3, label=f"t = {round(times[i])} s.")

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.set_xlabel("z, mm")
    ax.set_xlabel("temperature, C")
    ax.set_title("temperature distribution")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    plt.show()
