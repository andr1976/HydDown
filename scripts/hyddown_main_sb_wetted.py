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

    hyddown_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "src")
    sys.path.append(os.path.abspath(hyddown_path))

    from hyddown import HydDown


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)

    hdown = HydDown(input)

    hdown.run(disable_pbar=False)
    hdown.plot()
    
    hdown.analyze_rupture('/tmp/')
    
    from matplotlib import pyplot as plt
    plt.clf()
    import scienceplots 
    plt.figure(figsize=(4,3))
    plt.style.use(['science','nature'])
    plt.plot(hdown.peak_times, hdown.von_mises, label="von Mises stress")
    plt.plot(hdown.peak_times, hdown.ATS_wetted, label="Allowable Tensile Stength")
    plt.xlabel("Time (s)")
    plt.ylabel("von Mises stress /ATS (MPa)")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig('/home/anra/TML/rupture.png',dpi=600)
