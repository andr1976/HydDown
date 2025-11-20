import math
import numpy as np

# Propeties from Scandpower guideline
T_Cp = (
    np.array([20, 100, 200, 300, 400, 500, 600, 700, 750, 800, 900, 1000, 1100])
    + 273.15
)
SS316_Cp = np.array((472, 487, 503, 512, 520, 530, 541, 551, 555, 559, 565, 571, 577))
# Missing points added manually
Duplex_Cp = np.array([480, 500, 530, 560, 600, 635, 670, 710, 730, 750, 790, 840, 840])
SMo_Cp = np.array([500, 520, 540, 555, 570, 580, 590, 600, 610, 610, 610, 610, 610])
CS_LT_Cp = np.array([450, 480, 510, 550, 600, 660, 750, 900, 1450, 820, 540, 540, 540])

# Material UTS data from Scandpower guideline
T = (
    np.array(
        [
            20,
            50,
            100,
            150,
            200,
            250,
            300,
            350,
            400,
            450,
            500,
            550,
            600,
            650,
            700,
            750,
            800,
            900,
            1000,
            1100,
        ]
    )
    + 273.17
)
Duplex_UTS = (
    np.array(
        [
            730,
            701,
            668,
            650,
            640,
            630,
            621,
            606,
            591,
            540,
            482,
            423,
            358,
            299,
            234,
            164,
            124,
            88,
            69,
            58,
        ]
    )
    * 1e6
)
SS_UTS = (
    np.array(
        [
            575,
            549,
            523,
            503,
            489,
            477,
            472,
            466,
            463,
            460,
            449,
            431,
            397,
            357,
            299,
            242,
            184,
            98,
            70,
            60,
        ]
    )
    * 1e6
)
SMo_UTS = (
    np.array(
        [
            730,
            710,
            680,
            660,
            645,
            633,
            625,
            618,
            610,
            595,
            580,
            550,
            510,
            445,
            380,
            305,
            230,
            130,
            90,
            65,
        ]
    )
    * 1e6
)

CS_235LT_UTS = (
    np.array(
        [
            420,
            414,
            407,
            403,
            397,
            391,
            382,
            378,
            370,
            353,
            308,
            252,
            189,
            139,
            92,
            59,
            46,
            38,
            30,  # Manual hand-interpolation
            22,  # Manual hand-interpolation
        ]
    )
    * 1e6
)

CS_360LT_UTS = (
    np.array(
        [
            545,
            537,
            529,
            523,
            515,
            507,
            496,
            491,
            480,
            458,
            400,
            327,
            245,
            180,
            120,
            76,
            60,
            49,
            38,  # Manual hand-interpolation
            27,  # Manual hand-interpolation
        ]
    )
    * 1e6
)


def von_mises(p, d, wt, sigma_a=30e6):
    """
    von Mises stress calculated according to:
    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.

    As also applied in:
    Andreasen, A.; Borroni, F.; Zan Nieto, M.; Stegelmann, C.; P. Nielsen, R.
    On the Adequacy of API 521 Relief-Valve Sizing Method for Gas-Filled Pressure Vessels
    Exposed to Fire. Safety 2018, 4, 11. https://doi.org/10.3390/safety4010011

    Parameters
    ----------
    p : float
        Pressure (Pa)
    d : float
        Inner diameter (m)
    D : float
        Outer diameter (m)
    wt: float
        Wall thickness (m)
    sigma_a: float
        Default

    Returns
    ----------
    sigma_e : float
        von Mises stress (Pa)

    """

    D = d + 2 * wt

    sigma_e = math.sqrt(3 * ((p * D**2) / (D**2 - d**2)) ** 2 + sigma_a)
    return sigma_e


def ATS(temperature, material, k_s=0.85, k_y=1):
    """
    Calculation of Allowable Tensile Strength according to:
    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.

    Parameters
    ----------
    temperature : float
        Temperature (K)
    material : string
        Material type: 235LT, 360LT (ASTM A-333/A-671), 2205 (SA-790/ASTM A-790),
        316 (ASTM A-320, ASME A-358), 6Mo (ASTM B-677)
    k_s : float
        General safety factor. For typical materials 0.85 is used. If "guaranteed" minimum
        values a factor 1.0 can be used.
    k_y : float
        Additional factor used for materials with missing or uncertain material data.
        Normally 1.0.

    Returns
    ----------
    ATS : float
        Allowable Tensile Strength (Pa)

    """

    return UTS(temperature, material=material) * k_s * k_y


def UTS(temperature, material):
    """
    Tabulation look-up / interpolation to retrieve the Ultimate Tensile Strength
    as a function of temperature for various typical materials according to:

    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.


    Parameters
    ----------
    temperature : float
        Temperature (K)
    material : string
        Material type: 235LT, 360LT (ASTM A-333/A-671), Duplex (2205, SA-790/ASTM A-790),
        316 (ASTM A-320, ASME A-358), 6Mo (ASTM B-677)

    Return
    ----------
    UTS : float
        Ultimate Tensile Strength (Pa)

    """

    if material == "CS_235LT":
        return np.interp(temperature, T, CS_235LT_UTS)
    elif material == "CS_360LT":
        return np.interp(temperature, T, CS_360LT_UTS)
    elif material == "SS316":
        return np.interp(temperature, T, SS_UTS)
    elif material == "Duplex":
        return np.interp(temperature, T, Duplex_UTS)
    elif material == "6Mo":
        return np.interp(temperature, T, SMo_UTS)


def steel_Cp(temperature, material):
    """
    Tabulation look-up / interpolation to retrieve the heat capacity
    as a function of temperature for various typical materials according to:

    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.


    Parameters
    ----------
    temperature : float
        Temperature (K)
    material : string
        Material type: 235LT, 360LT (ASTM A-333/A-671), Duplex (2205, SA-790/ASTM A-790),
        316 (ASTM A-320, ASME A-358), 6Mo (ASTM B-677)

    Return
    ----------
    Cp : float
        Ultimate Tensile Strength (Pa)

    """
    if material == "CS_235LT":
        return np.interp(temperature, T_Cp, CS_LT_Cp)
    elif material == "CS_360LT":
        return np.interp(temperature, T_Cp, CS_LT_Cp)
    elif material == "SS316":
        return np.interp(temperature, T_Cp, SS316_Cp)
    elif material == "Duplex":
        return np.interp(temperature, T_Cp, Duplex_Cp)
    elif material == "6Mo":
        return np.interp(temperature, T_Cp, SMo_Cp)


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    plt.figure(1)
    plt.plot(T_Cp, Duplex_Cp, "--", label="22Cr Duplex")
    plt.plot(T_Cp, SS316_Cp, "-.", label="SS316")
    plt.plot(T_Cp, SMo_Cp, "-", label="6Mo")
    plt.plot(T_Cp, CS_LT_Cp, "k--", label="CS")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Steel heat capacity (J/kg K)")
    plt.savefig("heat_capacity.png", dpi=600)
    plt.legend(loc="best")

    plt.figure(2)
    plt.plot(T, Duplex_UTS, "--", label="22Cr Duplex")
    plt.plot(T, SS_UTS, "-.", label="SS316")
    plt.plot(T, SMo_UTS, "-", label="6Mo")
    plt.plot(T, CS_235LT_UTS, "k--", label="235LT")
    plt.plot(T, CS_360LT_UTS, "k-.", label="360LT")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Ultimate Tensile Strength (Pa)")
    plt.legend(loc="best")
    plt.savefig("UTS.png", dpi=600)
    plt.show()
