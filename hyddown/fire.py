# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license


def stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel):
    """
    Heat flux from flame to vessel surface from combined radiative
    and convective heat transfer. Returns the heat flux in W/m2
    See also API521 and Scandpower guidelines.
    """
    sigma = 5.67e-8  # Stefan-Botlzmann constant W/m^2 K^4
    return (
        alpha * e_flame * sigma * Tradiative ** 4
        + h * (Tflame - Tvessel)
        - e_surface * sigma * Tvessel ** 4
    )


def pool_fire_api521(Tvessel):
    """
    Incident heat flux of 60 kW/m2
    """
    alpha = 0.75
    e_flame = 0.75
    e_surface = 0.75
    h = 20
    Tflame = 600 + 273.15
    Tradiative = 750 + 273.15
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def pool_fire_scandpower(Tvessel):
    """
    Incident heat flux of 100 kW/m2
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 30
    Tflame = 804 + 273.15
    Tradiative = 804 + 273.15
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def jet_fire_api521(Tvessel):
    """
    Incident heat flux of 100 kW/m2
    """
    alpha = 0.75
    e_flame = 0.33
    e_surface = 0.75
    h = 40
    Tflame = 900 + 273.15
    Tradiative = 1100 + 273.15
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def jet_fire_scandpower(Tvessel):
    """
    Incident heat flux of 100 kW/m2
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


if __name__ == "__main__":
    print(jet_fire_scandpower(20 + 273.15))
    alpha = 1
    e_flame = 1
    e_surface = 0
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    print(
        stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, 20 + 273.15)
    )
