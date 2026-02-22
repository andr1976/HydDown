# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021 Anders Andreasen
# Published under an MIT license

"""
Fire heat load calculations using the Stefan-Boltzmann approach.

This module provides functions to calculate heat flux from fire scenarios (pool fire, jet fire)
to pressure vessel surfaces using combined radiative and convective heat transfer models.
The models implement industry-standard approaches from API 521 and Scandpower guidelines.

The Stefan-Boltzmann equation accounts for:
- Radiative heat transfer from flame to vessel
- Convective heat transfer from hot gases
- Radiative emission from vessel surface

Available fire scenarios:
- API 521 pool fire: 60 kW/m² incident heat flux
- API 521 jet fire: 100 kW/m² incident heat flux
- Scandpower pool fire: 100 kW/m² incident heat flux
- Scandpower jet fire: 100 kW/m² incident heat flux
- Scandpower peak fires: 150-350 kW/m² incident heat flux

References:
- API Standard 521, Pressure-relieving and Depressuring Systems
- Scandpower Risk Management AS guidelines for fire scenarios
"""


def stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel):
    """
    Calculate heat flux from flame to vessel surface using Stefan-Boltzmann equation.

    Combines radiative and convective heat transfer to determine the net heat flux
    to a vessel surface exposed to fire conditions.

    Parameters
    ----------
    alpha : float
        Absorptivity of the vessel surface (dimensionless, 0-1)
    e_flame : float
        Emissivity of the flame (dimensionless, 0-1)
    e_surface : float
        Emissivity of the vessel surface (dimensionless, 0-1)
    h : float
        Convective heat transfer coefficient [W/(m²·K)]
    Tflame : float
        Temperature of the flame/hot gases [K]
    Tradiative : float
        Effective radiative temperature of the flame [K]
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    The heat flux is calculated as:
    Q = α·ε_flame·σ·T_rad⁴ + h·(T_flame - T_vessel) - ε_surface·σ·T_vessel⁴

    where σ = 5.67×10⁻⁸ W/(m²·K⁴) is the Stefan-Boltzmann constant.

    See also API 521 and Scandpower guidelines for typical parameter values.
    """
    sigma = 5.67e-8  # Stefan-Botlzmann constant W/m^2 K^4
    return (
        alpha * e_flame * sigma * Tradiative**4
        + h * (Tflame - Tvessel)
        - e_surface * sigma * Tvessel**4
    )


def pool_fire_api521(Tvessel):
    """
    Calculate heat flux for API 521 pool fire scenario.

    Implements API 521 standard pool fire with typical incident heat flux of 60 kW/m².
    Uses conservative parameter values for absorptivity, emissivity, and heat transfer
    coefficient as specified in API Standard 521.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.75
    - Flame emissivity (ε_flame): 0.75
    - Surface emissivity (ε_surface): 0.75
    - Convective heat transfer coefficient (h): 20 W/(m²·K)
    - Flame temperature: 600°C (873.15 K)
    - Radiative temperature: 750°C (1023.15 K)

    References
    ----------
    API Standard 521, Pressure-relieving and Depressuring Systems
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
    Calculate heat flux for Scandpower pool fire scenario.

    Implements Scandpower pool fire model with typical incident heat flux of 100 kW/m².
    Uses higher absorptivity and emissivity values compared to API 521 for a more
    severe fire scenario.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.85
    - Flame emissivity (ε_flame): 1.0
    - Surface emissivity (ε_surface): 0.85
    - Convective heat transfer coefficient (h): 30 W/(m²·K)
    - Flame temperature: 804°C (1077.15 K)
    - Radiative temperature: 804°C (1077.15 K)

    References
    ----------
    Scandpower Risk Management AS guidelines for fire scenarios
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
    Calculate heat flux for API 521 jet fire scenario.

    Implements API 521 standard jet fire with typical incident heat flux of 100 kW/m².
    Jet fires have higher flame temperatures and convective heat transfer compared
    to pool fires due to forced convection from jet momentum.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.75
    - Flame emissivity (ε_flame): 0.33
    - Surface emissivity (ε_surface): 0.75
    - Convective heat transfer coefficient (h): 40 W/(m²·K)
    - Flame temperature: 900°C (1173.15 K)
    - Radiative temperature: 1100°C (1373.15 K)

    References
    ----------
    API Standard 521, Pressure-relieving and Depressuring Systems
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
    Calculate heat flux for Scandpower jet fire scenario.

    Implements Scandpower jet fire model with typical incident heat flux of 100 kW/m².
    Features high convective heat transfer coefficient due to forced convection from
    the jet impingement.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.85
    - Flame emissivity (ε_flame): 1.0
    - Surface emissivity (ε_surface): 0.85
    - Convective heat transfer coefficient (h): 100 W/(m²·K)
    - Flame temperature: 635°C (908.15 K)
    - Radiative temperature: 635°C (908.15 K)

    References
    ----------
    Scandpower Risk Management AS guidelines for fire scenarios
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def jet_fire_peak_large_scandpower(Tvessel):
    """
    Calculate heat flux for Scandpower large peak jet fire scenario.

    Implements severe jet fire scenario with very high incident heat flux of 350 kW/m².
    This represents peak impingement conditions from large high-pressure gas releases
    in close proximity to vessel surfaces.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.85
    - Flame emissivity (ε_flame): 1.0
    - Surface emissivity (ε_surface): 0.85
    - Convective heat transfer coefficient (h): 100 W/(m²·K)
    - Flame temperature: 1429.61 K
    - Radiative temperature: 1429.61 K

    References
    ----------
    Scandpower Risk Management AS guidelines for fire scenarios
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 100
    Tflame = 1429.61
    Tradiative = 1429.61
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def jet_fire_peak_small_scandpower(Tvessel):
    """
    Calculate heat flux for Scandpower small peak jet fire scenario.

    Implements severe jet fire scenario with high incident heat flux of 250 kW/m².
    This represents peak impingement conditions from moderate high-pressure gas
    releases or peripheral zones of large jet fires.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.85
    - Flame emissivity (ε_flame): 1.0
    - Surface emissivity (ε_surface): 0.85
    - Convective heat transfer coefficient (h): 100 W/(m²·K)
    - Flame temperature: 1279.29 K
    - Radiative temperature: 1279.29 K

    References
    ----------
    Scandpower Risk Management AS guidelines for fire scenarios
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 100
    Tflame = 1279.29
    Tradiative = 1279.29
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def pool_fire_peak_scandpower(Tvessel):
    """
    Calculate heat flux for Scandpower peak pool fire scenario.

    Implements severe pool fire scenario with high incident heat flux of 150 kW/m².
    This represents very large pool fires or locations close to the fire source
    with maximum radiative and convective heat transfer.

    Parameters
    ----------
    Tvessel : float
        Temperature of the vessel surface [K]

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Notes
    -----
    Fire parameters:
    - Absorptivity (α): 0.85
    - Flame emissivity (ε_flame): 1.0
    - Surface emissivity (ε_surface): 0.85
    - Convective heat transfer coefficient (h): 30 W/(m²·K)
    - Flame temperature: 1212.54 K
    - Radiative temperature: 1212.54 K

    References
    ----------
    Scandpower Risk Management AS guidelines for fire scenarios
    """
    alpha = 0.85
    e_flame = 1
    e_surface = 0.85
    h = 30
    Tflame = 1212.54
    Tradiative = 1212.54
    return stefan_boltzmann(alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel)


def sb_fire(T_vessel, fire_type):
    """
    Dispatcher function for Stefan-Boltzmann fire heat load calculations.

    Selects and calls the appropriate fire scenario function based on the specified
    fire type string. This is the main entry point for fire heat load calculations
    used by the HydDown energy balance solver.

    Parameters
    ----------
    T_vessel : float
        Temperature of the vessel surface [K]
    fire_type : str
        Fire scenario type. Valid options:
        - "api_jet": API 521 jet fire (100 kW/m²)
        - "api_pool": API 521 pool fire (60 kW/m²)
        - "scandpower_pool": Scandpower pool fire (100 kW/m²)
        - "scandpower_jet": Scandpower jet fire (100 kW/m²)
        - "scandpower_jet_peak_large": Scandpower large peak jet fire (350 kW/m²)
        - "scandpower_jet_peak_small": Scandpower small peak jet fire (250 kW/m²)
        - "scandpower_pool_peak": Scandpower peak pool fire (150 kW/m²)

    Returns
    -------
    float
        Net heat flux to vessel surface [W/m²]

    Raises
    ------
    ValueError
        If fire_type is not one of the recognized fire scenario strings

    Examples
    --------
    >>> T_vessel = 400  # K
    >>> Q = sb_fire(T_vessel, "api_pool")
    >>> print(f"Heat flux: {Q:.0f} W/m²")
    """
    if fire_type == "api_jet":
        Q = jet_fire_api521(T_vessel)
    elif fire_type == "api_pool":
        Q = pool_fire_api521(T_vessel)
    elif fire_type == "scandpower_pool":
        Q = pool_fire_scandpower(T_vessel)
    elif fire_type == "scandpower_jet":
        Q = jet_fire_scandpower(T_vessel)
    elif fire_type == "scandpower_jet_peak_large":
        Q = jet_fire_peak_large_scandpower(T_vessel)
    elif fire_type == "scandpower_jet_peak_small":
        Q = jet_fire_peak_small_scandpower(T_vessel)
    elif fire_type == "scandpower_pool_peak":
        Q = pool_fire_peak_scandpower(T_vessel)
    else:
        raise ValueError("Unknown Stefan-Bolzmann fire heat load")
    return Q
