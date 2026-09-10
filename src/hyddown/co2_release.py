# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
thermopack-based CO2 release model (HEM discharge rate + dry-ice atmospheric state).

Why this module exists
-----------------------
When liquefied/pressurised CO2 is released and expands towards atmospheric pressure it
crosses the **triple point** (P_triple = 5.18 bar, T_triple = -56.6 degC). Below the triple
point liquid CO2 cannot exist and the stream is a mixture of **vapour + solid CO2
(dry ice / "snow")** on the sublimation line (frost point ~ -79 degC / 194.14 K).

CoolProp - the main HydDown backend - cannot evaluate CO2 below the triple point, so it can
neither integrate the HEM isentrope through a sub-triple-point throat nor compute the
atmospheric end state. This module uses **thermopack (tcPR + solid CO2 EoS)** instead, whose
solid-aware ``two_phase_psflash`` returns vapour+solid states (phase index 7, solid fraction
in the ``betaL`` slot). All thermodynamics stay inside thermopack, so there is a single
consistent reference basis - no CoolProp/thermopack basis bridge (see DRY_ICE_HANDOVER.md
sec. 5.1/5.2 for the +415.8 kJ/mol offset trap this avoids).

Scope: tailored for **pure CO2**. Two release cases (DRY_ICE_HANDOVER.md sec. 2):
  * hole in the liquid space  -> saturated-liquid stagnation
  * hole in the vapour space  -> saturated-vapour stagnation

The class is deliberately self-contained: HydDown constructs it once and calls it per
timestep with the current tank pressure; no thermopack import leaks into the HydDown class.

Public API (:class:`CO2ReleaseModel`)
    stagnation(P, phase)            -> (h0[J/kg], s0[J/mol/K], rho0[kg/m3], T0[K])
    hem_rate(P0, phase, Cd, area)   -> dict with mdot, throat diagnostics, stagnation state
    atm_split(h0_mass)              -> {T, vapour_frac, solid_frac} isenthalpic flash to 1 atm
    atm_split_isentropic(s0_molar)  -> {T, vapour_frac, solid_frac} isentropic (blast bound)
    release_state(P0, phase, Cd, area) -> (rate_dict, atm_dict)  convenience wrapper
"""

import math

import numpy as np
from scipy import optimize

from thermopack.tcPR import tcPR


# --- CO2 constants -----------------------------------------------------------------
M_CO2 = 0.0440095          # kg/mol (molar mass of CO2)
P_TRIPLE = 5.18e5          # Pa    (CO2 triple-point pressure)
T_TRIPLE = 216.592         # K     (CO2 triple-point temperature)
P_CRIT = 7.3773e6          # Pa    (CO2 critical pressure)
P_ATM_DEFAULT = 101325.0   # Pa

# thermopack two_phase_psflash phase index for a vapour + solid equilibrium.
# In that regime the solid mole fraction is returned in the ``betaL`` slot
# (DRY_ICE_HANDOVER.md sec. 5.4).
PHASE_VAP_SOLID = 7


def _scal(x):
    """Coerce a thermopack return (scalar, 1-tuple or length-1 array) to a float."""
    try:
        return float(x[0])
    except (TypeError, IndexError, KeyError):
        return float(x)


class CO2ReleaseModel:
    """Homogeneous-equilibrium CO2 release rate and dry-ice atmospheric end state.

    Parameters
    ----------
    back_pressure : float
        Downstream/back pressure the HEM throat expands against [Pa]. For an
        atmospheric release this equals ``atm_pressure``.
    atm_pressure : float
        Pressure at which the atmospheric (dry-ice) state is evaluated [Pa].
    eos : str
        Equation of state. Only ``"tcPR"`` is supported (translated-consistent
        Peng-Robinson) because the solid CO2 model is calibrated against it.
    """

    def __init__(self, back_pressure=P_ATM_DEFAULT, atm_pressure=P_ATM_DEFAULT, eos="tcPR"):
        if str(eos).lower().replace("-", "") != "tcpr":
            raise ValueError(
                f"Unsupported eos '{eos}'. Only 'tcPR' is supported for CO2 dry-ice modelling."
            )
        self.eos = tcPR("CO2")
        self.eos.init_solid("CO2")
        self.z = np.array([1.0])
        self.LIQ = self.eos.LIQPH
        self.VAP = self.eos.VAPPH
        self.M = self.eos.compmoleweight(1) / 1000.0  # g/mol -> kg/mol
        self.p_back = back_pressure
        self.p_atm = atm_pressure
        self.P_TRIPLE = P_TRIPLE  # convenience: the solid-in-vessel validity floor
        self._init_atm_endpoints()

    # ------------------------------------------------------------------ setup
    def _init_atm_endpoints(self):
        """Pre-compute the (time-invariant) 1-atm endpoint constants.

        The frost point and the pure-phase gas/solid enthalpies at the atmospheric
        frost point depend only on the back pressure, so they are computed once and
        reused by the isenthalpic lever at every timestep (DRY_ICE_HANDOVER.md sec. 4).
        The frost point is obtained directly from a psflash whose entropy brackets the
        gas/solid values - this lands on the gas<->solid Gibbs-equality frost point
        (194.14 K), consistent with the lever it feeds (sec. 5.3).
        """
        s_gas = _scal(self.eos.entropy(195.0, self.p_atm, self.z, self.VAP))
        s_sol = _scal(self.eos.solid_entropy(195.0, self.p_atm, self.z))
        fr = self.eos.two_phase_psflash(self.p_atm, self.z, 0.5 * (s_gas + s_sol))
        self.T_frost = fr.T
        self.h_gas_atm = _scal(self.eos.enthalpy(self.T_frost, self.p_atm, self.z, self.VAP)) / self.M
        self.h_solid_atm = _scal(self.eos.solid_enthalpy(self.T_frost, self.p_atm, self.z)) / self.M

    # ------------------------------------------------------------- stagnation
    def stagnation(self, P, phase):
        """Saturated stagnation state feeding the hole at tank pressure ``P``.

        Parameters
        ----------
        P : float
            Tank (stagnation) pressure [Pa]. Must be below the CO2 critical pressure
            (the vessel is two-phase, so this always holds within scope).
        phase : str
            ``"liquid"`` -> saturated liquid (hole in the liquid space);
            anything else -> saturated vapour (hole in the vapour space).

        Returns
        -------
        h0 : float
            Mass-specific stagnation enthalpy [J/kg].
        s0 : float
            **Molar** stagnation entropy [J/mol/K] (the unit thermopack's flashes expect).
        rho0 : float
            Stagnation density [kg/m3].
        T0 : float
            Saturation temperature at ``P`` [K].
        """
        if P >= P_CRIT:
            raise ValueError(
                f"Stagnation pressure {P/1e5:.2f} bar >= CO2 critical pressure "
                f"{P_CRIT/1e5:.2f} bar - saturated stagnation is undefined."
            )
        Tsat = _scal(self.eos.bubble_temperature(P, self.z))
        ph = self.LIQ if phase == "liquid" else self.VAP
        h0 = _scal(self.eos.enthalpy(Tsat, P, self.z, ph)) / self.M
        s0 = _scal(self.eos.entropy(Tsat, P, self.z, ph))  # molar - fed straight to psflash
        v0 = _scal(self.eos.specific_volume(Tsat, P, self.z, ph))
        return h0, s0, self.M / v0, Tsat

    # ------------------------------------------------ isentrope mixture props
    def _iso_props(self, P, s0_molar):
        """Mass enthalpy, density, T and solid fraction on the isentrope at pressure ``P``.

        Uses the solid-aware psflash so that below the triple point the throat state is
        a vapour + dry-ice mixture and its density/enthalpy carry the solid fraction
        (the "dry-ice lever at the throat").
        """
        fr = self.eos.two_phase_psflash(P, self.z, s0_molar)
        T, bv, bl = fr.T, fr.betaV, fr.betaL
        if fr.phase == PHASE_VAP_SOLID:
            # vapour + solid: ``bl`` is the solid mole fraction
            hv = _scal(self.eos.enthalpy(T, P, self.z, self.VAP))
            vv = _scal(self.eos.specific_volume(T, P, self.z, self.VAP))
            hs = _scal(self.eos.solid_enthalpy(T, P, self.z))
            vs = _scal(self.eos.solid_volume(T, P, self.z))
            h = bv * hv + bl * hs
            v = bv * vv + bl * vs
            solid = bl
        elif 0.0 < bv < 1.0:
            # vapour + liquid (above the triple point)
            hv = _scal(self.eos.enthalpy(T, P, self.z, self.VAP))
            vv = _scal(self.eos.specific_volume(T, P, self.z, self.VAP))
            hl = _scal(self.eos.enthalpy(T, P, self.z, self.LIQ))
            vl = _scal(self.eos.specific_volume(T, P, self.z, self.LIQ))
            h = bv * hv + bl * hl
            v = bv * vv + bl * vl
            solid = 0.0
        else:
            # single phase
            ph = self.VAP if bv >= 0.5 else self.LIQ
            h = _scal(self.eos.enthalpy(T, P, self.z, ph))
            v = _scal(self.eos.specific_volume(T, P, self.z, ph))
            solid = 0.0
        return h / self.M, self.M / v, T, solid

    # --------------------------------------------------------------- HEM rate
    def hem_rate(self, P0, phase, Cd, area):
        """Homogeneous-equilibrium mass flow through the hole.

        Integrates the isentrope from the stagnation state down to the choked throat,
        maximising the mass flux ``G = rho * sqrt(2 * (h0 - h))``. The maximiser is the
        choked throat; if the flux keeps rising to the back pressure the flow is
        unchoked and the throat sits at ``back_pressure``.

        Parameters
        ----------
        P0 : float
            Stagnation/tank pressure [Pa].
        phase : str
            ``"liquid"`` or ``"gas"`` (see :meth:`stagnation`).
        Cd : float
            Discharge coefficient [-].
        area : float
            Hole area [m2].

        Returns
        -------
        dict
            ``mdot`` [kg/s], ``G`` [kg/m2/s], ``P_throat`` [Pa], ``T_throat`` [K],
            ``solid_frac_throat`` [-], ``choked`` [bool] and the stagnation state
            (``h0`` [J/kg], ``s0`` [J/mol/K], ``rho0`` [kg/m3], ``T0`` [K]).
        """
        h0, s0, rho0, T0 = self.stagnation(P0, phase)

        stag = {"h0": h0, "s0": s0, "rho0": rho0, "T0": T0}
        if P0 <= self.p_back:
            return {"mdot": 0.0, "G": 0.0, "P_throat": P0, "T_throat": T0,
                    "solid_frac_throat": 0.0, "choked": False, **stag}

        def neg_flux(P):
            hk, rho, _T, _solid = self._iso_props(P, s0)
            dh = h0 - hk
            return -rho * math.sqrt(2.0 * dh) if dh > 0.0 else 0.0

        # The flux curve G(P) can carry TWO local maxima across the triple-point kink:
        # one on the vapour+liquid branch (throat above 5.18 bar) and one on the
        # vapour+solid branch (throat below it). A bare bounded search can lock onto the
        # wrong one, so first locate the global maximum on a coarse grid, then refine
        # locally. The back-pressure boundary is included (unchoked case).
        lo, hi = self.p_back, P0
        n = 40
        p_grid = [lo + (hi - lo) * k / n for k in range(n + 1)]
        g_grid = [-neg_flux(P) for P in p_grid]
        kbest = max(range(len(p_grid)), key=lambda k: g_grid[k])
        Pth, Gmax = p_grid[kbest], g_grid[kbest]

        a = p_grid[max(kbest - 1, 0)]
        b = p_grid[min(kbest + 1, n)]
        if b > a:
            res = optimize.minimize_scalar(
                neg_flux, bounds=(a, b), method="bounded",
                options={"xatol": max((b - a) * 1e-3, 1.0)},
            )
            if -neg_flux(res.x) >= Gmax:
                Pth, Gmax = res.x, -neg_flux(res.x)
        Pth = float(Pth)
        choked = bool(Pth > lo * 1.001)

        _hk, _rho, Tth, solid = self._iso_props(Pth, s0)
        return {"mdot": Cd * area * Gmax, "G": Gmax, "P_throat": Pth, "T_throat": Tth,
                "solid_frac_throat": solid, "choked": choked, **stag}

    # ------------------------------------------------------- atmospheric state
    def atm_split(self, h0_mass):
        """Isenthalpic flash of the stagnation enthalpy to atmospheric pressure.

        The atmospheric state conserves the stagnation enthalpy ``h0`` (a free jet doing
        no external work). thermopack's ``two_phase_phflash`` is *not* solid-aware, so
        the vapour/solid split is computed by an explicit lever against the pre-computed
        1-atm endpoint enthalpies (DRY_ICE_HANDOVER.md sec. 2/5.4):

            beta_gas   = (h0 - h_solid) / (h_gas - h_solid)
            solid_frac = 1 - clip(beta_gas, 0, 1)

        For a pure component the mass and mole fractions coincide.

        Parameters
        ----------
        h0_mass : float
            Mass-specific stagnation enthalpy [J/kg].

        Returns
        -------
        dict
            ``T`` [K] (frost point when solid forms, else the superheated temperature),
            ``vapour_frac`` [-] and ``solid_frac`` [-] (mass fractions).
        """
        beta_gas = (h0_mass - self.h_solid_atm) / (self.h_gas_atm - self.h_solid_atm)
        if beta_gas >= 1.0:
            # Superheated: no dry ice - a plain (solid-free) phflash gives the temperature.
            fr = self.eos.two_phase_phflash(self.p_atm, self.z, h0_mass * self.M)
            return {"T": fr.T, "vapour_frac": 1.0, "solid_frac": 0.0}
        beta_gas = max(beta_gas, 0.0)
        return {"T": self.T_frost, "vapour_frac": beta_gas, "solid_frac": 1.0 - beta_gas}

    def atm_split_isentropic(self, s0_molar):
        """Isentropic flash to atmospheric pressure (available-work / blast-energy bound).

        Uses the solid-aware psflash directly (it *is* solid-aware for entropy), so the
        vapour/solid split comes straight from the flash result.
        """
        fr = self.eos.two_phase_psflash(self.p_atm, self.z, s0_molar)
        if fr.phase == PHASE_VAP_SOLID:
            return {"T": fr.T, "vapour_frac": fr.betaV, "solid_frac": fr.betaL}
        return {"T": fr.T, "vapour_frac": 1.0, "solid_frac": 0.0}

    # ------------------------------------------------------------- convenience
    def release_state(self, P0, phase, Cd, area):
        """Return both the HEM rate dict and the (isenthalpic) atmospheric split."""
        rate = self.hem_rate(P0, phase, Cd, area)
        atm = self.atm_split(rate["h0"])
        return rate, atm


# --------------------------------------------------------------------------- self-test
if __name__ == "__main__":
    # Reproduce the DRY_ICE_HANDOVER.md table (live-tcPR column) as a sanity check.
    model = CO2ReleaseModel()
    print(f"M = {model.M:.5f} kg/mol   frost point = {model.T_frost:.3f} K")
    print(f"{'release':>22} | {'P0[bar]':>8} | {'isenth':>7} | {'isentr':>7}")
    for label, TC, ph in (
        ("sat-liquid +17C", 17.0, "liquid"),
        ("sat-liquid -30C", -30.0, "liquid"),
        ("sat-vapour +17C", 17.0, "gas"),
        ("sat-vapour -30C", -30.0, "gas"),
    ):
        T = 273.15 + TC
        P0 = _scal(model.eos.bubble_pressure(T, model.z))
        h0, s0, _rho0, _T0 = model.stagnation(P0, ph)
        se = model.atm_split(h0)["solid_frac"]
        ss = model.atm_split_isentropic(s0)["solid_frac"]
        print(f"{label:>22} | {P0/1e5:8.2f} | {se:7.3f} | {ss:7.3f}")

    rate, atm = model.release_state(P0=53.39e5, phase="liquid", Cd=0.62, area=math.pi * 0.05 ** 2 / 4)
    print(f"\nHEM sat-liquid +17C, 50 mm hole: mdot = {rate['mdot']:.2f} kg/s, "
          f"throat {rate['P_throat']/1e5:.2f} bar / {rate['T_throat']:.1f} K, "
          f"solid_frac_throat = {rate['solid_frac_throat']:.3f}, choked = {rate['choked']}")
    print(f"atmospheric: T = {atm['T']:.1f} K, vapour = {atm['vapour_frac']:.3f}, "
          f"dry ice = {atm['solid_frac']:.3f}")
