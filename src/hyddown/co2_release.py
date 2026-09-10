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
        self.P_TRIPLE = P_TRIPLE  # convenience: the (literature) validity-floor constant
        self._init_atm_endpoints()
        self._init_triple_point()

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
        return self.hem_rate_from_stagnation(h0, s0, rho0, T0, P0, Cd, area)

    def hem_rate_from_stagnation(self, h0, s0, rho0, T0, P0, Cd, area):
        """HEM mass flow from an explicit stagnation state (see :meth:`hem_rate`).

        Used when the stagnation is not a simple saturated liquid/vapour at ``P0`` -
        e.g. a vapour leak on the sublimation line below the triple point, where
        ``bubble_temperature`` is undefined. ``h0`` [J/kg], ``s0`` [J/mol/K],
        ``rho0`` [kg/m3], ``T0`` [K], ``P0`` [Pa] describe the upstream state.
        """
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

    def gas_leak_rate(self, T, P, Cd, area):
        """HEM rate for a vapour leak from an explicit vessel state (T, P).

        Used in the solid-in-vessel regime (below the triple point), where the leaking
        gas sits on the sublimation line and ``stagnation()`` (which needs a bubble
        point) does not apply. Returns the usual hem_rate dict.
        """
        h0 = _scal(self.eos.enthalpy(T, P, self.z, self.VAP)) / self.M
        s0 = _scal(self.eos.entropy(T, P, self.z, self.VAP))  # molar
        v0 = _scal(self.eos.specific_volume(T, P, self.z, self.VAP))
        return self.hem_rate_from_stagnation(h0, s0, self.M / v0, T, P, Cd, area)

    # =====================================================================
    # Solid-in-vessel regime (below the triple point) - opt-in fallback
    # =====================================================================
    # Once the tank itself reaches the triple point, Gibbs' phase rule pins the
    # state: with three phases (solid+liquid+gas) the invariant point fixes T and
    # P, and the leak/heat only shift the phase amounts (regime B); with two
    # phases (solid+gas) the state rides the sublimation line (regime C). Both are
    # solved from the mass / volume / internal-energy balances against the
    # pure-phase properties - no flash iteration through the degenerate triple
    # point. All properties are in thermopack's own basis (no CoolProp mixing).

    def _g_fluid(self, T, P, ph):
        """Mass-specific Gibbs energy of a fluid phase [J/kg]."""
        h = _scal(self.eos.enthalpy(T, P, self.z, ph)) / self.M
        s = _scal(self.eos.entropy(T, P, self.z, ph)) / self.M
        return h - T * s

    def _g_solid(self, T, P):
        """Mass-specific Gibbs energy of solid CO2 [J/kg]."""
        h = _scal(self.eos.solid_enthalpy(T, P, self.z)) / self.M
        s = _scal(self.eos.solid_entropy(T, P, self.z)) / self.M
        return h - T * s

    def _init_triple_point(self):
        """Locate the tcPR self-consistent triple point and the three pure-phase
        (v, u, h) vertices used by the invariant lever.

        The triple point is where g_solid = g_liquid along the fluid saturation
        line (on which g_liquid = g_vapour already, for a pure component).
        """
        def gap(T):
            P = _scal(self.eos.bubble_pressure(T, self.z))
            return self._g_solid(T, P) - self._g_fluid(T, P, self.LIQ)

        self.T_TRIPLE_EOS = optimize.brentq(gap, 208.0, 220.0)
        self.P_TRIPLE_EOS = _scal(self.eos.bubble_pressure(self.T_TRIPLE_EOS, self.z))
        T, P = self.T_TRIPLE_EOS, self.P_TRIPLE_EOS
        # pure-phase mass-specific volume, enthalpy, internal energy at the triple point
        self.v_s = _scal(self.eos.solid_volume(T, P, self.z)) / self.M
        self.h_s = _scal(self.eos.solid_enthalpy(T, P, self.z)) / self.M
        self.v_l = _scal(self.eos.specific_volume(T, P, self.z, self.LIQ)) / self.M
        self.h_l = _scal(self.eos.enthalpy(T, P, self.z, self.LIQ)) / self.M
        self.v_g = _scal(self.eos.specific_volume(T, P, self.z, self.VAP)) / self.M
        self.h_g = _scal(self.eos.enthalpy(T, P, self.z, self.VAP)) / self.M
        self.u_s = self.h_s - P * self.v_s
        self.u_l = self.h_l - P * self.v_l
        self.u_g = self.h_g - P * self.v_g
        # Coldest physical vessel state on the sublimation line: the leak stops once the
        # vessel reaches the back pressure, so it cannot cool below the sublimation
        # temperature at p_back.
        self.T_subl_floor = optimize.brentq(
            lambda T: self._sublimation_pressure(T) - self.p_back,
            185.0, self.T_TRIPLE_EOS - 1e-4,
        )

    def triple_point_U(self, m_liquid, m_gas, m_solid=0.0):
        """Total internal energy [J] in thermopack basis for a phase split *at the
        triple point*. Used to re-base the vessel energy at the CoolProp->thermopack
        hand-off (masses are basis-independent; energy is not)."""
        return m_solid * self.u_s + m_liquid * self.u_l + m_gas * self.u_g

    def triple_LG_from_MV(self, M, V):
        """Liquid+vapour split at the triple point consistent with total mass ``M`` [kg]
        and vessel volume ``V`` [m3] (no solid yet). Returns (m_l, m_g, U).

        This is the physically-consistent way to place the vessel *on* the triple point at
        the CoolProp->thermopack hand-off (and to clamp back to it if it warms up): the
        volume constraint fixes the liquid/vapour split, and the internal energy follows.
        """
        m_g = (V - M * self.v_l) / (self.v_g - self.v_l)
        m_g = min(max(m_g, 0.0), M)
        m_l = M - m_g
        return m_l, m_g, m_l * self.u_l + m_g * self.u_g

    def triple_lever(self, M, U, V):
        """Invariant three-phase lever at the triple point (regime B).

        Solves the linear system (volume / internal-energy / mass) for the three
        phase masses at fixed T, P:

            [v_s v_l v_g] [m_s]   [V]
            [u_s u_l u_g] [m_l] = [U]
            [ 1   1   1 ] [m_g]   [M]

        Returns m_s, m_l, m_g [kg] (may be negative if the state is outside the
        three-phase triangle), plus the invariant T, P.
        """
        A = np.array([[self.v_s, self.v_l, self.v_g],
                      [self.u_s, self.u_l, self.u_g],
                      [1.0, 1.0, 1.0]])
        m = np.linalg.solve(A, np.array([V, U, M]))
        return {"m_s": m[0], "m_l": m[1], "m_g": m[2],
                "T": self.T_TRIPLE_EOS, "P": self.P_TRIPLE_EOS}

    def _sublimation_pressure(self, T):
        """Solid+gas equilibrium (sublimation) pressure at temperature T [Pa]."""
        return optimize.brentq(
            lambda P: self._g_fluid(T, P, self.VAP) - self._g_solid(T, P),
            1.0, self.P_TRIPLE_EOS * 1.0001,
        )

    def _sublimation_props(self, T):
        """(P, v_g, u_g, v_s, u_s) on the sublimation line at temperature T (mass basis)."""
        P = self._sublimation_pressure(T)
        vg = _scal(self.eos.specific_volume(T, P, self.z, self.VAP)) / self.M
        hg = _scal(self.eos.enthalpy(T, P, self.z, self.VAP)) / self.M
        vs = _scal(self.eos.solid_volume(T, P, self.z)) / self.M
        hs = _scal(self.eos.solid_enthalpy(T, P, self.z)) / self.M
        return P, vg, hg - P * vg, vs, hs - P * vs

    def sublimation_state(self, M, U, V, T_floor=None):
        """Solid+gas state on the sublimation line below the triple point (regime C).

        With liquid gone the system has one degree of freedom, so T (hence P and the
        pure-phase properties) is found by matching the internal energy; the solid/gas
        split follows from the mass and volume balances. T is clamped to the physical
        window [sublimation temperature at p_back, triple point].
        """
        lo = T_floor if T_floor is not None else self.T_subl_floor
        hi = self.T_TRIPLE_EOS - 1e-4

        def split(T):
            P, vg, ug, vs, us = self._sublimation_props(T)
            m_g = (V - M * vs) / (vg - vs)
            m_s = M - m_g
            return m_s, m_g, us, ug, P

        def energy_residual(T):
            m_s, m_g, us, ug, _P = split(T)
            return m_s * us + m_g * ug - U

        f_lo, f_hi = energy_residual(lo), energy_residual(hi)
        if f_lo * f_hi > 0:
            # U outside the bracketable window: clamp to the nearer physical bound
            T = lo if abs(f_lo) < abs(f_hi) else hi
        else:
            T = optimize.brentq(energy_residual, lo, hi)
        m_s, m_g, _us, _ug, P = split(T)
        return {"m_s": m_s, "m_g": m_g, "m_l": 0.0, "T": T, "P": P}

    def vessel_state_below_triple(self, M, U, V, tol=1e-9):
        """Dispatch a below-/at-triple-point vessel state (M [kg], U [J], V [m3]).

        Returns a dict with ``regime`` in {"above", "triple", "sublimation"} and the
        phase masses (m_s, m_l, m_g), temperature and pressure. ``"above"`` signals the
        state has warmed back onto the liquid+vapour saturation line above the triple
        point (hand control back to the CoolProp model).
        """
        tl = self.triple_lever(M, U, V)
        if tl["m_s"] < -tol:
            # warmer than the triple point: liquid+vapour above it
            return {"regime": "above", **tl}
        if tl["m_l"] < -tol:
            # colder than the triple point, liquid exhausted: solid+gas
            st = self.sublimation_state(M, U, V)
            st["regime"] = "sublimation"
            return st
        # genuine three-phase point (clamp tiny negatives from round-off)
        tl["m_s"] = max(tl["m_s"], 0.0)
        tl["m_l"] = max(tl["m_l"], 0.0)
        tl["regime"] = "triple"
        return tl


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
