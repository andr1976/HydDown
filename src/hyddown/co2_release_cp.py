# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license

"""
CoolProp-only CO2 release model (no thermopack runtime dependency).

This is a drop-in replacement for :class:`hyddown.co2_release.CO2ReleaseModel` that removes
the thermopack (Fortran) dependency. It reuses all of the base class's lever/table math and
overrides only the handful of methods that previously called thermopack, replacing them with:

  * **CoolProp** for the fluid - native above the triple point, and with an imposed gas/liquid
    phase (``AbstractState.specify_phase``) below it, which makes CoolProp extrapolate its
    Span-Wagner EoS past the triple-point guard (verified to match thermopack GERG-2008 to
    <0.01 % on density, ~0.2-0.3 % on reference-aligned h/s down to the frost point);
  * a small **solid-CO2 property table** (:mod:`hyddown.co2_solid`), generated offline and
    re-referenced into CoolProp's h/s basis via the latent heat of sublimation - the one thing
    CoolProp cannot supply (it has no solid model).

Why this works without a solid-aware flash: below the triple point pure CO2 is a one-component
system, so (Gibbs' phase rule) the vapour+solid region has a single degree of freedom and every
"flash" reduces to lever-rule arithmetic on the sublimation line - which reproduces thermopack's
solid-aware ``two_phase_psflash`` *exactly* (T to 0.01 K, phase fractions to 4 decimals). See the
project notes for the derivation. The public API is identical to the thermopack version.
"""

import math

import numpy as np
from scipy import optimize

import CoolProp as CP
from CoolProp.CoolProp import AbstractState, PropsSI

from hyddown import co2_solid
from hyddown.co2_release import CO2ReleaseModel, P_CRIT


class CO2ReleaseModelCP(CO2ReleaseModel):
    """CoolProp + solid-table implementation of :class:`CO2ReleaseModel` (thermopack-free)."""

    def __init__(self, back_pressure=101325.0, atm_pressure=101325.0, eos="CoolProp",
                 liquid_nonequilibrium=0.0, liquid_ne_pressure_scaled=False,
                 liquid_ne_pref=None):
        # NB: do NOT call super().__init__ (it would construct thermopack). Set up the same
        # attributes the inherited methods expect, backed by CoolProp + the solid table.
        self.p_back = back_pressure
        self.p_atm = atm_pressure
        self.eos = None  # no thermopack backend
        self.M = PropsSI("molar_mass", "CO2")
        self.liquid_nonequilibrium = float(liquid_nonequilibrium)
        self.liquid_ne_pressure_scaled = bool(liquid_ne_pressure_scaled)
        self.liquid_ne_pref = liquid_ne_pref
        self.P_TRIPLE = co2_solid.P_TRIPLE
        # forced-gas / forced-liquid states for sub-triple queries (CoolProp phase imposition)
        self._gas = AbstractState("HEOS", "CO2"); self._gas.specify_phase(CP.iphase_gas)
        self._init_atm_endpoints()
        self._init_triple_point()
        self._init_gas_table()

    # ------------------------------------------------------------- CoolProp helpers
    def _gasp(self, T, P):
        """Forced-gas CO2 properties at (T, P): (h [J/kg], s [J/kg/K], v [m3/kg])."""
        self._gas.update(CP.PT_INPUTS, P, T)
        return self._gas.hmass(), self._gas.smass(), 1.0 / self._gas.rhomass()

    def _gas_T_from_s(self, s, P, lo, hi=360.0):
        """Invert forced-gas entropy s(T,P)=s for T at fixed P (monotone in T)."""
        return optimize.brentq(lambda T: self._gasp(T, P)[1] - s, lo, hi)

    def _gas_T_from_h(self, h, P, lo, hi=420.0):
        """Invert forced-gas enthalpy h(T,P)=h for T at fixed P."""
        return optimize.brentq(lambda T: self._gasp(T, P)[0] - h, lo, hi)

    # ------------------------------------------------------------------ setup
    def _init_triple_point(self):
        T = self.T_TRIPLE_EOS = co2_solid.T_TRIPLE
        P = self.P_TRIPLE_EOS = co2_solid.P_TRIPLE
        # solid vertex from the table; liquid/vapour vertices from CoolProp saturation at Ptp
        self.v_s = co2_solid.v_solid(T); self.h_s = co2_solid.h_solid(T)
        self.v_l = 1.0 / PropsSI("Dmass", "P", P, "Q", 0, "CO2")
        self.h_l = PropsSI("Hmass", "P", P, "Q", 0, "CO2")
        self.v_g = 1.0 / PropsSI("Dmass", "P", P, "Q", 1, "CO2")
        self.h_g = PropsSI("Hmass", "P", P, "Q", 1, "CO2")
        self.u_s = self.h_s - P * self.v_s
        self.u_l = self.h_l - P * self.v_l
        self.u_g = self.h_g - P * self.v_g
        self.L_sub = self.h_g - self.h_s
        self.cp_solid = co2_solid.cp_solid(T)
        self.T_subl_floor = co2_solid.T_sub(self.p_back)

    def _init_atm_endpoints(self):
        self.T_frost = co2_solid.T_sub(self.p_atm)
        self.h_gas_atm = self._gasp(self.T_frost, self.p_atm)[0]
        self.h_solid_atm = co2_solid.h_solid(self.T_frost)

    def _init_gas_table(self):
        """Precompute the same superheated-gas tables as the base class, from CoolProp."""
        T = np.linspace(self.T_TRIPLE_EOS - 0.2, 345.0, 180)
        u = np.empty_like(T); rho = np.empty_like(T); h = np.empty_like(T); G = np.empty_like(T)
        P = self.P_TRIPLE_EOS
        for i, Ti in enumerate(T):
            hi, _si, v = self._gasp(Ti, P)
            u[i] = hi - P * v; rho[i] = 1.0 / v; h[i] = hi
            G[i] = self.gas_leak_rate(Ti, P, 1.0, 1.0)["mdot"]
        self._gt_T, self._gt_u, self._gt_rho, self._gt_h, self._gt_G = T, u, rho, h, G

        Pgrid = np.linspace(self.p_back * 0.85, self.P_TRIPLE_EOS * 1.001, 55)
        Tgrid = np.linspace(self.T_TRIPLE_EOS - 30.0, 345.0, 70)
        RHO = np.empty((Tgrid.size, Pgrid.size)); U2 = np.empty_like(RHO); H2 = np.empty_like(RHO)
        for a, Ti in enumerate(Tgrid):
            for b, Pi in enumerate(Pgrid):
                hi, _si, v = self._gasp(Ti, Pi)
                RHO[a, b] = 1.0 / v; U2[a, b] = hi - Pi * v; H2[a, b] = hi
        self._g2_T, self._g2_P, self._g2_rho, self._g2_u, self._g2_h = Tgrid, Pgrid, RHO, U2, H2

        self._sl_P = np.linspace(self.p_back * 0.9, self.P_TRIPLE_EOS * 0.999, 50)
        self._sl_T = np.array([co2_solid.T_sub(Pi) for Pi in self._sl_P])
        self._sl_us = np.array([co2_solid.u_solid(Ti, Pi)
                                for Ti, Pi in zip(self._sl_T, self._sl_P)])

    # ------------------------------------------------------------- stagnation
    def stagnation(self, P, phase):
        if P >= P_CRIT:
            raise ValueError(
                f"Stagnation pressure {P/1e5:.2f} bar >= CO2 critical pressure "
                f"{P_CRIT/1e5:.2f} bar - saturated stagnation is undefined.")
        Q = 0 if phase == "liquid" else 1
        Tsat = PropsSI("T", "P", P, "Q", Q, "CO2")
        h0 = PropsSI("Hmass", "P", P, "Q", Q, "CO2")
        s0 = PropsSI("Smass", "P", P, "Q", Q, "CO2")  # mass-specific (CoolProp basis)
        rho0 = PropsSI("Dmass", "P", P, "Q", Q, "CO2")
        return h0, s0, rho0, Tsat

    # ------------------------------------------------ isentrope mixture props
    def _iso_props(self, P, s0):
        """Mass enthalpy, density, T and solid fraction on the isentrope (s = s0) at P.

        Above the triple point: CoolProp PS-flash (vapour+liquid or single phase).
        Below it: a lever on the sublimation line (vapour+solid), the exact 1-DOF replacement
        for thermopack's solid-aware psflash. ``s0`` is mass-specific.
        """
        if P >= self.P_TRIPLE_EOS:
            h = PropsSI("Hmass", "P", P, "Smass", s0, "CO2")
            rho = PropsSI("Dmass", "P", P, "Smass", s0, "CO2")
            T = PropsSI("T", "P", P, "Smass", s0, "CO2")
            return h, rho, T, 0.0
        T = co2_solid.T_sub(P)
        hg, sg, vg = self._gasp(T, P)
        hs = co2_solid.h_solid(T); ss = co2_solid.s_solid(T); vs = co2_solid.v_solid(T)
        if s0 >= sg:  # superheated vapour below the triple point (no solid)
            Tsup = self._gas_T_from_s(s0, P, T)
            hh, _ss, vv = self._gasp(Tsup, P)
            return hh, 1.0 / vv, Tsup, 0.0
        bv = min(max((s0 - ss) / (sg - ss), 0.0), 1.0)
        bs = 1.0 - bv
        h = bv * hg + bs * hs
        v = bv * vg + bs * vs
        return h, 1.0 / v, T, bs

    def gas_leak_rate(self, T, P, Cd, area):
        h0, s0, v0 = self._gasp(T, P)
        return self.hem_rate_from_stagnation(h0, s0, 1.0 / v0, T, P, Cd, area)

    def dense_leak_rate(self, T, P, Cd, area, phase="liquid"):
        """HEM rate for a single-phase (subcooled/dense or supercritical) discharge at the
        actual vessel state (T, P) - not a saturated bubble/dew point. This is what feeds the
        outlet while the vessel is single-phase, before it flashes: there is only one phase,
        so the draw is phase-agnostic (a liquid-space or vapour-space outlet both draw it).
        Bypasses stagnation()'s P >= P_crit guard by taking the state directly."""
        h0 = PropsSI("Hmass", "T", T, "P", P, "CO2")
        s0 = PropsSI("Smass", "T", T, "P", P, "CO2")
        rho0 = PropsSI("Dmass", "T", T, "P", P, "CO2")
        return self.hem_rate_from_stagnation(h0, s0, rho0, T, P, Cd, area)

    # ------------------------------------------------------- atmospheric state
    def atm_split(self, h0_mass):
        beta_gas = (h0_mass - self.h_solid_atm) / (self.h_gas_atm - self.h_solid_atm)
        if beta_gas >= 1.0:  # superheated: no dry ice
            T = self._gas_T_from_h(h0_mass, self.p_atm, self.T_frost)
            return {"T": T, "vapour_frac": 1.0, "solid_frac": 0.0}
        beta_gas = max(beta_gas, 0.0)
        return {"T": self.T_frost, "vapour_frac": beta_gas, "solid_frac": 1.0 - beta_gas}

    def atm_split_isentropic(self, s0):
        T = self.T_frost
        hg, sg, vg = self._gasp(T, self.p_atm)
        ss = co2_solid.s_solid(T)
        if s0 >= sg:  # superheated vapour at 1 atm
            Tsup = self._gas_T_from_s(s0, self.p_atm, T)
            return {"T": Tsup, "vapour_frac": 1.0, "solid_frac": 0.0}
        bv = min(max((s0 - ss) / (sg - ss), 0.0), 1.0)
        return {"T": T, "vapour_frac": bv, "solid_frac": 1.0 - bv}

    # ---------------------------------------------- sublimation-line helpers
    def _sublimation_pressure(self, T):
        return co2_solid.P_sub(T)

    def _sublimation_props(self, T):
        P = co2_solid.P_sub(T)
        hg, _sg, vg = self._gasp(T, P)
        vs = co2_solid.v_solid(T)
        us = co2_solid.u_solid(T, P)
        return P, vg, hg - P * vg, vs, us


# --------------------------------------------------------------------------- self-test
if __name__ == "__main__":
    m = CO2ReleaseModelCP()
    print(f"M = {m.M:.5f} kg/mol   frost point = {m.T_frost:.3f} K")
    print(f"triple: T = {m.T_TRIPLE_EOS:.3f} K  P = {m.P_TRIPLE_EOS/1e5:.4f} bar")
    print(f"{'release':>22} | {'P0[bar]':>8} | {'isenth':>7} | {'isentr':>7}")
    for label, TC, ph in (("sat-liquid +17C", 17.0, "liquid"), ("sat-liquid -30C", -30.0, "liquid"),
                          ("sat-vapour +17C", 17.0, "gas"), ("sat-vapour -30C", -30.0, "gas")):
        T = 273.15 + TC
        P0 = PropsSI("P", "T", T, "Q", 0, "CO2")
        h0, s0, _r, _T = m.stagnation(P0, ph)
        se = m.atm_split(h0)["solid_frac"]
        ss = m.atm_split_isentropic(s0)["solid_frac"]
        print(f"{label:>22} | {P0/1e5:8.2f} | {se:7.3f} | {ss:7.3f}")
    rate, atm = m.release_state(P0=53.39e5, phase="liquid", Cd=0.62, area=math.pi * 0.05 ** 2 / 4)
    print(f"\nHEM sat-liquid +17C, 50 mm hole: mdot = {rate['mdot']:.2f} kg/s, "
          f"throat {rate['P_throat']/1e5:.2f} bar / {rate['T_throat']:.1f} K, "
          f"solid_frac_throat = {rate['solid_frac_throat']:.3f}, choked = {rate['choked']}")
    print(f"atmospheric: T = {atm['T']:.1f} K, vapour = {atm['vapour_frac']:.3f}, "
          f"dry ice = {atm['solid_frac']:.3f}")
