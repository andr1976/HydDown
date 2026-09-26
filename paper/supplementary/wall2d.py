# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license
"""
2-D axisymmetric (r-z) transient conjugate wall-conduction model
(:class:`WallConduction2D`).

Motivation
----------
The lumped and 1-D-radial wall models in :mod:`hyddown.wall` treat the vessel wall as two
independent columns (gas-contact and wetted) that never exchange heat, and they ignore the
extra thermal mass of the thick bottom plate and the wide/heavy top lid. For the SINTEF/NTNU
dense-phase CO2 tests (Hoydalsvik et al., 2026, Appl. Therm. Eng. 306, 133032) this makes the
modelled wetted wall run too cold at the liquid heel: in reality axial (lateral) conduction from
the warm gas-contact wall and the massive bottom plate re-warm the wetted wall.

This module reproduces the SINTEF approach: transient **axisymmetric** heat conduction through
the full steel domain (bottom plate + cylindrical shell + flange + lid), with an **adiabatic**
outer surface (the vessel is insulated) and a **Robin** boundary condition on the inner surface
that is **split at the moving liquid level** (pool-boiling coefficient below, free-convection
coefficient above). It is solved here with a cell-centred **finite-volume** method on a masked
structured r-z grid (equivalent to the SINTEF FEM for this flat geometry) using an implicit
backward-Euler step and a sparse direct solve. Zero new dependencies (numpy + scipy.sparse).

Cell types
----------
``STEEL``    - part of the wall, an unknown in the conduction problem;
``CAVITY``   - the internal fluid space; a steel face touching it carries the inner Robin BC;
``EXTERIOR`` - everything else (insulation/ambient side); steel faces touching it are adiabatic.

Sign convention: ``Q_gas`` / ``Q_liq`` returned by :meth:`step` are the heat flows FROM the wall
INTO the gas and liquid zones respectively (positive = wall heats the fluid), i.e. the same sense
as ``Q_inner`` in :mod:`hyddown.hdclass`.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

STEEL, CAVITY, EXTERIOR = 0, 1, 2


def k_ss316(T):
    """Approximate SS316 thermal conductivity [W/m/K] as a function of temperature [K].

    Austenitic stainless steel conductivity rises with temperature; over the -120..+50 C range
    of these tests a linear fit anchored at the room-temperature handbook value is adequate:
    ``k = 16.2 + 0.013*(T_C - 25)``  (approx. 16.2 at 25 C, 14.8 at -80 C, 14.0 at -140 C).
    Replace with the full k(T) correlation if higher fidelity is needed.
    """
    T_C = np.asarray(T) - 273.15
    return 16.2 + 0.013 * (T_C - 25.0)


def build_geometry(r_in, thickness, H, t_bot, t_flange, t_lid, r_lid):
    """Assemble the stepped-vessel steel rectangles + internal cavity from basic dimensions [m].

    z=0 is the inner floor (top of the bottom plate); r=0 is the axis. The cavity fills the bore
    up to the lid underside (z = H + t_flange). Set ``t_bot``/``t_flange``/``t_lid`` to 0 and
    ``r_lid = r_in + thickness`` to recover a plain flat-end cylinder.
    """
    r_out = r_in + thickness
    r_lid = max(r_lid, r_out)
    z_lidbot = H + t_flange
    steel = [
        # (r0, r1, z0, z1)
        (0.0,   r_out,  -t_bot,   0.0),                 # bottom plate
        (r_in,  r_out,   0.0,     H),                   # cylindrical shell
    ]
    if t_flange > 0:
        steel.append((r_in, r_lid, H, z_lidbot))        # flange ring
    if t_lid > 0:
        steel.append((0.0, r_lid, z_lidbot, z_lidbot + t_lid))  # lid disc
    cavity = (0.0, r_in, 0.0, z_lidbot)                 # internal fluid space (up to lid underside)
    return dict(steel=steel, cavity=cavity, r_in=r_in, r_out=r_out, r_lid=r_lid,
                H=H, t_bot=t_bot, t_flange=t_flange, t_lid=t_lid,
                z_top_cavity=z_lidbot)


def default_sintef_geometry():
    """SINTEF vessel geometry (Hoydalsvik et al., 2026): ID 273 mm (r_in 0.1365), wall 25.4 mm,
    internal height 1000 mm, bottom plate 50 mm, flange 83 mm, lid 80 mm x 580 mm dia."""
    return build_geometry(r_in=0.1365, thickness=0.0254, H=1.000,
                          t_bot=0.050, t_flange=0.083, t_lid=0.080, r_lid=0.290)


def _edges(breakpoints, target):
    """Structured 1-D edge array: subdivide each interval between sorted unique breakpoints so
    cells are ~<= target size, guaranteeing grid lines fall exactly on component boundaries."""
    bp = np.unique(np.asarray(breakpoints, float))
    edges = [bp[0]]
    for a, b in zip(bp[:-1], bp[1:]):
        n = max(1, int(np.ceil((b - a) / target - 1e-9)))
        edges.extend(np.linspace(a, b, n + 1)[1:])
    return np.array(edges)


class WallConduction2D:
    """Axisymmetric finite-volume conjugate wall solver for the stepped vessel geometry.

    Parameters
    ----------
    geometry : dict
        As returned by :func:`default_sintef_geometry` (steel rectangles + cavity rectangle).
    rho, cp : float
        Steel density [kg/m3] and specific heat [J/kg/K].
    T0 : float
        Uniform initial wall temperature [K] (the insulated wall equilibrates with the fill).
    k_func : callable or float
        ``k(T[K]) -> W/m/K``; a constant is wrapped automatically. Default :func:`k_ss316`.
    dr_wall, dr_bulk, dz_wall, dz_bulk : float
        Target cell sizes [m] in the thin shell band and elsewhere.
    """

    def __init__(self, geometry, rho, cp, T0, k_func=k_ss316,
                 dr_wall=0.004, dr_bulk=0.018, dz_wall=0.010, dz_bulk=0.028):
        self.geom = geometry
        self.rho = rho
        self.cp = cp
        self.k_func = k_func if callable(k_func) else (lambda T, _k=float(k_func): np.full_like(np.asarray(T, float), _k))
        self._build_grid(dr_wall, dr_bulk, dz_wall, dz_bulk)
        self._classify()
        self._build_faces()
        self.T = np.full(self.nunk, float(T0))

    # ---------------------------------------------------------------- grid
    def _build_grid(self, dr_wall, dr_bulk, dz_wall, dz_bulk):
        g = self.geom
        r_bp = [0.0, g["r_in"], g["r_out"], g["r_lid"]]
        z_bp = [-g["t_bot"], 0.0, g["H"], g["H"] + g["t_flange"], g["H"] + g["t_flange"] + g["t_lid"]]
        # fine radial spacing only across the shell band [r_in, r_out]
        re = _merge_targeted(r_bp, {(g["r_in"], g["r_out"]): dr_wall}, dr_bulk)
        ze = _edges(z_bp, dz_bulk)
        # refine the bottom plate and flange/lid thickness bands axially
        ze = _merge_targeted(z_bp, {(-g["t_bot"], 0.0): dz_wall,
                                    (g["H"], g["H"] + g["t_flange"]): dz_wall,
                                    (g["H"] + g["t_flange"], g["H"] + g["t_flange"] + g["t_lid"]): dz_wall},
                             dz_bulk)
        self.re, self.ze = re, ze
        self.nr, self.nz = len(re) - 1, len(ze) - 1
        self.rc = 0.5 * (re[:-1] + re[1:])
        self.zc = 0.5 * (ze[:-1] + ze[1:])
        self.dr = np.diff(re)
        self.dz = np.diff(ze)
        # cell volumes V = pi(r_e^2 - r_w^2) dz  (axisymmetric)
        ring = np.pi * (re[1:] ** 2 - re[:-1] ** 2)          # (nr,)
        self.vol = np.outer(ring, self.dz)                    # (nr, nz)
        self.ncell = self.nr * self.nz

    def _idx(self, ir, iz):
        return ir * self.nz + iz

    # ---------------------------------------------------------- classify
    def _classify(self):
        g = self.geom
        typ = np.full((self.nr, self.nz), EXTERIOR, dtype=np.int8)
        for ir in range(self.nr):
            for iz in range(self.nz):
                r, z = self.rc[ir], self.zc[iz]
                if _in_rect(r, z, g["cavity"]):
                    typ[ir, iz] = CAVITY
                elif any(_in_rect(r, z, rect) for rect in g["steel"]):
                    typ[ir, iz] = STEEL
        self.typ = typ
        self.steel_mask = typ == STEEL
        # unknown numbering over steel cells only
        self.uid = -np.ones((self.nr, self.nz), dtype=int)
        steel_cells = np.argwhere(self.steel_mask)
        for n, (ir, iz) in enumerate(steel_cells):
            self.uid[ir, iz] = n
        self.nunk = len(steel_cells)
        self.steel_cells = steel_cells

    # --------------------------------------------------------- face lists
    def _build_faces(self):
        """Pre-compute steel-steel conduction faces and steel-cavity (Robin) inner faces.

        Conduction face: (n_a, n_b, geomfac, ir_a, iz_a, ir_b, iz_b) where the conductance is
        ``k_face * geomfac`` and geomfac = A_face / distance. Inner face: (n, area, kind, ir, iz)
        with kind in {'side','floor','lid'} for classification against the liquid level; outer
        faces are adiabatic and simply omitted.
        """
        cond = []
        inner = []
        outer = []   # (n, area, kind) for reporting outer-face temperatures
        re, ze, rc, zc, dr, dz = self.re, self.ze, self.rc, self.zc, self.dr, self.dz
        for ir, iz in self.steel_cells:
            n = self.uid[ir, iz]
            # ---- radial (r-) faces ----
            for dirn, jr in ((-1, ir - 1), (+1, ir + 1)):
                rf = re[ir] if dirn < 0 else re[ir + 1]
                area = 2.0 * np.pi * rf * dz[iz]
                nb_type = self._neighbour_type(jr, iz, axis="r")
                if nb_type == STEEL:
                    if dirn > 0:  # add each internal face once (from the lower-index side)
                        dist = rc[ir + 1] - rc[ir]
                        cond.append((n, self.uid[ir + 1, iz], area / dist, ir, iz))
                elif nb_type == CAVITY:
                    # cavity is on the smaller-r side (inner cylindrical surface)
                    inner.append((n, area, "side", ir, iz))
                else:  # EXTERIOR or domain edge -> adiabatic
                    if (jr < 0 or jr >= self.nr):
                        outer.append((n, area, "r"))
                    elif nb_type == EXTERIOR:
                        outer.append((n, area, "r"))
            # ---- axial (z-) faces ----
            for dirn, jz in ((-1, iz - 1), (+1, iz + 1)):
                area = np.pi * (re[ir + 1] ** 2 - re[ir] ** 2)
                nb_type = self._neighbour_type(ir, jz, axis="z")
                if nb_type == STEEL:
                    if dirn > 0:
                        dist = zc[iz + 1] - zc[iz]
                        cond.append((n, self.uid[ir, iz + 1], area / dist, ir, iz))
                elif nb_type == CAVITY:
                    # cavity above steel (dirn +) -> floor ; cavity below steel (dirn -) -> lid underside
                    kind = "floor" if dirn > 0 else "lid"
                    inner.append((n, area, kind, ir, iz))
                else:
                    outer.append((n, area, "z"))
        self.cond = cond
        self.inner = inner
        self.outer = outer
        # static part of the conduction graph for assembly speed
        self._cond_np = (np.array([c[0] for c in cond]),
                         np.array([c[1] for c in cond]),
                         np.array([c[2] for c in cond]),
                         np.array([c[3] for c in cond]),
                         np.array([c[4] for c in cond])) if cond else None
        # flat volume per unknown
        self.vol_u = np.array([self.vol[ir, iz] for ir, iz in self.steel_cells])

        # Cylinder (shell) inner/outer face indices + their heights, for reporting a
        # cylinder-only wall temperature that matches the experimental thermocouples
        # (TT1x2 on the inner cylinder wall, TT1x1 on the outer, at six heights). The lid,
        # flange and bottom stay in the conduction solve but are excluded from the report.
        # Restrict the reported band to the instrumented span (the wall thermocouples sit at
        # 5-95% of the vessel height, not at the welded ends), so the model band matches the
        # sensor coverage rather than including the end cells.
        z_lo, z_hi = 0.05 * self.geom["H"], 0.95 * self.geom["H"]

        def _is_shell(ir, iz):
            r, z = self.rc[ir], self.zc[iz]
            return (self.geom["r_in"] - 1e-9) <= r <= (self.geom["r_out"] + 1e-9) \
                and z_lo <= z <= z_hi
        sin_n, sin_z, sout_n, sout_z = [], [], [], []
        for (nn, area, kind, ir, iz) in inner:
            if kind == "side" and _is_shell(ir, iz):
                sin_n.append(nn); sin_z.append(self.zc[iz])
        for (nn, area, orient) in outer:
            ir, iz = self.steel_cells[nn]
            if orient == "r" and _is_shell(ir, iz):
                sout_n.append(nn); sout_z.append(self.zc[iz])
        self.shell_inner_n = np.array(sin_n, dtype=int)
        self.shell_inner_z = np.array(sin_z, dtype=float)
        self.shell_outer_n = np.array(sout_n, dtype=int)
        self.shell_outer_z = np.array(sout_z, dtype=float)

    def _neighbour_type(self, ir, iz, axis):
        if ir < 0 or ir >= self.nr or iz < 0 or iz >= self.nz:
            return EXTERIOR
        return self.typ[ir, iz]

    # ------------------------------------------------------------- step
    def step(self, dt, liquid_level, liquid_present, h_gas, T_gas, h_liq, T_liq):
        """Advance the wall field by one implicit backward-Euler step.

        Inner faces below ``liquid_level`` (and the floor when ``liquid_present``) get the
        boiling Robin BC ``(h_liq, T_liq)``; all other inner faces get ``(h_gas, T_gas)``.
        Returns a summary dict with area-averaged inner/outer temperatures for the dry and wetted
        regions [K] and the heat flows into each fluid zone [W].
        """
        n = self.nunk
        # k at cell temperatures (lagged); face k = arithmetic mean of the two steel cells
        kcell = self.k_func(self.T)
        rows, cols, vals = [], [], []
        diag = self.rho * self.cp * self.vol_u / dt
        b = diag * self.T
        # conduction
        if self._cond_np is not None:
            na, nb, gf, ira, iza = self._cond_np
            kf = 0.5 * (kcell[na] + kcell[nb])
            g = kf * gf
            rows.extend(na); cols.extend(nb); vals.extend(-g)
            rows.extend(nb); cols.extend(na); vals.extend(-g)
            np.add.at(diag, na, g)
            np.add.at(diag, nb, g)
        # inner Robin BC
        for (nn, area, kind, ir, iz) in self.inner:
            wetted = (kind == "floor" and liquid_present) or \
                     (kind == "side" and self.zc[iz] < liquid_level)
            if wetted:
                h, Tf = h_liq, T_liq
            else:
                h, Tf = h_gas, T_gas
            diag[nn] += h * area
            b[nn] += h * area * Tf
        A = sp.coo_matrix((vals + list(diag), (rows + list(range(n)), cols + list(range(n)))),
                          shape=(n, n)).tocsr()
        self.T = spsolve(A, b)
        return self._summary(liquid_level, liquid_present, h_gas, T_gas, h_liq, T_liq)

    # ----------------------------------------------------------- reports
    def _summary(self, liquid_level, liquid_present, h_gas, T_gas, h_liq, T_liq):
        Tiw_dry_num = Tiw_dry_area = 0.0
        Tiw_wet_num = Tiw_wet_area = 0.0
        Qg = Ql = 0.0
        for (nn, area, kind, ir, iz) in self.inner:
            wetted = (kind == "floor" and liquid_present) or \
                     (kind == "side" and self.zc[iz] < liquid_level)
            Tw = self.T[nn]
            if wetted:
                Tiw_wet_num += Tw * area; Tiw_wet_area += area
                Ql += h_liq * area * (Tw - T_liq)
            else:
                Tiw_dry_num += Tw * area; Tiw_dry_area += area
                Qg += h_gas * area * (Tw - T_gas)
        # outer-face area-averaged temperatures, split wetted/dry by z vs level
        Tow_dry_num = Tow_dry_area = 0.0
        Tow_wet_num = Tow_wet_area = 0.0
        for (nn, area, orient) in self.outer:
            ir, iz = self.steel_cells[nn]
            below = self.zc[iz] < liquid_level or (self.zc[iz] < 0.0 and liquid_present)
            Tw = self.T[nn]
            if below:
                Tow_wet_num += Tw * area; Tow_wet_area += area
            else:
                Tow_dry_num += Tw * area; Tow_dry_area += area

        def avg(num, ar, fallback):
            return num / ar if ar > 0 else fallback
        Ti_dry = avg(Tiw_dry_num, Tiw_dry_area, np.nan)
        Ti_wet = avg(Tiw_wet_num, Tiw_wet_area, Ti_dry)

        # Cylinder-only report (matches the wall thermocouples): band statistics over the shell
        # faces (all heights, mixing gas and wetted like the six sensors) plus a gas/wetted split.
        def shell_stats(idx, zarr):
            out = dict(min=np.nan, max=np.nan, med=np.nan, dry=np.nan, wet=np.nan)
            if len(idx) == 0:
                return out
            T = self.T[idx]
            out["min"], out["max"], out["med"] = float(T.min()), float(T.max()), float(np.median(T))
            wet = zarr < liquid_level
            if np.any(wet):
                out["wet"] = float(T[wet].mean())
            if np.any(~wet):
                out["dry"] = float(T[~wet].mean())
            return out
        si = shell_stats(self.shell_inner_n, self.shell_inner_z)
        so = shell_stats(self.shell_outer_n, self.shell_outer_z)
        return dict(
            T_inner_dry=Ti_dry,
            T_inner_wet=Ti_wet,
            T_outer_dry=avg(Tow_dry_num, Tow_dry_area, np.nan),
            T_outer_wet=avg(Tow_wet_num, Tow_wet_area, avg(Tow_dry_num, Tow_dry_area, np.nan)),
            Q_gas=Qg, Q_liq=Ql,
            T_mean=float(self.T.mean()),
            # cylinder-only (thermocouple-comparable) inner/outer wall temperatures
            T_cyl_in_min=si["min"], T_cyl_in_max=si["max"], T_cyl_in_med=si["med"],
            T_cyl_in_dry=si["dry"], T_cyl_in_wet=si["wet"],
            T_cyl_out_min=so["min"], T_cyl_out_max=so["max"], T_cyl_out_med=so["med"],
            T_cyl_out_dry=so["dry"], T_cyl_out_wet=so["wet"],
        )

    def sample_heights(self, heights):
        """Inner- and outer-cylinder-wall temperatures [K] at the given axial heights [m].

        Returns ``(T_inner, T_outer)`` arrays, each the nearest shell inner/outer cell temperature
        to each requested height, for comparison with wall thermocouples at fixed positions.
        """
        heights = np.atleast_1d(np.asarray(heights, float))
        Ti = np.full(len(heights), np.nan)
        To = np.full(len(heights), np.nan)
        for k, z in enumerate(heights):
            if len(self.shell_inner_n):
                j = int(np.argmin(np.abs(self.shell_inner_z - z)))
                Ti[k] = self.T[self.shell_inner_n[j]]
            if len(self.shell_outer_n):
                j = int(np.argmin(np.abs(self.shell_outer_z - z)))
                To[k] = self.T[self.shell_outer_n[j]]
        return Ti, To

    def energy(self):
        """Total stored thermal energy [J] relative to 0 K (rho*cp*V*T summed over steel)."""
        return float(np.sum(self.rho * self.cp * self.vol_u * self.T))


def _in_rect(r, z, rect):
    r0, r1, z0, z1 = rect
    return (r0 - 1e-12) <= r <= (r1 + 1e-12) and (z0 - 1e-12) <= z <= (z1 + 1e-12)


def _merge_targeted(breakpoints, interval_targets, default_target):
    """Like :func:`_edges` but allows a per-interval target size (keyed by (a,b) breakpoint pair)."""
    bp = np.unique(np.asarray(breakpoints, float))
    edges = [bp[0]]
    for a, b in zip(bp[:-1], bp[1:]):
        tgt = default_target
        for (ka, kb), t in interval_targets.items():
            if abs(ka - a) < 1e-9 and abs(kb - b) < 1e-9:
                tgt = t
        nseg = max(1, int(np.ceil((b - a) / tgt - 1e-9)))
        edges.extend(np.linspace(a, b, nseg + 1)[1:])
    return np.array(edges)
