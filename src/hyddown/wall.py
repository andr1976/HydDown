# HydDown hydrogen/other gas depressurisation
# Copyright (c) 2021-2025 Anders Andreasen
# Published under an MIT license
"""
1-D transient wall-conduction model (:class:`WallConduction`), a thin object wrapper around the
:mod:`hyddown.thermesh` finite-element solver.

HydDown's ``detailed`` heat-transfer mode resolves the through-thickness temperature gradient in
the vessel wall (Type III/IV vessels, thick steel, composite liner+shell). The same mesh/domain/
solve sequence was previously inlined and duplicated eight times in the main loop - single-layer
and composite, each for the unwetted (gas-contact) and wetted (liquid-contact) faces, and again
in the fire path. This class owns one wall face's mesh, material model and temperature profile,
so the loop holds two instances (``wall``, ``wall_wetted``) and advances each with ``step``.

Node/BC conventions are preserved exactly from the original inline code:
  * single-layer: ``z = [0 .. thickness]`` -> node 0 is the OUTER face, node -1 the INNER face;
  * composite:    ``z = [-liner .. shell]`` -> node 0 is the INNER face, node -1 the OUTER face,
                  and the liner/shell bond is at node ``n_nodes - 1``.
The initial steady profile is solved with Dirichlet ``[{T: T0}, {T: Tamb}]`` on the two ends
(first node T0, last node Tamb), matching the original for both constructions. Each ``step`` then
applies flux boundary conditions (inner ``-q_inner``, outer ``+q_outer``) in the end order that
matches the node convention.
"""
import numpy as np

from hyddown import thermesh as tm


class WallConduction:
    """Transient 1-D conduction through one vessel-wall face (single-layer or composite).

    The node/BC convention is parametrized so one class reproduces every original inline block:
    the ``specified_h``/``detailed`` single-layer wall (inner face at the LAST node), and the
    composite and fire (s-b) walls (inner face at the FIRST node). ``inner_first`` selects the
    convention; ``init_ends``/``init_guess`` set the Dirichlet initial-profile solve (the detailed
    path uses a ``(T0, Tamb)`` gradient, the fire path a uniform ``(T0, T0)``).

    Parameters
    ----------
    thickness : float
        Shell (outer layer) thickness [m].
    k, rho, cp : float
        Shell thermal conductivity [W/m/K], density [kg/m3], heat capacity [J/kg/K].
    init_ends : (float, float)
        Dirichlet temperatures [K] on the first and last node for the initial steady profile.
    init_guess : float
        Uniform initial guess [K] for the initial-profile solve.
    n_nodes : int
        Nodes through each layer (default 11).
    theta : float
        Time-integration parameter (0.5 = Crank-Nicolson).
    liner : tuple or None
        ``(thickness, k, rho, cp)`` of an inner liner for a composite wall, or ``None``.
    inner_first : bool
        True if node 0 is the INNER (fluid-side) face - composite and fire walls; False if node 0
        is the OUTER face - the detailed single-layer wall.
    """

    def __init__(self, thickness, k, rho, cp, init_ends, init_guess,
                 n_nodes=11, theta=0.5, liner=None, inner_first=True):
        self.theta = theta
        self.n_nodes = n_nodes
        self.inner_first = inner_first
        self.composite = liner is not None
        if not self.composite:
            z = np.linspace(0.0, thickness, n_nodes)
            self.mesh = tm.Mesh(z, tm.LinearElement)
            self.models = [tm.isothermal_model(k, rho, cp)]
        else:
            lt, lk, lrho, lcp = liner
            z = np.hstack((np.linspace(-lt, 0.0, n_nodes),
                           np.linspace(0.0, thickness, n_nodes)[1:]))
            self.mesh = tm.Mesh(z, tm.LinearElement)
            for j, elem in enumerate(self.mesh.elem):
                if elem.nodes.mean() > 0.0:
                    self.mesh.subdomain[j] = 1
            self.models = [tm.isothermal_model(lk, lrho, lcp),
                           tm.isothermal_model(k, rho, cp)]
        self.inner_idx = 0 if inner_first else -1
        self.outer_idx = -1 if inner_first else 0
        self.bonded_idx = (n_nodes - 1) if self.composite else None
        self.z = z
        # initial steady profile (Dirichlet on the two ends, first node = init_ends[0])
        domain = tm.Domain(self.mesh, self.models,
                           [{"T": init_ends[0]}, {"T": init_ends[1]}])
        domain.set_T(init_guess * np.ones(len(self.mesh.nodes)))
        _t, prof = tm.solve_ht(domain, {"dt": 100, "t_end": 10000, "theta": theta})
        self.profile = prof[-1, :]

    def faces(self):
        """Current (T_inner, T_outer, T_bonded) [K]; T_bonded is None for a single layer."""
        T_bonded = self.profile[self.bonded_idx] if self.composite else None
        return self.profile[self.inner_idx], self.profile[self.outer_idx], T_bonded

    def step(self, q_inner, q_outer, dt, tstep):
        """Advance one macro timestep with inner/outer wall heat fluxes [W/m2].

        ``q_inner`` is the heat flux INTO the fluid at the inner face and ``q_outer`` the flux
        into the wall from the environment at the outer face, matching the original sign use
        (inner BC ``-q_inner``, outer BC ``+q_outer``). ``dt`` is the FE sub-step and ``tstep``
        the macro step. Returns ``(T_inner, T_outer, T_bonded)`` [K].
        """
        if self.inner_first:
            bc = [{"q": -q_inner}, {"q": q_outer}]      # node 0 = inner, node -1 = outer
        else:
            bc = [{"q": q_outer}, {"q": -q_inner}]      # node 0 = outer, node -1 = inner
        domain = tm.Domain(self.mesh, self.models, bc)
        domain.set_T(self.profile)
        _t, prof = tm.solve_ht(domain, {"dt": dt, "t_end": tstep, "theta": self.theta})
        self.profile = prof[-1, :]
        return self.faces()
