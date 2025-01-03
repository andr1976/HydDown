# Thermeshc code adapted verbatim from https://github.com/wjbg/thermesh
#
from __future__ import annotations
import numpy as np
from typing import Callable, Union
import numpy.typing as npt


vector = matrix = npt.NDArray[np.float64]


def solve_ht(domain: Domain, solver: dict[str, float]) -> tuple[vector, matrix]:
    """Solves heat transfer problem.

    Parameters
    ----------
    domain : Domain
        Domain with a mesh, material model and boundary conditions.
    solver : dict
        Solver information, dictionary with the following keys: dt,
        t_end, theta

    Returns
    -------
    t : nd.array(dtype=float, dim=1, len=int(t_end/dt))
        Times.
    T : nd.array(dtype=float, dim=2, shape=(int(t_end/dt), mesh.nn))
        Temperature data.

    """
    timesteps = int(solver["t_end"] / solver["dt"])

    # Data storage
    T = np.zeros((timesteps + 1, domain.mesh.nn))
    t = np.zeros(timesteps + 1)
    T[0] = domain.T
    t[0] = domain.t

    # Solve problem
    for i in range(1, timesteps + 1):
        t[i], T[i] = domain.timestep(solver["dt"], solver["theta"])
    return t, T


class Domain:
    """Class to represent a domain.

    Attributes
    ==========
    mesh : Mesh
        Object with mesh information.
    bc : two-item list of dicts
        The boundary conditions are provided in a two-item list of
        dictionaries. The first dictionary (or zeroth item in the list)
        applies to the start or left side of the domain, while the second
        item applies to the end or right side of the domain. The dictionaries
        can have the keys: "T" OR ( ("h" and "T_inf") AND/OR "q" ), with "T"
        an applied temperature, "h" and "T_inf" the convective heat transfer
        coefficient and far field temperature, respectively, while "q"
        represents a direct flux on the surface which is directed inwards.
    constitutive_model : list[function]
        Functions that takes temperature as an input and provides a
        dictionary with the following keys: k, cp, rho. The list has a
        length equal to the number of subdomains in the mesh. The i-th
        function in the list belongs to the i-th subdomain.
    t : float
        Time.
    T : nd.array(dim=1, len=mesh.nn)
        Temperature at time t.
    q : nd.array(dim=1, len=mesh.nn)
        Heat flux at time t.

    Methods
    =======
    timestep(dt, theta=0.5)
        Apply a timestep dt and update data.

    system_matrices()
        Returns domain stiffness and dampling matrix.

    check_bc()
        Checks if bc's are valid.

    set_T(T) and set_q(q)
        Set temperature and heat flux.

    clear()
        Clears data.

    """

    def __init__(
        self, mesh: Mesh, constitutive_model: list[Callable], bc: dict[str, float] = {}
    ):

        self.mesh = mesh
        self.constitutive_model = constitutive_model
        self.bc = bc

        self.t = 0
        self.T = np.zeros(self.mesh.nn)
        self.q = np.zeros(self.mesh.nn)

    def timestep(self, dt: float, theta: float = 0.5):
        """Apply timestep and update data.

        Parameters
        ----------
        dt : float
            Timestep size.
        theta: float (0 < theta <= 1)
            Timestepping type.

        NOTE: This is far from optimized as each matrix and vector is
              repeatedly constructed.

        NOTE: A more advanced solution strategy should be used when the
              problem becomes nonlinear, e.g. due to a radiative bc or
              due to temperature-dependent material properties.

        NOTE: The internal heat source Q is not implemented yet.

        """
        if self.check_bc():
            K, C = self.system_matrices()
            T_new, q_new = np.zeros(self.mesh.nn), np.zeros(self.mesh.nn)
            H = np.zeros(K.shape)
            j = np.arange(self.mesh.nn)

            for i in [0, -1]:
                if "T" in self.bc[i].keys():
                    T_new[i] = self.bc[i]["T"]
                    j = np.delete(j, i)
                if "q" in self.bc[i].keys():
                    q_new[i] += self.bc[i]["q"]
                if "h" in self.bc[i].keys() and "T_inf" in self.bc[i].keys():
                    H[i, i] = self.bc[i]["h"]
                    q_new[i] += self.bc[i]["h"] * self.bc[i]["T_inf"]

            A = C + dt * theta * (K + H)  # Matrix on LHS
            f = (
                1 - theta
            ) * dt * self.q + dt * theta * q_new  # misses contribution from Q
            b = np.dot(C - dt * (1 - theta) * (K + H), self.T) + f  # RHS
            b -= np.dot(A, T_new)  # take into account T boundary condition
            T_new[j] = np.linalg.solve(A[np.ix_(j, j)], b[j])
            self.T = T_new
            self.q = q_new
            self.t += dt
        return self.t, self.T

    def system_matrices(self) -> tuple[matrix, matrix]:
        """Returns domain stiffness and damping matrix."""
        K = np.zeros((self.mesh.nn, self.mesh.nn))
        C = np.zeros((self.mesh.nn, self.mesh.nn))
        for c, elem, sd in zip(self.mesh._conn, self.mesh.elem, self.mesh.subdomain):
            mat = self.constitutive_model[sd](self.T[c].mean())
            cols = np.tile(c, [len(c), 1])
            rows = np.tile(c[np.newaxis].T, [1, len(c)])
            K[rows, cols] += elem.K(mat)
            C[rows, cols] += elem.C(mat)
        return K, C

    def check_bc(self) -> bool:
        """Check if boundary conditions are valid."""
        for bc in self.bc:
            if "T" in bc.keys():
                if "q" in bc.keys() or "h" in bc.keys():
                    raise KeyError("invalid combination of bc's")
            if "q" in bc.keys():
                if "T" in bc.keys():
                    raise KeyError("invalid combination of bc's")
            if "h" in bc.keys():
                if "T" in bc.keys():
                    raise KeyError("invalid combination of bc's")
                if "T_inf" not in bc.keys():
                    raise KeyError("no temperature provided for convective bc")
        return True

    def set_T(self, T: Union[float, vector]):
        """Set temperature.

        Parameter
        ---------
        T : float OR np.ndarray(dim=1, dtype=float, len=nn)
            Temperature at nodes.

        """
        if type(T) == float:
            self.T = T * np.ones(self.mesh.nn)
        else:
            self.T = T * np.ones(self.mesh.nn)

    def set_q(self, q: Union[float, vector]):
        """Set heat flux.

        Parameter
        ---------
        q : float OR np.ndarray(dim=1, dtype=float, len=nn)
            Heat flux at nodes.

        """
        if type(q) == float:
            self.q = q * np.ones(self.mesh.nn)
        else:
            self.q = q * np.ones(self.mesh.nn)

    def clear(self):
        """Clears time and temperature data."""
        self.t = 0
        self.T = np.zeros(self.mesh.nn)
        self.q = np.zeros(self.mesh.nn)


def isothermal_model(k: float, rho: float, cp: float) -> Callable:
    """Returns a function that represents an isothermal material model.

    Parameter
    ---------
    k : float
        Thermal conductivity.
    rho : float
        Density
    cp : float
        Specific heat.

    Returns
    -------
    model : Callable
        Function that returns a dictionary with the provided
        constitutive properties.

    """

    def model(T: float) -> dict[str, float]:
        """A isothermal constitutive model.

        Parameter
        ---------
        T : float
            Temperature

        Returns
        -------
        mat : dict[str, float]
            Dictionary with the following keys: k, cp, rho.

        """
        return {"k": k, "rho": rho, "cp": cp}

    return model


def piecewise_linear_model(k: Matrix, rho: Matrix, cp: Matrix) -> Callable:
    """Returns a function that represents an isothermal material model.

    Parameter
    ---------
    k : np.ndarray(dim=2, dtype=float)
        Temperature vs. thermal conductivity.
    rho : np.ndarray(dim=2, dtype=float)
        Temperature vs. density
    cp : np.ndarray(dim=2, dtype=float)
        Temperature vs. specific heat.

    Returns
    -------
    model : Callable
        Function that returns a dictionary with the provided
        constitutive properties.

    """

    def model(T: float) -> dict[str, float]:
        """An piece-wise linear constitutive model.

        Parameter
        ---------
        T : float
            Temperature

        Returns
        -------
        mat : dict[str, float]
            Dictionary with the following keys: k, rho, cp.

        """
        k_ = np.interp(T, k[:, 0], k[:, 1])  # conductivity at T
        rho_ = np.interp(T, rho[:, 0], rho[:, 1])
        cp_ = np.interp(T, cp[:, 0], cp[:, 1])
        return {"k": k_, "rho": rho_, "cp": cp_}

    return model


class Mesh:
    """Class to represent a mesh.

    Attributes
    ==========
    nodes : nd.array()
        Node locations.
    elem : list[Element]
        List of elements.
    nn : int
        Number of nodes.
    nel : int
        Number of elements.
    subdomain : list[int] (defaults to a list with zeros)
        Number that indicates the subdomain in the mesh. Correlates
        with the constitutive model that will be used for the
        analysis.

    """

    def __init__(self, z: vector, element: LinearElement):
        """Initializes Mesh instance.

        Parameters
        ----------
        z : np.ndarray(dim=1, dtype=float)
            Node locations.
        element : Element
            Element type.

        """
        if (len(z) - 1) % (element.order) == 0:
            self.nn = len(z)
            self.nodes = z
            self._conn = np.array(
                [
                    np.arange(i, i + 1 + element.order)
                    for i in np.arange(0, self.nn - 1, element.order)
                ]
            )
            self.elem = [element(self.nodes[c]) for c in self._conn]
            self.nel = len(self.elem)
            self.subdomain = [0] * self.nel
        else:
            raise ValueError("node and element number mismatch")

    def __str__(self) -> str:
        s = "Mesh information\n"
        s += "----------------\n"
        s += f"Nodes: {self.nn:<3}\n"
        s += f"Elements: {self.nel:<3}\n"
        s += f"Element type: {type(self.elem[0])}\n"
        return s


class Element:
    """Class to represent an element.

    Attributes
    ==========
    order : int
        Order of the element.
    dim : int
        Dimension of the element.

    """

    dim = 1
    order = 0

    def __init__(self, nodes: vector):
        """Initializes Element instance.

        Parameters
        ----------
        nodes : np.ndarray(dim=1, dtype=float)
            Node locations.

        """
        self.nodes = nodes

    def length(self) -> float:
        return self.nodes[-1] - self.nodes[0]


class LinearElement(Element):
    """Class to represent a linear element (order = 1).

    Methods
    =======
    K()
        Returns element stiffness matrix.
    C()
        Returns element damping matrix.

    """

    order = 1

    def K(self, mat: dict[str, float]) -> matrix:
        """Returns element stiffness matrix."""
        K = (mat["k"] / self.length()) * np.array([[1, -1], [-1, 1]])
        return K

    def C(self, mat: dict[str, float]):
        """Returns element damping matrix."""
        C = (self.length() * mat["rho"] * mat["cp"] / 6) * np.array([[2, 1], [1, 2]])
        return C


class QuadraticElement(Element):
    """Class to represent a linear element (order = 2).

    Methods
    =======
    K()
        Returns element stiffness matrix.
    C()
        Returns element damping matrix.

    """

    order = 2

    def K(self, mat: dict[str, float]) -> matrix:
        """Returns element stiffness matrix."""
        K = (mat["k"] / self.length() / 3) * np.array(
            [[7, -8, 1], [-8, 16, -8], [1, -8, 7]]
        )
        return K

    def C(self, mat: dict[str, float]) -> matrix:
        """Returns element damping matrix."""
        C = (self.length() * mat["rho"] * mat["cp"] / 30) * np.array(
            [[4, 2, -1], [2, 16, 2], [-1, 2, 4]]
        )
        return C
