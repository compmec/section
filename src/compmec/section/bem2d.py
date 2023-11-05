"""
Solves the poisson problem using boundary element method

nabla^2_u = f(x, y)


"""
import math
from typing import Tuple, Union

import numpy as np


class Integration:
    """
    Nodes positions and weights to perform numerical integration
    On the interval [0, 1]
    """

    def log(npts: int) -> Tuple[Tuple[float]]:
        """Return nodes and weights to integrate:

        .. math::
            I = \\int_{0}^{1} f(x) \\cdot ln(x) dx

        .. math::
            I = \\sum_{i} f(x_i) \cdot w_i

        Values extracted from BEM book
        """
        if not isinstance(npts, int) or npts <= 0:
            raise ValueError("npts must be a positive integer")
        if npts > 8:
            raise NotImplementedError

        all_nodes = (
            ("1.78b56362cef38p-2",),
            ("1.cac9bef59e5f4p-4", "1.345da38f030f6p-1"),
            ("1.05b25a2d35842p-4", "1.79da5dc3ea182p-2", "1.88a48902ba847p-1"),
            (
                "1.538bc35d9baf2p-5",
                "1.f6bb8d47f68ddp-3",
                "1.1cc1b7e469ad6p-1",
                "1.b2add206cc41dp-1",
            ),
            (
                "1.dd56d5450d669p-6",
                "1.644e2a4bb324cp-3",
                "1.a595587137bcbp-2",
                "1.5ac8ec69e6a0ep-1",
                "1.ca1f78ca0d19ep-1",
            ),
            (
                "1.627398e53ec00p-6",
                "1.09630458ec1b0p-3",
                "1.418e93aaa2f0dp-2",
                "1.13cae0f88f8a8p-1",
                "1.838a6837c0e52p-1",
                "1.d8680d3b5cbe6p-1",
            ),
            (
                "1.11ee0f2c13284p-6",
                "1.9a5c4c22ce7e7p-4",
                "1.f8691e253f8d6p-3",
                "1.bbddda9e320bdp-2",
                "1.43c3823a84528p-1",
                "1.9f4af0ce0ce4bp-1",
                "1.e1b6d9d554160p-1",
            ),
            (
                "1.b47a4e85dbb85p-7",
                "1.46a862c74e6e4p-4",
                "1.953d67fe41326p-3",
                "1.6aa7583df4f58p-2",
                "1.0f1531c27102cp-1",
                "1.67543bebe45fcp-1",
                "1.b2e1d8a662f24p-1",
                "1.e81a678acec40p-1",
            ),
        )
        all_weights = (
            ("1.0000000000000p+0",),
            ("1.6fe462b840500p-1", "1.20373a8f7f600p-2"),
            ("1.06dcf622e93a3p-1", "1.91633746949acp-2", "1.838b71ce63c35p-4"),
            (
                "1.88d7b2183cef9p-2",
                "1.8c28724ee498fp-2",
                "1.85983ff68c584p-3",
                "1.419ddcecc25b0p-5",
            ),
            (
                "1.310afc7bff6b4p-2",
                "1.662bbd372b9b9p-2",
                "1.e03b658845f82p-3",
                "1.95381b033c2dcp-4",
                "1.35d8cc7e2f1abp-6",
            ),
            (
                "1.e8fcec51fbfacp-3",
                "1.3baf79b807afbp-2",
                "1.f668fba1d2b19p-3",
                "1.22d57ca996d0fp-3",
                "1.c648c5a98474fp-5",
                "1.4d376882a0614p-7",
            ),
            (
                "1.91c141c07636ap-3",
                "1.14ca376441e9dp-2",
                "1.eade547018af3p-3",
                "1.53823fda3d2f7p-3",
                "1.6c4fbbbc1c6b5p-4",
                "1.0fed8073ff8a0p-5",
                "1.84cfa6343fde8p-8",
            ),
            (
                "1.50b9a721cf1d1p-3",
                "1.e673d3b819b4bp-3",
                "1.d09287c3f569fp-3",
                "1.67f1c12bc1e40p-3",
                "1.ce896d8d7a3e3p-4",
                "1.da16d28c29bd2p-5",
                "1.57b89ce7dc83bp-6",
                "1.e32f4be730620p-9",
            ),
        )
        all_nodes = tuple(tuple(map(np.float64.fromhex, nodes)) for nodes in all_nodes)
        all_weights = tuple(
            tuple(map(np.float64.fromhex, weights)) for weights in all_weights
        )
        # https://dlmf.nist.gov/3.5#v

        nodes = np.array(all_nodes[npts - 1], dtype="float64")
        weights = -np.array(all_weights[npts - 1], dtype="float64")
        return nodes, weights

    def gauss(npts: int) -> Tuple[Tuple[float]]:
        """
        Returns nodes and weights to gaussian integration
        """
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        nodes, weights = np.polynomial.legendre.leggauss(npts)
        return (1 + nodes) / 2, weights / 2

    def closed(npts: int) -> Tuple[Tuple[float]]:
        """Closed newton cotes formula"""
        if not isinstance(npts, int) or npts < 2:
            raise ValueError(f"npts invalid: {npts}")
        if npts > 7:
            raise NotImplementedError
        nodes = np.linspace(0, 1, npts)
        weights = (
            (0.5, 0.5),
            (1 / 6, 2 / 3, 1 / 6),
            (1 / 8, 3 / 8, 3 / 8, 1 / 8),
            (7 / 90, 16 / 45, 2 / 15, 16 / 45, 7 / 90),
            (19 / 288, 25 / 96, 25 / 144, 25 / 144, 25 / 96, 19 / 288),
            (41 / 840, 9 / 35, 9 / 280, 34 / 105, 9 / 280, 9 / 35, 41 / 840),
        )
        weights = np.array(weights[npts - 2], dtype="float64")
        return nodes, weights

    def chebyshev(npts: int) -> Tuple[Tuple[float]]:
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        if npts > 6:
            raise NotImplementedError
        nums = range(1, 2 * npts, 2)
        nums = tuple(np.float64(num) / (2 * npts) for num in nums)
        nodes = tuple(math.sin(0.5 * math.pi * num) ** 2 for num in nums)
        root2 = np.sqrt(2)
        root3 = np.sqrt(3)
        root5 = np.sqrt(5)
        all_weights = (
            (1,),
            (0.5, 0.5),
            (2 / 9, 5 / 9, 2 / 9),
            (
                (3 - root2) / 12,
                (3 + root2) / 12,
                (3 + root2) / 12,
                (3 - root2) / 12,
            ),
            (
                (13 - 3 * root5) / 75,
                (13 + 3 * root5) / 75,
                23 / 75,
                (13 + 3 * root5) / 75,
                (13 - 3 * root5) / 75,
            ),
            (
                (14 - 5 * root3) / 90,
                17 / 90,
                (14 + 5 * root3) / 90,
                (14 + 5 * root3) / 90,
                17 / 90,
                (14 - 5 * root3) / 90,
            ),
        )
        nodes = np.array(nodes, dtype="float64")
        weights = np.array(all_weights[npts - 1], dtype="float64")
        return nodes, weights


def IntegralUVn(vertices):
    """
    Given vertices, this functions computes the matrix [H]

    H_ij = int_{Gamma} phi_j * dv/dn * dGamma

    phi_j are the basis function

    U(t) = sum_j U_j * phi_j
    v = ln r -> dv/dn = (r x p')/(<r, r> * abs(p'))
    H_{ij} = int_{tmin}^{tmax} phi_j * (r x p')/<r, r> dt

    For j != i, use standard gaussian integration

    For j == i,
        p(t) = p_{i} + z * vec
        r = (z-1/2) * vec
        r x p' = 0
        H_{ii} = 0

    """
    nverts = len(vertices)
    vectors = np.roll(vertices, -1, axis=0) - vertices
    sources = 0.5 * (np.roll(vertices, -1, axis=0) + vertices)
    matrix = np.zeros((nverts, nverts), dtype="float64")
    nodes, weights = Integration.gauss(5)
    phi1 = nodes  # Increase
    phi2 = 1 - nodes  # Decrease
    for i, source in enumerate(sources):
        for j, vector in enumerate(vectors):
            if j == i:
                continue
            j1 = (j + 1) % nverts
            vect0 = vertices[j] - source
            numera = np.cross(vect0, vector)  # r x p'
            radius = vect0 + np.tensordot(nodes, vector, axes=0)
            radsqu = np.einsum("ij,ij->i", radius, radius)  # <r, r>
            funcvals = numera / radsqu
            # Increase
            val = np.einsum("i,i,i", phi1, weights, funcvals)
            matrix[i, j1] += val
            # Decrease
            val = np.einsum("i,i,i", phi2, weights, funcvals)
            matrix[i, j] += val
    return matrix


def IntegralUnV(vertices):
    """
    Given vertices, this functions computes the matrix [H]

    G_i = int_{Gamma} du/dn * v * dGamma

    du/dn = <p, p'>/abs(p')
    v = ln r

    G_i = int_{tmin}^{tmax} <p, p'> * ln(r) dt
    G_i = sum_j int_{t_j}^{t_{j+1}} <p, p'> ln(r) dt

    For j != i, then the integral is computed by gaussian quadrature

    for j = i, then
        vec = p_{i+1} - p_{i}
        p(z) = p_{i} + z * vec
        r(z) = p(z) - (p_i + p_{i+1})/2
        r(z) = (z-1/2) * vec
        a = <p_{i}, vec>
        b = <vec, vec>
        <p, p'> = a + z * b
        ln |r| = ln |vec| + ln |z-1/2|

        int_{t_{i-1}}^{t_i} <p, p'> ln |r| dt =
        int_{0}^{1} (a + z * b) * (ln |vec| + ln |z-1/2|) dz =
        (ln |vec| - 1 - ln(2)) * (a + b/2)

    """
    vectors = np.roll(vertices, -1, axis=0) - vertices
    sources = 0.5 * (np.roll(vertices, -1, axis=0) + vertices)
    avals = np.einsum("ij,ij->i", vertices, vectors)
    bvals = np.einsum("ij,ij->i", vectors, vectors)
    result = (avals + bvals / 2) * (0.5 * np.log(bvals) - 1 - np.log(2))
    nodes, weights = Integration.gauss(5)
    for i, source in enumerate(sources):
        for j, vector in enumerate(vectors):
            if j == i:
                continue
            points = vertices[j] + np.tensordot(nodes, vector, axes=0)
            radius = points - source
            lograds = np.log(np.einsum("ij,ij->i", radius, radius)) / 2  # ln(r)
            funcvals = np.einsum("ij,j->i", points, vector)
            val = np.einsum("i,i,i", lograds, weights, funcvals)
            result[i] += val
    return result


def TorsionVector(vertices):
    """
    Computes the integral

    I = int w * <p, p'> dt

    In fact returns the vector

    V = [V_0, ..., V_{n-1}]

    with

    V_j = int phi_j * <p, p'> dt
    """
    nverts = len(vertices)
    vectors = np.roll(vertices, -1, axis=0) - vertices
    result = np.zeros(len(vectors), dtype="float64")
    nodes, weights = Integration.gauss(15)
    phi1 = nodes  # Increase
    phi2 = 1 - nodes  # Decrease
    for j, vector in enumerate(vectors):
        j1 = (j + 1) % nverts
        points = vertices[j] + np.tensordot(nodes, vector, axes=0)
        prods = np.einsum("ij,j->i", points, vector)
        result[j] += np.einsum("i,i,i", prods, phi2, weights)
        result[j1] += np.einsum("i,i,i", prods, phi1, weights)
    return result


class WarpingFunction:
    def __init__(self, all_vertices: Tuple[Tuple[Tuple[float]]]):
        self.__all_vertices = []
        for vertices in all_vertices:
            vertices = np.array(vertices, dtype="float64")
            self.__all_vertices.append(vertices)
        if len(self.__all_vertices) != 1:
            msg = "For now, only 1 group of vertices is allowed"
            raise ValueError(msg)
        self.__all_vertices = tuple(self.__all_vertices)
        self.__solution = None

    def __mount_matrix(self):
        vertices = self.__all_vertices[0]
        matrix = IntegralUVn(vertices)
        matrix -= np.pi * np.eye(len(matrix))
        return matrix

    def __warping_vector(self):
        vertices = self.__all_vertices[0]
        vector = IntegralUnV(vertices)
        return vector

    def __torsion_vector(self):
        vertices = self.__all_vertices[0]
        vector = TorsionVector(vertices)
        return vector

    def torsion_contribuition(self):
        vector = self.__torsion_vector()
        return np.inner(vector, self.__solution)

    def solve(self):
        matrix = self.__mount_matrix()
        matrix = np.pad(matrix, ((0, 1), (0, 1)), constant_values=1)
        matrix[-1, -1] = 0
        vector = self.__warping_vector()
        vector = np.pad(vector, (0, 1))
        solution = np.linalg.solve(matrix, vector)
        self.__solution = solution[:-1]

    def eval_boundary(self, param: float) -> float:
        pass

    def eval_interior(self, point: Tuple[float]) -> float:
        pass

    def eval(self, points: Tuple[Tuple[float]]) -> float:
        """
        Evaluates the warping function at given point
        """
        raise NotImplementedError
        points = np.array(points, dtype="float64")
        values = np.zeros(len(points))
        for i, point in enumerate(points):
            param = projection(point)
            projected_pt = curve(param)
            if np.linalg.norm(projected_pt - point) < 1e-9:
                values[i] = self.eval_boundary(param)
            elif point in shape:
                values[i] = self.eval_interior(point)
        return values

    def __call__(self, points: Union[Tuple[float], Tuple[Tuple[float]]]):
        points = np.array(points, dtype="float64")
        if points.ndim == 2:
            return self.eval(points)
        return self.eval([points])[0]
