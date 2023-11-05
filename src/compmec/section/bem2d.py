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


def IntegrateUVn(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    """
    Computes the integral

    2*pi * int_{Gamma} u * (dv/dn) ds

    with

    * u the objective function
    * v the green function with source at S_i
    S_i = closedcurve(source)

    It in fact returns the vector like

    B = [B_0, B_1, ..., B_{n-1}]

    (n) is the number of degree of freedom

    The integral is transformed to the equivalent

    B_{ij} = 2*pi * sum_k int_{a_k}^{b_k} phi_j (r x p')/(r^2) dt

    with r(t) = p(t) - S_i
    """
    vertices = np.array(vertices, dtype="float64")
    vectors = np.roll(vertices, -1, axis=0) - vertices
    nverts = len(vertices)
    result = np.zeros((nverts, nverts), dtype="float64")

    gauss_nodes, gauss_weights = Integration.gauss(8)
    for i, Si in enumerate(vertices):
        # Source point at Si
        for j0, dVj in enumerate(vectors):
            # Integral over the segment j
            # between points V[j] and V[j+1]
            j1 = (j0 + 1) % nverts
            Vj0 = vertices[j0] - Si
            radius = tuple(Vj0 + t * dVj for t in gauss_nodes)
            funcs = tuple(np.cross(rad, dVj) / np.inner(rad, rad) for rad in radius)
            phi1 = gauss_nodes
            phi2 = 1 - gauss_nodes
            if j0 != i:
                result[i, j0] += np.einsum("k,k,k", phi2, funcs, gauss_weights)
            if j1 != i:
                result[i, j1] += np.einsum("k,k,k", phi1, funcs, gauss_weights)
        vertices += Si
    return result


def IntegrateUnV(vertices: Tuple[Tuple[float]]) -> Tuple[float]:
    """
    Computes the integral

    2*pi * int_{Gamma} (du/dn) * v ds

    with

    * u the objective function
    * v the green function with source at S_i
    S_i = closedcurve(source)

    It in fact returns the vector like

    A = [A_0, A_1, ..., A_{n-1}]

    (n) is the number of degree of freedom

    The integral is transformed to the equivalent

    A_{ij} = 2*pi * sum_k int_{a_k}^{b_k} phi_j ln(r) * abs(p') dt

    with r(t) = p(t) - S_i
    """
    vertices = np.array(vertices, dtype="float64")
    nverts = len(vertices)
    vectors = np.roll(vertices, -1, axis=0) - vertices
    # normvectors = np.linalg.norm(vectors, axis=1)
    alphs = np.einsum("ij,ij->i", vertices, vectors)
    betas = np.einsum("ij,ij->i", vectors, vectors)
    result = np.zeros(nverts, dtype="float64")
    gauss_nodes, gauss_weights = Integration.gauss(5)
    for i, Si in enumerate(vertices):
        for j0, dVj in enumerate(vectors):
            alp, bet = alphs[j0], betas[j0]
            j1 = j0 % nverts
            Vj0 = vertices[j0] - Si
            if j0 == i:
                # Associated with basis function (t)
                result[i] -= alp + bet / 4
                absdV = np.linalg.norm(dVj)
                result[i] += (alp + bet / 2) * np.log(absdV)
            elif j1 == i:
                result[i] -= alp + 3 * bet / 4
                # Associated with basis function (1-t)
                absdV = np.linalg.norm(dVj)
                result[i] += (alp + bet / 2) * np.log(absdV)
            else:
                radius = np.array(
                    tuple(Vj0 + t * dVj for t in gauss_nodes), dtype="float64"
                )
                lnradius = np.log(np.linalg.norm(radius, axis=1))
                phis = alp + bet * gauss_nodes
                result[i] += np.einsum("i,i,i", phis, lnradius, gauss_weights)
    return result


def CornerAngles(vertices):
    """
    Returns the value of the angle of each vertice
    """
    vectors0 = np.roll(vertices, -1, axis=0) - vertices
    vectors1 = np.roll(vectors0, -1, axis=0)
    crosses = tuple(np.cross(v0, v1) for v0, v1 in zip(vectors0, vectors1))
    inners = tuple(np.inner(v0, v1) for v0, v1 in zip(vectors0, vectors1))
    angles = tuple(np.arctan2(cross, inner) for cross, inner in zip(crosses, inners))
    angles = np.array(angles, dtype="float64")
    return angles
