"""
This file contains integration methods to be used to compute geometric
properties and matrices M and F in BEM

The most useful are:

IntegratePolygon: Computes the polynomial integrals for polygons
Integration: Contains quadratures: nodes and weights for 1D integrals
Integration.log: Logarithm quadrature
Integration.gauss: Gauss-Legendre quadrature
Integration.closed: Closed Newton Cotes formulas
Integration.chebyshev: Polynomial integration at chebyshev nodes
"""

import math
from typing import Tuple

import numpy as np


# pylint: disable=invalid-name
def comb(a, b):
    """
    Computes the binomial coeffient of a, b:

    comb(a, b) = a! /((a-b)! * b!)

    Parameters
    ----------
    a: int
    b: int
    return: int

    Example
    -------
    >>> comb(1, 0)
    1
    >>> comb(5, 2)
    10

    """
    prod = 1
    for i in range(min(b, a - b)):
        prod *= a - i
        prod //= i + 1
    return prod


# pylint: disable=invalid-name
def integrate_polygon(
    xverts: Tuple[float], yverts: Tuple[float], amax: int = 3, bmax: int = 3
) -> Tuple[Tuple[float]]:
    """
    Computes integrals, returning a matrix II such

    II_{a, b} = int_D x^a * y^b dx dy

    for a = [0, ..., amax] and b = [0, ..., bmax]

    and D being a polygonal domain given by the vertices

    Parameters
    ----------
    xverts: tuple[float]
        Abssice values of the vertices
    yverts: tuple[float]
        Ordoned values of the vertices
    amax: int, default 3
        The maximum expoent multiplying 'x'
    bmax: int, default 3
        The maximum expoent multiplying 'y'
    return: tuple[tuple[float]]
        A matrix II of shape (amax+1, bmax+1) such
        II_{ij} = int_D x^a * y^b dx dy
    """
    xvan = np.vander(xverts, amax + 1, True)
    yvan = np.vander(yverts, bmax + 1, True)

    cross = xverts * np.roll(yverts, shift=-1)
    cross -= yverts * np.roll(xverts, shift=-1)

    geomprops = np.zeros((amax + 1, bmax + 1), dtype="float64")
    for a in range(amax + 1):
        for b in range(bmax + 1):
            M = np.zeros((a + 1, b + 1), dtype="float64")
            for i in range(a + 1):
                for j in range(b + 1):
                    M[i, j] = comb(i + j, i) * comb(a + b - i - j, b - j)
            X = np.roll(xvan[:, : a + 1], shift=-1, axis=0) * xvan[:, a::-1]
            Y = np.roll(yvan[:, : b + 1], shift=-1, axis=0) * yvan[:, b::-1]
            geomprops[a, b] = np.einsum("k,ki,ij,kj", cross, X, M, Y)
            geomprops[a, b] /= (a + b + 2) * (a + b + 1) * comb(a + b, a)
    return geomprops


class Integration:
    """
    This class contains nodes positions and weights to perform
    numerical integration on the interval [0, 1]

    It mostly gives nodes x_i and weights w_i such

    I = int_0^1 g(x) * f(x) dx = sum_i w_i * f(x_i)

    With f(x) being an arbitrary function
    and g(x) being a known function, like a weighted function

    Example
    -------
    >>> nodes, weights = Integration.log(4)
    >>> f = lambda x: 1 - 2*x + 4*x**2
    >>> I = sum(wi * f(xi) for xi, wi in zip(nodes, weights))
    """

    def log(npts: int) -> Tuple[Tuple[float]]:
        """Return nodes and weights to integrate:

        .. math::
            I = int_{0}^{1} f(x) * ln(x) dx

        .. math::
            I = sum_{i} f(x_i) * w_i

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
        all_nodes = tuple(
            tuple(map(np.float64.fromhex, nodes)) for nodes in all_nodes
        )
        all_weights = tuple(
            tuple(map(np.float64.fromhex, weights)) for weights in all_weights
        )
        # https://dlmf.nist.gov/3.5#v

        nodes = np.array(all_nodes[npts - 1], dtype="float64")
        weights = -np.array(all_weights[npts - 1], dtype="float64")
        return nodes, weights

    def gauss(npts: int) -> Tuple[Tuple[float]]:
        """
        Returns nodes and weights to perform gaussian integration
        in the interval [0, 1]

        Parameters
        ----------
        npts: int
            The number of points to use gauss integration.
            Must be at least 1
        return: tuple[tuple[float]]
            The pair (nodes, weights), with each node in the interval [0, 1]
        """
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        nodes, weights = np.polynomial.legendre.leggauss(npts)
        return (1 + nodes) / 2, weights / 2

    def closed(npts: int) -> Tuple[Tuple[float]]:
        """Closed newton cotes formula

        Parameters
        ----------
        npts: int
            The number of points to use gauss integration.
            Must be at least 2
        return: tuple[tuple[float]]
            The pair (nodes, weights), with each node in the interval [0, 1]
        """
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
        """Chebyshev integration

        Parameters
        ----------
        npts: int
            The number of points to use chebyshev integration.
            Must be at least 1
        return: tuple[tuple[float]]
            The pair (nodes, weights), with each node in the interval [0, 1]
        """
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        if npts > 6:
            raise NotImplementedError
        nums = range(1, 2 * npts, 2)
        nums = tuple(float(num) / (2 * npts) for num in nums)
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
