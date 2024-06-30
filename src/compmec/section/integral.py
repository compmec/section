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

from .abcs import ICurve


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


def winding_number_linear(
    pointa: Tuple[float], pointb: Tuple[float], center: Tuple[float]
) -> float:
    """
    Computes the winding number for a straight line AB with respect to C

    :param pointa: The pair (xa, ya)
    :type pointa: Tuple[float]
    :param pointb: The pair (xb, yb)
    :type pointb: Tuple[float]
    :param center: The pair (xc, yc)
    :type center: Tuple[float]
    :return: The winding number value in the interval [-0.5, 0.5]
    :rtype: float
    """
    x1 = float(pointa[0] - center[0])
    y1 = float(pointa[1] - center[1])
    x2 = float(pointb[0] - center[0])
    y2 = float(pointb[1] - center[1])
    cross = x1 * y2 - x2 * y1
    inner = x1 * x2 + y1 * y2
    wind = np.arctan2(cross, inner) / math.tau
    return wind


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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def opened(npts: int) -> Tuple[Tuple[float]]:
        """
        Open newton cotes formula

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
        if npts > 7:
            raise NotImplementedError
        nodes = tuple(num / (npts + 1) for num in range(1, npts + 1))
        nodes = np.array(nodes, dtype="float64")
        weights = (
            (1.0,),
            (0.5, 0.5),
            (2 / 3, -1 / 3, 2 / 3),
            (11 / 24, 1 / 24, 1 / 24, 11 / 24),
        )
        weights = np.array(weights[npts - 1], dtype="float64")
        return nodes, weights

    @staticmethod
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


class Polynomial:
    """
    Class responsible to compute the integral

    II(a, b) = int_D x^a * y^b dx dy

    Which D is the domain defined by the interior of a curve
    The (a, b) values are the expoents

    """

    expoents = [
        (0, 0),
        (0, 1),
        (1, 0),
        (0, 2),
        (1, 1),
        (2, 0),
        (0, 3),
        (2, 1),
        (1, 2),
        (3, 0),
    ]

    # pylint: disable=invalid-name
    @staticmethod
    def polygon(vertices: Tuple[Tuple[float]]) -> Tuple[float]:
        """
        Computes integrals, returning the values of the integral

        II_{a, b} = int_D x^a * y^b dx dy

        for (a, b) = [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0),
                      (0, 3), (2, 1), (1, 2), (3, 0)]

        and D being a polygonal domain given by the vertices

        Parameters
        ----------

        :param vertices: The vertices that defines the polygon
        :type vertices: tuple[tuple[float]]
        :return: A vector of lenght 10 containing the integral values
        :rtype: tuple[float]
        """
        vertices = np.array(vertices, dtype="float64")
        xvan = np.vander(vertices[:, 0], 4, True)
        yvan = np.vander(vertices[:, 1], 4, True)

        cross = vertices[:, 0] * np.roll(vertices[:, 1], shift=-1)
        cross -= vertices[:, 1] * np.roll(vertices[:, 0], shift=-1)

        expoents = Polynomial.expoents
        geomprops = np.zeros(len(expoents), dtype="float64")
        for k, (a, b) in enumerate(expoents):
            M = np.zeros((a + 1, b + 1), dtype="float64")
            for i in range(a + 1):
                for j in range(b + 1):
                    M[i, j] = comb(i + j, i) * comb(a + b - i - j, b - j)
            X = np.roll(xvan[:, : a + 1], shift=-1, axis=0) * xvan[:, a::-1]
            Y = np.roll(yvan[:, : b + 1], shift=-1, axis=0) * yvan[:, b::-1]
            geomprops[k] = np.einsum("k,ki,ij,kj", cross, X, M, Y)
            geomprops[k] /= (a + b + 2) * (a + b + 1) * comb(a + b, a)
        return geomprops

    @staticmethod
    def adaptative(curve: ICurve, tolerance: float = 1e-9) -> Tuple[float]:
        """
        Computes the polynomials integrals over the area defined by the curve.
        It's function is suitable for non-polygonal curves

        It uses an adaptative algorithm that allows computing the integrals
        over smooth curves using milne's (open newton quadrature 3 points)

        This function uses a recursive approach

        :param curve: The boundary curve around the area
        :type curve: Curve
        :return: The integral values
        :rtype: float

        """
        assert isinstance(curve, ICurve)
        expoents = Polynomial.expoents
        integrals = np.zeros(len(expoents), dtype="float64")
        integrator = AdaptativePolynomialIntegrator(curve, tolerance=tolerance)
        for k, (expx, expy) in enumerate(expoents):
            value = integrator.integrate(expx, expy)
            integrals[k] = value / (2 + expx + expy)
        return integrals


class AdaptativePolynomialIntegrator:
    """
    Adaptative Polynomial Integrator

    Receives a curve and computes the polynomials integrals recursivelly
    using open newton cotes quadrature with 3 points
    """

    integ_matrix = (
        np.array([[0, 4, -1], [-4, 0, 4], [1, -4, 0]], dtype="float64") / 3
    )

    def __init__(self, curve, *, max_depth: int = 9, tolerance=1e-9):
        self.curve = curve
        self.max_depth = max_depth

        mesh = set()
        for t0, t1 in zip(curve.knots, curve.knots[1:]):
            subtol = tolerance * (t1 - t0) / (curve.knots[-1] - curve.knots[0])
            mesh |= self.find_mesh(t0, t1, tolerance=subtol)
        self.mesh = tuple(sorted(mesh))

    def find_mesh(
        self, t0: float, t1: float, tolerance: float = 1e-9
    ) -> Tuple[float]:
        """
        Finds the curve's subdivisions such the area computed
        by the adaptative subdivision is lower than the tolerance
        """
        subknots = np.linspace(t0, t1, 5)
        xvals, yvals = np.transpose(self.curve.eval(subknots))
        area_mid = xvals[::2] @ self.integ_matrix @ yvals[::2]
        area_lef = xvals[:3] @ self.integ_matrix @ yvals[:3]
        area_rig = xvals[2:] @ self.integ_matrix @ yvals[2:]
        diff = abs(area_lef + area_rig - area_mid)
        mesh = {t0, t1}
        if diff > tolerance:
            mesh |= self.find_mesh(t0, subknots[2], tolerance / 2)
            mesh |= self.find_mesh(subknots[2], t1, tolerance / 2)
        return mesh

    def milne_formula(self, t0: float, t1: float, a: int, b: int) -> float:
        """
        Computes the integral using Milne formula

        int_{t0}^{t1} x^a * y^b * (x * dy - y * dx)

        with (x(t), y(t)) being the curve

        Parameters
        ----------

        :param t0: The lower interval's value
        :type t0: float
        :param t1: The upper interval's value
        :type t1: float
        :param a: The x expoent
        :type a: int
        :param b: The y expoent
        :type b: int
        :return: The interval's integral value
        :rtype: float

        """
        ts = np.linspace(t0, t1, 5)
        points = self.curve.eval(ts)
        xs = points[:, 0]
        ys = points[:, 1]
        f1 = xs[1] ** a * ys[1] ** b
        f2 = xs[2] ** a * ys[2] ** b
        f3 = xs[3] ** a * ys[3] ** b
        f1 *= xs[1] * (ys[2] - ys[0]) - ys[1] * (xs[2] - xs[0])
        f2 *= xs[2] * (ys[3] - ys[1]) - ys[2] * (xs[3] - xs[1])
        f3 *= xs[3] * (ys[4] - ys[2]) - ys[3] * (xs[4] - xs[2])
        return (4 * f1 - 2 * f2 + 4 * f3) / 3

    def integrate(self, expx: int, expy: int) -> float:
        """
        Computes the integral over the entire curve

        int_{Gamma} x^expx * y^expy * (x * dy - y * dx)

        with (x(t), y(t)) being the curve

        Parameters
        ----------

        :param tolerance: The tolerance to how when stop adaptative quadrature
        :type tolerance: float
        :return: The interval's integral value
        :rtype: float
        """
        value = 0
        for t0, t1 in zip(self.mesh, self.mesh[1:]):
            tm = 0.5 * (t0 + t1)
            value += self.milne_formula(t0, tm, expx, expy)
            value += self.milne_formula(tm, t1, expx, expy)
        return value
