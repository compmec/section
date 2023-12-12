"""
Solves the poisson problem using boundary element method

nabla^2_u = f(x, y)


"""
import math
from typing import Tuple

import numpy as np

from compmec import nurbs


def comb(a, b):
    prod = 1
    for i in range(min(b, a - b)):
        prod *= a - i
        prod /= i + 1
    return prod


def IntegratePolygon(
    xverts: Tuple[float], yverts: Tuple[float], amax: int, bmax: int
) -> Tuple[Tuple[float]]:
    """
    Computes integrals, returning a matrix II such

    II_{a, b} = int_D x^a * y^b dx dy

    for a = [0, ..., amax] and b = [0, ..., bmax]

    and D being a polygonal domain given by the vertices
    """
    xvan = np.vander(xverts, amax + 1, True)
    yvan = np.vander(yverts, bmax + 1, True)

    cross = xverts * np.roll(yverts, shift=-1)
    cross -= yverts * np.roll(xverts, shift=-1)

    II = np.zeros((amax + 1, bmax + 1), dtype="float64")

    for a in range(amax + 1):
        for b in range(bmax + 1):
            M = np.zeros((a + 1, b + 1), dtype="float64")
            for i in range(a + 1):
                for j in range(b + 1):
                    M[i, j] = comb(i + j, i) * comb(a + b - i - j, b - j)
            X = np.roll(xvan[:, : a + 1], shift=-1, axis=0) * xvan[:, a::-1]
            Y = np.roll(yvan[:, : b + 1], shift=-1, axis=0) * yvan[:, b::-1]
            II[a, b] = np.einsum("k,ki,ij,kj", cross, X, M, Y)
            II[a, b] /= (a + b + 2) * (a + b + 1) * comb(a + b, a)
    return II


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


def IntegralUVn(curve, basis, tsources, tmesh):
    """Computes the matrix [M] which

    M_ij = int phi_j * (r x p')/<r, r> dt

    Which represents the integral

    int u * (dv/dn) ds

    Parameters
    ----------
    curve: A jordan curve of degree 1, that contains a knotvector and
        control points, which are the vertices
    basis: A basis function of the given basis
    tsources: The parameters which the source points are applied
        tuple[float]
    tmesh: The nodes separation to easier integrating
        tuple[float]

    The number of tsources must be equal to the number of dofs of basis
    """
    tknots = curve.knotvector.knots
    tknots = np.array(tknots, dtype="float64")
    uknots = sorted(set(tmesh) | set(tknots) | set(tsources))
    uknots = np.array(uknots, dtype="float64")
    tsources = np.array(tsources, dtype="float64")
    sources = curve(tsources)
    nbasis = basis.npts
    nsources = len(tsources)
    matrix = np.zeros((nsources, nbasis), dtype="float64")
    nodes, weights = Integration.gauss(10)

    if curve.degree != 1:
        raise NotImplementedError
    for k0, tk0 in enumerate(uknots[:-1]):
        tk1 = uknots[k0 + 1]
        vertex0 = curve(tk0)
        vertex1 = curve(tk1)
        vector = vertex1 - vertex0
        tvals = tk0 + nodes * (tk1 - tk0)
        phis = basis(tvals)
        points = vertex0 + np.tensordot(nodes, vector, axes=0)
        for i, (tsi, source) in enumerate(zip(tsources, sources)):
            if tk0 == tsi or tsi == tk1:
                continue
            radius = tuple(point - source for point in points)
            radius = np.array(radius, dtype="float64")
            rcrossdp = vector[1] * radius[:, 0] - vector[0] * radius[:, 1]  # r x p'
            rinnerr = np.einsum("ij,ij->i", radius, radius)  # <r, r>
            funcvals = rcrossdp / rinnerr
            matrix[i] += np.einsum("i,i,ji->j", weights, funcvals, phis)
    return matrix


def IntegralUnV(curve, tsources):
    """Computes the matrix [F] which

    F_i = int ln|r| * <p, p'> dt

    Which represents the integral

    int (du/dn) * v ds

    Parameters
    ----------
    curve: A jordan curve of degree 1, that contains a knotvector and
        control points, which are the vertices
    tsources: The parameters which the source points are applied

    """
    tknots = curve.knotvector.knots
    tknots = sorted(set(tknots) | set(tsources))
    tknots = np.array(tknots, dtype="float64")
    if np.any(tsources < tknots[0]) or np.any(tknots[-1] < tsources):
        raise ValueError
    if len(set(tsources)) != len(tsources):
        raise ValueError

    if curve.degree != 1:
        raise NotImplementedError
    nsources = len(tsources)
    sources = curve(tsources)
    result = np.zeros(nsources, dtype="float64")
    nodes, weights = Integration.gauss(10)
    for k0, tk0 in enumerate(tknots[:-1]):
        tk1 = tknots[k0 + 1]
        vertex0 = curve(tk0)
        vertex1 = curve(tk1)
        vector = vertex1 - vertex0
        p0innerv = np.inner(vertex0, vector)
        vinnerv = np.inner(vector, vector)
        pinnerv = p0innerv + nodes * vinnerv
        points = vertex0 + np.tensordot(nodes, vector, axes=0)
        for i, (tsi, source) in enumerate(zip(tsources, sources)):
            if tsi != tk0 and tsi != tk1:
                radius = points - source
                rinnerr = np.einsum("ij,ij->i", radius, radius)
                lograds = np.log(rinnerr) / 2
                result[i] += np.einsum("i,i,i", weights, lograds, pinnerv)
            else:
                result[i] += np.log(vinnerv) * (p0innerv + vinnerv / 2) / 2
                result[i] -= p0innerv
                if tsi == tk0:
                    result[i] -= vinnerv / 4
                else:
                    result[i] -= 3 * vinnerv / 4

    return result


def TorsionVector(vertices):
    """Computes the vector V

    V_j = int phi_j * <p, p'> dt

    to help computing

    J_w = int w * <p, p'> dt
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


def WeightedTorsionVector(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    """Returns a matrix A of shape (n, 3) such

    A_{0j} = int phi_j * (p x p') dt
    A_{1j} = int phi_j * (x^2) dy
    A_{2j} = - int phi_j * (y^2) dx

    B_0 = (1/2) * int <p, p> * <p, p'> dt
    B_1 = (1/3) * int x^3 * <p, p'> dt
    B_2 = (1/3) * int y^3 * <p, p'> dt
    """
    nverts = len(vertices)
    vectors = np.roll(vertices, -1, axis=0) - vertices
    amatrix = np.zeros((3, nverts), dtype="float64")
    bvector = np.zeros(3, dtype="float64")
    nodes, weights = Integration.gauss(5)
    phi1 = nodes  # Increase
    phi2 = 1 - nodes  # Decrease
    for j0, vector in enumerate(vectors):
        vertex = vertices[j0]
        j1 = (j0 + 1) % nverts
        amatrix[0, j0] += np.cross(vertex, vector) * np.einsum("i,i", weights, phi2)
        amatrix[0, j1] += np.cross(vertex, vector) * np.einsum("i,i", weights, phi1)
        xvals = vertex[0] + nodes * vector[0]
        yvals = vertex[1] + nodes * vector[1]
        xvals2 = xvals**2
        yvals2 = yvals**2
        xvals3 = xvals2 * xvals
        yvals3 = yvals2 * yvals
        amatrix[1, j0] += vector[1] * np.einsum("i,i,i", weights, phi2, xvals2)
        amatrix[1, j1] += vector[1] * np.einsum("i,i,i", weights, phi1, xvals2)
        amatrix[2, j0] -= vector[0] * np.einsum("i,i,i", weights, phi2, yvals2)
        amatrix[2, j1] -= vector[0] * np.einsum("i,i,i", weights, phi1, yvals2)

        pinnerdp = xvals * vector[0] + yvals * vector[1]
        bvector[0] += np.einsum("i,i,i", weights, pinnerdp, xvals2 + yvals2)
        bvector[1] += np.einsum("i,i,i", weights, pinnerdp, xvals3)
        bvector[2] += np.einsum("i,i,i", weights, pinnerdp, yvals3)

    amatrix /= 2
    bvector[0] *= 1 / 4
    bvector[1:] *= 1 / 6
    return amatrix, bvector


def ShearVector(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    """
    Computes the matrix [X] of shape (n, 6) such
                      [x^2   * dx]
                      [2*x*y * dx]
    X_i = int ln(r) * [y^2   * dx]  dt
                      [x^2   * dy]
                      [2*x*y * dy]
                      [y^2   * dy]
    """
    pass


def AreaProperties(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    """
    Receives n vertices, a matrix of shape (n, 3)
    Returns a matrix of shape C = (3, 3) such

            [1    x    y]
    C = int [x  x^2   xy] dx dy
            [y   xy  y^2]

    These values are computed suposing we use a polygon
    """
    nverts = len(vertices)
    area = 0
    mqx = 0
    mqy = 0
    ixx = 0
    ixy = 0
    iyy = 0
    for j0, (x0, y0) in enumerate(vertices):
        j1 = (j0 + 1) % nverts
        x1, y1 = vertices[j1]
        dx = x1 - x0
        dy = y1 - y0
        area += (x0 + x1) * dy / 2
        mqx += (x0**2 + x0 * x1 + x1**2) * dy / 6
        mqy -= dx * (y0**2 + y0 * y1 + y1**2) / 6
        ixx += dy * (x1**3 + x1**2 * x0 + x1 * x0**2 + x0**3) / 12
        val = 2 * x0 * y0 + x1 * y0 + x0 * y1 + 2 * x1 * y1
        ixy += val * (x0 * y1 - x1 * y0) / 24
        iyy -= dx * (y1**3 + y1**2 * y0 + y1 * y0**2 + y0**3) / 12

    matrix = np.array([[area, mqx, mqy], [mqx, ixx, ixy], [mqy, ixy, iyy]])
    return matrix


def WarpingValues(vertices: Tuple[Tuple[float]]) -> Tuple[float]:
    """Finds the vector W that satisfies the laplace's equation

    nabla^2 w = 0
    <grad w, n> = <p, p'>/|p'|
    w = sum_j phi_j * W_j

    (M-A)*W = F

    """
    matrix = IntegralUVn(vertices)
    matrix -= np.pi * np.eye(len(matrix))
    vector = IntegralUnV(vertices)
    matrix = np.pad(matrix, ((0, 1), (0, 1)), constant_values=1)
    matrix[-1, -1] = 0
    vector = np.pad(vector, (0, 1), constant_values=0)
    warpvalues = np.linalg.solve(matrix, vector)[:-1]
    return warpvalues


def ShearValues(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    raise NotImplementedError


class WarpingEvaluator:
    def interior(
        vertices: Tuple[Tuple[float]],
        sources: Tuple[Tuple[float]],
        warpvals: Tuple[float],
    ) -> Tuple[float]:
        """
        Computes w(s), when s is a interior point by using the
        equation

        alpha(s)*w(s) = int w * dv/dn * ds - int v * dw/dn * ds

        Rewrite it as

        alpha(s)*w(s) = int w * (r x p')/<r, r> dt - int ln|r| * <p, p'> dt

        When w is a linear combination
        """
        vertices = tuple(tuple(map(np.float64, point)) for point in vertices)
        vertices = np.array(vertices, dtype="float64")
        vectors = np.roll(vertices, -1, axis=0) - vertices
        nverts = vectors.shape[0]
        nodes, weights = Integration.gauss(5)
        avals = np.einsum("ij,ij->i", vertices, vectors)
        bvals = np.einsum("ij,ij->i", vectors, vectors)

        phi1 = nodes  # Increase
        phi2 = 1 - nodes  # Decrease
        result = np.zeros(len(sources))
        for k0, vector in enumerate(vectors):
            k1 = (k0 + 1) % nverts
            wvals = warpvals[k0] * phi2 + warpvals[k1] * phi1
            pinnerdp = avals[k0] + nodes * bvals[k0]
            tempmat = np.tensordot(nodes, vector, axes=0)
            for i, source in enumerate(sources):
                radius = vertices[k0] - source + tempmat
                rcrossdp = radius[:, 0] * vector[1] - radius[:, 1] * vector[0]
                rinnerr = np.einsum("ij,ij->i", radius, radius)
                funcvals = rcrossdp / rinnerr
                lnradius = np.log(rinnerr) / 2
                result[i] += np.einsum("i,i,i", weights, funcvals, wvals)
                result[i] -= np.einsum("i,i,i", weights, lnradius, pinnerdp)
        return result / (2 * np.pi)
