"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

We can divide it in two parts:

* Numerical Integration
* BEM object



"""
from typing import Tuple

import numpy as np


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
