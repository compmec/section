"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

"""

import math
from typing import Tuple

import numpy as np

from .abcs import IBasisFunc, ICurve, ISection
from .integral import Integration


class ComputeStiffness:
    """
    This class is resposible to compute the matrix [M]
    made by the following integral

    M_ij = int phi_j * (r x p')/<r, r> * dt

    Which represents the integral

    int u * (dv/dn) ds
    """

    # pylint: disable=too-many-locals
    @staticmethod
    def incurve(curve: ICurve, basis: IBasisFunc, tsources: Tuple[float]):
        """
        Computes the integral when the sources are placed at the curve.
        The emplacement of these sources are given by parameter 'tsources'

        We suppose the curve is a polygon.

        Parameters
        ----------

        :param tsources: The parametric emplacement of sources
        :type tsources: Tuple[float]
        :return: The output matrix, integral of UVn
        :rtype: Tuple[Tuple[float]]
        """
        if not isinstance(curve, ICurve):
            raise NotImplementedError
        if not isinstance(basis, IBasisFunc):
            raise NotImplementedError
        if curve.degree != 1:
            raise NotImplementedError
        cknots = np.array(curve.knots, dtype="float64")
        tknots = set(curve.knots) | set(basis.knots) | set(tsources)
        tknots = np.array(sorted(tknots))
        matrix = np.zeros((len(tsources), basis.ndofs), dtype="float64")
        nodes, weights = Integration.gauss(10)

        vertices = curve.eval(cknots[:-1])
        vectors = np.roll(vertices, shift=-1, axis=0) - vertices
        sources = curve.eval(tsources)

        for i, vector in enumerate(vectors):
            tva, tvb = cknots[i], cknots[i + 1]
            vertex = vertices[i]
            mask = (tva <= tknots) * (tknots <= tvb)
            tmesh = tknots[mask]
            for tk0, tk1 in zip(tmesh, tmesh[1:]):
                tvals = tk0 + nodes * (tk1 - tk0)
                taus = (tvals - tva) / (tvb - tva)
                phis = basis.eval(tvals)
                points = vertex + np.tensordot(taus, vector, axes=0)
                for i, (tsi, source) in enumerate(zip(tsources, sources)):
                    if tsi in (tk0, tk1):
                        continue
                    radius = tuple(point - source for point in points)
                    radius = np.array(radius, dtype="float64")
                    rcrossdp = vector[1] * radius[:, 0]
                    rcrossdp -= vector[0] * radius[:, 1]
                    rcrossdp *= (tk1 - tk0) / (tvb - tva)
                    rinnerr = np.einsum("ij,ij->i", radius, radius)  # <r, r>
                    funcvals = rcrossdp / rinnerr
                    matrix[i] += np.einsum(
                        "ij,j,j->i", phis, weights, funcvals
                    )
        matrix[np.abs(matrix) < 1e-9] = 0
        return matrix / math.tau

    @staticmethod
    def outcurve(
        curve: ICurve, basis: IBasisFunc, sources: Tuple[Tuple[float]]
    ):
        """
        Computes the integral when the sources are placed outside (or not)
        of the curve.

        The emplacement of these sources can be any, including on the curve.

        Not yet implemented cases which source lies on the curve

        Parameters
        ----------

        :return: The output matrix, integral of UVn
        :rtype: Tuple[Tuple[float]]
        """
        matrix = np.zeros((len(sources), basis.ndofs), dtype="float64")

        tknots = set(curve.knots) | set(basis.knots)
        for source in sources:
            tknots |= set(curve.projection(source))
        tknots = np.array(sorted(tknots), dtype="float64")

        nodes, weights = Integration.gauss(10)
        vertices = curve.eval(curve.knots[:-1])
        vectors = np.roll(vertices, shift=-1, axis=0) - vertices

        for i, vector in enumerate(vectors):
            tva, tvb = curve.knots[i], curve.knots[i + 1]
            vertex = vertices[i]
            mask = (tva <= tknots) * (tknots <= tvb)
            tmesh = tknots[mask]
            for tk0, tk1 in zip(tmesh, tmesh[1:]):
                tvals = tk0 + nodes * (tk1 - tk0)
                taus = (tvals - tva) / (tvb - tva)
                phis = basis.eval(tvals)
                points = vertex + np.tensordot(taus, vector, axes=0)
                for i, source in enumerate(sources):
                    radius = tuple(point - source for point in points)
                    radius = np.array(radius, dtype="float64")
                    rcrossdp = vector[1] * radius[:, 0]
                    rcrossdp -= vector[0] * radius[:, 1]
                    rcrossdp *= (tk1 - tk0) / (tvb - tva)
                    rinnerr = np.einsum("ij,ij->i", radius, radius)  # <r, r>
                    funcvals = rcrossdp / rinnerr
                    matrix[i] += np.einsum(
                        "ij,j,j->i", phis, weights, funcvals
                    )
        matrix[np.abs(matrix) < 1e-9] = 0
        return matrix


class TorsionEvaluator:

    @staticmethod
    def torsion_center_matrix(
        curve: ICurve, basis: IBasisFunc
    ) -> Tuple[Tuple[float]]:
        """
        Computes the matrix used to compute the torsion center

                            [1]
        B = int_{Omega} w * [x] dOmega
                            [y]

        Where w is the warping function.

        This function returns a matrix [M], of shape (n, 3), such

        M_{i0} = int_{Omega} w*1 dOmega
        M_{i1} = int_{Omega} w*x dOmega
        M_{i2} = int_{Omega} w*y dOmega

        """
        raise NotImplementedError

    @staticmethod
    def constant_vector(curve: ICurve, basis: IBasisFunc) -> Tuple[float]:
        """
        Computes the vector used to compute the torsion constant

        It's interested to compute J:

        J = int_{tmin}^{tmax} w * <p, p'> dt

        Where w(t) is the warping function on the boundary of the curve.
        This function, on the boundary is defined as

        w(t) = sum_{i} F_{i}(t) * W_{i}

        This function returns a vector [V], of lenght 'n', such

        V_i = int_{tmin}^{tmax} F_{i}(t) * <p, p'> dt

        .. note:
            We suppose the curve is a polygon

        """
        if curve.degree != 1:
            raise NotImplementedError
        result = np.zeros(basis.ndofs, dtype="float64")
        vertices = curve.eval(curve.knots[:-1])
        vectors = np.roll(vertices, shift=-1, axis=0) - vertices
        alphas = np.einsum("ij,ij->i", vertices, vectors)
        betas = np.einsum("ij,ij->i", vectors, vectors)

        cknots = np.array(curve.knots, dtype="float64")
        tknots = np.array(sorted(set(curve.knots) | set(basis.knots)))
        nodes, weights = Integration.closed(2)

        for i, (alpha, beta) in enumerate(zip(alphas, betas)):
            tva, tvb = cknots[i], cknots[i + 1]
            mask = (tva <= tknots) * (tknots <= tvb)
            tmesh = tknots[mask]
            for tk0, tk1 in zip(tmesh, tmesh[1:]):
                diff = tk1 - tk0
                tvals = tk0 + nodes * diff
                zvals = (tvals - tva) / (tvb - tva)
                phis = basis.eval(tvals)
                result += diff * alpha * np.einsum("ij,j->i", phis, weights)
                result += (
                    diff * beta * np.einsum("ij,j,j->i", phis, weights, zvals)
                )
        return result

    @staticmethod
    def warping_source(
        curve: ICurve, sources: Tuple[Tuple[float]]
    ) -> Tuple[float]:
        """
        Computes the integral

        I = int_{tmin}^{tmax} <p, p'> * ln |r(t)| dt

        Where

        r(t) = p(t) - source

        We suppose the curve is a polygon

        I = sum_{k} dt_k * int_{0}^{1} (alpha + z * beta) * ln |r| dz
        """
        if curve.degree != 1:
            raise NotImplementedError
        sources = np.array(sources, dtype="float64")
        vertices = curve.eval(curve.knots[:-1])
        vectors = np.roll(vertices, shift=-1, axis=0) - vertices
        alphas = np.einsum("ij,ij->i", vertices, vectors)
        betas = np.einsum("ij,ij->i", vectors, vectors)

        nodes, weights = Integration.chebyshev(4)
        results = np.zeros(len(sources), dtype="float64")
        for k, source in enumerate(sources):
            for i, (alpha, beta) in enumerate(zip(alphas, betas)):
                tva, tvb = curve.knots[i], curve.knots[i + 1]
                projection = np.inner(vertices[i] - source, vectors[i]) / beta
                projection = min(1, max(0, projection))

                gamma = np.inner(vertices[i] - source, vertices[i] - source)
                delta = np.inner(vertices[i] - source, vectors[i])
                if projection != 0:
                    znodes = projection * nodes
                    logvals = np.log(
                        gamma + 2 * delta * znodes + beta * znodes**2
                    )
                    results[k] += (
                        (tvb - tva)
                        * alpha
                        * projection
                        * np.inner(weights, logvals)
                    )
                    results[k] += (
                        (tvb - tva)
                        * beta
                        * projection
                        * np.einsum("i,i,i", znodes, weights, logvals)
                    )
                if projection != 1:
                    znodes = projection + (1 - projection) * nodes
                    logvals = np.log(
                        gamma + 2 * delta * znodes + beta * znodes**2
                    )
                    results[k] += (
                        (tvb - tva)
                        * alpha
                        * (1 - projection)
                        * np.inner(weights, logvals)
                    )
                    results[k] += (
                        (tvb - tva)
                        * beta
                        * (1 - projection)
                        * np.einsum("i,i,i", znodes, weights, logvals)
                    )
        return results


class BEMModel:
    """
    A BEM2D Model to solve laplace's equation
    """

    def __check_model(self):
        if set(self.basis.keys()) ^ set(self.curves.keys()):
            raise NotImplementedError

    def __init__(self, section: ISection):
        self.curves = {}
        self.basis = {}
        self.solution = None
        for geometry in section.geometries:
            for curve in geometry.curves:
                if curve.label not in self.curves:
                    self.curves[curve.label] = curve

    def solve(self, sources: Tuple[Tuple[float]]):
        """
        Solves the BEM problem, computing
        """
        self.__check_model()
        sources = np.array(sources, dtype="float64")
        if sources.ndim != 2 or sources.shape[1] != 2:
            raise NotImplementedError
        total_ndofs = sum(base.ndofs for base in self.basis.values())
        nsources = len(set(tuple(map(float, source)) for source in sources))
        if nsources != total_ndofs:
            raise NotImplementedError
        matrix = np.zeros((total_ndofs, total_ndofs), dtype="float64")
        vector = np.zeros((total_ndofs, 1), dtype="float64")

        basis_labels = tuple(sorted(self.basis.keys()))
        index_base = 0
        for label in basis_labels:
            base = self.basis[label]
            curve = self.curves[label]
            slicej = slice(index_base, index_base + base.ndofs)
            winds = map(curve.winding, sources)
            mask = np.array(tuple(0 < wind < 1 for wind in winds))

            tsources = tuple(
                curve.projection(source)[0] for source in sources[mask]
            )
            submatrix = ComputeStiffness.incurve(curve, base, tsources)
            matrix[mask, slicej] += submatrix
            submatrix = ComputeStiffness.outcurve(curve, base, sources[~mask])
            matrix[~mask, slicej] += submatrix
            vector[:, 0] += TorsionEvaluator.warping_source(curve, sources)
            index_base += base.ndofs

        result = np.linalg.solve(matrix, vector)

        self.solution = {}
        index_base = 0
        for label in basis_labels:
            base = self.basis[label]
            curve = self.curves[label]
            slicej = slice(index_base, index_base + base.ndofs)
            self.solution[label] = result[slicej]
            index_base += base.ndofs
