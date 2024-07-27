"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

"""

import math
from typing import Tuple, Union

import numpy as np

from .abcs import IBasisFunc, ICurve, IPoissonEvaluator, ISection
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
                jacobin = (tk1 - tk0) / (tvb - tva)
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
                    over_rinr = 1 / np.einsum("ij,ij->i", radius, radius)
                    matrix[i] += jacobin * np.einsum(
                        "ij,j,j,j->i", phis, weights, rcrossdp, over_rinr
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
                jacobin = (tk1 - tk0) / (tvb - tva)
                tvals = tk0 + nodes * (tk1 - tk0)
                taus = (tvals - tva) / (tvb - tva)
                phis = basis.eval(tvals)
                points = vertex + np.tensordot(taus, vector, axes=0)
                for i, source in enumerate(sources):
                    radius = tuple(point - source for point in points)
                    radius = np.array(radius, dtype="float64")
                    rcrossdp = vector[1] * radius[:, 0]
                    rcrossdp -= vector[0] * radius[:, 1]
                    over_rinr = 1 / np.einsum("ij,ij->i", radius, radius)
                    matrix[i] += jacobin * np.einsum(
                        "ij,j,j,j->i", phis, weights, rcrossdp, over_rinr
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

    BASIS_DEGREE = 1
    NDOFS_BY_CURVE = 20

    def __check_model(self):
        if set(self.basis.keys()) ^ set(self.curves.keys()):
            raise NotImplementedError

    def __init__(self, section: ISection):
        if not isinstance(section, ISection):
            raise TypeError
        self.section = section
        self.curves = {}
        self.basis = {}
        for homosection in section:
            for curve in homosection.geometry.curves:
                if curve.label not in self.curves:
                    self.curves[curve.label] = curve

    def add_basis(self, curve: Union[int, ICurve], basis: IBasisFunc):
        curve_label = curve.label if isinstance(curve, ICurve) else curve
        if curve_label not in self.curves:
            raise ValueError
        if not isinstance(basis, IBasisFunc):
            raise TypeError
        self.basis[curve_label] = basis

    def solve(self, sources: Tuple[Tuple[float]]):
        """
        Solves the BEM problem, computing
        """
        self.__check_model()
        sources = np.array(sources, dtype="float64")
        if sources.ndim != 2 or sources.shape[1] != 2:
            raise NotImplementedError
        total_ndofs = sum(base.ndofs for base in self.basis.values())
        nsources = len(sources)
        if nsources != total_ndofs:
            msg = f"The number of sources ({nsources}) must be"
            msg += f" equal to the total ndofs ({total_ndofs})"
            raise ValueError(msg)
        matrix = np.zeros((total_ndofs, total_ndofs), dtype="float64")
        vector = np.zeros((total_ndofs, 1), dtype="float64")

        all_labels = tuple(sorted(self.basis.keys()))
        all_basis = tuple(self.basis[key] for key in all_labels)
        all_curves = tuple(self.curves[key] for key in all_labels)
        all_winds = tuple(
            tuple(map(curve.winding, sources)) for curve in all_curves
        )
        all_winds = np.array(all_winds)
        if not np.all(np.any((0 < all_winds) * (all_winds < 1), axis=0)):
            msg = str(all_winds)
            raise ValueError(msg)

        index_basis = 0
        for label, basis, curve, winds in zip(
            all_labels, all_basis, all_curves, all_winds
        ):
            mask = np.array(tuple(0 < wind < 1 for wind in winds))
            slicej = slice(index_basis, index_basis + basis.ndofs)
            if any(mask):
                subsources = sources[mask]
                tsources = tuple(
                    curve.projection(source)[0] for source in subsources
                )
                submatrix = ComputeStiffness.incurve(curve, basis, tsources)
                matrix[mask, slicej] += submatrix
            if not all(mask):
                subsources = sources[~mask]
                submatrix = ComputeStiffness.outcurve(curve, basis, subsources)
                matrix[~mask, slicej] += submatrix
            vector[:, 0] += TorsionEvaluator.warping_source(curve, sources)
            index_basis += basis.ndofs

        # Constraint solution
        matrix = np.pad(matrix, ((0, 1), (0, 1)), constant_values=1)
        matrix[-1, -1] = 0
        vector = np.pad(vector, ((0, 1), (0, 0)), constant_values=0)
        result = np.linalg.solve(matrix, vector)

        solution = {}
        index_base = 0
        for label in all_labels:
            base = self.basis[label]
            curve = self.curves[label]
            slicej = slice(index_base, index_base + base.ndofs)
            solution[label] = result[slicej]
            index_base += base.ndofs

        for homosection in self.section:
            geometry = homosection.geometry
            evaluators = []
            for curve in geometry.curves:
                base = self.basis[curve.label]
                ctrlpoints = solution[curve.label]
                warping = ScalarFunction(base, ctrlpoints)
                ctrlpoints = [0 for _ in range(base.ndofs)]
                normal = ScalarFunction(base, ctrlpoints)
                poisson_evalcurve = PoissonEvaluatorCurve(
                    curve, warping, normal
                )
                evaluators.append(poisson_evalcurve)
            homosection.warping = PoissonEvaluatorGeometry(evaluators)


class ScalarFunction:

    def __init__(self, basisfunc: IBasisFunc, ctrlpoints: Tuple[float]):
        self.basis = basisfunc
        self.ctrlpoints = ctrlpoints

    @property
    def knots(self) -> Tuple[float]:
        return self.basis.knots

    def eval(self, parameters: Tuple[float]) -> Tuple[float]:
        matrix = self.basis.eval(parameters)
        return np.dot(np.transpose(matrix), self.ctrlpoints)

    def deval(self, parameters: Tuple[float]) -> Tuple[float]:
        matrix = self.basis.deval(parameters)
        return np.dot(np.transpose(matrix), self.ctrlpoints)


class PoissonEvaluatorCurve(IPoissonEvaluator):

    def __init__(
        self,
        curve: ICurve,
        boundary: ScalarFunction,
        normal: ScalarFunction,
    ):
        self.curve = curve
        self.bound = boundary
        self.normal = normal

    def __integrate_eval(
        self, source: Tuple[float], ta: float, tb: float
    ) -> float:
        """
        Direct integral for

            int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        nodes, weights = Integration.chebyshev(5)
        tvals = ta + (tb - ta) * nodes
        boundvals = self.bound.eval(tvals)
        normavals = self.normal.eval(tvals)
        points = self.curve.eval(tvals)
        dpoints = self.curve.deval(tvals)
        radius = tuple(point - source for point in points)
        rinnerr = np.einsum("ij,ij->i", radius, radius)
        rcrossdp = tuple(np.cross(r, dp) for r, dp in zip(radius, dpoints))
        logradius = tuple(np.log(rinr) / 2 for rinr in rinnerr)
        absdpts = tuple(map(np.linalg.norm, dpoints))
        result = np.einsum(
            "i,i,i,i", weights, boundvals, rcrossdp, 1 / rinnerr
        )
        result -= np.einsum("i,i,i,i", weights, normavals, logradius, absdpts)
        return (tb - ta) * result

    def __integrate_grad(
        self, source: Tuple[float], ta: float, tb: float
    ) -> float:
        """
        Direct integral for

            int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        nodes, weights = Integration.chebyshev(5)
        tvals = ta + (tb - ta) * nodes
        boundvals = self.bound.eval(tvals)
        normavals = self.normal.eval(tvals)
        points = self.curve.eval(tvals)
        dpoints = self.curve.deval(tvals)
        radius = tuple(point - source for point in points)
        rinnerr = np.einsum("ij,ij->i", radius, radius)
        over_rinr = 1 / rinnerr
        rcrossdp = tuple(np.cross(r, dp) for r, dp in zip(radius, dpoints))
        absdpts = tuple(map(np.linalg.norm, dpoints))
        result = 2 * np.einsum(
            "i,i,i,i,i,ij->j",
            weights,
            boundvals,
            rcrossdp,
            over_rinr,
            over_rinr,
            radius,
        )
        result[0] -= np.einsum(
            "i,i,i,i", weights, boundvals, dpoints[:, 1], over_rinr
        )
        result[1] += np.einsum(
            "i,i,i,i", weights, boundvals, dpoints[:, 0], over_rinr
        )
        result += np.einsum(
            "i,i,i,i,ij->j", weights, normavals, over_rinr, absdpts, radius
        )
        return (tb - ta) * result

    def __eval_adapt(
        self,
        source: Tuple[float],
        ta: float,
        tb: float,
        tolerance: float = 1e-9,
    ) -> float:
        """
        Adapdative integral for

            int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        tm = (ta + tb) / 2
        midd = self.__integrate_eval(source, ta, tb)
        left = self.__integrate_eval(source, ta, tm)
        righ = self.__integrate_eval(source, tm, tb)
        if abs(left + righ - midd) > tolerance:
            left = self.__eval_adapt(source, ta, tm, tolerance / 2)
            righ = self.__eval_adapt(source, tm, tb, tolerance / 2)
        return left + righ

    def __grad_adapt(
        self,
        source: Tuple[float],
        ta: float,
        tb: float,
        tolerance: float = 1e-9,
    ) -> float:
        """
        Adapdative integral for

            int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        tm = (ta + tb) / 2
        middx, middy = self.__integrate_grad(source, ta, tb)
        leftx, lefty = self.__integrate_grad(source, ta, tm)
        righx, righy = self.__integrate_grad(source, tm, tb)
        diffx = abs(leftx + righx - middx)
        diffy = abs(lefty + righy - middy)
        if diffx > tolerance or diffy > tolerance:
            leftx, lefty = self.__grad_adapt(source, ta, tm, tolerance / 2)
            righx, righy = self.__grad_adapt(source, tm, tb, tolerance / 2)
        return (leftx + righx, lefty + righy)

    def eval(self, source: Tuple[float]) -> float:
        wind = self.curve.winding(source)
        if abs(wind) < 1e-9:  # Outside
            return 0
        if wind < 1:  # At boundary
            param = self.curve.projection(source)[0]
            return self.bound.eval(param)

        tknots = set(self.curve.projection(source))
        tknots |= set(self.curve.knots)
        tknots |= set(self.bound.knots)
        tknots |= set(self.normal.knots)
        tknots = sorted(tknots)
        tknots += [(ta + tb) / 2 for ta, tb in zip(tknots, tknots[1:])]
        tknots = sorted(tknots)

        result = 0
        for ta, tb in zip(tknots, tknots[1:]):
            result += self.__eval_adapt(source, ta, tb)
        return result / (2 * np.pi)

    def grad(self, source: Tuple[float]) -> Tuple[float]:
        wind = self.curve.winding(source)
        if abs(wind) < 1e-9:  # Outside
            return (0, 0)
        if wind < 1:  # At boundary
            param = self.curve.projection(source)[0]
            dpdt = self.curve.deval(param)
            norm = np.linalg.norm(dpdt)
            tx, ty = dpdt / norm  # tangent
            dudt = self.bound.deval(param) / norm
            dudn = self.normal.eval(param)
            matrix = ((tx, ty), (ty, -tx))
            vector = (dudt, dudn)
            return tuple(np.dot(matrix, vector))

        tknots = set(self.curve.projection(source))
        tknots |= set(self.curve.knots)
        tknots |= set(self.bound.knots)
        tknots |= set(self.normal.knots)
        tknots = tuple(sorted(tknots))

        result = np.zeros(2, dtype="float64")
        for ta, tb in zip(tknots, tknots[1:]):
            result += self.__grad_adapt(source, ta, tb)
        return result / (2 * np.pi)


class PoissonEvaluatorGeometry(IPoissonEvaluator):

    def __init__(self, evaluators: Tuple[PoissonEvaluatorCurve]):
        self.evaluators = evaluators

    def eval(self, source: Tuple[float]) -> float:
        wind_tolerance = 1e-6
        winds = np.zeros(len(self.evaluators), dtype="float64")
        for i, evaluator in enumerate(self.evaluators):
            winds[i] = evaluator.curve.winding(source)
            if wind_tolerance < winds[i] < 1 - wind_tolerance:
                return evaluator.eval(source)
        if np.any(np.abs(winds) < wind_tolerance):
            return 0
        result = 0
        for evaluator in self.evaluators:
            result += evaluator.eval(source)
        return result

    def grad(self, source: Tuple[float]) -> Tuple[float]:
        wind_tolerance = 1e-6
        winds = np.zeros(len(self.evaluators), dtype="float64")
        for i, evaluator in enumerate(self.evaluators):
            winds[i] = evaluator.curve.winding(source)
            if wind_tolerance < winds[i] < 1 - wind_tolerance:
                return evaluator.grad(source)
        if np.any(np.abs(winds) < wind_tolerance):
            return (0, 0)
        result = np.zeros(2, dtype="float64")
        for evaluator in self.evaluators:
            result += evaluator.grad(source)
        return result
