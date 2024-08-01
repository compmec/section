"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

"""

import math
from typing import Optional, Tuple, Union

import numpy as np

from .abcs import IBasisFunction, ICurve, IModel, IPoissonEvaluator, ISection
from .basisfunc import SplineBasisFunction, distributed_knots
from .integral import Integration


class ScalarFunction:

    def __init__(self, basisfunc: IBasisFunction, ctrlpoints: Tuple[float]):
        self.basis = basisfunc
        self.ctrlpoints = ctrlpoints

    @property
    def knots(self) -> Tuple[float]:
        return self.basis.knots

    def eval(self, parameters: Tuple[float]) -> Tuple[float]:
        matrix = self.basis.eval(parameters)
        return np.dot(matrix, self.ctrlpoints)

    def deval(self, parameters: Tuple[float]) -> Tuple[float]:
        matrix = self.basis.deval(parameters)
        return np.dot(matrix, self.ctrlpoints)


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
    def incurve(curve: ICurve, basis: IBasisFunction, tsources: Tuple[float]):
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
        if not isinstance(basis, IBasisFunction):
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
                        "ji,j,j,j->i", phis, weights, rcrossdp, over_rinr
                    )
        matrix[np.abs(matrix) < 1e-9] = 0
        return matrix / math.tau

    @staticmethod
    def outcurve(
        curve: ICurve, basis: IBasisFunction, sources: Tuple[Tuple[float]]
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
                        "ji,j,j,j->i", phis, weights, rcrossdp, over_rinr
                    )
        matrix[np.abs(matrix) < 1e-9] = 0
        return matrix


class TorsionEvaluator:

    @staticmethod
    def torsion_center_matrix(
        curve: ICurve, basis: IBasisFunction
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
    def constant_torsion(curve: ICurve, warping: ScalarFunction) -> float:
        """
        Returns the value of

        J = int_{tmin}^{tmax} w * <p, p'> dt

        Where w(t) is the warping function on the boundary of the curve.

        """
        if not isinstance(curve, ICurve):
            raise TypeError
        if not isinstance(warping, ScalarFunction):
            raise TypeError
        result = 0
        tknots = tuple(sorted(set(curve.knots) | set(warping.knots)))
        nodes, weights = Integration.chebyshev(6)

        def direct_integral(ta: float, tb: float):
            tvalues = ta + (tb - ta) * nodes
            points = curve.eval(tvalues)
            dpoints = curve.deval(tvalues)
            wvalues = warping.eval(tvalues)
            pinnerdp = np.einsum("ij,ij->i", points, dpoints)
            return (tb - ta) * np.einsum("i,i,i", weights, wvalues, pinnerdp)

        def adaptative_integral(ta: float, tb: float, tolerance: float = 1e-9):
            tm = (ta + tb) / 2
            integ_mid = direct_integral(ta, tb)
            integ_lef = direct_integral(ta, tm)
            integ_rig = direct_integral(tm, tb)
            if abs(integ_lef + integ_rig - integ_mid) > tolerance:
                integ_lef = adaptative_integral(ta, tm, tolerance / 2)
                integ_rig = adaptative_integral(tm, tb, tolerance / 2)
            return integ_lef + integ_rig

        for ta, tb in zip(tknots, tknots[1:]):
            result += adaptative_integral(ta, tb, 1e-6)

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


class BEMModel(IModel):
    """
    A BEM2D Model to solve laplace's equation
    """

    BASIS_DEGREE = 2
    NDOFS_BY_CURVE = 20

    def __check_model(self):
        if set(self.basis.keys()) ^ set(self.curves.keys()):
            msg = "There are curves that has no basis functions:\n"
            msg += f"    curves: {sorted(self.curves.keys())}\n"
            msg += f"     basis: {sorted(self.basis.keys())}"
            raise ValueError(msg)
        if self.sources is None:
            raise NotImplementedError
        total_ndofs = sum(basis.ndofs for basis in self.basis.values())
        if total_ndofs != len(self.sources):
            raise NotImplementedError

    def __init__(self, section: ISection):
        if not isinstance(section, ISection):
            raise TypeError
        self.section = section
        self.curves = {}
        self.basis = {}
        self.__sources = None
        for homosection in section:
            for curve in homosection.geometry.curves:
                if curve.label not in self.curves:
                    self.curves[curve.label] = curve
        self.__within_context = False

    def __enter__(self):
        self.__within_context = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__within_context = False

    @property
    def sources(self) -> Tuple[Tuple[float]]:
        return self.__sources

    @sources.setter
    def sources(self, points: Tuple[Tuple[float]]):
        points = np.array(points, dtype="float64")
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError(f"Invalid points: {points.shape} != (n, 2)")
        self.__sources = points

    def add_basis(self, curve: Union[int, ICurve], basis: IBasisFunction):
        if not self.__within_context:
            raise NotImplementedError
        curve_label = curve.label if isinstance(curve, ICurve) else curve
        if curve_label not in self.curves:
            raise ValueError
        if not isinstance(basis, IBasisFunction):
            raise TypeError
        self.basis[curve_label] = basis

    def generate_mesh(self, mesh_size: Optional[float] = None):
        if not self.__within_context:
            raise NotImplementedError

        for curve in self.curves.values():
            print("   -")
            self.__generate_mesh_curve(curve, mesh_size)

        sources = []
        for label, curve in self.curves.items():
            basis = self.basis[label]
            tsources = distributed_knots(basis)
            sources += list(curve.eval(tsources))
        self.sources = sources

    def __generate_mesh_curve(
        self, curve: ICurve, mesh_size: Optional[float] = None
    ):

        nodes, weights = Integration.chebyshev(3)

        def direct_lenght(curve: ICurve, ta: float, tb: float) -> float:
            tvals = ta + (tb - ta) * nodes
            absdpoints = np.linalg.norm(curve.deval(tvals), axis=1)
            return (tb - ta) * np.inner(weights, absdpoints)

        def adaptive_lenght(
            curve: ICurve, ta: float, tb: float, tolerance: float
        ) -> float:
            tm = (ta + tb) / 2
            len_mid = direct_lenght(curve, ta, tb)
            len_lef = direct_lenght(curve, ta, tm)
            len_rig = direct_lenght(curve, tm, tb)
            if abs(len_lef + len_rig - len_mid) > tolerance:
                len_lef = adaptive_lenght(curve, ta, tm, tolerance / 2)
                len_rig = adaptive_lenght(curve, tm, tb, tolerance / 2)
            return len_lef + len_rig

        sublenghts = []
        for ta, tb in zip(curve.knots, curve.knots[1:]):
            sublenght = adaptive_lenght(curve, ta, tb, 1e-9)
            sublenghts.append(sublenght)
        if mesh_size is None:
            mesh_size = sum(sublenghts) / self.NDOFS_BY_CURVE

        basis_knots = set(curve.knots)
        for i, (ta, tb) in enumerate(zip(curve.knots, curve.knots[1:])):
            ndiv = int(np.ceil(sublenghts[i] / mesh_size))
            basis_knots |= set(
                ((ndiv - i) * ta + i * tb) / ndiv for i in range(ndiv)
            )
        basis_knots = tuple(sorted(basis_knots))
        basis = SplineBasisFunction.cyclic(basis_knots, degree=1)
        self.add_basis(curve, basis)

    def solve(self):
        """
        Solves the BEM problem, computing
        """
        if not self.__within_context:
            raise NotImplementedError
        self.__check_model()
        nsources = len(self.sources)
        matrix = np.zeros((nsources, nsources), dtype="float64")
        vector = np.zeros((nsources, 1), dtype="float64")

        all_labels = tuple(sorted(self.basis.keys()))
        all_basis = tuple(self.basis[key] for key in all_labels)
        all_curves = tuple(self.curves[key] for key in all_labels)
        all_winds = tuple(
            tuple(map(curve.winding, self.sources)) for curve in all_curves
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
                subsources = self.sources[mask]
                tsources = tuple(
                    curve.projection(source)[0] for source in subsources
                )
                submatrix = ComputeStiffness.incurve(curve, basis, tsources)
                matrix[mask, slicej] += submatrix
            if not all(mask):
                subsources = self.sources[~mask]
                submatrix = ComputeStiffness.outcurve(curve, basis, subsources)
                matrix[~mask, slicej] += submatrix
            vector[:, 0] += TorsionEvaluator.warping_source(
                curve, self.sources
            )
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
                warping = ScalarFunction(base, ctrlpoints[:, 0])
                ctrlpoints = [0 for _ in range(base.ndofs)]
                normal = ScalarFunction(base, ctrlpoints)
                poisson_evalcurve = PoissonEvaluatorCurve(
                    curve, warping, normal
                )
                evaluators.append(poisson_evalcurve)
            homosection.warping = PoissonEvaluatorGeometry(evaluators)


class PoissonEvaluatorCurve(IPoissonEvaluator):

    def __init__(
        self,
        curve: ICurve,
        boundary: ScalarFunction = None,
        normal: ScalarFunction = None,
    ):
        self.__curve = curve
        self.bound = boundary
        self.normal = normal

    @property
    def curve(self) -> ICurve:
        return self.__curve

    @property
    def bound(self) -> ScalarFunction:
        return self.__bound

    @property
    def normal(self) -> ScalarFunction:
        return self.__normal

    @bound.setter
    def bound(self, new_func: ScalarFunction):
        if new_func is None:
            pass
        elif not isinstance(new_func, ScalarFunction):
            raise TypeError
        self.__bound = new_func

    @normal.setter
    def normal(self, new_func: ScalarFunction):
        if new_func is None:
            pass
        elif not isinstance(new_func, ScalarFunction):
            raise TypeError
        self.__normal = new_func

    def eval(self, source: Tuple[float], tolerance: float = 1e-9) -> float:
        """
        Evaluates the integral of

        int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        wind = self.curve.winding(source)
        if abs(wind) < 1e-9:  # Outside
            return 0
        if self.bound is None:
            raise NotADirectoryError
        if wind < 1:  # At boundary
            param = self.curve.projection(source)[0]
            return self.bound.eval(param)
        if self.normal is None:
            raise NotADirectoryError

        tknots = set(self.curve.projection(source))
        tknots |= set(self.curve.knots)
        tknots |= set(self.bound.knots)
        tknots |= set(self.normal.knots)
        tknots = tuple(sorted(tknots))

        nodes, weights = Integration.chebyshev(5)

        def direct_integral(ta: float, tb: float):
            tvals = ta + (tb - ta) * nodes
            boundvals = self.bound.eval(tvals)
            normavals = self.normal.eval(tvals)
            points = self.curve.eval(tvals)
            dpoints = self.curve.deval(tvals)
            radius = tuple(point - source for point in points)
            rinnerr = np.einsum("ij,ij->i", radius, radius)
            rcrossdp = tuple(np.cross(r, dp) for r, dp in zip(radius, dpoints))
            logradius = tuple(np.log(rinr) / 2 for rinr in rinnerr)
            absdpts = np.linalg.norm(dpoints, axis=1)
            result = np.einsum(
                "i,i,i,i", weights, boundvals, rcrossdp, 1 / rinnerr
            )
            result -= np.einsum(
                "i,i,i,i", weights, normavals, logradius, absdpts
            )
            return (tb - ta) * result

        def adaptative_integral(ta: float, tb: float, tolerance: float = 1e-9):
            tm = (ta + tb) / 2
            midd = direct_integral(ta, tb)
            left = direct_integral(ta, tm)
            righ = direct_integral(tm, tb)
            if abs(left + righ - midd) > tolerance:
                left = adaptative_integral(ta, tm, tolerance / 2)
                righ = adaptative_integral(tm, tb, tolerance / 2)
            return left + righ

        tolerance /= tknots[-1] - tknots[0]
        result = 0
        for ta, tb in zip(tknots, tknots[1:]):
            tol = tolerance * (tb - ta)
            result += adaptative_integral(ta, tb, tol)
        return result / (2 * np.pi)

    def grad(
        self, source: Tuple[float], tolerance: float = 1e-9
    ) -> Tuple[float]:
        """
        Evaluates the integral of

        d/ds int_{ta}^{tb} u * (r x p')/<r, r> - du/dn * |p'| * ln |r| * dt
        """
        wind = self.curve.winding(source)
        if abs(wind) < 1e-9:  # Outside
            return (0, 0)
        if self.bound is None or self.normal is None:
            raise NotADirectoryError
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

        nodes, weights = Integration.chebyshev(5)

        def direct_integral(ta: float, tb: float):
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

        def adaptative_integral(ta: float, tb: float, tolerance: float):
            tm = (ta + tb) / 2
            middx, middy = direct_integral(ta, tb)
            leftx, lefty = direct_integral(ta, tm)
            righx, righy = direct_integral(tm, tb)
            diffx = abs(leftx + righx - middx)
            diffy = abs(lefty + righy - middy)
            if diffx > tolerance or diffy > tolerance:
                leftx, lefty = adaptative_integral(ta, tm, tolerance / 2)
                righx, righy = adaptative_integral(tm, tb, tolerance / 2)
            return (leftx + righx, lefty + righy)

        tolerance /= tknots[-1] - tknots[0]
        result = np.zeros(2, dtype="float64")
        for ta, tb in zip(tknots, tknots[1:]):
            tol = (tb - ta) * tolerance
            result += adaptative_integral(ta, tb, tol)
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
