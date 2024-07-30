"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest

from compmec.section.basisfunc import SplineBasisFunction
from compmec.section.bem2d import (
    ComputeStiffness,
    PoissonEvaluatorCurve,
    PoissonEvaluatorGeometry,
    ScalarFunction,
    TorsionEvaluator,
)
from compmec.section.curve import NurbsCurve


@pytest.mark.order(8)
@pytest.mark.dependency(
    depends=[
        "tests/test_integral.py::test_end",
        "tests/test_basisfunc.py::test_end",
        "tests/test_material.py::test_end",
        "tests/test_curve.py::test_end",
        "tests/test_geometry.py::test_end",
        "tests/test_geomprop.py::test_end",
        "tests/test_axial.py::test_end",
        "tests/test_bending.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


class TestComputeMatrixIncurve:
    """
    Computes the integrals

    int_{tmin}^{tmax} phi_j(t) * (r x p')/<r, r> dt

    when

    r = p(t) - s = p(t) - p(u)

    s is a point on the curve
    """

    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_eqtriangle_corner(self):
        good_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        good_matrix = np.array(good_matrix) / 12

        angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0, 1, 2)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_eqtriangle_middle(self):
        pi, ln3, sq3 = np.pi, np.log(3), np.sqrt(3)
        good_matrix = [[-1, -1, 2], [2, -1, -1], [-1, 2, -1]]
        good_matrix = sq3 * ln3 * np.array(good_matrix) / (16 * pi)
        good_matrix += np.array([[3, 3, 2], [2, 3, 3], [3, 2, 3]]) / 16

        angles = np.linspace(0, 2 * pi, 3, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0.5, 1.5, 2.5)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_square_corner(self):
        a = 1 / 8 - np.log(2) / (4 * np.pi)
        b = np.log(2) / (2 * np.pi)
        good_matrix = [[0, a, b, a], [a, 0, a, b], [b, a, 0, a], [a, b, a, 0]]
        good_matrix = np.array(good_matrix)

        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        tsources = (0, 1, 2, 3)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_square_middle(self):
        a = np.arctan(2) / (2 * np.pi) - np.log(5) / (8 * np.pi)
        b = 1 / 4 + np.log(5) / (8 * np.pi) - np.arctan(2) / (2 * np.pi)
        good_matrix = [[a, a, b, b], [b, a, a, b], [b, b, a, a], [a, b, b, a]]
        good_matrix = np.array(good_matrix)

        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        tsources = (0.5, 1.5, 2.5, 3.5)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_hexagon_corner(self):
        pi, ln2, ln3, sq3 = np.pi, np.log(2), np.log(3), np.sqrt(3)
        a = 1 / 8 - sq3 * ln3 / (8 * pi)
        b = 1 / 24 + sq3 * (3 * ln3 - 4 * ln2) / (8 * pi)
        c = sq3 * (2 * ln2 - ln3) / (2 * pi)
        good_matrix = [
            [0, a, b, c, b, a],
            [a, 0, a, b, c, b],
            [b, a, 0, a, b, c],
            [c, b, a, 0, a, b],
            [b, c, b, a, 0, a],
            [a, b, c, b, a, 0],
        ]
        good_matrix = np.array(good_matrix)

        angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0, 1, 2, 3, 4, 5)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixIncurve::test_begin"])
    def test_hexagon_middle(self):
        pi, ln7, ln13, sq3 = np.pi, np.log(7), np.log(13), np.sqrt(3)
        a = -5 / 48 - sq3 * ln7 / (16 * pi) + 5 * np.arctan(5 / sq3) / (8 * pi)
        b = 1 / 48 + sq3 * (4 * ln7 - 3 * ln13) / (16 * pi)
        b += 5 * np.arctan(3 * sq3 / 8) / (8 * pi) - np.arctan(5 / sq3) / (
            8 * pi
        )
        c = 3 * sq3 * (ln13 - ln7) / (16 * pi) - np.arctan(5 * sq3 / 9) / (
            8 * pi
        )
        c += np.arctan(sq3 / 9) / (8 * pi) + np.arctan(sq3 / 6) / (2 * pi)
        good_matrix = [
            [a, a, b, c, c, b],
            [b, a, a, b, c, c],
            [c, b, a, a, b, c],
            [c, c, b, a, a, b],
            [b, c, c, b, a, a],
            [a, b, c, c, b, a],
        ]
        good_matrix = np.array(good_matrix)

        angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.incurve(curve, basis, tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestComputeMatrixIncurve::test_eqtriangle_corner",
            "TestComputeMatrixIncurve::test_eqtriangle_middle",
            "TestComputeMatrixIncurve::test_square_corner",
            "TestComputeMatrixIncurve::test_square_middle",
            "TestComputeMatrixIncurve::test_hexagon_corner",
            "TestComputeMatrixIncurve::test_hexagon_middle",
        ]
    )
    def test_end(self):
        pass


class TestComputeMatrixOutcurve:
    """
    Test the values of the integrals:

    int_{Gamma}
    """

    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrixOutcurve::test_begin"])
    def test_square(self):
        ln3, ln5, ln13 = np.log(3), np.log(5), np.log(13)
        atan = np.arctan
        atan2, atan3 = np.arctan(2), np.arctan(3)
        good_matrix = np.zeros((3, 4), dtype="float64")
        good_matrix[0, :] = np.pi / 2
        good_matrix[1, 0] = (
            -ln5 / 4
            + atan(1 / 2) / 4
            + atan(3 / 2) / 4
            + atan(2 / 3)
            + ln13 / 4
        )
        good_matrix[1, 1] = (
            -ln13 / 4
            + 3 * atan(1 / 2) / 4
            + ln5 / 4
            + 3 * atan(3 / 2) / 4
            + atan(2)
        )
        good_matrix[1, 2] = (
            -ln13 / 4
            + 3 * atan(1 / 2) / 4
            + ln5 / 4
            + 3 * atan(3 / 2) / 4
            + atan2
        )
        good_matrix[1, 3] = (
            -ln5 / 4
            + atan(1 / 2) / 4
            + atan(3 / 2) / 4
            + atan(2 / 3)
            + ln13 / 4
        )
        good_matrix[2, 0] = (
            -3 * ln5 / 4 + atan(1 / 3) / 2 + np.pi / 8 + 3 * ln3 / 2
        )
        good_matrix[2, 1] = (
            -3 * ln3 / 4
            + 3 * atan(1 / 3) / 4
            + atan(3) / 4
            + np.pi / 4
            + ln5 / 2
        )
        good_matrix[2, 2] = -ln5 / 4 + 3 * np.pi / 8 + 3 * atan(3) / 2
        good_matrix[2, 3] = (
            -3 * ln3 / 4
            + 3 * atan(1 / 3) / 4
            + atan3 / 4
            + np.pi / 4
            + ln5 / 2
        )

        vertices = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
        sources = [(0, 0), (1, 0), (1, 1)]
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots)
        test_matrix = ComputeStiffness.outcurve(curve, basis, sources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestComputeMatrixOutcurve::test_begin",
            "TestComputeMatrixOutcurve::test_square",
        ]
    )
    def test_end(self):
        pass


class TestPoissonEvaluator:
    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestPoissonEvaluator::test_begin"])
    def test_solid_square(self):
        a, b, c = np.random.uniform(-1, 1, 3)
        a, b, c = 1, 0, 0
        a, b, c = 0, 1, 0
        a, b, c = 0.5, 2, 3
        funct = lambda x, y: a * x + b * y + c
        gradi = lambda x, y: (a, b)

        vertices = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        vertices = np.array(vertices, dtype="float64")
        curve = NurbsCurve.from_vertices(vertices)
        basis = SplineBasisFunction.cyclic(curve.knots, degree=1)
        ctrlpoints = tuple(funct(x, y) for x, y in vertices)
        bound = ScalarFunction(basis, ctrlpoints)
        basis = SplineBasisFunction.cyclic(curve.knots, degree=0)
        ctrlpoints = (-b, a, b, -a)
        normal = ScalarFunction(basis, ctrlpoints)
        evaluator = PoissonEvaluatorCurve(curve, bound, normal)

        # At boundary
        sources = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sources += [(0, -1), (1, 0), (0, 1), (-1, 0)]
        for source in sources:
            x, y = source
            assert evaluator.eval(source) == funct(x, y)
            assert evaluator.grad(source) == gradi(x, y)

        # Exterior
        sources = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
        sources += [(0, -2), (2, 0), (0, 2), (-2, 0)]
        for source in sources:
            x, y = source
            assert evaluator.eval(source) == 0
            assert evaluator.grad(source) == (0, 0)

        # Interior
        sources = [(0, 0), (0, -0.5), (0.5, 0), (0, 0.5), (-0.5, 0)]
        sources += [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
        for source in sources:
            x, y = source
            assert abs(evaluator.eval(source) - funct(x, y)) < 1e-9
            grad_test = evaluator.grad(source)
            grad_good = gradi(x, y)
            assert abs(grad_test[0] - grad_good[0]) < 1e-9
            assert abs(grad_test[1] - grad_good[1]) < 1e-9

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestPoissonEvaluator::test_begin"])
    def test_hollow_square(self):
        a, b, c = np.random.uniform(-1, 1, 3)
        a, b, c = 1, 0, 0
        a, b, c = 0, 1, 0
        a, b, c = 0.5, 2, 3
        funct = lambda x, y: a * x + b * y + c
        gradi = lambda x, y: (a, b)

        vertices = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
        vertices = np.array(vertices, dtype="float64")
        curve_in = NurbsCurve.from_vertices(vertices)
        basis_in = SplineBasisFunction.cyclic(curve_in.knots, degree=1)
        ctrlpoints = tuple(funct(x, y) for x, y in vertices)
        bound_in = ScalarFunction(basis_in, ctrlpoints)
        basis_in = SplineBasisFunction.cyclic(curve_in.knots, degree=0)
        ctrlpoints = (a, -b, -a, b)
        normal = ScalarFunction(basis_in, ctrlpoints)
        evaluator_in = PoissonEvaluatorCurve(curve_in, bound_in, normal)

        vertices = [(-2, -2), (2, -2), (2, 2), (-2, 2)]
        vertices = np.array(vertices, dtype="float64")
        curve_out = NurbsCurve.from_vertices(vertices)
        basis_out = SplineBasisFunction.cyclic(curve_out.knots, degree=1)
        ctrlpoints = tuple(funct(x, y) for x, y in vertices)
        bound_out = ScalarFunction(basis_out, ctrlpoints)
        basis_out = SplineBasisFunction.cyclic(curve_out.knots, degree=0)
        ctrlpoints = (-b, a, b, -a)
        normal = ScalarFunction(basis_out, ctrlpoints)
        evaluator_out = PoissonEvaluatorCurve(curve_out, bound_out, normal)

        evaluator = PoissonEvaluatorGeometry([evaluator_in, evaluator_out])

        # At boundary
        sources = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sources += [(0, -1), (1, 0), (0, 1), (-1, 0)]
        sources += [(0, -2), (2, 0), (0, 2), (-2, 0)]
        sources += [(-2, -2), (2, -2), (2, 2), (-2, 2)]
        for source in sources:
            x, y = source
            assert evaluator.eval(source) == funct(x, y)
            assert evaluator.grad(source) == gradi(x, y)

        # Exterior
        sources = [(0, 0)]
        sources += [(-0.5, 0), (0.5, 0), (0, -0.5), (0, 0.5)]
        sources += [(-3, -3), (3, -3), (3, 3), (-3, 3)]
        sources += [(0, -3), (3, 0), (0, 3), (-3, 0)]
        for source in sources:
            x, y = source
            assert evaluator.eval(source) == 0
            assert evaluator.grad(source) == (0, 0)

        # Interior
        sources = [(0, -1.5), (1.5, 0), (0, 1.5), (-1.5, 0)]
        sources += [(-1.5, -1.5), (1.5, -1.5), (1.5, 1.5), (-1.5, 1.5)]
        for source in sources:
            x, y = source
            assert abs(evaluator.eval(source) - funct(x, y)) < 1e-9
            grad_test = evaluator.grad(source)
            grad_good = gradi(x, y)
            assert abs(grad_test[0] - grad_good[0]) < 1e-9
            assert abs(grad_test[1] - grad_good[1]) < 1e-9

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestPoissonEvaluator::test_solid_square",
            "TestPoissonEvaluator::test_hollow_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(8)
@pytest.mark.dependency(
    depends=[
        "TestComputeMatrixIncurve::test_end",
        "TestComputeMatrixOutcurve::test_end",
        "TestPoissonEvaluator::test_end",
    ]
)
def test_end():
    pass
