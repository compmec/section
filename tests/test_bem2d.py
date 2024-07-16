"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest

from compmec.section.basisfunc import BasisFunc
from compmec.section.bem2d import ComputeStiffness, TorsionEvaluator
from compmec.section.curve import Curve


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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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
        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
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


class TestTorsionVectors:
    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestTorsionVectors::test_begin"])
    def test_square1_tensor_const_vector(self):
        vertices = [(0, 0), (1, 0), (1, 1), (0, 1)]
        vertices = np.array(vertices, dtype="float64")
        good_vector = np.array([0, 1 / 2, 0, -1 / 2])

        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        test_vector = TorsionEvaluator.constant_vector(curve, basis)

        assert len(test_vector) == len(good_vector)
        np.testing.assert_allclose(test_vector, good_vector)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestTorsionVectors::test_begin"])
    def test_square2_tensor_const_vector(self):
        vertices = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        vertices = np.array(vertices, dtype="float64")
        good_vector = np.array([0, 0, 0, 0])

        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        test_vector = TorsionEvaluator.constant_vector(curve, basis)

        assert len(test_vector) == len(good_vector)
        np.testing.assert_allclose(test_vector, good_vector)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestTorsionVectors::test_begin"])
    def test_square3_tensor_const_vector(self):
        vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)]
        vertices = np.array(vertices, dtype="float64")
        good_vector = np.array([0, 0, 0, 0])

        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        test_vector = TorsionEvaluator.constant_vector(curve, basis)

        assert len(test_vector) == len(good_vector)
        np.testing.assert_allclose(test_vector, good_vector)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestTorsionVectors::test_begin"])
    def test_square4_tensor_const_vector(self):
        a, b = 7, 13
        vertices = [
            (a - 1, b - 1),
            (a + 1, b - 1),
            (a + 1, b + 1),
            (a - 1, b + 1),
        ]
        vertices = np.array(vertices, dtype="float64")
        good_vector = np.array([a - b, a + b, b - a, -b - a])

        curve = Curve.from_vertices(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        test_vector = TorsionEvaluator.constant_vector(curve, basis)

        assert len(test_vector) == len(good_vector)
        np.testing.assert_allclose(test_vector, good_vector)

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestTorsionVectors::test_square1_tensor_const_vector",
            "TestTorsionVectors::test_square2_tensor_const_vector",
            "TestTorsionVectors::test_square3_tensor_const_vector",
            "TestTorsionVectors::test_square4_tensor_const_vector",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(8)
@pytest.mark.dependency(
    depends=[
        "TestComputeMatrixIncurve::test_end",
        "TestComputeMatrixOutcurve::test_end",
        "TestTorsionVectors::test_end",
    ]
)
def test_end():
    pass
