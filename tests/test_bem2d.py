"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest

from compmec.section.basisfunc import BasisFunc
from compmec.section.bem2d import ComputeMatrix
from compmec.section.curve import PolygonCurve


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


class TestComputeMatrix:
    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrix::test_begin"])
    def test_eqtriangle_corner(self):
        good_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        good_matrix = np.array(good_matrix) / 12

        angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0, 1, 2)
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrix::test_begin"])
    def test_eqtriangle_middle(self):
        pi, ln3, sq3 = np.pi, np.log(3), np.sqrt(3)
        good_matrix = [[-1, -1, 2], [2, -1, -1], [-1, 2, -1]]
        good_matrix = sq3 * ln3 * np.array(good_matrix) / (16 * pi)
        good_matrix += np.array([[3, 3, 2], [2, 3, 3], [3, 2, 3]]) / 16

        angles = np.linspace(0, 2 * pi, 3, endpoint=False)
        vertices = tuple((np.cos(angle), np.sin(angle)) for angle in angles)
        tsources = (0.5, 1.5, 2.5)
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrix::test_begin"])
    def test_square_corner(self):
        a = 1 / 8 - np.log(2) / (4 * np.pi)
        b = np.log(2) / (2 * np.pi)
        good_matrix = [[0, a, b, a], [a, 0, a, b], [b, a, 0, a], [a, b, a, 0]]
        good_matrix = np.array(good_matrix)

        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        tsources = (0, 1, 2, 3)
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrix::test_begin"])
    def test_square_middle(self):
        a = np.arctan(2) / (2 * np.pi) - np.log(5) / (8 * np.pi)
        b = 1 / 4 + np.log(5) / (8 * np.pi) - np.arctan(2) / (2 * np.pi)
        good_matrix = [[a, a, b, b], [b, a, a, b], [b, b, a, a], [a, b, b, a]]
        good_matrix = np.array(good_matrix)

        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        tsources = (0.5, 1.5, 2.5, 3.5)
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)

        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestComputeMatrix::test_eqtriangle_corner",
            "TestComputeMatrix::test_eqtriangle_middle",
            "TestComputeMatrix::test_square_corner",
            "TestComputeMatrix::test_square_middle",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(8)
@pytest.mark.dependency(depends=["TestComputeMatrix::test_end"])
def test_end():
    pass
