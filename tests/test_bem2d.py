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
    def test_eqtriangle(self):

        vertices = [[1, 0], [-0.5, np.sqrt(3) / 2], [-0.5, -np.sqrt(3) / 2]]
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        tsources = basis.knots[: basis.ndofs]
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)
        good_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        good_matrix = np.array(good_matrix) / 12
        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestComputeMatrix::test_begin"])
    def test_square(self):

        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        curve = PolygonCurve(vertices)
        basis = BasisFunc.cyclic(curve.knots)
        tsources = basis.knots[: basis.ndofs]
        computer = ComputeMatrix(curve, basis)
        test_matrix = computer.inpolygon(tsources)
        good_matrix = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
        good_matrix = np.array(good_matrix) / 12
        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_allclose(test_matrix, good_matrix)

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestComputeMatrix::test_eqtriangle",
            "TestComputeMatrix::test_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(8)
@pytest.mark.dependency(depends=["TestComputeMatrix::test_end"])
def test_end():
    pass
