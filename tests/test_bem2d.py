"""
File to test solving the Poisson's equation by using Boundary Element Method
"""
import numpy as np
import pytest

from compmec.section.bem2d import IntegralUnV, IntegralUVn


@pytest.mark.order(6)
@pytest.mark.skip()
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
        "tests/test_axial.py::test_end",
        "tests/test_bending.py::test_end",
        "tests/test_integral.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass



class TestMatrixPolygonal:
    """
    Verifies the values of the matrix [UVn]
    when the shape is polygonal and the source points lie
    on the vertices

    U(t) = sum_j phi_j(t) * U_j
    V = ln(r)
    dV/dn = <p, p'>/abs(p')

    """

    @pytest.mark.order(6)
    @pytest.mark.dependency(depends=["test_begin", "TestIntegrals::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(6)
    @pytest.mark.skip(reason="Source points must be at vertices")
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_equilateral_triangle_corner(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        good_matrix = np.array(good_matrix) * np.pi / 6
        test_matrix = IntegralUVn(vertices)
        np.testing.assert_almost_equal(test_matrix, good_matrix)

    @pytest.mark.order(6)
    @pytest.mark.skip(reason="Source points must be at vertices")
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_square_corner(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_matrix0 = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]
        good_matrix1 = [[0, -1, 2, -1], [-1, 0, -1, 2], [2, -1, 0, -1], [-1, 2, -1, 0]]
        good_matrix = np.array(good_matrix0, dtype="float64") * np.pi / 4
        good_matrix += np.array(good_matrix1, dtype="float64") * np.log(2) / 2
        test_matrix = IntegralUVn(vertices)
        np.testing.assert_almost_equal(test_matrix, good_matrix)

    @pytest.mark.order(6)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_equilateral_triangle_middle(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

    @pytest.mark.order(6)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_square_middle(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

    @pytest.mark.order(6)
    @pytest.mark.dependency(
        depends=[
            "TestMatrixPolygonal::test_begin",
            "TestMatrixPolygonal::test_equilateral_triangle_middle",
            "TestMatrixPolygonal::test_square_middle",
        ]
    )
    def test_end(self):
        pass


class TestWarpingVector:
    """
    Verifies the values of the vector [UVn]

    """

    @pytest.mark.order(6)
    @pytest.mark.dependency(depends=["test_begin", "TestIntegrals::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(6)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestWarpingVector::test_begin"])
    def test_equilateral_triangle(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_vector = np.zeros(len(vertices))
        test_vector = IntegralUnV(vertices)
        np.testing.assert_almost_equal(test_vector, good_vector)

    @pytest.mark.order(6)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestWarpingVector::test_begin"])
    def test_square(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_vector = np.zeros(len(vertices))
        test_vector = IntegralUnV(vertices)
        np.testing.assert_almost_equal(test_vector, good_vector)

    @pytest.mark.order(6)
    @pytest.mark.dependency(
        depends=[
            "TestWarpingVector::test_begin",
            "TestWarpingVector::test_equilateral_triangle",
            "TestWarpingVector::test_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(6)
@pytest.mark.dependency(
    depends=[
        "TestIntegrals::test_end",
        "TestMatrixPolygonal::test_end",
        "TestWarpingVector::test_end",
    ]
)
def test_end():
    pass
