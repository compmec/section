"""
This file tests the basic geometry properties, such as
* Area
* First moment of inertia
* Second moment of inertia
"""
import numpy as np
import pytest

from compmec.section.bem2d import CornerAngles, IntegrateUnV, IntegrateUVn, Integration


@pytest.mark.order(4)
@pytest.mark.dependency(
    depends=["tests/test_material.py::test_end"],
    scope="session",
)
@pytest.mark.dependency()
def test_begin():
    pass


class TestIntegrals:
    """
    Verifies the values of the matrix [UVn]

    The
    """

    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestIntegrals::test_begin"])
    def test_singular_logarithm(self):
        """
        Tests if the given nodes and weights integrates exactly
        a polynomial of degree p, with p+1 evaluation points

        We start with P(x) = x^{p}
        And after we evaluate P(x) = sum_i c_i * x^i
        for random c_i
        """
        for nptsinteg in range(1, 9):
            nodes, weights = Integration.log(nptsinteg)
            for degree in range(nptsinteg):  # Integrate each basis
                good_integral = -1 / ((1 + degree) ** 2)
                fvalues = nodes**degree
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert abs(diff) < 1e-5

        for nptsinteg in range(1, 9):
            nodes, weights = Integration.log(nptsinteg)
            fvalues = np.zeros(len(nodes))
            for degree in range(nptsinteg):
                coefs = np.random.uniform(-1, 1, degree + 1)
                good_integral = -sum(ci / ((1 + i) ** 2) for i, ci in enumerate(coefs))
                fvalues.fill(0)
                for i, node in enumerate(nodes):
                    for ck in coefs[::-1]:  # Horner's method
                        fvalues[i] *= node
                        fvalues[i] += ck
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert diff < 1e-15

    @pytest.mark.order(4)
    @pytest.mark.dependency(
        depends=[
            "TestIntegrals::test_singular_logarithm",
        ]
    )
    def test_end(self):
        pass


class TestMatrixPolygonal:
    """
    Verifies the values of the matrix [UVn]

    """

    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin", "TestIntegrals::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_equilateral_triangle(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_matrix = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
        good_matrix = np.array(good_matrix) * np.pi / 6
        test_matrix = IntegrateUVn(vertices)
        np.testing.assert_almost_equal(test_matrix, good_matrix)

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_square(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_matrix0 = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]
        good_matrix1 = [[0, -1, 2, -1], [-1, 0, -1, 2], [2, -1, 0, -1], [-1, 2, -1, 0]]
        good_matrix = np.array(good_matrix0, dtype="float64") * np.pi / 4
        good_matrix += np.array(good_matrix1, dtype="float64") * np.log(2) / 2
        test_matrix = IntegrateUVn(vertices)
        np.testing.assert_almost_equal(test_matrix, good_matrix)

    @pytest.mark.order(4)
    @pytest.mark.dependency(
        depends=[
            "TestMatrixPolygonal::test_begin",
            "TestMatrixPolygonal::test_equilateral_triangle",
            "TestMatrixPolygonal::test_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(4)
@pytest.mark.dependency(
    depends=["TestIntegrals::test_end", "TestMatrixPolygonal::test_end"]
)
def test_end():
    pass
