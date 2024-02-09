"""
File to test solving the Poisson's equation by using Boundary Element Method
"""
import numpy as np
import pytest

from compmec.section.bem2d import IntegralUnV, IntegralUVn, Integration


@pytest.mark.order(4)
@pytest.mark.skip()
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
        "tests/test_axial.py::test_end",
        "tests/test_bending.py::test_end",
    ],
    scope="session",
)
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
    def test_closed_newton(self):
        """
        Tests if the given nodes and weights integrates exactly
        a polynomial of degree p, with p+1 evaluation points

        We start with P(x) = x^{p}
        And after we evaluate P(x) = sum_i c_i * x^i
        for random c_i
        """
        for nptsinteg in range(2, 8):
            nodes, weights = Integration.closed(nptsinteg)
            print("npts = ", nptsinteg)
            for degree in range(nptsinteg):  # Integrate each basis
                good_integral = 1 / (1 + degree)
                fvalues = nodes**degree
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert abs(diff) < 1e-5

        for nptsinteg in range(2, 8):
            nodes, weights = Integration.closed(nptsinteg)
            fvalues = np.zeros(len(nodes))
            for degree in range(nptsinteg):
                coefs = np.random.uniform(-1, 1, degree + 1)
                good_integral = sum(ci / (1 + i) for i, ci in enumerate(coefs))
                fvalues.fill(0)
                for i, node in enumerate(nodes):
                    for ck in coefs[::-1]:  # Horner's method
                        fvalues[i] *= node
                        fvalues[i] += ck
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert diff < 1e-15

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestIntegrals::test_begin"])
    def test_chebyshev(self):
        """
        Tests if the given nodes and weights integrates exactly
        a polynomial of degree p, with p+1 evaluation points

        We start with P(x) = x^{p}
        And after we evaluate P(x) = sum_i c_i * x^i
        for random c_i
        """
        for nptsinteg in range(1, 7):
            nodes, weights = Integration.chebyshev(nptsinteg)
            for degree in range(nptsinteg):  # Integrate each basis
                good_integral = 1 / (1 + degree)
                fvalues = nodes**degree
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert abs(diff) < 1e-5

        for nptsinteg in range(1, 7):
            nodes, weights = Integration.chebyshev(nptsinteg)
            fvalues = np.zeros(len(nodes))
            for degree in range(nptsinteg):
                coefs = np.random.uniform(-1, 1, degree + 1)
                good_integral = sum(ci / (1 + i) for i, ci in enumerate(coefs))
                fvalues.fill(0)
                for i, node in enumerate(nodes):
                    for ck in coefs[::-1]:  # Horner's method
                        fvalues[i] *= node
                        fvalues[i] += ck
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert diff < 1e-15

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestIntegrals::test_begin"])
    def test_gauss(self):
        """
        Tests if the given nodes and weights integrates exactly
        a polynomial of degree p, with p+1 evaluation points

        We start with P(x) = x^{p}
        And after we evaluate P(x) = sum_i c_i * x^i
        for random c_i
        """
        for nptsinteg in range(1, 9):
            nodes, weights = Integration.gauss(nptsinteg)
            for degree in range(2 * nptsinteg):  # Integrate each basis
                good_integral = 1 / (1 + degree)
                fvalues = nodes**degree
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert abs(diff) < 1e-5

        for nptsinteg in range(1, 9):
            nodes, weights = Integration.gauss(nptsinteg)
            fvalues = np.zeros(len(nodes))
            for degree in range(2 * nptsinteg):
                coefs = np.random.uniform(-1, 1, degree + 1)
                good_integral = sum(ci / (1 + i) for i, ci in enumerate(coefs))
                fvalues.fill(0)
                for i, node in enumerate(nodes):
                    for ck in coefs[::-1]:  # Horner's method
                        fvalues[i] *= node
                        fvalues[i] += ck
                test_integral = np.inner(fvalues, weights)
                diff = abs(test_integral - good_integral)
                assert diff < 1e-9

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
            "TestIntegrals::test_closed_newton",
            "TestIntegrals::test_chebyshev",
            "TestIntegrals::test_gauss",
            "TestIntegrals::test_singular_logarithm",
        ]
    )
    def test_end(self):
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

    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin", "TestIntegrals::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
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

    @pytest.mark.order(4)
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

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_equilateral_triangle_middle(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestMatrixPolygonal::test_begin"])
    def test_square_middle(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

    @pytest.mark.order(4)
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

    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin", "TestIntegrals::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestWarpingVector::test_begin"])
    def test_equilateral_triangle(self):
        theta = np.linspace(0, 2 * np.pi, 3, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_vector = np.zeros(len(vertices))
        test_vector = IntegralUnV(vertices)
        np.testing.assert_almost_equal(test_vector, good_vector)

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestWarpingVector::test_begin"])
    def test_square(self):
        theta = np.linspace(0, 2 * np.pi, 4, False)
        xvals, yvals = np.cos(theta), np.sin(theta)
        vertices = np.transpose([xvals, yvals])

        good_vector = np.zeros(len(vertices))
        test_vector = IntegralUnV(vertices)
        np.testing.assert_almost_equal(test_vector, good_vector)

    @pytest.mark.order(4)
    @pytest.mark.dependency(
        depends=[
            "TestWarpingVector::test_begin",
            "TestWarpingVector::test_equilateral_triangle",
            "TestWarpingVector::test_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(4)
@pytest.mark.dependency(
    depends=[
        "TestIntegrals::test_end",
        "TestMatrixPolygonal::test_end",
        "TestWarpingVector::test_end",
    ]
)
def test_end():
    pass
