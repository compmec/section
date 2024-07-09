import numpy as np
import pytest

from compmec.section.curve import Curve
from compmec.section.integral import Integration, Polynomial, comb


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_comb():
    """
    Tests the binomial function
    """
    assert comb(1, 0) == 1
    assert comb(1, 1) == 1

    assert comb(2, 0) == 1
    assert comb(2, 1) == 2
    assert comb(2, 2) == 1

    assert comb(3, 0) == 1
    assert comb(3, 1) == 3
    assert comb(3, 2) == 3
    assert comb(3, 3) == 1

    assert comb(4, 0) == 1
    assert comb(4, 1) == 4
    assert comb(4, 2) == 6
    assert comb(4, 3) == 4
    assert comb(4, 4) == 1

    assert comb(5, 0) == 1
    assert comb(5, 1) == 5
    assert comb(5, 2) == 10
    assert comb(5, 3) == 10
    assert comb(5, 4) == 5
    assert comb(5, 5) == 1

    assert comb(6, 0) == 1
    assert comb(6, 1) == 6
    assert comb(6, 2) == 15
    assert comb(6, 3) == 20
    assert comb(6, 4) == 15
    assert comb(6, 5) == 6
    assert comb(6, 6) == 1


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_closed_newton():
    """
    Tests if the given nodes and weights integrates exactly
    a polynomial of degree p, with p+1 evaluation points

    We start with P(x) = x^{p}
    And after we evaluate P(x) = sum_i c_i * x^i
    for random c_i
    """
    for nptsinteg in range(2, 8):
        nodes, weights = Integration.closed(nptsinteg)
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


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_opened_newton():
    """
    Tests if the given nodes and weights integrates exactly
    a polynomial of degree p, with p+1 evaluation points

    We start with P(x) = x^{p}
    And after we evaluate P(x) = sum_i c_i * x^i
    for random c_i
    """
    for nptsinteg in range(1, 4):
        nodes, weights = Integration.opened(nptsinteg)
        for degree in range(nptsinteg):  # Integrate each basis
            good_integral = 1 / (1 + degree)
            fvalues = nodes**degree
            test_integral = np.inner(fvalues, weights)
            diff = abs(test_integral - good_integral)
            assert abs(diff) < 1e-5

    for nptsinteg in range(1, 4):
        nodes, weights = Integration.opened(nptsinteg)
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


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_chebyshev():
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


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_gauss():
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


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_singular_logarithm():
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
            good_integral = -sum(
                ci / ((1 + i) ** 2) for i, ci in enumerate(coefs)
            )
            fvalues.fill(0)
            for i, node in enumerate(nodes):
                for ck in coefs[::-1]:  # Horner's method
                    fvalues[i] *= node
                    fvalues[i] += ck
            test_integral = np.inner(fvalues, weights)
            diff = abs(test_integral - good_integral)
            assert diff < 1e-15


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_integral_polygon():
    """
    Tests the polynomial integrals of polygons
    """
    vertices = ((1, 1), (-1, 1), (-1, -1), (1, -1))
    geomprops = Polynomial.polygon(vertices)
    assert geomprops[0] == 4  # area
    assert geomprops[1] == 0  # Qx
    assert geomprops[2] == 0  # Qy
    assert geomprops[3] == 4 / 3  # Ixx
    assert geomprops[4] == 0  # Ixy
    assert geomprops[5] == 4 / 3  # Iyy
    assert geomprops[6] == 0  # Ixxx
    assert geomprops[7] == 0  # Ixxy
    assert geomprops[8] == 0  # Ixyy
    assert geomprops[9] == 0  # Iyyy


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_integral_polygon_adaptative():
    """
    Tests the polynomial integrals of polygons
    """
    vertices = ((1, 1), (-1, 1), (-1, -1), (1, -1))
    curve = Curve.from_vertices(vertices)
    geomprops = Polynomial.adaptative(curve)
    assert geomprops[0] == 4  # area
    assert geomprops[1] == 0  # Qx
    assert geomprops[2] == 0  # Qy
    assert abs(geomprops[3] - 4 / 3) < 1e-9  # Ixx
    assert abs(geomprops[4]) < 1e-9  # Ixy
    assert abs(geomprops[5] - 4 / 3) < 1e-9  # Iyy
    assert abs(geomprops[6]) < 1e-9  # Ixxx
    assert abs(geomprops[7]) < 1e-9  # Ixxy
    assert abs(geomprops[8]) < 1e-9  # Ixyy
    assert abs(geomprops[9]) < 1e-9  # Iyyy


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_fail():
    """
    Tests the fail cases
    """
    with pytest.raises(ValueError):
        Integration.opened(0)
    with pytest.raises(ValueError):
        Integration.closed(1)
    with pytest.raises(ValueError):
        Integration.chebyshev(0)
    with pytest.raises(ValueError):
        Integration.gauss(0)
    with pytest.raises(ValueError):
        Integration.log(0)


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[
        "test_comb",
        "test_closed_newton",
        "test_opened_newton",
        "test_chebyshev",
        "test_gauss",
        "test_singular_logarithm",
        "test_integral_polygon",
        "test_integral_polygon_adaptative",
        "test_fail",
    ]
)
def test_end():
    pass
