import pytest

import numpy as np
from compmec.section.integral import Integration

class TestIntegrals:
    """
    Verifies the values of the matrix [UVn]

    The
    """

    @pytest.mark.order(5)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(5)
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

    @pytest.mark.order(5)
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

    @pytest.mark.order(5)
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

    @pytest.mark.order(5)
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

    @pytest.mark.order(5)
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
