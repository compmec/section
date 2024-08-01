"""
File to tests the basis functions
"""

import numpy as np
import pytest

from compmec.section.basisfunc import KnotVector, SplineBasisFunction


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[],
    scope="session",
)
def test_begin():
    pass


class TestCyclicKnotVector:

    @pytest.mark.order(1)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_begin"])
    def test_standard(self):

        for degree in range(5):
            knots = [0] * (degree + 1) + [1] * (degree + 1)
            knotvector = KnotVector.cyclic(knots, degree)
            np.testing.assert_equal(knotvector, knots)

        for degree in range(5):
            for npts in range(degree + 1, degree + 11):
                ninters = npts - degree
                knots = [0] * degree
                knots += [i / ninters for i in range(ninters + 1)]
                knots += [1] * degree
                knotvector = KnotVector.cyclic(knots, degree)
                np.testing.assert_equal(knotvector, knots)

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_standard"])
    def test_constant(self):
        knots = (0, 1)
        knotvector = KnotVector.cyclic(knots, 0)
        np.testing.assert_equal(knotvector, (0, 1))

        knots = (0, 1, 3)
        knotvector = KnotVector.cyclic(knots, 0)
        np.testing.assert_equal(knotvector, (0, 1, 3))

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_constant"])
    def test_linear(self):

        knots = (0, 1)
        knotvector = KnotVector.cyclic(knots, 1)
        np.testing.assert_equal(knotvector, (-1, 0, 1, 2))

        knots = (0, 1, 3)
        knotvector = KnotVector.cyclic(knots, 1)
        np.testing.assert_equal(knotvector, (-2, 0, 1, 3, 4))

        knots = (0, 0, 1, 1)
        knotvector = KnotVector.cyclic(knots, 1)
        np.testing.assert_equal(knotvector, (0, 0, 1, 1))

        knots = (0, 1, 2)
        knotvector = KnotVector.cyclic(knots, 1)
        np.testing.assert_equal(knotvector, (-1, 0, 1, 2, 3))

        knots = (0, 0, 1, 2, 2)
        knotvector = KnotVector.cyclic(knots, 1)
        np.testing.assert_equal(knotvector, (0, 0, 1, 2, 2))

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_linear"])
    def test_quadratic(self):

        knots = (0, 1)
        knotvector = KnotVector.cyclic(knots, 2)
        np.testing.assert_equal(knotvector, (-2, -1, 0, 1, 2, 3))

        knots = (0, 1, 3)
        knotvector = KnotVector.cyclic(knots, 2)
        np.testing.assert_equal(knotvector, (-3, -2, 0, 1, 3, 4, 6))

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_quadratic"])
    def test_cubic(self):

        knots = (0, 1)
        knotvector = KnotVector.cyclic(knots, 3)
        np.testing.assert_equal(knotvector, (-3, -2, -1, 0, 1, 2, 3, 4))

    @pytest.mark.order(1)
    @pytest.mark.dependency(
        depends=[
            "TestCyclicKnotVector::test_standard",
            "TestCyclicKnotVector::test_constant",
            "TestCyclicKnotVector::test_linear",
            "TestCyclicKnotVector::test_quadratic",
            "TestCyclicKnotVector::test_cubic",
        ]
    )
    def test_end(self):
        pass


class TestCyclic:
    @pytest.mark.order(1)
    @pytest.mark.dependency(depends=["TestCyclicKnotVector::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclic::test_begin"])
    def test_linear(self):
        knots = (0, 0.5, 1)
        basis = SplineBasisFunction(knots, degree=1)
        assert basis.ndofs == 2

        umesh = (0, 0.25, 0.5, 0.75, 1)
        valus = basis.eval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (1, 0.5, 0, 0.5, 1))
        np.testing.assert_equal(valus[:, 1], (0, 0.5, 1, 0.5, 0))

        umesh = (0.1, 0.49, 0.51, 0.99)
        valus = basis.deval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (-2.0, -2.0, 2.0, 2.0))
        np.testing.assert_equal(valus[:, 1], (2.0, 2.0, -2.0, -2.0))

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclic::test_begin"])
    def test_linear2(self):
        knots = (0, 1, 3)
        basis = SplineBasisFunction(knots, degree=1)
        assert basis.ndofs == 2

        umesh = (0, 0.5, 1, 2, 3)
        valus = basis.eval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (1, 0.5, 0, 0.5, 1))
        np.testing.assert_equal(valus[:, 1], (0, 0.5, 1, 0.5, 0))

        umesh = (0.1, 0.99, 1.01, 2.99)
        valus = basis.deval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (-1.0, -1.0, 0.5, 0.5))
        np.testing.assert_equal(valus[:, 1], (1.0, 1.0, -0.5, -0.5))

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclic::test_begin"])
    def test_quadratic(self):
        knots = (0, 0.5, 1)
        basis = SplineBasisFunction(knots, degree=2)
        assert basis.ndofs == 2

        umesh = (0, 0.25, 0.5, 0.75, 1)
        valus = basis.eval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (0.5, 0.25, 0.5, 0.75, 0.5))
        np.testing.assert_equal(valus[:, 1], (0.5, 0.75, 0.5, 0.25, 0.5))

        umesh = (0, 0.5, 1)
        valus = basis.deval(umesh)
        assert valus.shape == (len(umesh), basis.ndofs)
        np.testing.assert_equal(valus[:, 0], (-2.0, 2.0, -2.0))
        np.testing.assert_equal(valus[:, 1], (2.0, -2.0, 2.0))

    @pytest.mark.order(1)
    @pytest.mark.dependency(
        depends=[
            "TestCyclic::test_linear",
            "TestCyclic::test_linear2",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(1)
@pytest.mark.dependency(depends=["TestCyclic::test_end"])
def test_end():
    pass
