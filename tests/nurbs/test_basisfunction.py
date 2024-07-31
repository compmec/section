from fractions import Fraction

import numpy as np
import pytest
from matplotlib import pyplot as plt
from pynurbs import Function, GeneratorKnotVector

from compmec.section.nurbs import CyclicSplineBasisFunction, KnotVector


@pytest.mark.order(2)
@pytest.mark.dependency(
    depends=[
        "tests/nurbs/test_knotvector.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


@pytest.mark.order(2)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_bezier():
    usample = np.linspace(0, 1, 128, endpoint=False)
    for degree in range(2, 6):
        good_knotvector = GeneratorKnotVector.bezier(degree, Fraction)
        test_knotvector = KnotVector(tuple(good_knotvector))
        good_basis = Function(good_knotvector)
        test_basis = CyclicSplineBasisFunction(test_knotvector)
        assert test_basis.npts == degree + 1

        good_matrix = good_basis.eval(usample)
        good_matrix = np.array(good_matrix)
        test_matrix = test_basis.eval(usample)
        assert test_matrix.shape == good_matrix.shape
        np.testing.assert_equal(test_matrix, good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_bezier"])
def test_uniform():
    usample = np.linspace(0, 1, 128, endpoint=False)
    for degree in range(6):
        for npts in range(degree + 1, degree + 11):
            good_knotvector = GeneratorKnotVector.uniform(
                degree, npts, Fraction
            )
            test_knotvector = KnotVector(tuple(good_knotvector))

            good_basis = Function(good_knotvector)
            test_basis = CyclicSplineBasisFunction(test_knotvector)
            assert test_basis.npts == npts

            good_matrix = good_basis.eval(usample)
            good_matrix = np.array(good_matrix)
            test_matrix = test_basis.eval(usample)
            assert test_matrix.shape == good_matrix.shape
            np.testing.assert_equal(test_matrix, good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(30)
@pytest.mark.dependency(depends=["test_uniform"])
def test_custom1():
    nodes = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
    knotvector = (-1, 0, 1, 2, 3)
    knotvector = KnotVector(knotvector, degree=1)
    test_basis = CyclicSplineBasisFunction(knotvector)
    assert test_basis.npts == 2

    good_matrix = [
        [1.0, 0.75, 0.5, 0.25, 0.0, 0.25, 0.5, 0.75, 1.0],
        [0.0, 0.25, 0.5, 0.75, 1.0, 0.75, 0.5, 0.25, 0.0],
    ]
    test_matrix = test_basis.eval(nodes)
    assert np.all(test_matrix == good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(30)
@pytest.mark.dependency(depends=["test_uniform"])
def test_custom2():
    nodes = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)
    knotvector = (-1, 0, 1, 2, 3, 4)
    knotvector = KnotVector(knotvector, degree=1)
    test_basis = CyclicSplineBasisFunction(knotvector)
    assert test_basis.npts == 3

    good_matrix = [
        [1, 0.75, 0.5, 0.25, 0, 0, 0, 0, 0, 0.25, 0.5, 0.75, 1],
        [0, 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25, 0],
    ]
    test_matrix = test_basis.eval(nodes)
    assert np.all(test_matrix == good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(30)
@pytest.mark.dependency(depends=["test_uniform"])
def test_custom3():
    nodes = (0, 0.25, 0.5, 0.75, 1)
    knotvector = (-2, -1, 0, 1, 2, 3)
    knotvector = KnotVector(knotvector, degree=2)
    test_basis = CyclicSplineBasisFunction(knotvector)
    assert test_basis.npts == 1
    good_matrix = [
        [1, 1, 1, 1, 1],
    ]
    test_matrix = test_basis.eval(nodes)
    assert np.all(test_matrix == good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(30)
@pytest.mark.dependency(depends=["test_uniform"])
def test_custom4():
    nodes = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
    knotvector = (-2, -1, 0, 1, 2, 3, 4)
    knotvector = KnotVector(knotvector, degree=2)
    test_basis = CyclicSplineBasisFunction(knotvector)
    assert test_basis.npts == 2
    good_matrix = [
        [0.5, 0.3125, 0.25, 0.3125, 0.5, 0.6875, 0.75, 0.6875, 0.5],
        [0.5, 0.6875, 0.75, 0.6875, 0.5, 0.3125, 0.25, 0.3125, 0.5],
    ]
    test_matrix = test_basis.eval(nodes)
    assert np.all(test_matrix == good_matrix)


@pytest.mark.order(2)
@pytest.mark.timeout(30)
@pytest.mark.dependency(depends=["test_uniform"])
def test_custom5():
    nodes = (0, 0.5, 1, 1.5, 2, 2.5, 3)
    knotvector = (-2, -1, 0, 1, 2, 3, 4, 5)
    knotvector = KnotVector(knotvector, degree=2)
    test_basis = CyclicSplineBasisFunction(knotvector)
    assert test_basis.npts == 3
    good_matrix = [
        [0.5, 0.125, 0, 0.125, 0.5, 0.75, 0.5],
        [0.5, 0.75, 0.5, 0.125, 0, 0.125, 0.5],
        [0, 0.125, 0.5, 0.75, 0.5, 0.125, 0],
    ]
    test_matrix = test_basis.eval(nodes)
    assert np.all(test_matrix == good_matrix)


@pytest.mark.order(2)
@pytest.mark.dependency(
    depends=[
        "test_begin",
        "test_bezier",
        "test_uniform",
        "test_custom1",
        "test_custom2",
        "test_custom3",
        "test_custom4",
        "test_custom5",
    ]
)
def test_end():
    pass


if __name__ == "__main__":
    test_custom()
