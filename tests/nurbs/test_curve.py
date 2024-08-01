import numpy as np
import pytest

from compmec.section.nurbs import CyclicScalarSpline, KnotVector


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=[
        "tests/nurbs/test_knotvector.py::test_end",
        "tests/nurbs/test_basisfunction.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


@pytest.mark.order(3)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_scalar_linear1():
    nodes = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
    knotvector = (-1, 0, 1, 2, 3)
    knotvector = KnotVector(knotvector, degree=1)
    ctrlpoints = (12, 56)
    function = CyclicScalarSpline(knotvector, ctrlpoints)

    good_values = (12, 23, 34, 45, 56, 45, 34, 23, 12)
    test_values = function.eval(nodes)
    assert np.all(test_values == good_values)


@pytest.mark.order(3)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_scalar_linear2():
    nodes = (0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)
    knotvector = (-1, 0, 1, 2, 3, 4)
    knotvector = KnotVector(knotvector, degree=1)
    ctrlpoints = (12, 56, 36)
    function = CyclicScalarSpline(knotvector, ctrlpoints)

    good_values = (12, 23, 34, 45, 56, 51, 46, 41, 36, 30, 24, 18, 12)
    test_values = function.eval(nodes)
    assert np.all(test_values == good_values)


@pytest.mark.order(3)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_scalar_quadratic1():
    nodes = (0, 0.5, 1, 1.5, 2, 2.5, 3)
    knotvector = (-2, -1, 0, 1, 2, 3, 4, 5)
    knotvector = KnotVector(knotvector, degree=2)
    ctrlpoints = (96, 448, 288)
    function = CyclicScalarSpline(knotvector, ctrlpoints)
    assert function.npts == 3

    good_values = (272, 384, 368, 284, 192, 164, 272)
    test_values = function.eval(nodes)
    assert np.all(test_values == good_values)


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=[
        "test_begin",
        "test_scalar_linear1",
        "test_scalar_linear2",
        "test_scalar_quadratic1",
    ]
)
def test_end():
    pass
