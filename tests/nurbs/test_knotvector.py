import numpy as np
import pytest
from pynurbs import GeneratorKnotVector

from compmec.section.nurbs import KnotVector


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_bezier():
    usample = np.linspace(0, 1, 129)
    for degree in range(6):
        good_knotvector = GeneratorKnotVector.bezier(degree)
        test_knotvector = KnotVector(tuple(good_knotvector))

        assert test_knotvector.degree == good_knotvector.degree
        assert test_knotvector.npts == good_knotvector.npts
        assert test_knotvector.knots == good_knotvector.knots
        for u in test_knotvector.knots:
            assert test_knotvector.span(u) == good_knotvector.span(u)
            assert test_knotvector.mult(u) == good_knotvector.mult(u)
        for test, good in zip(test_knotvector, good_knotvector):
            assert test == good
        for u in usample:
            assert test_knotvector.span(u) == good_knotvector.span(u)
            assert test_knotvector.mult(u) == good_knotvector.mult(u)


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_bezier"])
def test_uniform():
    usample = np.linspace(0, 1, 129)
    for degree in range(6):
        for npts in range(degree + 1, degree + 11):
            good_knotvector = GeneratorKnotVector.uniform(degree, npts)
            test_knotvector = KnotVector(tuple(good_knotvector))

            assert test_knotvector.degree == good_knotvector.degree
            assert test_knotvector.npts == good_knotvector.npts
            assert test_knotvector.knots == good_knotvector.knots
            for u in test_knotvector.knots:
                assert test_knotvector.span(u) == good_knotvector.span(u)
                assert test_knotvector.mult(u) == good_knotvector.mult(u)
            for test, good in zip(test_knotvector, good_knotvector):
                assert test == good
            for u in usample:
                assert test_knotvector.span(u) == good_knotvector.span(u)
                assert test_knotvector.mult(u) == good_knotvector.mult(u)


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[
        "test_begin",
        "test_bezier",
        "test_uniform",
    ]
)
def test_end():
    pass
