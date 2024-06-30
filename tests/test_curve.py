import numpy as np
import pynurbs
import pytest

from compmec.section.curve import NurbsCurve, PolygonCurve


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_winding_square():
    knotvector = pynurbs.GeneratorKnotVector.uniform(1, 5)
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1], [1, 1]]
    vertices = np.array(vertices, dtype="float64")
    curve = pynurbs.Curve(knotvector, vertices)
    curve = NurbsCurve(curve)

    assert curve.winding((0, 0)) == 1
    assert curve.winding((1, 0)) == 0.5
    assert curve.winding((1, 1)) == 0.25
    assert curve.winding((0, 1)) == 0.5
    assert curve.winding((-1, 1)) == 0.25
    assert curve.winding((-1, 0)) == 0.5
    assert curve.winding((-1, -1)) == 0.25
    assert curve.winding((0, -1)) == 0.5
    assert curve.winding((1, -1)) == 0.25
    assert curve.winding((1, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_square_polygon():
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    vertices = np.array(vertices, dtype="float64")
    curve = PolygonCurve(vertices)

    parameters = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4)
    points = curve.eval(parameters)

    assert np.all(points[0] == (1.0, 1.0))
    assert np.all(points[1] == (0.0, 1.0))
    assert np.all(points[2] == (-1.0, 1.0))
    assert np.all(points[3] == (-1.0, 0.0))
    assert np.all(points[4] == (-1.0, -1.0))
    assert np.all(points[5] == (0.0, -1.0))
    assert np.all(points[6] == (1.0, -1.0))
    assert np.all(points[7] == (1.0, 0.0))
    assert np.all(points[8] == (1.0, 1.0))

    assert curve.winding((0, 0)) == 1
    assert curve.winding((1, 0)) == 0.5
    assert curve.winding((1, 1)) == 0.25
    assert curve.winding((0, 1)) == 0.5
    assert curve.winding((-1, 1)) == 0.25
    assert curve.winding((-1, 0)) == 0.5
    assert curve.winding((-1, -1)) == 0.25
    assert curve.winding((0, -1)) == 0.5
    assert curve.winding((1, -1)) == 0.25
    assert curve.winding((1, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.dependency(depends=["test_winding_square", "test_square_polygon"])
def test_end():
    pass
