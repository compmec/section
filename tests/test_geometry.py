"""
File to test src/compmec/section/geometry.py module
"""

import numpy as np
import pytest

from compmec import nurbs
from compmec.section.curve import NurbsCurve
from compmec.section.geometry import Geometry


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_full_square():
    knotvector = nurbs.GeneratorKnotVector.uniform(1, 5)
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1], [1, 1]]
    vertices = np.array(vertices, dtype="float64")
    curve = nurbs.Curve(knotvector, vertices)
    curve = NurbsCurve(curve)

    geometry = Geometry([curve.label])

    assert geometry.winding((0, 0)) == 1
    assert geometry.winding((1, 0)) == 0.5
    assert geometry.winding((1, 1)) == 0.25
    assert geometry.winding((0, 1)) == 0.5
    assert geometry.winding((-1, 1)) == 0.25
    assert geometry.winding((-1, 0)) == 0.5
    assert geometry.winding((-1, -1)) == 0.25
    assert geometry.winding((0, -1)) == 0.5
    assert geometry.winding((1, -1)) == 0.25
    assert geometry.winding((1, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_hollow_square():
    knotvector = nurbs.GeneratorKnotVector.uniform(1, 5)
    vertices = [[3, 3], [-3, 3], [-3, -3], [3, -3], [3, 3]]
    vertices = np.array(vertices, dtype="float64")
    curve_ext = nurbs.Curve(knotvector, vertices)
    curve_ext = NurbsCurve(curve_ext)

    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1], [1, 1]]
    vertices = np.array(vertices, dtype="float64")
    curve_int = nurbs.Curve(knotvector, vertices)
    curve_int = NurbsCurve(curve_int)

    geometry = Geometry([curve_ext.label, -curve_int.label])

    assert geometry.winding((0, 0)) == 0
    assert geometry.winding((1, 0)) == 0.5
    assert geometry.winding((1, 1)) == 0.25
    assert geometry.winding((0, 1)) == 0.5
    assert geometry.winding((-1, 1)) == 0.25
    assert geometry.winding((-1, 0)) == 0.5
    assert geometry.winding((-1, -1)) == 0.25
    assert geometry.winding((0, -1)) == 0.5
    assert geometry.winding((1, -1)) == 0.25
    assert geometry.winding((1, 0)) == 0.5

    assert geometry.winding((2, 0)) == 1
    assert geometry.winding((2, 2)) == 1
    assert geometry.winding((0, 2)) == 1
    assert geometry.winding((-2, 2)) == 1
    assert geometry.winding((-2, 0)) == 1
    assert geometry.winding((-2, -2)) == 1
    assert geometry.winding((0, -2)) == 1
    assert geometry.winding((2, -2)) == 1
    assert geometry.winding((2, 0)) == 1

    assert geometry.winding((3, 0)) == 0.5
    assert geometry.winding((3, 3)) == 0.25
    assert geometry.winding((0, 3)) == 0.5
    assert geometry.winding((-3, 3)) == 0.25
    assert geometry.winding((-3, 0)) == 0.5
    assert geometry.winding((-3, -3)) == 0.25
    assert geometry.winding((0, -3)) == 0.5
    assert geometry.winding((3, -3)) == 0.25
    assert geometry.winding((3, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[
        "test_full_square",
        "test_hollow_square",
    ]
)
def test_end():
    pass
