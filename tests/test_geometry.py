"""
File to test src/compmec/section/geometry.py module
"""

import pytest

from compmec.section.curve import Curve
from compmec.section.geometry import Geometry


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[
        "tests/test_curve.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_full_square():
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    curve = Curve.from_vertices(vertices)

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
    vertices = [[3, 3], [-3, 3], [-3, -3], [3, -3]]
    curve_ext = Curve.from_vertices(vertices)

    vertices = [[-1, -1], [-1, 1], [1, 1], [1, -1]]
    curve_int = Curve.from_vertices(vertices)

    geometry = Geometry([curve_ext, curve_int])

    assert geometry.winding((0, 0)) == 0
    assert geometry.winding((1, 0)) == 0.5
    assert geometry.winding((1, 1)) == 0.75
    assert geometry.winding((0, 1)) == 0.5
    assert geometry.winding((-1, 1)) == 0.75
    assert geometry.winding((-1, 0)) == 0.5
    assert geometry.winding((-1, -1)) == 0.75
    assert geometry.winding((0, -1)) == 0.5
    assert geometry.winding((1, -1)) == 0.75
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
