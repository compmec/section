import pytest

from compmec.section.curve import NurbsCurve, PolygonCurve


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_winding_square_counterclock():
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    curve = NurbsCurve.from_vertices(vertices)

    assert float(curve) == 4

    assert curve.projection((1, 0)) == (3.5,)
    assert curve.projection((1, 1)) == (0, 4)
    assert curve.projection((0, 1)) == (0.5,)
    assert curve.projection((-1, 1)) == (1,)
    assert curve.projection((-1, 0)) == (1.5,)
    assert curve.projection((-1, -1)) == (2,)
    assert curve.projection((0, -1)) == (2.5,)
    assert curve.projection((1, -1)) == (3,)
    assert curve.projection((1, 0)) == (3.5,)

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
def test_winding_square_clockwise():
    vertices = [[1, 1], [1, -1], [-1, -1], [-1, 1]]
    curve = NurbsCurve.from_vertices(vertices)

    assert float(curve) == -4

    assert curve.projection((0, 0)) == (0.5, 1.5, 2.5, 3.5)
    assert curve.projection((1, 0)) == (0.5,)
    assert curve.projection((1, 1)) == (0, 4)
    assert curve.projection((0, 1)) == (3.5,)
    assert curve.projection((-1, 1)) == (3,)
    assert curve.projection((-1, 0)) == (2.5,)
    assert curve.projection((-1, -1)) == (2,)
    assert curve.projection((0, -1)) == (1.5,)
    assert curve.projection((1, -1)) == (1,)
    assert curve.projection((1, 0)) == (0.5,)

    assert curve.winding((0, 0)) == 0
    assert curve.winding((1, 0)) == 0.5
    assert curve.winding((1, 1)) == 0.75
    assert curve.winding((0, 1)) == 0.5
    assert curve.winding((-1, 1)) == 0.75
    assert curve.winding((-1, 0)) == 0.5
    assert curve.winding((-1, -1)) == 0.75
    assert curve.winding((0, -1)) == 0.5
    assert curve.winding((1, -1)) == 0.75
    assert curve.winding((1, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_polygon_winding_square_counterclock():
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
    curve = PolygonCurve(vertices)

    assert curve.projection((0, 0)) == (0.5, 1.5, 2.5, 3.5)
    assert curve.projection((1, 0)) == (3.5,)
    assert curve.projection((1, 1)) == (0, 4)
    assert curve.projection((0, 1)) == (0.5,)
    assert curve.projection((-1, 1)) == (1,)
    assert curve.projection((-1, 0)) == (1.5,)
    assert curve.projection((-1, -1)) == (2,)
    assert curve.projection((0, -1)) == (2.5,)
    assert curve.projection((1, -1)) == (3,)
    assert curve.projection((1, 0)) == (3.5,)

    assert float(curve) == 4
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
def test_polygon_winding_square_clockwise():
    vertices = [[1, 1], [1, -1], [-1, -1], [-1, 1]]
    curve = PolygonCurve(vertices)

    assert float(curve) == -4

    assert curve.projection((1, 0)) == (0.5,)
    assert curve.projection((1, 1)) == (0, 4)
    assert curve.projection((0, 1)) == (3.5,)
    assert curve.projection((-1, 1)) == (3,)
    assert curve.projection((-1, 0)) == (2.5,)
    assert curve.projection((-1, -1)) == (2,)
    assert curve.projection((0, -1)) == (1.5,)
    assert curve.projection((1, -1)) == (1,)
    assert curve.projection((1, 0)) == (0.5,)

    assert curve.winding((0, 0)) == 0
    assert curve.winding((1, 0)) == 0.5
    assert curve.winding((1, 1)) == 0.75
    assert curve.winding((0, 1)) == 0.5
    assert curve.winding((-1, 1)) == 0.75
    assert curve.winding((-1, 0)) == 0.5
    assert curve.winding((-1, -1)) == 0.75
    assert curve.winding((0, -1)) == 0.5
    assert curve.winding((1, -1)) == 0.75
    assert curve.winding((1, 0)) == 0.5


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[
        "test_winding_square_counterclock",
        "test_winding_square_clockwise",
        "test_polygon_winding_square_counterclock",
        "test_polygon_winding_square_clockwise",
    ]
)
def test_end():
    pass
