import numpy as np
import pytest

from compmec import nurbs
from compmec.section.curve import NurbsCurve


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.timeout(10)
@pytest.mark.dependency(depends=["test_begin"])
def test_winding_square():
    knotvector = nurbs.GeneratorKnotVector.uniform(1, 5)
    vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1], [1, 1]]
    vertices = np.array(vertices, dtype="float64")
    curve = nurbs.Curve(knotvector, vertices)
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
@pytest.mark.dependency(
    depends=[
        "test_winding_square",
    ]
)
def test_end():
    pass


if __name__ == "__main__":
    test_winding_square()
