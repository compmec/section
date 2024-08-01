from typing import Union

from shapepy import ConnectedShape, JordanCurve, SimpleShape

from .abcs import ICurve, IGeometry
from .curve import PolygonCurve
from .geometry import ConnectedGeometry


def jordan2curve(jordan: JordanCurve) -> ICurve:
    if not isinstance(jordan, JordanCurve):
        raise TypeError
    if all(seg.degree == 1 for seg in jordan.segments):
        vertices = [tuple(seg.ctrlpoints[0]) for seg in jordan.segments]
        return PolygonCurve(vertices)
    raise NotImplementedError


def shape2geometry(shape: Union[SimpleShape, ConnectedShape]) -> IGeometry:
    if not isinstance(shape, (SimpleShape, ConnectedShape)):
        raise TypeError
    curves = tuple(map(jordan2curve, shape.jordans))
    return ConnectedGeometry(curves)
