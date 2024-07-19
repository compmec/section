"""
File that contains Geometry class, which defines a planar region
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

from shapepy import ConnectedShape, SimpleShape

from .abcs import IGeometry, NamedTracker
from .curve import Curve
from .integral import AdaptativePolynomialIntegrator


class ConnectedGeometry(IGeometry, NamedTracker):
    """
    Connected Geometry class that represents a bidimensional region
    on the plane and has some functions such as to decide if points
    are inside the region
    """

    @classmethod
    def from_shape(cls, shape: Union[SimpleShape, ConnectedShape]):
        """
        Creates a Geometry instance from a shape of shapepy package
        """
        assert isinstance(shape, (SimpleShape, ConnectedShape))
        curves = []
        for jordan in shape.jordans:
            curve = Curve.from_jordan(jordan)
            curves.append(curve)
        return cls(curves)

    def __init__(
        self, curves: Tuple[Union[int, Curve]], *, name: Optional[str] = None
    ):
        curves = list(curves)
        for i, curve in enumerate(curves):
            if isinstance(curve, Curve):
                continue
            if not isinstance(curve, int):
                raise NotImplementedError
            if abs(curve) not in Curve.instances:
                raise NotImplementedError
            if curve > 0:
                curves[i] = Curve.instances[curve]
            else:
                curves[i] = ~Curve.instances[-curve]
        self.__curves = tuple(curves)
        self.name = name

    def __iter__(self):
        yield from self.__curves

    def integrate(
        self, expx: int, expy: int, tolerance: Optional[float] = 1e-9
    ) -> float:
        assert isinstance(expx, int) and expx >= 0
        assert isinstance(expy, int) and expy >= 0

        result = 0
        for curve in self:
            integrator = AdaptativePolynomialIntegrator(
                curve, tolerance=tolerance
            )
            value = integrator.integrate(expx, expy)
            result += value / (2 + expx + expy)
        return result

    def winding(self, point: Tuple[float]) -> float:
        wind_tolerance = 1e-9
        for curve in self:
            wind = curve.winding(point)
            if wind < 1 - wind_tolerance:
                return wind
        return 1
