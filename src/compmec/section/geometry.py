"""
File that contains Geometry class, which defines a planar region
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

from .abcs import ICurve, IGeometry, NamedTracker
from .integral import Bidimensional


class ConnectedGeometry(IGeometry, NamedTracker):
    """
    Connected Geometry class that represents a bidimensional region
    on the plane and has some functions such as to decide if points
    are inside the region
    """

    def __init__(
        self, curves: Tuple[Union[int, ICurve]], *, name: Optional[str] = None
    ):
        curves = list(curves)
        for i, curve in enumerate(curves):
            if isinstance(curve, ICurve):
                continue
            if not isinstance(curve, int):
                raise NotImplementedError
            if curve not in ICurve.instances:
                raise NotImplementedError
            curves[i] = ICurve.instances[curve]
        self.__curves = tuple(curves)
        self.name = name

    @property
    def curves(self) -> Tuple[ICurve]:
        return self.__curves

    def integrate(
        self, expx: int, expy: int, tolerance: Optional[float] = 1e-9
    ) -> float:
        assert isinstance(expx, int) and expx >= 0
        assert isinstance(expy, int) and expy >= 0

        result = 0
        for curve in self.curves:
            result += Bidimensional.general(curve, expx, expy)
        return result

    def winding(self, point: Tuple[float]) -> float:
        wind_tolerance = 1e-9
        for curve in self.curves:
            wind = curve.winding(point)
            if wind < 1 - wind_tolerance:
                return wind
        return 1
