"""
File that contains Geometry class, which defines a planar region
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

import numpy as np
from shapepy import ConnectedShape, SimpleShape

from .abcs import NamedTracker
from .curve import Curve


class Geometry(NamedTracker):
    """
    Connected Geometry class that represents a bidimensional region
    on the plane and has some functions such as to decide if points
    are inside the region
    """

    @classmethod
    def from_shape(cls, shape: Union[SimpleShape, ConnectedShape]):
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
            if curve in Curve.instances:
                curves[i] = Curve.instances[curve]
            elif (-curve) in Curve.instances:
                curves[i] = ~Curve.instances[-curve]
            else:
                raise NotImplementedError
        self.__curves = tuple(curves)
        self.name = name

    @property
    def labels(self) -> Tuple[int]:
        """
        Gives the curve labels that defines the geometry
        """
        return tuple(curve.label for curve in self.curves)

    @property
    def curves(self) -> Tuple[Curve]:
        """
        Gives the curves that defines the geometry
        """
        return self.__curves

    def winding(self, point: Tuple[float]) -> float:
        """
        Computes the winding number of the giving point

        Normally winding number is one for each curve, but this method gives
        the overall winding number
        """
        wind_tolerance = 1e-9
        labels = np.array(sorted(self.labels), dtype="int32")
        winds = []
        for label in labels:
            curve = Curve.instances[abs(label)]
            wind = curve.winding(point)
            if wind < 1 - wind_tolerance:
                return wind
            winds.append(wind)
        return 1
