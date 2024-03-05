"""
File that contains Geometry class, which defines a planar region
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from compmec.shape.shape import ConnectedShape, DefinedShape, SimpleShape

from .abcs import NamedTracker
from .curve import Curve, NurbsCurve


class Geometry(NamedTracker):
    """
    Connected Geometry class that represents a bidimensional region
    on the plane and has some functions such as to decide if points
    are inside the region
    """

    def __init__(self, curve_labels: Tuple[int], name: Optional[str] = None):
        curve_labels = tuple(map(int, curve_labels))
        for label in curve_labels:
            assert abs(label) in Curve.instances
        self.__curve_labels = curve_labels
        self.name = name

    @property
    def labels(self) -> Tuple[int]:
        """
        Gives the curve labels that defines the geometry
        """
        return self.__curve_labels

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
            if wind_tolerance < wind < 1 - wind_tolerance:
                return wind
            winds.append(wind)
        winds = np.array(winds, dtype="int64")
        mask = (winds < 0.5) & (labels < 0)
        mask |= (winds > 0.5) & (labels > 0)
        return float(np.all(mask))


def shapes2geometries(shapes: Tuple[DefinedShape]) -> Tuple[Geometry]:
    """
    Transform shapes instances into geometry instances

    :param shapes: The shape to be converted
    :type shapes: Union[SimpleShape, ConnectedShape]
    :return: The geometry used
    :rtype: Tuple[]
    """
    geometries = []
    for shape in shapes:
        assert isinstance(shape, (SimpleShape, ConnectedShape))
        curve_labels = []
        for jordan in shape.jordans:
            signal = 1 if float(jordan) > 0 else -1
            if signal < 0:
                jordan = ~jordan
            curve = NurbsCurve.from_jordan(jordan)
            curve_labels.append(signal * curve.label)
        geometry = Geometry(curve_labels)
        geometries.append(geometry)
    return geometries
