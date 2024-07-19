"""
This file contains the principal class of this packaged, called "Section"
It uses all the other files to compute geometric properties, torsional
constant, torsion and shear center and others.


"""

from __future__ import annotations

from collections import OrderedDict
from typing import Optional, Tuple, Union

import numpy as np
from shapepy.shape import DefinedShape

from .abcs import IGeometry, IMaterial, ISection, NamedTracker
from .field import ChargedField
from .geometry import ConnectedGeometry
from .material import Material


class HomogeneousSection(ISection, NamedTracker):
    """
    HomogeneousSection's class

    """

    instances = OrderedDict()

    @classmethod
    def from_shape(cls, shape: DefinedShape, material: Material):
        geometry = ConnectedGeometry.from_shape(shape)
        return cls(geometry, material)

    def __init__(
        self,
        geometry: Union[str, IGeometry],
        material: Union[str, IMaterial],
        *,
        name: Optional[str] = None,
    ):
        if not isinstance(geometry, IGeometry):
            if not isinstance(geometry, str):
                raise NotImplementedError
            geometry = ConnectedGeometry.instances[geometry]
        if not isinstance(material, IMaterial):
            if not isinstance(material, str):
                raise NotImplementedError
            material = Material.instances[material]
        self.geometry = geometry
        self.material = material
        self.__geomintegs = None
        self.name = name

    def __compute_geomintegs(self):
        """
        Compute the geometric integrals over the domain
        creating the object __geomintegs
        """
        geomintegs = np.zeros(10, dtype="float64")
        geomintegs[0] = self.geometry.integrate(0, 0)
        geomintegs[1] = self.geometry.integrate(0, 1)
        geomintegs[2] = self.geometry.integrate(1, 0)
        geomintegs[3] = self.geometry.integrate(0, 2)
        geomintegs[4] = self.geometry.integrate(1, 1)
        geomintegs[5] = self.geometry.integrate(2, 0)
        geomintegs[6] = self.geometry.integrate(0, 3)
        geomintegs[7] = self.geometry.integrate(1, 2)
        geomintegs[8] = self.geometry.integrate(2, 1)
        geomintegs[9] = self.geometry.integrate(3, 0)
        self.__geomintegs = geomintegs

    def area(self) -> float:
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        return self.__geomintegs[0]

    def first_moment(self) -> Tuple[float]:
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        iqx, iqy = self.__geomintegs[1:3]
        return iqx, iqy

    def second_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        area = self.area()
        ixx, ixy, iyy = self.__geomintegs[3:6]
        ixx -= area * center[1] ** 2
        ixy -= area * center[0] * center[1]
        iyy -= area * center[0] ** 2
        return ixx, ixy, iyy

    def third_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        area = self.area()
        ixxx, ixxy, ixyy, iyyy = self.__geomintegs[6:10]
        ixxx -= area * center[1] ** 3
        ixxy -= area * center[0] * center[1] ** 2
        ixyy -= area * center[0] ** 2 * center[1]
        iyyy -= area * center[0] ** 3
        return ixxx, ixxy, ixyy, iyyy

    def geometric_center(self) -> Tuple[float]:
        iqx, iqy = self.first_moment()
        area = self.area()
        return iqx / area, iqy / area

    def bending_center(self) -> Tuple[float]:
        return self.geometric_center()

    def gyradius(self) -> Tuple[float]:
        area = self.area()
        ixx, _, iyy = self.second_moment()
        return np.sqrt(iyy / area), np.sqrt(ixx / area)

    def charged_field(self) -> ChargedField:
        return ChargedField(self)

    def torsion_center(self) -> Tuple[float]:
        return (0, 0)

    def torsion_constant(self) -> float:
        center = self.geometric_center()
        ixx, _, iyy = self.second_moment(center)
        return ixx + iyy

    def shear_center(self) -> Tuple[float]:
        return (0, 0)

    def __iter__(self):
        yield self
