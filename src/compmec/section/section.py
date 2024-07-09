"""
This file contains the principal class of this packaged, called "Section"
It uses all the other files to compute geometric properties, torsional
constant, torsion and shear center and others.


"""

from __future__ import annotations

from collections import OrderedDict
from typing import Iterable, Optional, Tuple, Union

import numpy as np
from shapepy.shape import DefinedShape, DisjointShape

from .abcs import ISection, NamedTracker
from .curve import Curve
from .field import ChargedField
from .geometry import Geometry
from .integral import Polynomial
from .material import Material


class BaseSection(ISection, NamedTracker):
    """
    BaseSection class that is the base for others section classes

    Basicaly this class is responsible to construct the initial section,
    verifying if the inputs are correct, and generate getter only properties
    such as 'curves', 'geomlabels' and 'materials'.
    """

    instances = OrderedDict()

    def __init__(
        self,
        geometries: Union[str, Geometry, Tuple[Union[str, Geometry]]],
        materials: Union[str, Material, Tuple[Union[str, Material]]],
        *,
        name: Optional[str] = None,
    ):
        if isinstance(geometries, (str, Geometry)):
            geometries = [geometries]
        else:
            geometries = list(geometries)
        if isinstance(materials, (str, Material)):
            materials = [materials]
        else:
            materials = list(materials)
        for i, geometry in enumerate(geometries):
            if isinstance(geometry, Geometry):
                continue
            if not isinstance(geometry, str):
                raise NotImplementedError
            if geometry not in Geometry.instances:
                raise NotImplementedError
            geometries[i] = Geometry.instances[geometry]
        for i, material in enumerate(materials):
            if isinstance(material, Material):
                continue
            if not isinstance(material, str):
                raise NotImplementedError
            if material not in Material.instances:
                raise NotImplementedError
            materials[i] = Material.instances[material]
        self.name = name
        self.__geometries = tuple(geometries)
        self.__materials = tuple(materials)

    @property
    def geometries(self) -> Iterable[Geometry]:
        """
        Geometric curves labels that defines shapes

        :getter: Returns the curve labels
        :type: Tuple[Tuple[int]]
        """
        return self.__geometries

    @property
    def materials(self) -> Iterable[Material]:
        """
        Used materials for every shape

        :getter: Returns the used materials, in the shapes' order
        :type: Tuple[Material]
        """
        return self.__materials

    @classmethod
    def from_shapes(
        cls,
        shapes: Union[DefinedShape, Tuple[DefinedShape]],
        materials: Union[Material, Tuple[Material]],
    ) -> BaseSection:
        """
        Creates an Section instance based on given shapes and materials
        It's used along the shapepy packaged

        :param name: The section's name, defaults to "custom-section"
        :type name: str, optional
        :return: The dictionary ready to be saved in JSON file
        :rtype: Dict[str, Any]
        """
        if isinstance(shapes, DefinedShape):
            shapes = [shapes]
        if isinstance(materials, Material):
            materials = [materials]
        geome_names = []
        mater_names = []
        for shape, material in zip(shapes, materials):
            if not isinstance(shape, DisjointShape):
                geometrie = Geometry.from_shape(shape)
                geome_names.append(geometrie.name)
                mater_names.append(material.name)
                continue
            for subshape in shape.subshapes:
                geometrie = Geometry.from_shape(subshape)
                geome_names.append(geometrie.name)
                mater_names.append(material.name)
        return cls(geome_names, mater_names)


class GeometricSection(BaseSection):
    """
    GeometricSection's class

    """

    def __init__(
        self,
        geome_names: Tuple[Geometry],
        mater_names: Tuple[Material],
    ):
        super().__init__(geome_names, mater_names)
        self.__geomintegs = None

    def __compute_geomintegs(self):
        """
        Compute the geometric integrals over the domain
        creating the object __geomintegs
        """
        integrals = {}
        all_labels = set()
        for geometry in self.geometries:
            all_labels |= set(map(abs, geometry.labels))
        for label in all_labels:
            curve = Curve.instances[label]
            integrals[label] = Polynomial.adaptative(curve)

        geomintegs = np.zeros(10, dtype="float64")
        for geometry in self.geometries:
            for label in geometry.labels:
                signal = 1 if label > 0 else -1
                geomintegs += signal * integrals[abs(label)]
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


class Section(GeometricSection):
    """
    Section's class

    """

    def torsion_center(self) -> Tuple[float]:
        return (0, 0)

    def torsion_constant(self) -> float:
        return 0

    def shear_center(self) -> Tuple[float]:
        return (0, 0)
