"""
This file contains the principal class of this packaged, called "Section"
It uses all the other files to compute geometric properties, torsional
constant, torsion and shear center and others.


"""

from __future__ import annotations

from collections import OrderedDict
from typing import Dict, Optional, Tuple, Union

import numpy as np
from compmec.shape.shape import DefinedShape

from .curve import Curve, shapes_to_curves
from .field import ChargedField
from .integral import Polynomial
from .material import Material


class BaseSection:
    """
    BaseSection class that is the base for others section classes

    Basicaly this class is responsible to construct the initial section,
    verifying if the inputs are correct, and generate getter only properties
    such as 'curves', 'geomlabels' and 'materials'.
    """

    instances = OrderedDict()

    @staticmethod
    def clear(names: Optional[Tuple[str]] = None):
        """
        Removes all given instances of Curve
        """
        if names is None:
            BaseSection.instances.clear()
            return
        for name in names:
            if name in BaseSection.instances:
                BaseSection.instances.pop(name)

    @staticmethod
    def __next_available_name() -> str:
        index = 1
        while True:
            name = f"custom-section-{index}"
            if name not in Material.instances:
                return name
            index += 1

    def __init__(
        self,
        geom_labels: Tuple[Tuple[int]],
        mater_names: Tuple[str],
        name: Optional[str] = None,
    ):
        for labels in geom_labels:
            for label in labels:
                assert abs(label) in Curve.instances
        for mat_name in mater_names:
            assert mat_name in Material.instances
        if name is None:
            name = BaseSection.__next_available_name()
        elif name in BaseSection.instances:
            raise ValueError
        self.__geom_labels = tuple(
            tuple(map(int, labels)) for labels in geom_labels
        )
        self.__mater_names = tuple(mater_names)
        self.__name = name
        self.instances[name] = self

    @property
    def name(self) -> str:
        """
        Gives the material name

        :getter: Returns the material's name
        :setter: Attribuates a new name for material
        :type: str

        """
        return self.__name

    @name.setter
    def name(self, new_name: str):
        if self.name == new_name:
            return
        if new_name in self.instances:
            msg = f"Section name '{new_name}' is already used"
            raise ValueError(msg)
        self.instances[new_name] = self.instances.pop(self.name)
        self.__name = new_name

    @property
    def geom_labels(self) -> Tuple[Tuple[int]]:
        """
        Geometric curves labels that defines shapes

        :getter: Returns the curve labels
        :type: Tuple[Tuple[int]]
        """
        return self.__geom_labels

    @property
    def mater_names(self) -> Tuple[str]:
        """
        Material names in each shape

        :getter: Returns the material names
        :type: Tuple[str]
        """
        return self.__mater_names

    @property
    def materials(self) -> Tuple[Material]:
        """
        Used materials for every shape

        :getter: Returns the used materials, in the shapes' order
        :type: Tuple[Material]
        """
        return tuple(Material.instances[name] for name in self.mater_names)

    @classmethod
    def from_shapes(
        cls,
        shapes: Union[DefinedShape, Tuple[DefinedShape]],
        materials=Union[Material, Tuple[Material]],
    ) -> BaseSection:
        """
        Creates an Section instance based on given shapes and materials
        It's used along the compmec-shape packaged

        :param name: The section's name, defaults to "custom-section"
        :type name: str, optional
        :return: The dictionary ready to be saved in JSON file
        :rtype: Dict[str, Any]
        """
        if isinstance(shapes, DefinedShape):
            shapes = [shapes]
        if isinstance(materials, Material):
            materials = [materials]
        geom_labels = shapes_to_curves(shapes)
        mater_names = tuple(material.name for material in materials)
        return cls(geom_labels, mater_names)

    @classmethod
    def from_dict(cls, dictionary: Dict) -> BaseSection:
        """
        Transforms the dictionary, read from json, into a section instance

        :param dictionary: The informations that defines the section
        :type dictionary: Dict
        :return: The created section instance
        :rtype: BaseSection
        """
        geom_labels = dictionary["geom_labels"]
        mater_names = dictionary["materials"]
        return cls(geom_labels, mater_names)


class GeometricSection(BaseSection):
    """
    GeometricSection's class

    """

    def __init__(
        self,
        shapes: Tuple[Tuple[int]],
        materials: Tuple[Material],
    ):
        super().__init__(shapes, materials)
        self.__geomintegs = None

    def __compute_geomintegs(self):
        """
        Compute the geometric integrals over the domain
        creating the object __geomintegs
        """
        integrals = {}
        all_labels = {val for labels in self.geom_labels for val in labels}
        for label in all_labels:
            curve = Curve.instances[label]
            integrals[label] = Polynomial.adaptative(curve)

        geomintegs = np.zeros(10, dtype="float64")
        for labels in self.geom_labels:
            for label in labels:
                signal = 1 if label > 0 else -1
                geomintegs += signal * integrals[abs(label)]
        self.__geomintegs = geomintegs

    def area(self) -> float:
        """
        Gives the cross-section area

        A = int 1 dx dy

        :return: The value of cross-section area
        :rtype: float

        Example use
        -----------

        >>> section.area()
        1.0

        """
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        return self.__geomintegs[0]

    def first_moment(self) -> Tuple[float]:
        """Gives the first moments of area

        Qx = int y dx dy
        Qy = int x dx dy

        :return: The first moment of inertia (Qx, Qy)
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section.first_moment()
        (0., 0.)

        """
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        iqx, iqy = self.__geomintegs[1:3]
        return iqx, iqy

    def second_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        """Gives the second moment of inertia with respect to ``center``

        Ixx = int (y-cy)^2 dx dy
        Ixy = int (x-cx)*(y-cy) dx dy
        Iyy = int (x-cx)^2 dx dy

        If no ``center`` is given, it assumes the origin (0, 0) and
        returns the global second moment of inertia

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The values of Ixx, Ixy, Iyy
        :rtype: tuple[float, float, float]

        Example use
        -----------

        >>> section.second_moment()
        (1., 0., 1.)

        """
        if self.__geomintegs is None:
            self.__compute_geomintegs()
        area = self.area()
        ixx, ixy, iyy = self.__geomintegs[3:6]
        ixx -= area * center[1] ** 2
        ixy -= area * center[0] * center[1]
        iyy -= area * center[0] ** 2
        return ixx, ixy, iyy

    def third_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        """Gives the third moment of inertia with respect to ``center``

        Ixxx = int (y-cy)^3 dx dy
        Ixxy = int (x-cx)*(y-cy)^2 dx dy
        Ixyy = int (x-cx)^2*(y-cy) dx dy
        Iyyy = int (x-cx)^3 dx dy

        If no ``center`` is given, it assumes the origin (0, 0)

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The values of Ixxx, Ixxy, Ixyy, Iyyy
        :rtype: tuple[float, float, float, float]

        Example use
        -----------

        >>> section.second_moment()
        (0., 0., 0., 0.)

        """
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
        """Gives the geometric center G

        G = (x_gc, y_gc)
        x_gc = (1/A) * Qy
        y_gc = (1/A) * Qx

        This center depends only on the geometry,
        not on the material

        :return: The value of geometric center G
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.geometric_center()
        (0., 0.)

        """
        iqx, iqy = self.first_moment()
        area = self.area()
        return iqx / area, iqy / area

    def bending_center(self) -> Tuple[float]:
        """Gives the bendin center B

        The bending center is the point of the
        intersection of two neutral lines, where
        the stress and strain are always zero

        :return: The value of bending center B
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.bending_center()
        (0., 0.)

        """
        return self.geometric_center()

    def gyradius(self) -> Tuple[float]:
        """Gives the gyradius (radii of gyration)

        R = (sqrt(Ixx/A), sqrt(Iyy/A))

        """
        area = self.area()
        ixx, _, iyy = self.second_moment()
        return np.sqrt(iyy / area), np.sqrt(ixx / area)

    def charged_field(self) -> ChargedField:
        """
        Gives the charged field instance to evaluate stresses

        :return: The field evaluator
        :rtype: ChargedField

        """
        return ChargedField(self)


class Section(GeometricSection):
    """
    Section's class

    """

    def torsion_center(self) -> Tuple[float]:
        """Gives the torsion center T

        :return: The value of torsion center T
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.torsion_center()
        (0., 0.)

        """
        raise NotImplementedError

    def torsion_constant(self) -> float:
        """Gives the torsion constant J

        J = Ixx + Iyy - Jw

        Careful: This function solves a linear system

        :return: The value of torsion constant J
        :rtype: float

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.torsion_constant()
        1.

        """
        raise NotImplementedError

    def shear_center(self) -> Tuple[float]:
        """Gives the shear center S

        The shear center is the point which,
        when applied a transverse force, this
        force doesn't cause torsion

        :return: The value of shear center S
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section.shear_center()
        (0., 0.)

        """
        raise NotImplementedError
