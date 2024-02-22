"""
This file contains the principal class of this packaged, called "Section"
It uses all the other files to compute geometric properties, torsional
constant, torsion and shear center and others.


"""
from __future__ import annotations

from collections import OrderedDict
from typing import Tuple

import numpy as np

from .material import Material


class BaseSection:
    """
    BaseSection class that is the base for others section classes

    Basicaly this class is responsible to construct the initial section,
    verifying if the inputs are correct, and generate getter only properties
    such as 'curves', 'geomlabels' and 'materials'.
    """

    instances = OrderedDict()

    def __init__(
        self,
        geomlabels: Tuple[Tuple[int]],
        materials: Tuple[str],
    ):
        self.__geomlabels = tuple(
            tuple(map(int, labels)) for labels in geomlabels
        )
        self.__materials = tuple(materials)

    @property
    def geomlabels(self) -> Tuple[Tuple[int]]:
        """
        Geometric curves labels that defines shapes

        :getter: Returns the curve labels
        :type: Tuple[Tuple[int]]
        """
        return self.__geomlabels

    @property
    def materials(self) -> Tuple[Material]:
        """
        Used materials for every shape

        :getter: Returns the used materials, in the shapes' order
        :type: Tuple[Tuple[int]]
        """
        return self.__materials


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
        raise NotImplementedError

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
        iqx, iqy = self.__geomintegs[1]
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
        ixx, ixy, iyy = self.__geomintegs[2]
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
        ixxx, ixxy, ixyy, iyyy = self.__geomintegs[3]
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
