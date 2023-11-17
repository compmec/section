from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from compmec.shape.shape import DefinedShape, IntegrateShape

from .bem2d import TorsionVector, WarpingEvaluator, WarpingValues
from .material import Material


class Section:
    """
    Section's class


    """

    AUTO_SOLVE = True

    def __init__(self, shapes: Tuple[DefinedShape], materials: Tuple[Material]):
        for shape in shapes:
            if not isinstance(shape, DefinedShape):
                raise TypeError
        for material in materials:
            if not isinstance(material, Material):
                raise TypeError
        self.__shapes = tuple(shapes)
        self.__materials = tuple(materials)
        if len(self.__shapes) != 1 or len(self.__materials) != 1:
            raise ValueError
        self.__area = None
        self.__first = None
        self.__second = None
        self.__warping = None
        self.__strain = StrainField(self)
        self.__stress = StressField(self)

    def __iter__(self) -> Tuple[DefinedShape, Material]:
        for shape, material in zip(self.__shapes, self.__materials):
            yield (shape, material)

    def __getitem__(self, index) -> Tuple[DefinedShape, Material]:
        index = int(index) % len(self.__shapes)
        return (self.__shapes[index], self.__materials[index])

    def area(self) -> float:
        """Gives the cross-section area

        A = int 1 dx dy

        :return: The value of cross-section area
        :rtype: float

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__area is None:
            areas = tuple(map(IntegrateShape.area, self.__shapes))
            self.__area = sum(areas)
        return self.__area

    def first_moment(self) -> Tuple[float]:
        """Gives the first moments of area

        Qx = int y dx dy
        Qy = int x dx dy

        :return: The first moment of inertia Qx, Qy
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__first is None:
            values = np.zeros(2, dtype="float64")
            for shape in self.__shapes:
                values[0] = IntegrateShape.polynomial(shape, 0, 1)
                values[1] = IntegrateShape.polynomial(shape, 1, 0)
            self.__first = values
        return tuple(self.__first)

    def second_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[Tuple[float]]:
        """Gives the second moment of inertia with respect to center

        Ixx = int (y-cy)^2 dx dy
        Ixy = int (x-cx)*(y-cy) dx dy
        Iyy = int (x-cx)^2 dx dy

        If no center is given, it assumes the origin (0, 0) and
        returns the global second moment of inertia

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The values of Ixx, Ixy, Iyy
        :rtype: tuple[float, float, float]

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__second is None:
            values = np.zeros(3, dtype="float64")
            for shape in self.__shapes:
                values[0] += IntegrateShape.polynomial(shape, 0, 2)
                values[1] += IntegrateShape.polynomial(shape, 1, 1)
                values[2] += IntegrateShape.polynomial(shape, 2, 0)
            self.__second = values
        area = self.area()
        Ixx, Ixy, Iyy = self.__second
        Ixx -= area * center[1] ** 2
        Ixy -= area * center[0] * center[1]
        Iyy -= area * center[0] ** 2
        return Ixx, Ixy, Iyy

    def torsion_constant(self) -> float:
        """Gives the torsion constant J

        J = Ixx + Iyy - Jw

        Careful: This function solves a linear system

        :return: The value of torsion constant J
        :rtype: float

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__warping is None:
            self.solve()
            shape = self.__shapes[0]
            curve = shape.jordans[0]
            vertices = curve.vertices
            vector = TorsionVector(vertices)
            second = self.second_moment()
            polar = second[0, 0] + second[1, 1]
            self.__torsion = polar - np.inner(vector, self.__warping)
        return self.__torsion

    def elastic_modulus(self) -> Tuple[Tuple[float]]:
        raise NotImplementedError

    def plastic_modulus(self) -> Tuple[float]:
        raise NotImplementedError

    def geometric_center(self) -> Tuple[float]:
        """Gives the geometric center G

        G = (x_gc, y_gc)
        x_gc = (1/A) * Qy
        y_gc = (1/A) * Qx

        :return: The value of geometric center G
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        Qx, Qy = self.first_moment()
        area = self.area()
        return Qy / area, Qx / area

    def bending_center(self) -> Tuple[float]:
        """Gives the bendin center B

        The bending center is the point of the
        intersection of two neutral lines, where
        the stress and strain are always zero

        :return: The value of bending center B
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        return self.geometric_center()

    def torsion_center(self) -> Tuple[float]:
        """Gives the torsion center T

        :return: The value of torsion center T
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

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

        >>> from compmec.shape import JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        raise NotImplementedError

    def plastic_center(self) -> Tuple[float]:
        raise NotImplementedError

    def gyradius(self) -> Tuple[float]:
        """Gives the gyradius, or radii of gyration

        R = (sqrt(Ixx/A), sqrt(Iyy/A))

        """
        area = self.area()
        second = self.second_moment()
        return np.sqrt(np.diag(second) / area)

    def warping_constant(self) -> float:
        """Gives the warping constant"""
        raise NotImplementedError

    def monosymmetry_constants(self) -> Tuple[float]:
        raise NotImplementedError

    def solve(self):
        if self.__warping is None:
            shape = self.__shapes[0]
            curve = shape.jordans[0]
            vertices = curve.vertices
            vertices = tuple(tuple(map(np.float64, pt)) for pt in vertices)
            warpvalues = WarpingValues(vertices)
            self.__warping = WarpingField(self, warpvalues)

    def warping(self) -> WarpingField:
        if self.__warping is None:
            if not self.AUTO_SOLVE:
                msg = "You must solve the system before getting warping"
                raise RuntimeError(msg)
            self.solve()
        return self.__warping

    def strain(self) -> StrainField:
        return self.__strain

    def stress(self) -> StressField:
        return self.__stress


class FieldEvaluator(object):
    def __init__(self, section: Section):
        self.section = section

    @property
    def section(self) -> Section:
        return self.__section

    @section.setter
    def section(self, section: Section):
        if not isinstance(section, Section):
            msg = f"section must be a Section instance, type is {type(section)}"
            raise TypeError(msg)
        self.__section = section

    def __call__(self, points: Tuple[Tuple[float]]) -> Tuple[float]:
        values = []
        shape, _ = self.section[0]
        curve = shape.jordans[0]
        for point in points:
            param = projection(point, curve)
            qoint = curve(param)
            vector = point - qoint
            dist2 = vector[0] ** 2 + vector[1] ** 2
            if dist2 < 1e-9:
                value = self.eval_boundary(param)
            elif point in shape:
                value = self.eval_interior(point)
            else:
                value = 0
            values.append(value)
        return np.array(values, dtype="float64")


class WarpingField(FieldEvaluator):
    def __init__(self, section: Section, warpvalues: Tuple[float]):
        super().__init__(section)
        self.__ctrlpoints = np.array(warpvalues, dtype="float64")

    def eval_boundary(self, param: float) -> float:
        raise NotImplementedError

    def eval_interior(self, points: Tuple[Tuple[float]]) -> Tuple[float]:
        shape, _ = self.section[0]
        curve = shape.jordans[0]
        vertices = curve.vertices
        warpvals = self.ctrlpoints
        values = WarpingEvaluator.interior(vertices, points, warpvals)
        return values

    @property
    def ctrlpoints(self) -> Tuple[float]:
        return self.__ctrlpoints


class ChargedField(FieldEvaluator):
    def __init__(self, section: Section):
        super().__init__(section)
        self.__charges = np.zeros(6, dtype="float64")

    @property
    def forces(self) -> Tuple[float]:
        return self.__charges[:3]

    @property
    def moments(self) -> Tuple[float]:
        return self.__charges[3:]

    @forces.setter
    def forces(self, forces: Tuple[float]):
        forces = np.array(forces, dtype="float64")
        if forces.shape != (3,):
            msg = "Forces must be a vector of lenght 3"
            raise ValueError(msg)
        if forces[0] != 0 or forces[1] != 0:
            msg = "Shear is not yet supported"
            raise NotImplementedError(msg)
        self.__charges[:3] = forces

    @moments.setter
    def moments(self, moments: Tuple[float]):
        moments = np.array(moments, dtype="float64")
        if moments.shape != (3,):
            msg = "Moments must be a vector of lenght 3"
            raise ValueError(msg)
        if moments[2] != 0:
            msg = "Torsion is not yet supported"
            raise NotImplementedError(msg)
        self.__charges[3:] = moments


class StrainField(ChargedField):
    def eval_boundary(self, param: float):
        """
        Returns the values of
        [exx, eyy, ezz, exz, eyz]
        The coordinate exy is always zero
        """
        pass

    def eval_interior(self, point: Tuple[float]):
        """
        Returns the strain values of
        [exx, eyy, ezz, exz, eyz]
        The coordinate exy is always zero
        """
        pass


class StressField(ChargedField):
    def eval_boundary(self, param: float):
        """
        Returns the strain values of
        [sxz, syz, szz]
        The coordinate exy is always zero
        """
        value = 0
        if self.forces[2] != 0:
            area = self.section.area()
            value += self.forces[2] / area

    def eval_interior(self, point: Tuple[float]):
        """
        Returns the strain values of
        [sxz, syz, szz]
        The coordinate exy is always zero
        """
        pass
