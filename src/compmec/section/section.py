from __future__ import annotations

from collections import OrderedDict
from typing import Dict, Tuple, Union

import numpy as np
from compmec.shape.shape import DefinedShape, IntegrateShape

from compmec import nurbs

from . import dataio
from .material import Material


class BaseSection(object):
    @classmethod
    def from_json(cls, filepath: str) -> Dict[str, BaseSection]:
        raise NotImplementedError

    @classmethod
    def from_dict(cls, infos: Dict) -> BaseSection:
        raise NotImplementedError

    def __init__(
        self,
        shapes: Union[DefinedShape, Tuple[DefinedShape]],
        materials: Union[Material, Tuple[Material]],
    ):
        if isinstance(shapes, DefinedShape):
            shapes = [shapes]
        else:
            for shape in shapes:
                if not isinstance(shape, DefinedShape):
                    msg = "shape is not a DefinedShape, but %s"
                    raise TypeError(msg % type(shape))
        if isinstance(materials, Material):
            materials = [materials]
        else:
            for material in materials:
                if not isinstance(material, Material):
                    msg = "material is not a DefinedShape, but %s"
                    raise TypeError(msg % type(material))
        if len(shapes) != len(materials):
            msg = "len of shapes (%d) != (%d) len of materials"
            msg %= len(shapes), len(materials)
            raise ValueError(msg)
        self.__shapes = tuple(shapes)
        self.__materials = tuple(materials)

    def __iter__(self):
        for pair in zip(self.shapes, self.materials):
            yield pair

    def __getitem__(self, key: int):
        return self.shapes[key], self.materials[key]

    @property
    def shapes(self):
        return self.__shapes

    @property
    def materials(self):
        return self.__materials

    def to_dict(self, name: str = "custom-section") -> Dict:
        dicionary = OrderedDict()
        nodesdict = {}
        curves = {}
        shapes = OrderedDict()
        materials = OrderedDict()
        sections = OrderedDict()
        sections[name] = OrderedDict()
        sections[name]["shapes"] = []
        sections[name]["materials"] = []

        for i, (shape, material) in enumerate(self):
            matname = "material-%d" % i
            shaname = "shape-%d" % i
            sections[name]["shapes"].append(shaname)
            sections[name]["materials"].append(matname)

            materials[matname] = OrderedDict()
            materials[matname]["young_modulus"] = material.young_modulus
            materials[matname]["poissons_ratio"] = material.poissons_ratio

            for jordan in shape.jordans:
                curve = jordan.full_curve()
                curvedict = OrderedDict()
                curvedict["degree"] = curve.degree
                curvedict["knotvector"] = tuple(map(float, curve.knotvector))
                curvedict["ctrlpoints"] = []
                for point in curvedict.ctrlpoints:
                    curvedict["ctrlpoints"].append(id(point))
                    nodesdict[id(point)] = tuple(map(float, point))
                curves[id(curve)] = curvedict

        nodes = []
        for label in sorted(nodesdict.keys()):
            x, y = nodesdict[label]
            nodes.append((label, x, y))

        dicionary["nodes"] = nodes
        dicionary["curves"] = curves
        dicionary["shapes"] = shapes
        dicionary["materials"] = materials
        dicionary["sections"] = sections
        return dicionary

    def to_json(self, filepath: str, name: str = "custom-section"):
        dicionary = self.to_dict(name)
        dataio.save_json(dicionary, filepath)


class Section(BaseSection):
    """
    Section's class

    Which is defined by only one DefinedShape and one material

    That means, it's possible only to use homogeneous sections

    For composite sections, use CompositeSection class
    """

    AUTO_SOLVE = True

    def __init__(self, shapes: Tuple[DefinedShape], materials: Tuple[Material]):
        super().__init__(shapes, materials)
        self.__area = None
        self.__first = None
        self.__second = None
        self.__charged_field = ChargedField(self)

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
            self.__area = sum(map(IntegrateShape.area, self.shapes))
        return self.__area

    def first_moment(self) -> Tuple[float]:
        """Gives the first moments of area

        Qx = int y dx dy
        Qy = int x dx dy

        :return: The first moment of inertia (Qx, Qy)
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> from compmec.shape import shapelib.JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = shapelib.JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__first is None:
            Qx = sum(IntegrateShape.polynomial(shape, 0, 1) for shape in self.shapes)
            Qy = sum(IntegrateShape.polynomial(shape, 1, 0) for shape in self.shapes)
            self.__first = (Qx, Qy)
        return tuple(self.__first)

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

        >>> from compmec.shape import shapelib.JordanCurve
        >>> vertices = [(0, 0), (4, 0), (0, 3)]
        >>> jordan = shapelib.JordanCurve.from_vertices(vertices)
        >>> jordan.move((2, 3))
        Jordan Curve of degree 1 and vertices
        ((2, 3), (6, 3), (2, 6))

        """
        if self.__second is None:
            Ixx = sum(IntegrateShape.polynomial(shape, 0, 2) for shape in self.shapes)
            Ixy = sum(IntegrateShape.polynomial(shape, 1, 1) for shape in self.shapes)
            Iyy = sum(IntegrateShape.polynomial(shape, 2, 0) for shape in self.shapes)
            self.__second = (Ixx, Ixy, Iyy)
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

        >>> section = Section(shapes, materials)
        >>> section.torsion_constant()
        1.

        """
        raise NotImplementedError

    def elastic_modulus(self) -> Tuple[Tuple[float]]:
        raise NotImplementedError

    def plastic_modulus(self) -> Tuple[float]:
        raise NotImplementedError

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

        >>> section = Section(shapes, materials)
        >>> section.bending_center()
        (0., 0.)

        """
        return self.geometric_center()

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

    def shear_center(self) -> Tuple[float]:
        """Gives the shear center S

        The shear center is the point which,
        when applied a transverse force, this
        force doesn't cause torsion

        :return: The value of shear center S
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.shear_center()
        (0., 0.)

        """
        raise NotImplementedError

    def plastic_center(self) -> Tuple[float]:
        raise NotImplementedError

    def gyradius(self) -> Tuple[float]:
        """Gives the gyradius (radii of gyration)

        R = (sqrt(Ixx/A), sqrt(Iyy/A))

        """
        area = self.area()
        Ixx, _, Iyy = self.second_moment()
        return np.sqrt(Iyy / area), np.sqrt(Ixx / area)

    def warping_constant(self) -> float:
        """Gives the warping constant"""
        raise NotImplementedError

    def monosymmetry_constants(self) -> Tuple[float]:
        raise NotImplementedError

    def solve(self, meshsize: float = None, degree: int = 1):
        raise NotImplementedError

    def warping(self) -> WarpingField:
        if self.__warping is None:
            if not self.AUTO_SOLVE:
                msg = "You must solve the system before getting warping"
                raise RuntimeError(msg)
            self.solve()
        return self.__warping

    def charged_field(self) -> ChargedField:
        return self.__charged_field


class Field(object):
    def __init__(self, section: Section):
        if not isinstance(section, Section):
            msg = f"section is not a Section, but {type(section)}"
            raise TypeError(msg)
        self.__section = section

    @property
    def section(self) -> Section:
        return self.__section


class WarpingField(Field):
    def __init__(self, section: Section):
        super().__init__(section)
        self.__values = None

    @property
    def values(self):
        return self.__values

    @values.setter
    def values(self, vals: Tuple[Tuple[float]]):
        self.__values = np.array(vals, dtype="float64")


class ChargedField(Field):
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

    def eval(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        """
        Receives a matrix P (list of points) of shape (n, 2)
        and returns a matrix M of shape (n, 8) that contains
        the strain and stress

        M_i = [s_zz, s_xz, s_yz, e_xx, e_yy, e_zz, e_xz, e_yz]

        """
        if not isinstance(points, np.ndarray):
            points = np.array(points, dtype="float64")
        if len(self.section.shapes) != 1:
            raise NotImplementedError
        values = np.zeros((len(points), 8), dtype="float64")

        contains = self.section.shapes[0].contains_point
        mask = tuple(contains(point, True) for point in points)
        mask = np.array(mask, dtype="bool")

        if np.any(mask):
            values[mask, :3] = self.__eval_stresses(points[mask])
        material = self.section.materials[0]
        values[:, 5] = values[:, 0] / material.young_modulus  # e_zz
        values[:, 3] = -material.poissons_ratio * values[:, 1]  # e_xx
        values[:, 4] = values[:, 3]  # e_yy
        values[:, 6] = 0.5 * values[:, 1] / material.shear_modulus  # e_xz
        values[:, 7] = 0.5 * values[:, 2] / material.shear_modulus  # e_yz
        return values

    def __eval_stresses(self, points: np.ndarray) -> np.ndarray:
        """
        Returns the strain values of
        [szz, sxz, syz]

        We suppose that all points are inside the section.
        No verification are made here
        """
        assert len(self.section.shapes) == 1  # For now
        stresses = np.zeros((len(points), 3), dtype="float64")
        stresses[:, 0] = self.__axial_stresses(points)

        contains = self.section.shapes[0].contains_point
        mask = tuple(contains(point, False) for point in points)
        mask = np.array(mask, dtype="bool")

        if np.any(mask):
            inside_points = points[mask]
            inside_shear_stresses = self.__inside_shear_stress(inside_points)
            stresses[mask, 1:] = inside_shear_stresses

        if not np.all(mask):
            bound_points = points[~mask]
            bound_shear_stresses = self.__boundary_shear_stress(bound_points)
            stresses[~mask, 1:] = bound_shear_stresses

        return stresses

    def __axial_stresses(self, points: np.ndarray) -> Tuple[float]:
        """
        Computes the values of s_zz due to normal force and bending moments
        """
        values = np.zeros(len(points), dtype="float64")

        _, _, forcez = self.forces
        momentx, momenty, _ = self.moments

        if forcez:
            values.fill(forcez / self.section.area())
        if momentx or momenty:
            center = self.section.bending_center()
            Ixx, Ixy, Iyy = self.section.second_moment(center)
            xvals, yvals = np.transpose(points)
            xvals -= center[0]
            yvals -= center[1]

            values += (Iyy * momentx + Ixy * momenty) * yvals
            values -= (Ixy * momentx + Ixx * momenty) * xvals
            values /= Ixx * Iyy - Ixy**2
        return values

    def __boundary_shear_stress(
        self, params: Tuple[float], curve: nurbs.Curve
    ) -> np.ndarray:
        """
        Computes the values of (s_xz, s_yz), the shear stresses

        It's caused by two phenomenums: torsion and shear loading


        """
        values = np.zeros((len(params), 2), dtype="float64")
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion moment
            raise NotImplementedError
        return values

    def __inside_shear_stress(self, points: np.ndarray) -> np.ndarray:
        """
        Computes the values of (s_xz, s_yz), the shear stresses

        It's caused by two phenomenums: torsion and shear loading

        We suppose that all points are inside the shape
        """
        values = np.zeros(points.shape, dtype="float64")
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion moment
            raise NotImplementedError
        return values
