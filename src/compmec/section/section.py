from __future__ import annotations

import json
from importlib import resources
from typing import Optional, Tuple, Union

import jsonschema
import numpy as np
from compmec.shape import JordanCurve, Point2D
from compmec.shape.shape import DefinedShape, IntegrateShape, ShapeFromJordans

from compmec import nurbs

from . import bem2d
from .material import Material


class BaseSection:
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

    def __iter__(self) -> Tuple[DefinedShape, Material]:
        for shape, material in zip(self.__shapes, self.__materials):
            yield (shape, material)

    def __getitem__(self, index) -> Tuple[DefinedShape, Material]:
        index = int(index) % len(self.__shapes)
        return (self.__shapes[index], self.__materials[index])

    def external_jordans(self) -> Tuple[JordanCurve]:
        pass

    def jordans(self) -> Tuple[JordanCurve]:
        pass


class Section:
    """
    Section's class


    """

    AUTO_SOLVE = True

    @classmethod
    def __validate_json(cls, filepath: str):
        schema_path = resources.files("compmec.section")
        schema_path = schema_path.joinpath("schema/section.json")
        with schema_path.open() as file:
            schema = json.load(file)
        with open(filepath, "r") as file:
            data = json.load(file)
        jsonschema.validate(data, schema)

        node_labels = tuple(line[0] for line in data["nodes"])
        if len(node_labels) != len(set(node_labels)):
            msg = "There are nodes with same label number!"
            raise ValueError(msg)
        node_labels = set(node_labels)
        curve_labels = tuple(info["label"] for info in data["curves"])
        if len(curve_labels) != len(set(curve_labels)):
            msg = "There are curves with same label number!"
            raise ValueError(msg)
        curve_labels = set(curve_labels)
        for info in data["curves"]:
            diff = set(info["ctrlpoints"]) - node_labels
            if diff:
                msg = f"The curve {info['label']} refer to ctrlpoints "
                msg += f"{sorted(diff)} which don't exist in the header "
                msg += f"of all node labels: {sorted(node_labels)}"
                raise ValueError(msg)
        for shape_name, labels in data["shapes"].items():
            diff = set(labels) - curve_labels
            if diff:
                msg = f"The shape {shape_name} refer to curves "
                msg += f"{sorted(diff)} which don't exist in the header "
                msg += f"of all curve labels: {sorted(curve_labels)}"
                raise ValueError(msg)
        all_shape_names = set(data["shapes"].keys())
        all_material_names = set(data["materials"].keys())
        for sec_name, info in data["sections"].items():
            shape_names = set(info["shapes"])
            material_names = set(info["materials"])
            diff = shape_names - all_shape_names
            if diff:
                msg = f"The section {sec_name} refer to shapes "
                msg += f"{sorted(diff)} which don't exist in the header "
                msg += f"of all curve labels: {sorted(all_shape_names)}"
                raise ValueError(msg)
            diff = material_names - all_material_names
            if diff:
                msg = f"The section {sec_name} refer to shapes "
                msg += f"{sorted(diff)} which don't exist in the header "
                msg += f"of all curve labels: {sorted(all_shape_names)}"
                raise ValueError(msg)

    @classmethod
    def from_json(cls, filepath: str) -> Union[Section, Tuple[Section]]:
        """
        Creates instances of section reading the data from json file.
        """
        cls.__validate_json(filepath)
        with open(filepath, "r") as file:
            data = json.load(file)

        nodes = {}
        for line in data["nodes"]:
            nodes[line[0]] = Point2D(line[1:])

        all_jordans = {}
        for info in data["curves"]:
            curve_label = info["label"]
            degree = info["degree"] if "degree" in info else None
            knotvector = info["knotvector"]
            knotvector = nurbs.KnotVector(knotvector, degree=degree)
            curve = nurbs.Curve(knotvector)
            points = tuple(nodes[lab] for lab in info["ctrlpoints"])
            curve.ctrlpoints = points
            if "weights" in info:
                curve.weight = info["weights"]
            jordan = JordanCurve.from_full_curve(curve)
            all_jordans[curve_label] = jordan

        all_shapes = {}
        for name, curve_labels in data["shapes"].items():
            jordans = tuple(all_jordans[lab] for lab in curve_labels)
            shape = ShapeFromJordans(jordans)
            all_shapes[name] = shape

        all_materials = {}
        for name, info in data["materials"].items():
            material = Material.from_dict(info)
            all_materials[name] = material

        sections = {}
        for name, info in data.items():
            shapes = tuple(all_shapes[sh] for sh in info["shapes"])
            materials = tuple(all_materials)
            section = Section(shapes, materials)
            sections[name] = section
        return sections

    def __init__(self, shapes: Tuple[DefinedShape], materials: Tuple[Material]):
        self.__area = None
        self.__first = None
        self.__second = None
        self.__warping = None
        self.__strain = StrainField(self)
        self.__stress = StressField(self)
        self.__base = BaseSection(shapes, materials)
        self.__shapes = shapes
        self.__materials = materials

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

        :return: The first moment of inertia (Qx, Qy)
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
            polar = second[0] + second[2]
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
        """Gives the gyradius (radii of gyration)

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

    def solve(self, meshsize: float = None, degree: int = 1):
        shape = self.__shapes[0]
        jordan = shape.jordans[0]
        mesh = bem2d.CurveMesh(jordan)
        vertices = jordan.vertices
        vertices = tuple(tuple(map(np.float64, pt)) for pt in vertices)
        warpvalues = WarpingValues(curve, basis, tsources, tmesh)
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

    def eval_boundary(self, params: Tuple[float]) -> Tuple[float]:
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
    def __eval_bend_moment(self, points: Tuple[Tuple[float]]) -> Tuple[float]:
        """
        Evaluates
        """

    def eval_boundary(self, params: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Returns the values of
        [exx, eyy, ezz, exz, eyz]
        The coordinate exy is always zero
        """
        nvals = len(params)
        result = np.zeros((nvals, 5), dtype="float64")
        if self.forces[2]:  # Axial force
            pass
        if self.moments[0] or self.moments[1]:  # Bending moments
            pass
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion force
            raise NotImplementedError
        return result

    def eval_interior(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        """
        Returns the strain values of
        [exx, eyy, ezz, exz, eyz]
        The coordinate exy is always zero
        """
        nvals = len(points)
        result = np.zeros((nvals, 5), dtype="float64")
        if self.forces[2]:  # Axial force
            pass
        if self.moments[0] or self.moments[1]:  # Bending moments
            pass
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion force
            raise NotImplementedError
        return result


class StressField(ChargedField):
    def __axial_stresses(self, points: Tuple[Tuple[float]]) -> Tuple[float]:
        area = self.section.area()
        szz = self.forces[2] / area
        values = szz * np.ones(len(points), dtype="float64")
        return values

    def __bending_stresses(self, points: Tuple[Tuple[float]]) -> Tuple[float]:
        values = np.zeros(len(points), dtype="float64")
        return values

    def eval_boundary(
        self, params: Tuple[float], material: Material
    ) -> Tuple[Tuple[float]]:
        """
        Returns the strain values of
        [s_xz, s_yz, s_zz]
        The coordinate s_xy is always zero
        """
        nvals = len(params)
        result = np.zeros((nvals, 3), dtype="float64")
        if self.forces[2]:  # Axial force
            result[:, 2] += self.__axial_stresses(points)
        if self.moments[0] or self.moments[1]:  # Bending moments
            result[:, 2] += self.__bending_stresses(points)
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion force
            raise NotImplementedError

    def eval_interior(
        self, points: Tuple[Tuple[float]], subshape: DefinedShape, material: Material
    ) -> Tuple[Tuple[float]]:
        """
        Returns the strain values of
        [sxz, syz, szz]
        The coordinate exy is always zero
        """
        nvals = len(points)
        result = np.zeros((nvals, 3), dtype="float64")
        if self.forces[2]:  # Axial force
            area = self.section.area()
            result[:, 2] += self.forces[2] / area
        if self.moments[0] or self.moments[1]:  # Bending moments
            result[:, 2] += self.__b
        if self.forces[0] or self.forces[1]:  # Shear force
            raise NotImplementedError
        if self.moments[2]:  # Torsion force
            raise NotImplementedError
