"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""

from __future__ import annotations

from collections import OrderedDict
from typing import Dict, Optional, Tuple

import numpy as np
from shapepy import JordanCurve

import pynurbs

from . import integral
from .abcs import ICurve, ISegment, LabeledTracker


class Node(LabeledTracker):
    """
    Class that stores all nodes
    """

    instances = {}

    def __new__(cls, label: int) -> Tuple[float]:
        return Node.instances[label]

    @staticmethod
    def insert_matrix(matrix: Tuple[Tuple[int, float, float]]):
        """
        Inserts the values of a matrix inside the
        'labels', 'xcoords' and 'ycoords'

        :param matrix: The matrix with nodes coordinates
        :type matrix: Tuple[Tuple[int, float, float]]

        Example
        -------
        >>> matrix = [[1, 0.0, 0.0],
                      [2, 1.0, 0.0],
                      [3, 0.0, 1.0]]
        >>> Node.insert_matrix(matrix)

        """
        for line in matrix:
            label = int(line[0])
            assert label not in Node.instances
            point = tuple(map(float, line[1:3]))
            Node.instances[label] = point

    @staticmethod
    def from_labels(labels: Tuple[int]) -> Tuple[Tuple[float]]:
        """
        Gives the coordinates of the points

        :param labels: The desired node labels
        :type labels: Tuple[int]
        :return: A matrix of shape (n, 2)
        :rtype: Tuple[Tuple[float]]

        Example
        -------
        >>> Node.from_labels([1, 2, 3])
        ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

        """
        return tuple(Node.instances[label] for label in labels)


class LinearSegment(ISegment):
    """
    Defines a Linear segment, that connects two points.

    The parametric space is defined as [0, 1]

    """

    def __init__(self, pointa: Tuple[float], pointb: Tuple[float]):
        self.pointa = np.array(pointa, dtype="float64")
        self.pointb = np.array(pointb, dtype="float64")

    @property
    def vector(self) -> Tuple[float]:
        """
        Gives the directional vector V = B - A

        :getter: Returns the vector V = B - A
        :type: Tuple[float]
        """
        return self.pointb - self.pointa

    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        return self.pointa + np.tensordot(parameters, self.vector, axes=0)

    def projection(self, point: Tuple[float]) -> float:
        vector = self.vector
        parameter = np.inner(point - self.pointa, vector)
        parameter /= np.inner(vector, vector)
        return max(0, min(1, parameter))


class Curve(LabeledTracker, ICurve):
    """
    Base class that tracks the instances
    """

    instances = OrderedDict()

    @staticmethod
    def new_instance(tipo: str, dictionary: Dict) -> Curve:
        """
        Creates a new instance of Curve depending on
        the 'tipo' and the informations from the dictionary

        :param tipo: The curve-subclass to be called, in ["nurbs"]
        :type tipo: str
        :return: The created curve instance
        :rtype: Curve
        """
        tipo2class = {"nurbs": NurbsCurve}
        if tipo not in tipo2class:
            raise NotImplementedError
        return tipo2class[tipo].from_dict(dictionary)


class NurbsCurve(Curve):
    """
    Nurbs Curve instance, a child of "Curve" that serves as interface
    and uses the functions/methods from compmec-nurbs package
    """

    @classmethod
    def from_jordan(cls, jordan: JordanCurve) -> NurbsCurve:
        """
        Converts a jordan curve into a NurbsCurve instance

        :param jordan: A jordan curve from shapepy packaged
        :type jordan: shapepy.JordanCurve
        :return: A NurbsCurve instance
        :rtype: NurbsCurve
        """
        bezier_curves = []
        for i, segment in enumerate(jordan.segments):
            knotvector = pynurbs.GeneratorKnotVector.bezier(segment.degree)
            knotvector.shift(i)
            new_bezier = pynurbs.Curve(knotvector)
            new_bezier.ctrlpoints = segment.ctrlpoints
            bezier_curves.append(new_bezier)

        curve = bezier_curves[0]
        for i in range(1, len(bezier_curves)):
            bezier = bezier_curves[i]
            curve |= bezier
        ctrlpoints = tuple(
            np.array(tuple(point), dtype="float64")
            for point in curve.ctrlpoints
        )
        curve.ctrlpoints = ctrlpoints
        return cls(curve)

    @classmethod
    def from_dict(cls, dictionary: Dict) -> NurbsCurve:
        degree = dictionary["degree"] if "degree" in dictionary else None
        knotvector = pynurbs.KnotVector(dictionary["knotvector"], degree)
        nurbs_curve = pynurbs.Curve(knotvector)
        if "ctrllabels" in dictionary:
            labels = tuple(dictionary["ctrllabels"])
            ctrlpoints = Node.from_labels(labels)
        else:
            ctrlpoints = dictionary["ctrlpoints"]
        nurbs_curve.ctrlpoints = np.array(ctrlpoints, dtype="float64")
        if "weights" in dictionary:
            nurbs_curve.weights = dictionary["weights"]
        return cls(nurbs_curve)

    def to_dict(self) -> Dict:
        dictionary = OrderedDict()
        dictionary["degree"] = self.internal.degree
        dictionary["knotvector"] = self.internal.knotvector
        dictionary["ctrlpoints"] = self.internal.ctrlpoints
        weights = self.internal.weights
        if weights is not None:
            dictionary["weights"] = weights
        return dictionary

    def __new__(cls, nurbs_curve: pynurbs.Curve, label: Optional[int] = None):
        if not isinstance(nurbs_curve, pynurbs.Curve):
            msg = "Invalid internal curve"
            raise TypeError(msg)
        instance = super().__new__(cls)
        instance.label = label
        instance.internal = nurbs_curve
        return instance

    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        values = self.internal.eval(parameters)
        return np.array(values, dtype="float64")

    @property
    def knots(self) -> Tuple[float]:
        return self.internal.knotvector.knots

    @property
    def internal(self) -> pynurbs.Curve:
        """
        Gives the internal pynurbs.Curve object

        :getter: Returns the curve's label
        :setter: Sets the new curve instance
        :type: compmec.pynurbs.Curve
        """
        return self.__internal

    @internal.setter
    def internal(self, new_curve: pynurbs.Curve):
        if not isinstance(new_curve, pynurbs.Curve):
            msg = "Invalid internal curve"
            raise TypeError(msg)
        self.__internal = new_curve

    def winding(self, point: Tuple[float]) -> float:
        # Verify if the point is at any vertex
        vertices = self.eval(self.knots[:-1])
        for i, vertex in enumerate(vertices):
            if np.all(point == vertex):
                vec_left = vertices[(i - 1) % len(vertices)] - point
                vec_righ = vertices[(i + 1) % len(vertices)] - point
                wind = 0.5 * np.arccos(np.inner(vec_left, vec_righ)) / np.pi
                return wind

        wind = 0
        for vertexa, vertexb in zip(vertices, np.roll(vertices, -1, axis=0)):
            sub_wind = integral.winding_number_linear(vertexa, vertexb, point)
            if abs(sub_wind) == 0.5:
                wind = 0.5
                break
            wind += sub_wind
        return wind

    def projection(self, point: Tuple[float]) -> Tuple[float]:
        return pynurbs.advanced.Projection.point_on_curve(point, self.internal)


class PolygonCurve(Curve):
    """
    A Curve that defines a polygon
    """

    def __init__(
        self, vertices: Tuple[Tuple[float]], label: Optional[int] = None
    ):
        self.label = label
        self.vertices = np.array(vertices)
        if np.all(self.vertices[0] == self.vertices[-1]):
            raise ValueError

    @property
    def knots(self) -> Tuple[float]:
        return np.arange(0, len(self.vertices) + 1, 1.0)

    @classmethod
    def from_dict(cls, dictionary: Dict) -> ICurve:
        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        return cls(vertices)

    def to_dict(self) -> Dict:
        return {}

    def eval(self, parameters):
        nvertices = len(self.vertices)
        points = []
        for parameter in parameters:
            lower = int(np.floor(parameter))
            middle = parameter - lower
            new_point = (1 - middle) * self.vertices[lower % nvertices]
            new_point += middle * self.vertices[(lower + 1) % nvertices]
            points.append(new_point)
        points = np.array(points, dtype="float64")
        return points

    def winding(self, point: Tuple[float]) -> float:
        # Verify if the point is at any vertex
        vertices = self.vertices
        nverts = len(vertices)
        for i, vertex in enumerate(vertices):
            if np.all(point == vertex):
                verta = vertices[(i - 1) % nverts]
                vertb = vertices[(i + 1) % nverts]
                wind = integral.winding_number_linear(vertb, verta, point)
                return wind

        wind = 0
        for i, vertexa in enumerate(vertices):
            vertexb = vertices[(i + 1) % nverts]
            sub_wind = integral.winding_number_linear(vertexa, vertexb, point)
            if abs(sub_wind) == 0.5:
                wind = 0.5
                break
            wind += sub_wind
        return wind

    def projection(self, point: Tuple[float]) -> Tuple[float]:
        point = np.array(point, dtype="float64")
        proj_params = []
        distances = []
        nverts = len(self.vertices)
        for i, vertexa in enumerate(self.vertices):
            vertexb = self.vertices[(i + 1) % nverts]
            segment = LinearSegment(vertexa, vertexb)
            temp_params = segment.projection(point)
            proj_point = segment.eval(temp_params)[0]
            distance = np.linalg.norm(point - proj_point)
            proj_params.append(i + temp_params[0])
            distances.append(distance)

        min_distance = min(distances)
        params = tuple(
            param
            for param, dist in zip(proj_params, distances)
            if dist == min_distance
        )
        return params
