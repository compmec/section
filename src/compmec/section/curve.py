"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Dict, Optional, Tuple

import numpy as np
from compmec.shape import JordanCurve
from compmec.shape.shape import DefinedShape

from compmec import nurbs

from . import integral
from .abcs import LabeledTracker


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


class ICurve(ABC):
    """
    Interface abstract curve to be parent of others.

    This class serves as interface between the curves from others packaged
    like nurbs.Curve, to this package, expliciting the minimum requirements
    of a curve must have.
    With this, it's possible to implement your own type of parametric curve
    """

    @classmethod
    @abstractmethod
    def from_dict(cls, dictionary: Dict) -> Curve:
        """
        Gives a curve instance based on the given parameters
        """
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> Dict:
        """
        Transforms a curve instance into a dictionary
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def limits(self) -> Tuple[float]:
        """
        Gives the curve's parametric interval

        :getter: Returns the pair [a, b] in which curve is parametric defined
        :type: Tuple[float]

        """
        raise NotImplementedError

    @property
    @abstractmethod
    def knots(self) -> Tuple[float]:
        """
        Gives the curve's knots, in which the parametric interval is divided

        :getter: Returns the knots that divides the curve's interval
        :type: Tuple[float]

        """
        raise NotImplementedError

    @abstractmethod
    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Evaluates the curves at given parameters.

        :param parameters: A vector-like of lenght n
        :type parameters: Tuple[float]
        :return: A matrix of shape (n, 2)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError

    @abstractmethod
    def winding(self, point: Tuple[float]) -> float:
        """
        Computes the winding number of the given curve

        The possible results are one in the interval [0, 1]

        * 0 means the point is outside the internal curve
        * 1 means the point is inside the internal curve
        * in (0, 1) means it's on the boundary

        :param point: A 2D-point
        :type point: Tuple[float]
        :return: The winding value with respect to the curve
        :rtype: float

        """
        raise NotImplementedError


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

        :param jordan: A jordan curve from compmec-shape packaged
        :type jordan: compmec.shape.JordanCurve
        :return: A NurbsCurve instance
        :rtype: NurbsCurve
        """
        bezier_curves = []
        for i, segment in enumerate(jordan.segments):
            knotvector = nurbs.GeneratorKnotVector.bezier(segment.degree)
            knotvector.shift(i)
            new_bezier = nurbs.Curve(knotvector)
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
        knotvector = nurbs.KnotVector(dictionary["knotvector"], degree)
        nurbs_curve = nurbs.Curve(knotvector)
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

    def __new__(cls, nurbs_curve: nurbs.Curve, label: Optional[int] = None):
        if not isinstance(nurbs_curve, nurbs.Curve):
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
    def limits(self) -> Tuple[float]:
        return self.internal.knotvector.limits

    @property
    def internal(self) -> nurbs.Curve:
        """
        Gives the internal nurbs.Curve object

        :getter: Returns the curve's label
        :setter: Sets the new curve instance
        :type: compmec.nurbs.Curve
        """
        return self.__internal

    @internal.setter
    def internal(self, new_curve: nurbs.Curve):
        if not isinstance(new_curve, nurbs.Curve):
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


def shapes_to_curves(shapes: Tuple[DefinedShape]) -> Tuple[Tuple[int]]:
    """
    Transform shapes instances into the pair (curves, geom_labels)
    curves contains all the parametric curves used to define the shapes
    geom_labels relates which curves are used to describe each shape

    :param shapes: The group of shapes used in section
    :type shapes: Tuple[DefinedShape]
    :return: The pair (curves, geom_labels)
    :rtype: Tuple[Tuple[int]]
    """
    geom_labels = []
    for shape in shapes:
        assert isinstance(shape, DefinedShape)
        new_geom_labels = []
        for jordan in shape.jordans:
            signal = 1 if float(jordan) > 0 else -1
            if signal < 0:
                jordan = ~jordan
            curve = NurbsCurve.from_jordan(jordan)
            new_geom_labels.append(signal * curve.label)
        geom_labels.append(new_geom_labels)
    return geom_labels
