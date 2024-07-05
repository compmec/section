"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""

from __future__ import annotations

from collections import OrderedDict
from typing import Optional, Tuple

import numpy as np
import pynurbs
from compmec.shape import JordanCurve

from . import integral
from .abcs import ICurve, LabeledTracker


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


class Curve(LabeledTracker, ICurve):
    """
    Base class that tracks the instances
    """

    instances = OrderedDict()

    @classmethod
    def from_jordan(cls, jordan: JordanCurve) -> Curve:
        """
        Converts a jordan curve into a Curve instance

        :param jordan: A jordan curve from compmec-shape packaged
        :type jordan: compmec.shape.JordanCurve
        :return: A Curve instance
        :rtype: Curve
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
        return cls(curve.knotvector, ctrlpoints, weights = curve.weights)

    @classmethod
    def from_vertices(cls, vertices: Tuple[Tuple[float]]) -> Curve:
        """
        Creates a Curve instance based on given vertices

        :param vertices: The polygon vertices
        :type vertices: tuple[tuple[float]]
        :return: A Curve instance
        :rtype: Curve
        """
        npts = len(vertices)
        ctrlpoints = np.array(vertices + [vertices[0]], dtype="float64")
        knotvector = [0] + list(range(npts + 1)) + [npts]
        knotvector = pynurbs.KnotVector(knotvector)
        return cls(knotvector, ctrlpoints)

    def __init__(
        self,
        knotvector: Tuple[float],
        ctrlpoints: Tuple[Tuple[float]],
        *,
        weights: Tuple[float] = None,
        label: Optional[int] = None,
    ):
        self.__internal = pynurbs.Curve(knotvector, ctrlpoints, weights)
        self.__dinternal = pynurbs.Derivate(self.__internal)
        self.label = label

    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        values = self.__internal.eval(parameters)
        return np.array(values, dtype="float64")

    def deval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        values = self.__dinternal.eval(parameters)
        return np.array(values, dtype="float64")

    @property
    def knots(self) -> Tuple[float]:
        return self.__internal.knotvector.knots

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
