"""
This file contains classes responsibles to evaluate some
specific values within the section, at any point.

For example, ChargedField is responsible to compute the
stress and strain of the section for every
"""

from typing import Tuple

import numpy as np

from .abcs import IField


class Field(IField):
    """
    This is a base abstract class parent of others

    It's responsible to decide if a given point is
    inside/outside of given section.
    """

    def eval_interior(
        self, points: Tuple[Tuple[float]]
    ) -> Tuple[Tuple[float]]:
        """
        Evaluate the field on points in the interior

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The field value for every point
        :rtype: Tuple[bool]
        """
        raise NotImplementedError

    def eval_boundary(
        self, points: Tuple[Tuple[float]]
    ) -> Tuple[Tuple[float]]:
        """
        Evaluate the field on points at the boundary (internal or external)

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The field value for every point
        :rtype: Tuple[bool]
        """
        raise NotImplementedError

    def on_boundary(self, points: Tuple[Tuple[float]]) -> Tuple[bool]:
        """
        Tells if given points are at any boundary (internal or external)

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The result for every point
        :rtype: Tuple[bool]
        """
        raise NotImplementedError

    def is_inside(self, points: Tuple[Tuple[float]]) -> Tuple[bool]:
        """
        Tells if given points are inside the section or not

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The result for every point
        :rtype: Tuple[bool]
        """
        raise NotImplementedError

    def eval(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        """
        Evaluate the field at given points

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The results in a matrix of shape (n, ndata)
        :rtype: Tuple[Tuple[float]]
        """
        points = np.array(points, dtype="float64")
        results = np.zeros((len(points), self.ndata), dtype="float64")
        mask_inside = self.is_inside(points)
        results[mask_inside] = self.eval_interior(points[mask_inside])
        mask_boundary = self.on_boundary(points)
        results[mask_boundary] = self.eval_boundary(points[mask_boundary])
        return results

    def __call__(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        try:
            return self.eval(points)
        except TypeError:
            return self.eval([points])[0]


class ChargedField(Field):
    """
    A stress and strain field evaluator when
    given forces and momentums are applied

    This evaluator returns the values
    (S13, S23, S33, E33, E13, E23, E33, E11, E22)
    """

    def __init__(self, section):
        self.section = section

    @property
    def ndata(self):
        return 8

    def eval_interior(self, points):
        raise NotImplementedError

    def eval_boundary(self, points):
        raise NotImplementedError
