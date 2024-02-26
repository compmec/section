"""
This file contains classes responsibles to evaluate some
specific values within the section, at any point.

For example, ChargedField is responsible to compute the
stress and strain of the section for every
"""

from abc import ABC, abstractmethod
from typing import Tuple


class Field(ABC):
    """
    This is a base abstract class parent of others

    It's responsible to decide if a given point is
    inside/outside of given section.
    """

    @abstractmethod
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

    @abstractmethod
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

    @property
    @abstractmethod
    def ndata(self) -> int:
        """
        Gives the quantity of output numbers for one point

        :getter: Returns the lenght of output data
        :type: int
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
        raise NotImplementedError

    def __call__(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        try:
            return self.eval(points)
        except TypeError:
            return self.eval([points])[0]
