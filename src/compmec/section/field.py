"""
This file contains classes responsibles to evaluate some
specific values within the section, at any point.

For example, ChargedField is responsible to compute the
stress and strain of the section for every
"""

from typing import Tuple

import numpy as np

from .abcs import IField, ISection


class Field(IField):
    """
    This is a base abstract class parent of others

    It's responsible to decide if a given point is
    inside/outside of given section.
    """

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

    def __init__(self, section: ISection):
        self.section = section
        self.__charges = np.zeros(6, dtype="float64")

    @property
    def ndata(self):
        return 8

    @property
    def forces(self):
        """
        Forces used in the charged field.

        * Fx: shear force
        * Fy: shear force
        * Fz: axial force

        :getter: Returns the (Fx, Fy, Fz)
        :setter: Sets the new forces
        :type: Tuple[float]
        """
        return self.__charges[:3]

    @property
    def momentums(self):
        """
        Forces used in the charged field.

        * Mx: bending momentum
        * My: bending momentum
        * Mz: torsion momentum

        :getter: Returns the (Mx, My, Mz)
        :setter: Sets the new momentums
        :type: Tuple[float]
        """
        return self.__charges[3:]

    @forces.setter
    def forces(self, new_forces: tuple[float]):
        self.__charges[:3] = new_forces

    @momentums.setter
    def momentums(self, new_momentums: tuple[float]):
        self.__charges[3:] = new_momentums

    def eval(self, points):
        raise NotImplementedError
