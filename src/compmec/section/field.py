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

    def __stress_axial(self, winds):
        sigmazz = self.forces[2] / self.section.area()
        return np.any(winds, axis=1) * sigmazz

    def __stress_bending(self, points, winds):
        bend_center = self.section.bending_center()
        ixx, ixy, iyy = self.section.second_moment(bend_center)
        detii = ixx * iyy - ixy**2
        matrix = np.array([[iyy, -ixy], [-ixy, ixx]])
        momx, momy, _ = self.momentums
        vector = np.dot(matrix, [-momy, momx]) / detii
        return np.dot(points, vector) * np.any(winds, axis=1)

    def eval(self, points):
        points = np.array(points, dtype="float64")
        results = np.zeros((len(points), self.ndata), dtype="float64")
        geometries = self.section.geometries
        winds = tuple(
            tuple(map(geome.winding, points)) for geome in geometries
        )
        winds = np.array(winds, dtype="float64")
        forx, fory, forz, momx, momy, momz = self.__charges
        if forz:  # Axial force
            results[:, 2] += self.__stress_axial(winds)
        if momx or momy:  # Bending moments
            results[:, 2] += self.__stress_bending(points, winds)
        if forx or fory or momz:
            raise NotImplementedError
        return results
