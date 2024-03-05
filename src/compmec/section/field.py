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
            stress, strain = self.eval([points])
            return stress[0], strain[0]


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
    def forces(self, new_forces: Tuple[float]):
        self.__charges[:3] = new_forces

    @momentums.setter
    def momentums(self, new_momentums: Tuple[float]):
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

    def __winding_numbers(
        self, points: Tuple[Tuple[float]]
    ) -> Tuple[Tuple[float]]:
        """
        Computes the winding number of every point,
        for every geometry

        """
        geometries = tuple(self.section.geometries)
        winds = np.zeros((len(points), len(geometries)), dtype="float64")
        for i, point in enumerate(points):
            for j, geome in enumerate(geometries):
                winds[i, j] = geome.winding(point)
        return winds

    def __stress_eval(self, points, winds):
        """
        Evaluates the stress values

        The stress tensor is given by
            [  0    0  E13]
        S = [  0    0  E23]
            [E13  E23  E33]

        Returned values are a matrix of shape (n, 3)
        each line are the stress components: E13, E23, E33

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The strain matrix of shape (n, 3)
        :rtype: Tuple[Tuple[float]]

        """
        results = np.zeros((len(points), 3), dtype="float64")
        if self.forces[2]:  # Axial force
            results[:, 2] += self.__stress_axial(winds)
        if np.any(self.momentums[:2]):  # Bending moments
            results[:, 2] += self.__stress_bending(points, winds)
        if np.any(self.forces[:2]):  # Shear force
            raise NotImplementedError
        if self.momentums[2]:  # Torsion
            raise NotImplementedError
        return results

    def __strain_eval(self, winds, stresses):
        """
        Evaluates the strain values from stress values by
        using Hook's law for isotropic materials

        The winds serves to know the position of the points,
        to decide which material will be used

        The strain tensor is given by
            [E11    0  E13]
        E = [  0  E22  E23]
            [E13  E23  E33]
        The values E22 and E11 are the same

        Returned values are a matrix of shape (n, 4)
        each line are the strain components: E11, E33, E13, E23

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The strain matrix of shape (n, 4)
        :rtype: Tuple[Tuple[float]]

        """
        strain = np.zeros((len(winds), 4), dtype="float64")
        for i, material in enumerate(self.section.materials):
            subwinds = winds[:, i]
            mask = subwinds != 0
            young = material.young_modulus
            poiss = material.poissons_ratio
            fract11 = -poiss / young
            fract13 = (1 + poiss) / young
            strain[mask, 1] += subwinds[mask] * stresses[mask, 2] / young
            strain[mask, 0] += fract11 * subwinds[mask] * stresses[mask, 2]
            strain[mask, 2] += fract13 * subwinds[mask] * stresses[mask, 0]
            strain[mask, 3] += fract13 * subwinds[mask] * stresses[mask, 1]
        return strain

    def eval(self, points):
        """
        Evaluate the stress and strain at given points

        The inputs/outputs are matrices:
        * points is a (n, 2) matrix
        * stress is a (n, 3) matrix
        * strain is a (n, 4) matrix

        The stress components are S13, S23, S33, meaning
        two shear stresses and one normal stress
        The strain components are E11, E33, E13, E23

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The pair (stress, strain)
        :rtype: Tuple[Tuple[Tuple[float]]]
        """
        points = np.array(points, dtype="float64")
        winds = self.__winding_numbers(points)
        stress = self.__stress_eval(points, winds)
        strain = self.__strain_eval(winds, stress)
        return stress, strain
