"""
This file contains classes responsibles to evaluate some
specific values within the section, at any point.

For example, ChargedField is responsible to compute the
stress and strain of the section for every
"""

from typing import Optional, Tuple, Union

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from .abcs import IField, ISection


class Field(IField):
    """
    This is a base abstract class parent of others

    It's responsible to decide if a given point is
    inside/outside of given section.
    """

    def __init__(self, section: ISection):
        self.section = section


class ChargedField(Field):
    """
    A stress and strain field evaluator when
    given forces and momentums are applied

    This evaluator returns the values
    (S13, S23, S33, E33, E13, E23, E33, E11, E22)
    """

    def __init__(self, section: ISection):
        super().__init__(section)
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
        homosections = tuple(self.section)
        winds = np.zeros((len(points), len(homosections)), dtype="float64")
        for i, point in enumerate(points):
            for j, homosection in enumerate(homosections):
                winds[i, j] = homosection.geometry.winding(point)
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
        for i, homosection in enumerate(self.section):
            subwinds = winds[:, i]
            mask = subwinds != 0
            young = homosection.material.young_modulus
            poiss = homosection.material.poissons_ratio
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

        Example
        -------

        >>> points = [(0, 0), (1, 0)]
        >>> stress, strain = field.eval(points)

        """
        points = np.array(points, dtype="float64")
        winds = self.__winding_numbers(points)
        stress = self.__stress_eval(points, winds)
        strain = self.__strain_eval(winds, stress)
        return stress, strain


def field_scalar_values(
    field: IField, points: Tuple[Tuple[float]], subfield_name: Union[str, None]
):
    if not isinstance(field, ChargedField):
        raise NotImplementedError
    if subfield_name is None:
        raise NotImplementedError
    if subfield_name in ("Sxx", "Syy", "Sxy", "Exy"):
        return (0,) * len(points)

    stress, strain = field.eval(points)
    if subfield_name == "Sxz":
        return stress[:, 0]
    if subfield_name == "Syz":
        return stress[:, 1]
    if subfield_name == "Szz":
        return stress[:, 2]
    if subfield_name in ("Exx", "Eyy"):
        return strain[:, 0]
    if subfield_name == "Ezz":
        return strain[:, 1]
    if subfield_name == "Exz":
        return strain[:, 2]
    if subfield_name == "Eyz":
        return strain[:, 3]
    Sxz = stress[:, 0]
    Syz = stress[:, 1]
    Szz = stress[:, 2]
    if subfield_name == "VM":
        return np.sqrt(Szz**2 + 3 * (Sxz**2 + Syz**2))
    if subfield_name == "TR":
        return np.sqrt((0.5 * Szz) ** 2 + Sxz**2 + Syz**2)
    raise NotImplementedError(f"{subfield_name} is unknown")


def plot_section(section: ISection, *, axes: Optional[plt.Axes] = None):
    if not isinstance(section, ISection):
        raise NotImplementedError
    if axes is None:
        axes = plt.gca()
    bending_center = section.bending_center()
    geometric_center = section.geometric_center()
    torsion_center = section.torsion_center()
    shear_center = section.shear_center()
    axes.scatter(
        geometric_center[0], geometric_center[1], label="G", marker="1"
    )
    axes.scatter(bending_center[0], bending_center[1], label="B", marker="2")
    axes.scatter(torsion_center[0], torsion_center[1], label="T", marker="3")
    axes.scatter(shear_center[0], shear_center[1], label="S", marker="4")

    usample = np.linspace(0, 1, 17)
    for geometry in section.geometries:
        for curve in geometry:
            knots = curve.knots
            tsample = set(knots)
            for ta, tb in zip(knots, knots[1:]):
                tsample |= set(ta + (tb - ta) * usample)
            tsample = tuple(sorted(tsample))
            points = curve.eval(tsample)
            axes.plot(points[:, 0], points[:, 1], color="k")

    axes.legend()


def plot_field(
    field: IField,
    *,
    subfield_name: Optional[str] = None,
    axes: Optional[plt.Axes] = None,
):
    if not isinstance(field, IField):
        raise NotImplementedError
    if axes is None:
        axes = plt.gca()
    valid_subfield_names = ["Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]
    valid_subfield_names += ["Exx", "Eyy", "Ezz", "Exy", "Exz", "Eyz"]
    valid_subfield_names += ["VM", "TR"]
    if subfield_name is None:
        subfield_name = "VM"
    elif subfield_name not in valid_subfield_names:
        raise NotImplementedError
    if axes is None:
        axes = plt.gca()

    section = field.section

    # First plot internal nodes

    if np.any(field.momentums[:2] != 0):  # Plot bending neutral line
        bending_center = section.bending_center()
        Ixx, Ixy, Iyy = section.second_moment(bending_center)
        vector = np.dot([[Ixx, Ixy], [Ixy, Iyy]], field.momentums[:2])
        axes.axline(
            bending_center, bending_center + vector, color="k", ls="dotted"
        )
        axes.scatter(
            bending_center[0], bending_center[1], label="B", marker="x"
        )

    # Second, plot the countour values
    curves = {}
    for geometry in section.geometries:
        for curve in geometry:
            if curve.label not in curves:
                curves[curve.label] = curve

    xmin, xmax = float("inf"), -float("inf")
    ymin, ymax = xmin, xmax
    usample = np.linspace(0, 1, 17)
    for _, curve in curves.items():
        knots = curve.knots
        tsample = set(knots)
        for ta, tb in zip(knots, knots[1:]):
            tsample |= set(ta + (tb - ta) * usample)
        tsample = tuple(sorted(tsample))
        points = curve.eval(tsample)
        scalars = field_scalar_values(field, points, subfield_name)

        xvals = points[:, 0]
        yvals = points[:, 1]
        xmin = min(xmin, min(xvals))
        xmax = max(xmax, max(xvals))
        ymin = min(ymin, min(yvals))
        ymax = max(ymax, max(yvals))
        x_midpts = np.hstack(
            (xvals[0], 0.5 * (xvals[1:] + xvals[:-1]), xvals[-1])
        )
        y_midpts = np.hstack(
            (yvals[0], 0.5 * (yvals[1:] + yvals[:-1]), yvals[-1])
        )
        coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[
            :, np.newaxis, :
        ]
        coord_mid = np.column_stack((xvals, yvals))[:, np.newaxis, :]
        coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[
            :, np.newaxis, :
        ]
        segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

        default_kwargs = {"capstyle": "butt"}
        lc = LineCollection(segments, **default_kwargs)
        lc.set_array(scalars)  # set the colors of each segment
        lines = axes.add_collection(lc)

    dx = xmax - xmin
    dy = ymax - ymin
    xmin -= 0.1 * dx
    xmax += 0.1 * dx
    ymin -= 0.1 * dy
    ymax += 0.1 * dy
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax)
    fig = plt.gcf()
    fig.colorbar(lines)  # add a color legend
    axes.legend()
