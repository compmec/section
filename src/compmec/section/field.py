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

    def __stress_axial(self):
        return self.forces[2] / self.section.area()

    def __stress_bending(self, point):
        bend_center = self.section.bending_center()
        ixx, ixy, iyy = self.section.second_moment(bend_center)
        detii = ixx * iyy - ixy**2
        matrix = np.array([[iyy, -ixy], [-ixy, ixx]])
        momx, momy, _ = self.momentums
        return point @ matrix @ [-momy, momx] / detii

    def __stress_shear(self, points, homosection):
        raise NotImplementedError

    def __stress_torsion(self, point, homosection):
        result = np.array([-point[1], point[0]], dtype="float64")
        result += homosection.warping.grad(point)
        return result * self.momentums[2] / self.section.torsion_constant()

    def __winding_numbers(self, point: Tuple[float]) -> Tuple[float]:
        """
        Computes the winding number of every point,
        for every geometry

        """
        wind_tolerance = 1e-6
        homosections = tuple(self.section)
        winds = np.zeros(len(homosections), dtype="float64")
        for i, homosection in enumerate(homosections):
            winds[i] = homosection.geometry.winding(point)
        winds[np.abs(winds) < wind_tolerance] = 0
        return winds

    def __stress_eval(self, point, winds):
        """
        Evaluates the stress values

        The stress tensor is given by
            [  0    0  S13]
        S = [  0    0  S23]
            [S13  S23  S33]

        Returned values are a matrix of shape (m, 3)
        each line are the stress components: S33, S13, S23

        :param point: The wanted point, a vector of two components
        :type point: Tuple[float]
        :return: The strain matrix of shape (m, 3)
        :rtype: Tuple[Tuple[float]]

        """
        results = np.zeros((len(winds), 3), dtype="float64")
        fx, fy, fz = self.forces
        mx, my, mz = self.momentums
        for i, homosection in enumerate(self.section):
            if winds[i] == 0:
                continue
            if fz:  # Axial force
                results[i, 0] += self.__stress_axial()
            if mx or my:  # Bending moments
                results[i, 0] += self.__stress_bending(point)
            if mz:  # Torsion
                results[i, 1:] += self.__stress_torsion(point, homosection)
            if fx or fy:  # Poisson's evaluation
                results[i, 1:] += self.__stress_shear(point, homosection)
        return results

    def __strain_eval(self, stress):
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

        Returned values are a matrix of shape (4, )
        each line are the strain components: E11, E33, E13, E23

        :param stress: The values of S33, S13, S23
        :type points: Tuple[float]
        :return: The values of E11, E33, E13, E23
        :rtype: Tuple[float]

        """
        strain = np.zeros((len(stress), 4), dtype="float64")
        for i, homosection in enumerate(self.section):
            young = homosection.material.young_modulus
            poiss = homosection.material.poissons_ratio
            fract11 = -poiss / young
            fract13 = (1 + poiss) / young
            strain[i, 1] += stress[i, 0] / young
            strain[i, 0] += fract11 * stress[i, 0]
            strain[i, 2] += fract13 * stress[i, 1]
            strain[i, 3] += fract13 * stress[i, 2]
        return strain

    def __eval_one(self, point: Tuple[float]):
        point = np.array(point, dtype="float64")
        if point.ndim != 1 or len(point) != 2:
            raise ValueError(f"point is not a 2D vector! {point.shape}")
        winds = self.__winding_numbers(point)
        if np.all(winds == 0):
            return (0,) * 7
        stress = self.__stress_eval(point, winds)
        strain = self.__strain_eval(stress)
        stress = np.average(stress, axis=0, weights=winds)
        strain = np.average(strain, axis=0, weights=winds)
        return tuple(stress) + tuple(strain)

    def eval(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        """
        Evaluate the stress and strain at given point

        The inputs/outputs are matrices:
        * points is a 2D vector
        * stress is a vector of (S13, S23, S33)
        * strain is a vector with (E11, E33, E13, E23)

        The stress components are S13, S23, S33, meaning
        two shear stresses and one normal stress
        The strain components are E11, E33, E13, E23

        :param points: The wanted point, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The pair (stress, strain)
        :rtype: Tuple[Tuple[Tuple[float]]]

        Example
        -------

        >>> points = [(0, 0), (1, 0)]
        >>> S33, S13, S23, E11, E33, E13, E23 = field.eval(points)

        """
        values = tuple(map(self.__eval_one, points))
        values = np.transpose(np.array(values))
        return values


def field_scalar_values(
    field: IField, points: Tuple[Tuple[float]], subfield_name: Union[str, None]
):
    if not isinstance(field, ChargedField):
        raise NotImplementedError
    if subfield_name is None:
        raise NotImplementedError
    if subfield_name in ("Sxx", "Syy", "Sxy", "Exy"):
        return (0,) * len(points)

    values = field.eval(points)
    Szz, Sxz, Syz = np.transpose(values[:, :3])
    Exx, Ezz, Exz, Eyz = np.transpose(values[:, 3:])
    if subfield_name == "Sxz":
        return Sxz
    if subfield_name == "Syz":
        return Syz
    if subfield_name == "Szz":
        return Szz
    if subfield_name in ("Exx", "Eyy"):
        return Exx
    if subfield_name == "Ezz":
        return Ezz
    if subfield_name == "Exz":
        return Exz
    if subfield_name == "Eyz":
        return Eyz
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
