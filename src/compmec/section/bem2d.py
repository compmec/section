"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

"""

import math
from typing import Tuple

import numpy as np

from .abcs import IBasisFunc, ISection
from .curve import Curve
from .integral import Integration


class ComputeMatrix:
    """
    This class is resposible to compute the matrix [M]
    made by the following integral

    M_ij = int phi_j * (r x p')/<r, r> * dt

    Which represents the integral

    int u * (dv/dn) ds
    """

    def __init__(self, curve: Curve, basis: IBasisFunc):
        self.curve = curve
        self.basis = basis

    @property
    def tmesh(self) -> Tuple[float]:
        """
        The subdivisions of parametric space

        :getter: Returns the union of curve and base knots
        :type: Tuple[float]
        """
        tmesh = set(self.curve.knots) | set(self.basis.knots)
        return np.array(sorted(tmesh))

    # pylint: disable=too-many-locals
    def inpolygon(self, tsources: Tuple[float]):
        """
        Computes the integral when the sources are placed at the curve.
        The emplacement of these sources are given by parameter 'tsources'

        We suppose the curve is a polygon.

        Parameters
        ----------

        :param tsources: The parametric emplacement of sources
        :type tsources: Tuple[float]
        :return: The output matrix, integral of UVn
        :rtype: Tuple[Tuple[float]]
        """
        ndofs = self.basis.ndofs
        assert ndofs == len(tsources)
        cknots = np.array(self.curve.knots, dtype="float64")
        tknots = np.array(sorted(set(self.tmesh) | set(tsources)))
        matrix = np.zeros((ndofs, ndofs), dtype="float64")
        nodes, weights = Integration.gauss(10)

        vertices = self.curve.eval(cknots[:-1])
        vectors = np.roll(vertices, shift=-1, axis=0) - vertices
        sources = self.curve.eval(tsources)

        for i, vector in enumerate(vectors):
            tva, tvb = cknots[i], cknots[i + 1]
            vertex = vertices[i]
            mask = (tva <= tknots) * (tknots <= tvb)
            tmesh = tknots[mask]
            for tk0, tk1 in zip(tmesh, tmesh[1:]):
                tvals = tk0 + nodes * (tk1 - tk0)
                taus = (tvals - tva) / (tvb - tva)
                phis = self.basis.eval(tvals)
                points = vertex + np.tensordot(taus, vector, axes=0)
                for i, (tsi, source) in enumerate(zip(tsources, sources)):
                    if tsi in (tk0, tk1):
                        continue
                    radius = tuple(point - source for point in points)
                    radius = np.array(radius, dtype="float64")
                    rcrossdp = vector[1] * radius[:, 0]
                    rcrossdp -= vector[0] * radius[:, 1]
                    rcrossdp *= (tk1 - tk0) / (tvb - tva)
                    rinnerr = np.einsum("ij,ij->i", radius, radius)  # <r, r>
                    funcvals = rcrossdp / rinnerr
                    matrix[i] += np.einsum(
                        "ij,j,j->i", phis, weights, funcvals
                    )
        matrix[np.abs(matrix) < 1e-9] = 0
        return matrix / math.tau

    def outcurve(self, sources: Tuple[Tuple[float]]):
        """
        Computes the integral when the sources are placed outside (or not)
        of the curve.

        The emplacement of these sources can be any, including on the curve.

        Not yet implemented cases which source lies on the curve

        Parameters
        ----------

        :return: The output matrix, integral of UVn
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError


class BEMModel:
    """
    A BEM2D Model to solve laplace's equation
    """

    def __init__(self, section: ISection):
        self.section = section
        labels = set()
        for geometry in self.section.geometries:
            labels |= set(map(abs, geometry.labels))
        self.all_labels = tuple(sorted(labels))
        self.__meshes = {}

    def make_mesh(self, distance: float):
        """
        Create the mesh on the boundary for every curve

        :param distance: The maximum distance to compute mesh
        :type distance: float
        """
        assert distance > 0
        for label in self.all_labels:
            curve = Curve.instances[label]
            knots = curve.knots
            new_mesh = set(knots)
            vertices = curve.eval(knots)
            vectors = vertices[1:] - vertices[:-1]
            for i, vector in enumerate(vectors):
                ndiv = np.linalg.norm(vector) / distance
                ndiv = max(2, int(np.ceil(ndiv)))
                new_mesh |= set(np.linspace(knots[i], knots[i + 1], ndiv))
            self[label] = tuple(sorted(new_mesh))

    def solve(self):
        """
        Solves the BEM problem, computing
        """
        raise NotImplementedError

    def __getitem__(self, key: int):
        return self.__meshes[key]

    def __setitem__(self, label: int, value: Tuple[float]):
        if label not in self.all_labels:
            msg = f"Given label {label} is not in {self.all_labels}"
            raise ValueError(msg)
        self.__meshes[label] = value
