"""
This file contains the data structure and functions to solve
the poisson problem using boundary element method

nabla^2_u = h(x, y)

subject only to neumman's boundary condition

"""

from typing import Tuple

import numpy as np

from .abcs import ISection
from .curve import Curve


class BEMModel:
    """
    A BEM2D Model to solve laplace's equation
    """

    def __init__(self, section: ISection):
        self.section = section
        self.__meshes = {}

    def make_mesh(self, distance: float):
        """
        Create the mesh on the boundary for every curve

        :param distance: The maximum distance to compute mesh
        :type distance: float
        """
        assert distance > 0
        labels = set()
        for geometry in self.section.geometries:
            labels |= set(map(abs, geometry.labels))
        for label in labels:
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

    def __setitem__(self, key: int, value: Tuple[float]):
        self.__meshes[key] = value
