"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""

from __future__ import annotations

import math
from collections import OrderedDict
from typing import Tuple

import numpy as np

from . import integral
from .abcs import ICurve
from .nurbs import CyclicScalarSpline, KnotVector, vectorize


class NurbsCurve(ICurve):
    """
    Base class that tracks the instances
    """

    instances = OrderedDict()

    def __init__(
        self,
        knotvector: Tuple[float],
        ctrlpoints: Tuple[Tuple[float]],
        *,
        weights: Tuple[float] = None,
    ):
        knotvector = KnotVector(knotvector)
        self.knotvector = knotvector
        xvalues, yvalues = np.transpose(ctrlpoints)
        if weights is not None:
            raise NotImplementedError
        self.__xfunction = CyclicScalarSpline(knotvector, xvalues)
        self.__yfunction = CyclicScalarSpline(knotvector, yvalues)

        area = integral.Bidimensional.general(self, 0, 0)
        self.__area = float(area)

    @vectorize
    def eval(self, parameter: float) -> Tuple[float]:
        xvalue = self.__xfunction.eval(parameter, 0)
        yvalue = self.__yfunction.eval(parameter, 0)
        return (xvalue, yvalue)

    @vectorize
    def deval(self, parameter: float) -> Tuple[float]:
        dxvalue = self.__xfunction.eval(parameter, 1)
        dyvalue = self.__yfunction.eval(parameter, 1)
        return (dxvalue, dyvalue)

    def __float__(self) -> float:
        return self.__area

    @property
    def knots(self) -> Tuple[float]:
        return self.knotvector.knots

    @property
    def degree(self) -> int:
        return self.knotvector.degree

    def winding(self, point: Tuple[float]) -> float:
        # Verify if the point is at any vertex
        vertices = self.eval(self.knots[:-1])
        for i, vertex in enumerate(vertices):
            if np.all(point == vertex):
                vec_left = vertices[(i - 1) % len(vertices)] - point
                vec_righ = vertices[(i + 1) % len(vertices)] - point
                wind = 0.5 * np.arccos(np.inner(vec_left, vec_righ)) / np.pi
                return wind if float(self) > 0 else 1 - wind

        wind = 0
        for vertexa, vertexb in zip(vertices, np.roll(vertices, -1, axis=0)):
            sub_wind = integral.winding_number_linear(vertexa, vertexb, point)
            if abs(sub_wind) == 0.5:
                wind = 0.5
                break
            wind += sub_wind
        if float(self) > 0:
            return wind
        return 0 if wind == -1 else 1 - wind

    def projection(self, point: Tuple[float]) -> Tuple[float]:
        x0, y0 = point
        ndiv = 4
        knots = self.knotvector.knots
        projvals = [float("inf")]
        distsquares = [float("inf")]
        for ta, tb in zip(knots, knots[1:]):
            tvs = (((ndiv - i) * ta + i * tb) / ndiv for i in range(1, ndiv))
            for tv in tvs:
                for _ in range(3):  # Newton iteration
                    deltax = self.__xfunction.eval(tv, 0) - x0
                    derix = self.__xfunction.eval(tv, 1)
                    deltay = self.__yfunction.eval(tv, 0) - y0
                    deriy = self.__yfunction.eval(tv, 1)
                    ftv = deltax * derix + deltay * deriy
                    if ftv == 0:
                        break
                    deri2x = self.__xfunction.eval(tv, 2)
                    deri2y = self.__yfunction.eval(tv, 2)
                    dftv = deltax * deri2x + deltay * deri2y
                    dftv += derix**2 + deriy**2
                    tv -= ftv / dftv
                    tv = max(ta, min(tb, tv))
                deltax = self.__xfunction.eval(tv, 0) - x0
                deltay = self.__yfunction.eval(tv, 0) - y0
                dist_square = deltax**2 + deltay**2
                min_dist_square = min(distsquares)
                if dist_square <= min_dist_square:
                    projvals.append(tv)
                    distsquares.append(dist_square)
        min_dist_square = min(distsquares)
        projvals = set(
            tv
            for tv, ds in zip(projvals, distsquares)
            if ds == min_dist_square
        )
        return tuple(sorted(projvals))


class PolygonCurve(ICurve):

    def __init__(self, vertices: Tuple[Tuple[float]]):
        vertices = np.array(vertices, dtype="float64")
        if vertices.ndim != 2 or vertices.shape[1] != 2:
            raise ValueError
        cross = vertices[:, 0] * np.roll(vertices[:, 1], shift=-1)
        cross -= vertices[:, 1] * np.roll(vertices[:, 0], shift=-1)
        self.__vertices = vertices
        self.__vectors = np.roll(vertices, shift=-1, axis=0) - vertices
        self.__area = float(sum(cross) / 2)

    @property
    def knots(self) -> Tuple[float]:
        return tuple(range(len(self.vertices) + 1))

    @property
    def vertices(self) -> Tuple[Tuple[float]]:
        return self.__vertices

    @property
    def vectors(self) -> Tuple[Tuple[float]]:
        return self.__vectors

    def __float__(self):
        return self.__area

    @vectorize
    def eval(self, parameter: float) -> Tuple[float]:
        index = math.floor(parameter) % len(self.vertices)
        return self.vertices[index] + (parameter % 1) * self.vectors[index]

    @vectorize
    def deval(self, parameter: float) -> Tuple[float]:
        index = math.floor(parameter) % len(self.vertices)
        return self.vectors[index]

    def projection(self, point: Tuple[float]) -> Tuple[float]:
        vertices = (point - vertex for vertex in self.vertices)
        projects = []
        dist_squares = []
        for i, vertex in enumerate(vertices):
            vector = self.vectors[i]
            param = np.inner(vertex, vector) / np.inner(vector, vector)
            param = max(0, min(1, param))
            vectdist = param * vector - vertex
            projects.append(i + param)
            dist_squares.append(np.inner(vectdist, vectdist))
        min_distsquare = min(dist_squares)
        params = set(
            p for p, d2 in zip(projects, dist_squares) if d2 == min_distsquare
        )
        return tuple(sorted(params))

    def winding(self, point: Tuple[float]) -> float:
        nverts = len(self.vertices)
        proj = self.projection(point)[0]
        proj_vec = self.eval([proj])[0] - point
        if np.inner(proj_vec, proj_vec) < 1e-6:
            if not isinstance(proj, int):
                return 0.5
            v0, v1 = self.vectors[proj], self.vectors[(proj + 1) % nverts]
            inner = np.inner(v0, v1)
            cross = np.cross(v0, v1)
            wind = np.arctan2(cross, inner) / math.tau
            return wind % 1

        vertices = tuple(vertex - point for vertex in self.vertices)
        wind = 0
        for i, vertex0 in enumerate(vertices):
            vertex1 = vertices[(i + 1) % nverts]
            cross = np.cross(vertex0, vertex1)
            inner = np.inner(vertex0, vertex1)
            wind += np.arctan2(cross, inner) / math.tau
        if float(self) > 0:
            return wind
        return 0 if wind == -1 else 1 - wind
