"""
File that contains the Basis Functions to compute BEM2D
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from pynurbs import Function, KnotVector
from pynurbs.heavy import Calculus, ImmutableKnotVector

from .abcs import IBasisFunction


def cyclic_knotvector(
    knots: Tuple[float], degree: int = 1, limits: Optional[Tuple[float]] = None
) -> ImmutableKnotVector:
    knots = tuple(sorted(knots))
    knots = np.array(knots, dtype="object")
    for knot in set(knots):
        assert sum(knots == knot) <= degree + 1
    if limits is None:
        limits = (knots[0], knots[-1])
    mult_left = sum(knots == limits[0])
    mult_righ = sum(knots == limits[1])
    assert knots[0] == limits[0]
    assert knots[-1] == limits[-1]
    assert mult_left == mult_righ
    knots = list(knots)
    for _ in range(mult_righ):
        knots.pop(-1)
    knots = np.array(knots, dtype="object")
    differences = list(knots[1:] - knots[:-1])
    differences.append(limits[1] - knots[-1])
    mult_left = sum(1 for knot in knots if knot == knots[0])

    npts = len(knots) + degree - mult_left + 1
    knotvector = [None] * (npts + degree + 1)
    for i, knot in enumerate(knots):
        knotvector[i + degree - mult_left + 1] = knot
    index = degree
    while index >= 0 and knotvector[index] is not None:
        index -= 1
    j = -1
    while index >= 0:
        diff = differences[j % len(differences)]
        knotvector[index] = knotvector[index + 1] - diff
        j -= 1
        index -= 1
    index = npts
    while index <= npts + degree and knotvector[index] is not None:
        index += 1
    index -= 1
    j = -1
    while index < npts + degree:
        diff = differences[j % len(differences)]
        knotvector[index + 1] = knotvector[index] + diff
        index += 1
        j += 1

    knotvector = np.array(knotvector, dtype="object")
    for knot in set(knots):
        assert sum(knots == knot) == sum(knotvector == knot)
    return tuple(knotvector)


class SplineBasisFunction(IBasisFunction):
    """
    Basis function that uses BSplines
    """

    @classmethod
    def cyclic(cls, knots: Tuple[float], degree: int = 1) -> SplineBasisFunction:
        """
        Creates a cyclic basis function

        :param knots: The division emplacements
        :type knots: Tuple[float]
        :param degree: The maximum curve degree, default 1
        :type degree: int
        :return: The SplineBasisFunction instance
        :rtype: SplineBasisFunction
        """
        knotvector = cyclic_knotvector(knots, degree=degree)
        knotvector = ImmutableKnotVector(knotvector, degree=degree)
        return cls(knotvector)

    def __init__(self, knotvector: Tuple[float]):
        if isinstance(knotvector, (KnotVector, ImmutableKnotVector)):
            degree = knotvector.degree
        else:
            degree = None
        knotvector = ImmutableKnotVector(knotvector, degree)
        degree = knotvector.degree
        npts = knotvector.npts
        self.basis = Function(knotvector)
        if degree != 0:
            delta = Calculus.difference_matrix(knotvector)
        else:
            delta = np.zeros((npts, npts), dtype="float64")
        self.derivate_matrix = delta
        mult = knotvector.mult(knotvector[knotvector.degree])
        self.__ndofs = knotvector.npts + mult - knotvector.degree - 1

    @property
    def knots(self):
        return self.basis.knots

    @property
    def ndofs(self):
        return self.__ndofs

    def eval(self, parameters):
        ndofs = self.ndofs
        start = self.basis.degree
        values = self.basis.eval(parameters)
        values = np.array(values)
        values[:start] += values[ndofs:]
        return values[:ndofs]

    def deval(self, parameters):
        knotvector = self.basis.knotvector
        degree = knotvector.degree
        ndofs = self.ndofs
        values = self.basis[:, degree - 1](parameters)
        values = np.dot(self.derivate_matrix, values)
        values[:degree] += values[ndofs:]
        return values[:ndofs]


def distributed_knots(basis: SplineBasisFunction):
    """
    Find 'n' distributed knots on a cyclic basis functions
    such 'n' is the number of degrees of freedom of basis
    and these are knots are not repeted
    """
    ndofs = basis.ndofs
    if len(basis.knots) != ndofs + 1:
        raise NotImplementedError
    return basis.knots[:-1]
