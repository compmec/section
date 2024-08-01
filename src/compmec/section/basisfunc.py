"""
File that contains the Basis Functions to compute BEM2D
"""

from __future__ import annotations

from fractions import Fraction
from typing import Tuple

import numpy as np

from .abcs import IBasisFunction
from .nurbs import CyclicSplineBasisFunction, KnotVector


class SplineBasisFunction(IBasisFunction):
    """
    Creates a cyclic basis function

    :param knots: The division emplacements
    :type knots: Tuple[float]
    :param degree: The maximum curve degree, default 1
    :type degree: int
    :return: The SplineBasisFunction instance
    :rtype: SplineBasisFunction
    """

    def __init__(self, knots: Tuple[float], degree: int = 1):
        knotvector = KnotVector.cyclic(knots, degree)
        self.__internal = CyclicSplineBasisFunction(knotvector)

    @property
    def degree(self) -> int:
        return self.__internal.degree

    @property
    def ndofs(self):
        return self.__internal.npts

    @property
    def knots(self):
        return self.__internal.knotvector.knots

    def eval(self, parameters):
        return self.__internal.eval(parameters, 0)

    def deval(self, parameters):
        return self.__internal.eval(parameters, 1)


def distributed_knots(basis: SplineBasisFunction):
    """
    Find 'n' distributed knots on a cyclic basis functions
    such 'n' is the number of degrees of freedom of basis
    and these are knots are not repeted
    """
    ndofs = basis.ndofs
    degree = basis.degree
    knotvector = basis.knotvector
    knots = knotvector.knots
    outknots = [None] * ndofs

    for knot in set(knots[:-1]):
        if knotvector.mult(knot) == degree:
            values = basis.eval(knot)
            index = np.where(values == max(values))[0][0]
            outknots[index] = knot

    if all(knot is not None for knot in outknots):
        return tuple(outknots)

    nsubdiv = degree + 1
    uniform_nodes = tuple(
        Fraction(i, 2 * nsubdiv) for i in range(1, 2 * nsubdiv, 2)
    )
    for ta, tb in zip(knots, knots[1:]):
        for uni0, uni1 in zip(uniform_nodes, uniform_nodes[1:]):
            tc0 = (1 - uni0) * ta + uni0 * tb
            tc1 = (1 - uni1) * ta + uni1 * tb
            mask = (basis.eval(tc0) != 0) * (basis.eval(tc1) != 0)
            deri0vals = basis.deval(tc0)
            deri1vals = basis.deval(tc1)
            for j, outknot in enumerate(outknots):
                if outknot is not None:
                    continue
                if not mask[j]:
                    continue
                df0, df1 = deri0vals[j], deri1vals[j]
                if df0 * df1 > 0:
                    continue
                outknots[j] = (tc0 * df1 - tc1 * df0) / (df1 - df0)
    return tuple(outknots)
