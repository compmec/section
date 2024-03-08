"""
File that contains the Basis Functions to compute BEM2D
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from compmec.nurbs import Function, KnotVector

from .abcs import IBasisFunc


class BasisFunc(IBasisFunc):
    """
    Basis function that uses BSplines
    """

    @classmethod
    def cyclic(cls, knots: Tuple[float], degree: int = 1) -> BasisFunc:
        """
        Creates a cyclic basis function

        :param knots: The division emplacements
        :type knots: Tuple[float]
        :param degree: The maximum curve degree, default 1
        :type degree: int
        :return: The BasisFunc instance
        :rtype: BasisFunc
        """
        knots = np.array(sorted(knots), dtype="float64")
        start = len(knots) - degree - 1
        end = 1 + degree
        left = knots[start:-1] + knots[0] - knots[-1]
        righ = knots[1:end] + knots[-1] - knots[0]
        knots = np.concatenate([left, knots, righ])
        knotvector = KnotVector(knots, degree=degree)
        return cls(knotvector)

    def __init__(self, knotvector: Tuple[float]):
        knotvector = KnotVector(knotvector)
        self.basis = Function(knotvector)

    @property
    def knots(self):
        return self.basis.knots

    @property
    def ndofs(self):
        maxknot = self.basis.knotvector.limits[1]
        ndofs = self.basis.npts - 1 - self.basis.degree
        ndofs += self.basis.knotvector.mult(maxknot)
        return ndofs

    def eval(self, parameters):
        ndofs = self.ndofs
        start = self.basis.degree
        values = self.basis.eval(parameters)
        values = np.transpose(values)
        values[:, :start] += values[:, ndofs:]
        return values[:, :ndofs]
