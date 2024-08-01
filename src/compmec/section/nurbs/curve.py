from typing import Tuple

import numpy as np

from .basisfunction import global_speval_matrix, horner_method
from .knotvector import KnotVector
from .util import vectorize


def derivate_matrix(matrix: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    degree = matrix.shape[-1] - 1
    shape = tuple(matrix.shape[:-1]) + (degree,)
    result = np.zeros(shape, dtype=matrix.dtype)
    if degree == 0:
        return result
    for k in range(degree):
        result[:, k] = (k + 1) * matrix[:, k + 1]
    return result


class CyclicScalarSpline:

    def __init__(self, knotvector: KnotVector, ctrlpoints: Tuple[float]):
        knotvector = KnotVector(knotvector)
        self.knotvector = knotvector

        matrix3d = global_speval_matrix(knotvector)
        mult = knotvector.mult(knotvector[knotvector.degree])
        ndofs = knotvector.npts + mult - knotvector.degree - 1
        self.__ndofs = ndofs

        matrix2d = []
        for i, span in enumerate(knotvector.spans):
            line = 0 * matrix3d[i, 0, :]
            for y in range(self.degree + 1):
                z = (y + span - self.degree) % self.npts
                line += ctrlpoints[z] * matrix3d[i, y, :]
            matrix2d.append(line)
        matrices = [np.array(matrix2d)]
        for i in range(self.degree):
            matrices.append(derivate_matrix(matrices[i]))
        self.matrices = tuple(matrices)

    @property
    def degree(self) -> int:
        return self.knotvector.degree

    @property
    def npts(self) -> int:
        return self.__ndofs

    @property
    def knots(self) -> Tuple[float]:
        return self.knotvector.knots

    @vectorize
    def eval(self, node: float, derivate: int = 0) -> float:
        float(node)
        matrix = self.matrices[min(derivate, self.degree)]
        span = self.knotvector.span(node)
        ind = self.knotvector.spans.index(span)
        return horner_method(matrix[ind], node)
