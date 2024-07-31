from typing import Tuple

import numpy as np

from .basisfunction import horner_method, speval_matrix
from .knotvector import KnotVector


class CyclicScalarSpline:

    def __init__(self, knotvector: KnotVector, ctrlpoints: Tuple[float]):
        knotvector = KnotVector(knotvector)
        self.knotvector = knotvector

        matrix3d = speval_matrix(knotvector, knotvector.degree)
        mult = knotvector.mult(knotvector[knotvector.degree])
        ndofs = knotvector.npts + mult - knotvector.degree - 1
        self.__ndofs = ndofs

        matrix2d = []
        for i, span in enumerate(knotvector.spans):
            line = 0 * matrix3d[i, 0, :]
            for y in range(self.knotvector.degree + 1):
                z = (y + span - self.degree) % self.npts
                line += ctrlpoints[z] * matrix3d[i, y, :]
            matrix2d.append(line)
        self.matrix = np.array(matrix2d)

    @property
    def npts(self) -> int:
        return self.__ndofs

    @property
    def degree(self) -> int:
        return self.knotvector.degree

    def eval(self, nodes: Tuple[float]) -> Tuple[Tuple[float]]:
        knots = self.knotvector.knots
        spans = self.knotvector.spans
        result = [0] * len(nodes)
        for j, node in enumerate(nodes):
            span = self.knotvector.span(node)
            ind = spans.index(span)
            shifnode = node - knots[ind]
            shifnode /= knots[ind + 1] - knots[ind]
            result[j] = horner_method(self.matrix[ind], shifnode)
        return np.array(result)
