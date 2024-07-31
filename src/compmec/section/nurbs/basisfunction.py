from typing import Tuple

import numpy as np

from .knotvector import KnotVector


def speval_matrix(
    knotvector: KnotVector, reqdegree: int
) -> Tuple[Tuple[Tuple[float]]]:
    """
    Given a knotvector, it has properties like
        - number of points: npts
        - polynomial degree: degree
        - knots: A list of non-repeted knots
        - spans: The span of each knot
    This function returns a matrix of size
        (m) x (j+1) x (j+1)
    which
        - m is the number of segments: len(knots)-1
        - j is the requested degree
    """
    if reqdegree < 0:
        raise ValueError(f"reqdegree must be in [0, {knotvector.degree}]")
    knots = knotvector.knots
    spans = knotvector.spans
    j = reqdegree

    matrix = np.zeros((len(knots) - 1, j + 1, j + 1), dtype="object")
    if j == 0:
        one = knotvector.knots[-1] - knotvector.knots[0]
        matrix.fill(one / one)
        return matrix
    matrix_less1 = speval_matrix(knotvector, j - 1)
    for y in range(j):
        for z, sz in enumerate(spans):
            i = y + sz - j + 1
            denom = knotvector[i + j] - knotvector[i]
            matrix_less1[z, y, :] /= denom

            a0 = knots[z] - knotvector[i]
            a1 = knots[z + 1] - knots[z]
            b0 = knotvector[i + j] - knots[z]
            b1 = knots[z] - knots[z + 1]

            matrix[z, y, :-1] += b0 * matrix_less1[z, y]
            matrix[z, y, 1:] += b1 * matrix_less1[z, y]
            matrix[z, y + 1, :-1] += a0 * matrix_less1[z, y]
            matrix[z, y + 1, 1:] += a1 * matrix_less1[z, y]
    return matrix


def derivate_matrix(matrix: Tuple[Tuple[Tuple[float]]]) -> float:
    m, p, q = matrix.shape
    result = np.zeros((m, p, q - 1), dtype=matrix.dtype)
    if q == 1:
        return result
    for i in range(m):
        for j in range(p):
            for k in range(1, q):
                value = matrix[i, j, k]
                result[i, j, k - 1] = k * value
    return result


def horner_method(coefs: Tuple[float], node: float):
    """
    Horner method is a efficient method of computing polynomials
    Let's say you have a polynomial
        P(x) = a_0 + a_1 * x + ... + a_n * x^n
    A way to compute P(x_0) is
        P(x_0) = a_0 + a_1 * x_0 + ... + a_n * x_0^n
    But a more efficient way is to use
        P(x) = ((...((x * a_n + a_{n-1})*x)...)*x + a_1)*x + a_0

    Input:
        coefs : Tuple[float] = (a_0, a_1, ..., a_n)
        value : float = x_0
    """
    soma = 0
    for ck in coefs[::-1]:
        soma *= node
        soma += ck
    return soma



class CyclicSplineBasisFunction:

    def __init__(self, knotvector: KnotVector):
        knotvector = KnotVector(knotvector)
        self.knotvector = knotvector
        self.matrix = speval_matrix(knotvector, knotvector.degree)
        mult = knotvector.mult(knotvector[knotvector.degree])
        self.__ndofs = knotvector.npts + mult - knotvector.degree - 1

    @property
    def npts(self) -> int:
        return self.__ndofs

    @property
    def degree(self) -> int:
        return self.knotvector.degree

    def eval(self, nodes: Tuple[float]) -> Tuple[Tuple[float]]:
        knots = self.knotvector.knots
        spans = self.knotvector.spans
        result = np.zeros((self.npts, len(nodes)), dtype="object")
        lima, limb = self.knotvector.knots[0], self.knotvector.knots[-1]
        diff = limb - lima
        nodes = (
            node if lima <= node < limb else lima + ((node - lima) % diff)
            for node in nodes
        )
        for j, node in enumerate(nodes):
            span = self.knotvector.span(node)
            ind = spans.index(span)
            shifnode = node - knots[ind]
            shifnode /= knots[ind + 1] - knots[ind]
            for y in range(self.degree + 1):
                i = (y + span - self.degree) % self.npts
                coefs = self.matrix[ind, y]
                result[i, j] += horner_method(coefs, shifnode)
        return result
