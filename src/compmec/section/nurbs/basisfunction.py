from typing import Optional, Tuple

import numpy as np

from .knotvector import KnotVector
from .util import vectorize


# pylint: disable=invalid-name
def comb(a, b):
    """
    Computes the binomial coeffient of a, b:

    comb(a, b) = a! /((a-b)! * b!)

    Parameters
    ----------
    a: int
    b: int
    return: int

    Example
    -------
    >>> comb(1, 0)
    1
    >>> comb(5, 2)
    10

    """
    prod = 1
    for i in range(min(b, a - b)):
        prod *= a - i
        prod //= i + 1
    return prod


def local_speval_matrix(
    knotvector: KnotVector, reqdegree: Optional[int] = None
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
    knotvector = KnotVector(knotvector)
    if reqdegree is None:
        reqdegree = knotvector.degree
    elif reqdegree < 0:
        raise ValueError(f"reqdegree must be in [0, {knotvector.degree}]")
    knots = knotvector.knots
    spans = knotvector.spans
    j = reqdegree

    matrix = np.zeros((len(knots) - 1, j + 1, j + 1), dtype="object")
    if j == 0:
        one = knotvector.knots[-1] - knotvector.knots[0]
        matrix.fill(one / one)
        return matrix
    matrix_less1 = local_speval_matrix(knotvector, j - 1)
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


def local_to_global_matrix(
    knots: Tuple[float], matrix: Tuple[Tuple[Tuple[float]]]
):
    matrix3d = 0 * np.copy(matrix)
    degree = matrix.shape[2] - 1
    for i, (ta, tb) in enumerate(zip(knots, knots[1:])):
        delta = tb - ta
        for j in range(degree + 1):
            line = matrix[i, :, j] / delta**j
            for k in range(j + 1):
                val = comb(j, k) * line * ta ** (j - k)
                if (j + k) % 2:
                    val *= -1
                matrix3d[i, :, k] += val
    return matrix3d


def global_speval_matrix(
    knotvector: KnotVector, reqdegree: Optional[int] = None
):
    knotvector = KnotVector(knotvector)
    local_matrix = local_speval_matrix(knotvector, reqdegree)
    global_matrix = local_to_global_matrix(knotvector.knots, local_matrix)
    return global_matrix


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


def derivate_matrix(
    matrix: Tuple[Tuple[Tuple[float]]],
) -> Tuple[Tuple[Tuple[float]]]:
    degree = matrix.shape[-1] - 1
    if degree == 0:
        return np.zeros(matrix.shape, dtype=matrix.dtype)
    shape = tuple(matrix.shape[:-1]) + (degree,)
    result = np.zeros(shape, dtype=matrix.dtype)
    for k in range(degree):
        result[:, :, k] = (k + 1) * matrix[:, :, k + 1]
    return result


class CyclicSplineBasisFunction:

    def __init__(self, knotvector: KnotVector):
        knotvector = KnotVector(knotvector)
        self.knotvector = knotvector
        mult = knotvector.mult(knotvector[knotvector.degree])
        self.__ndofs = knotvector.npts + mult - knotvector.degree - 1

        matrices = [global_speval_matrix(knotvector)]
        for i in range(self.degree + 1):
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
    def eval(self, node: float, derivate: int = 0) -> Tuple[float]:
        float(node)
        matrix = self.matrices[min(derivate, self.degree + 1)]
        degree = self.knotvector.degree
        spans = self.knotvector.spans
        result = [0] * self.npts
        if node < self.knots[0] or self.knots[-1] <= node:
            node -= self.knots[0]
            node %= self.knots[-1] - self.knots[0]
            node += self.knots[0]

        span = self.knotvector.span(node)
        ind = spans.index(span)
        for y in range(degree + 1):
            i = (y + span - degree) % self.npts
            coefs = matrix[ind, y]
            result[i] += horner_method(coefs, node)
        return np.array(result)
