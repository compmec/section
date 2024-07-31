from typing import Tuple

import numpy as np

from .basisfunction import CyclicSplineBasisFunction
from .curve import CyclicScalarSpline


def derivate_matrix(matrix: Tuple[Tuple[Tuple[float]]]) -> float:
    m, p, q = matrix.shape
    result = np.zeros((m, p, q - 1), dtype=matrix.dtype)
    if q == 1:
        return result
    for k in range(1, q):
        result[:, :, k - 1] = k * matrix[:, :, k]
    return result


def derivate_basis(
    basis: CyclicSplineBasisFunction,
) -> CyclicSplineBasisFunction:
    deri_matrix = derivate_matrix(basis.matrix)
    deri_basis = basis.__class__(basis.knotvector)
    deri_basis.matrix = deri_matrix
    deri_basis.degree = max(0, basis.degree - 1)
    return deri_basis


def derivate_scalar(scalar: CyclicScalarSpline) -> CyclicScalarSpline:
    deri_matrix = derivate_matrix(scalar.matrix)
    deri_scalar = scalar.__class__(scalar.knotvector, [0] * scalar.npts)
    deri_scalar.matrix = deri_matrix
    deri_scalar.degree = max(0, scalar.degree - 1)
    return deri_scalar
