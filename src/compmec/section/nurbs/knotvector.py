from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


class KnotVector(tuple):

    def __new__(cls, knotvector: Tuple[float], degree: Optional[int] = None):
        if isinstance(knotvector, cls):
            return knotvector
        return super(KnotVector, cls).__new__(cls, tuple(knotvector))

    def __init__(self, knotvector: Tuple[float], degree: Optional[int] = None):
        if isinstance(knotvector, KnotVector):
            return
        if degree is None:
            degree = 0
            while knotvector[degree] == knotvector[degree + 1]:
                degree += 1
        knotvector = tuple(knotvector)
        npts = len(knotvector) - degree - 1
        if npts <= degree:
            raise ValueError
        self.degree = degree
        self.npts = npts
        slic = slice(degree, npts + 1)
        self.knots = tuple(sorted(set(knotvector[slic])))
        spans = []
        span = 0
        for knot in self.knots[:-1]:
            while knotvector[span] != knot:
                span += 1
            while knotvector[span + 1] == knot:
                span += 1
            spans.append(span)
        self.spans = tuple(spans)

    def mult(self, node: float) -> int:
        return sum(knot == node for knot in self)

    def span(self, node: float) -> int:
        if not self.knots[0] <= node < self.knots[-1]:
            node -= self.knots[0]
            node % self.knots[-1] - self.knots[0]
            node += self.knots[0]
        for i, span in enumerate(self.spans):
            if self.knots[i] <= node < self.knots[i + 1]:
                return span
        return span


def cyclic_knotvector(
    knots: Tuple[float], degree: int = 1, limits: Optional[Tuple[float]] = None
) -> KnotVector:
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
    return KnotVector(knotvector, degree)
