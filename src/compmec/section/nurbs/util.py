import types
from collections.abc import Sequence
from typing import Any

import numpy as np


def vectorize(func):
    uncastable = (types.GeneratorType, range)

    def wrapper(self, nodes: Sequence[Any]) -> Sequence[Any]:
        try:
            iter(nodes)
        except TypeError:
            return func(self, nodes)
        values_generator = (func(self, node) for node in nodes)
        if isinstance(nodes, np.ndarray):
            return np.array(tuple(values_generator))
        elif not isinstance(nodes, uncastable):
            values_generator = np.array(tuple(values_generator))
        return values_generator

    return wrapper
