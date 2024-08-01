import types
from collections.abc import Sequence
from typing import Any

import numpy as np


def vectorize(func):
    uncastable = (types.GeneratorType, range)

    def wrapper(self, nodes: Sequence[Any], *args, **kwargs) -> Sequence[Any]:
        try:
            iter(nodes)
        except TypeError:
            return func(self, nodes, *args, **kwargs)
        generator = (func(self, node, *args, **kwargs) for node in nodes)
        if isinstance(nodes, np.ndarray):
            return np.array(tuple(generator))
        elif not isinstance(nodes, uncastable):
            return np.array(tuple(generator))
        return generator

    return wrapper
