"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict


class Curve(ABC):
    """
    Abstract curve to be parent of others.

    This class serves as interface between the curves from others packaged
    like nurbs.Curve, to this package, expliciting the minimum requirements
    of a curve must have.
    With this, it's possible to implement your only type of parametric curve
    """

    @abstractmethod
    @classmethod
    def from_dict(cls, dictionary: Dict) -> Curve:
        """
        Gives a curve instance based on the given parameters
        """
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> Dict:
        """
        Transforms a curve instance into a dictionary
        """
        raise NotImplementedError
