"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""
from __future__ import annotations

from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Dict, Optional, Tuple

import numpy as np

from compmec import nurbs


class Curve(ABC):
    """
    Abstract curve to be parent of others.

    This class serves as interface between the curves from others packaged
    like nurbs.Curve, to this package, expliciting the minimum requirements
    of a curve must have.
    With this, it's possible to implement your only type of parametric curve
    """

    instances = OrderedDict()

    @staticmethod
    def __next_available_label() -> str:
        index = 1
        while index in Curve.instances:
            index += 1
        return index

    def __new__(cls, label: Optional[int] = None) -> Curve:
        if label is None:
            label = Curve.__next_available_label()
        elif label in Curve.instances:
            msg = f"Cannot create curve '{label}'! There's already one!"
            raise ValueError(msg)
        instance = super().__new__(cls)
        instance.label = label
        Curve.instances[label] = instance
        return instance

    @classmethod
    @abstractmethod
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

    @property
    def label(self) -> int:
        """
        Gives the curve label

        :getter: Returns the curve's label
        :setter: Attribuates a new label for curve
        :type: int

        """
        return self.__label

    @label.setter
    def label(self, new_label: int):
        try:
            cur_label = self.label
        except AttributeError:
            self.__label = new_label
            return
        if cur_label == new_label:
            return
        if new_label in self.instances:
            msg = f"Cannot set the label '{new_label}' for the curve "
            msg += f"'{cur_label}' cause there's already a curve with "
            msg += "the same label!\n"
            msg += f"Cur: '{self.instances[new_label]}'\n"
            msg += f"New: '{self}'"
            raise ValueError(msg)
        self.instances[new_label] = self.instances.pop(cur_label)

    @property
    @abstractmethod
    def limits(self) -> Tuple[float]:
        """
        Gives the curve's parametric interval

        :getter: Returns the pair [a, b] in which curve is parametric defined
        :type: Tuple[float]

        """
        raise NotImplementedError

    @property
    @abstractmethod
    def knots(self) -> Tuple[float]:
        """
        Gives the curve's knots, in which the parametric interval is divided

        :getter: Returns the knots that divides the curve's interval
        :type: Tuple[float]

        """
        raise NotImplementedError

    @abstractmethod
    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Evaluates the curves at given parameters.

        :param parameters: A vector-like of lenght n
        :type parameters: Tuple[float]
        :return: A matrix of shape (n, 2)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError


class NurbsCurve(Curve):
    """
    Nurbs Curve instance, a child of "Curve" that serves as interface
    and uses the functions/methods from compmec-nurbs package
    """

    @classmethod
    def from_dict(cls, dictionary: Dict) -> NurbsCurve:
        degree = dictionary["degree"] if "degree" in dictionary else None
        knotvector = nurbs.KnotVector(dictionary["knotvector"], degree)
        nurbs_curve = nurbs.Curve(knotvector)
        ctrlpoints = tuple(dictionary["ctrlpoints"])
        nurbs_curve.ctrlpoints = np.array(ctrlpoints, dtype="float64")
        if "weights" in dictionary:
            nurbs_curve.weights = dictionary["weights"]
        return cls(nurbs_curve)

    def to_dict(self) -> Dict:
        dictionary = OrderedDict()
        dictionary["degree"] = self.internal.degree
        dictionary["knotvector"] = self.internal.knotvector
        dictionary["ctrlpoints"] = self.internal.ctrlpoints
        weights = self.internal.weights
        if weights is not None:
            dictionary["weights"] = weights
        return dictionary

    def __new__(cls, nurbs_curve: nurbs.Curve, label: Optional[int] = None):
        if not isinstance(nurbs_curve, nurbs.Curve):
            msg = "Invalid internal curve"
            raise TypeError(msg)
        instance = super().__new__(cls, label)
        instance.internal = nurbs_curve
        return instance

    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        return self.internal.eval(parameters)

    @property
    def knots(self) -> Tuple[float]:
        return self.internal.knotvector.knots

    @property
    def limits(self) -> Tuple[float]:
        return self.internal.knotvector.limits

    @property
    def internal(self) -> nurbs.Curve:
        """
        Gives the internal nurbs.Curve object

        :getter: Returns the curve's label
        :setter: Sets the new curve instance
        :type: compmec.nurbs.Curve
        """
        return self.__internal

    @internal.setter
    def internal(self, new_curve: nurbs.Curve):
        if not isinstance(new_curve, nurbs.Curve):
            msg = "Invalid internal curve"
            raise TypeError(msg)
        self.__internal = new_curve
