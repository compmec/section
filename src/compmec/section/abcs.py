"""
File that contains base classes
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Tuple, Union


class Tracker(ABC):  # pylint: disable=too-few-public-methods
    """
    Parent class to track all the instances
    """

    instances = {}

    @classmethod
    def clear(cls, keys: Optional[Tuple[int]] = None):
        """
        Remove node labels
        """
        if keys is None:
            cls.instances.clear()
            return
        for key in keys:
            if key in cls.instances:
                cls.instances.pop(key)

    @classmethod
    @abstractmethod
    def _next_available_key(cls) -> Any:
        """
        Protected method that gives next key
        """
        raise NotImplementedError

    def _get_key(self) -> Any:
        """
        Protected method that gives the key of the instance by iterating
        on the instances dictionary

        Parameters
        ----------

        :return: The instances key from the dictionary
        :rtype: Any

        """
        for key, instance in self.instances.items():
            if id(instance) == id(self):
                return key
        return None

    def _set_key(self, new_key: Any):
        """
        Protected method that sets a new key for the instance.
        Iterates over the dictionary and changes the name.

        If `None` is given, it will attribute the next available key

        Parameters
        ----------

        :param new_key: The new instance key"
        :type new_key: Any

        """
        cur_key = self._get_key()
        if new_key is None:
            new_key = self._next_available_key()
        elif cur_key == new_key:
            return
        elif new_key in self.instances:
            msg = f"Key {new_key} is already used"
            raise ValueError(msg)
        elif cur_key is not None:
            self.instances.pop(cur_key)
        self.instances[new_key] = self


class LabeledTracker(Tracker):
    """
    Labeled Base Tracker, to track instances of Curve and Node
    """

    @classmethod
    def _next_available_key(cls) -> str:
        """
        Protected method that gives next key
        """
        index = 1
        while index in cls.instances:
            index += 1
        return index

    @property
    def label(self) -> int:
        """
        Gives the instance label

        :getter: Returns the instance's label
        :setter: Attribuates a new label for instance
        :type: str

        """
        return int(self._get_key())

    @label.setter
    def label(self, new_label: Union[int, None]):
        if new_label is not None:
            new_label = int(new_label)
        self._set_key(new_label)


class NamedTracker(Tracker):
    """
    Named Base Tracker, to track instances of Material, Geometry and Section
    """

    @classmethod
    def _next_available_key(cls) -> str:
        """
        Protected method that gives next key
        """
        index = 1
        while f"instance-{index}" in cls.instances:
            index += 1
        return f"instance-{index}"

    @property
    def name(self) -> str:
        """
        Gives the instance name

        :getter: Returns the instance's name
        :setter: Attribuates a new name for instance
        :type: str

        """
        return str(self._get_key())

    @name.setter
    def name(self, new_name: Union[str, None]):
        if new_name is not None:
            new_name = str(new_name)
        self._set_key(new_name)


class ICurve(ABC):
    """
    Interface abstract curve to be parent of others.

    This class serves as interface between the curves from others packaged
    like pynurbs.Curve, to this package, expliciting the minimum requirements
    of a curve must have.
    With this, it's possible to implement your own type of parametric curve
    """

    @property
    @abstractmethod
    def knots(self) -> Tuple[float]:
        """
        Gives the curve's knots, in which the parametric interval is divided

        :getter: Returns the knots that divides the curve's interval
        :type: Tuple[float]

        """
        raise NotImplementedError

    @property
    @abstractmethod
    def degree(self) -> int:
        """
        Gives the curve's polynomial degree

        Degree 1 means the curve is linear/polygonal

        :getter: Returns the curve's polynomial degree
        :type: int

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

    @abstractmethod
    def deval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Evaluates the derivative of curve at given parameters.

        :param parameters: A vector-like of lenght n
        :type parameters: Tuple[float]
        :return: A matrix of shape (n, 2)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError

    @abstractmethod
    def winding(self, point: Tuple[float]) -> float:
        """
        Computes the winding number of the given curve

        The possible results are one in the interval [0, 1]

        * 0 means the point is outside the internal curve
        * 1 means the point is inside the internal curve
        * in (0, 1) means it's on the boundary

        :param point: A 2D-point
        :type point: Tuple[float]
        :return: The winding value with respect to the curve
        :rtype: float

        """
        raise NotImplementedError

    @abstractmethod
    def projection(self, point: Tuple[float]) -> Tuple[float]:
        """
        Finds the parameter t* such the distance between the
        `curve(t*)` and `point` is the nearest possible.

        :param parameters: The point to project
        :type parameters: Tuple[float]
        :return: A parameters that minimizes the distance
        :rtype: Tuple[float]
        """
        raise NotImplementedError


class IGeometry(ABC):
    """
    Connected Geometry class that represents a bidimensional region
    on the plane and has some functions such as to decide if points
    are inside the region
    """

    def integrate(
        self, expx: int, expy: int, tolerance: Optional[float] = 1e-9
    ) -> float:
        """
        Evaluates the integral

        I = int_{Omega} x^expx * y^expy dOmega


        """
        raise NotImplementedError

    def winding(self, point: Tuple[float]) -> float:
        """
        Computes the winding number of the giving point

        Normally winding number is one for each curve, but this method gives
        the overall winding number
        """
        raise NotImplementedError


class IMaterial(ABC):
    """
    Material abstract class, the parent of other more specific materials
    """

    @abstractmethod
    def to_dict(self) -> Dict:
        """
        Converts the material to a dictionary
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_dict(cls, dictionary: Dict) -> IMaterial:
        """
        Converts the dictionary in a material instance
        """
        raise NotImplementedError


class IFileIO(ABC):
    """
    Abstract Reader class that serves as basic interface to read/save
    informations from/to given file.
    The file format is defined by child class.
    """

    @abstractmethod
    def is_open(self) -> bool:
        """
        Tells if the reader is open
        """
        raise NotImplementedError

    @abstractmethod
    def open(self, mode: str = "r"):
        """
        Opens the file with given mode
        """
        raise NotImplementedError

    @abstractmethod
    def close(self):
        """
        Closes the file
        """
        raise NotImplementedError

    @abstractmethod
    def load_nodes(self):
        """
        Saves all the nodes from file into Node class
        """
        raise NotImplementedError

    @abstractmethod
    def load_curves(self):
        """
        Creates all the curves instances from file
        """
        raise NotImplementedError

    @abstractmethod
    def load_materials(self):
        """
        Creates all the materials instances from file
        """
        raise NotImplementedError

    @abstractmethod
    def load_sections(self):
        """
        Creates all the sections instances from file
        """
        raise NotImplementedError


class IPoissonEvaluator(ABC):
    """
    This class helps evaluating the function 'u' which satisfies
    the poissons equation:

    nabla^2 u = 0

    For this class, it's required the tangent and normal functions
    We know if it's on the boundary, we use the directly the function
    on the boundary.
    If it's in the interior, we use the equation:

    alpha(s) * u(s) = int_{Gamma} u * dv/dn * dGamma
                    - int_{Gamma} v * du/dn dGamma

    dv/dn * ds = r x p' / <r, r> * dt
    v * ds = ln |r| * |p'| * dt

    """

    @classmethod
    def eval(self, source: Tuple[float]) -> float:
        """
        Evaluates the function 'u' at given point

        :param source: The point where to evaluate
        :type source: Tuple[float]
        :return: The value of the function at given point
        :rtype: float

        """
        raise NotImplementedError

    @classmethod
    def grad(self, source: Tuple[float]) -> Tuple[float]:
        """
        Evaluates the gradient of 'u' at given point

        :param source: The point where to evaluate
        :type source: Tuple[float]
        :return: The gradient of the function at given point
        :rtype: Tuple[float]
        """
        raise NotImplementedError


class IField(ABC):
    """
    This is a base abstract class parent of others

    It's responsible to decide if a given point is
    inside/outside of given section.
    """

    @property
    @abstractmethod
    def ndata(self) -> int:
        """
        Gives the quantity of output numbers for one point

        :getter: Returns the lenght of output data
        :type: int
        """
        raise NotImplementedError

    @abstractmethod
    def eval(self, points: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
        """
        Evaluate the field at given points

        :param points: The wanted points, a matrix of shape (n, 2)
        :type points: Tuple[Tuple[float]]
        :return: The results in a matrix of shape (n, ndata)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError


class ISection(ABC):  # pylint: disable=too-few-public-methods
    """
    Section abstract class to serve as interface
    """

    @abstractmethod
    def area(self) -> float:
        """
        Gives the cross-section area

        A = int 1 dx dy

        :return: The value of cross-section area
        :rtype: float

        Example use
        -----------

        >>> section.area()
        1.0

        """
        raise NotImplementedError

    @abstractmethod
    def first_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        """Gives the first moments of area

        Qx = int (y-cy) dx dy
        Qy = int (x-cx) dx dy

        If no ``center`` is given, it assumes the origin (0, 0)

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The first moment of inertia (Qx, Qy)
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section.first_moment()
        (0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def second_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        """Gives the second moment of inertia with respect to ``center``

        Ixx = int (y-cy)^2 dx dy
        Ixy = int (x-cx)*(y-cy) dx dy
        Iyy = int (x-cx)^2 dx dy

        If no ``center`` is given, it assumes the origin (0, 0) and
        returns the global second moment of inertia

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The values of Ixx, Ixy, Iyy
        :rtype: tuple[float, float, float]

        Example use
        -----------

        >>> section.second_moment()
        (1., 0., 1.)

        """
        raise NotImplementedError

    @abstractmethod
    def third_moment(self, center: Tuple[float] = (0, 0)) -> Tuple[float]:
        """Gives the third moment of inertia with respect to ``center``

        Ixxx = int (y-cy)^3 dx dy
        Ixxy = int (x-cx)*(y-cy)^2 dx dy
        Ixyy = int (x-cx)^2*(y-cy) dx dy
        Iyyy = int (x-cx)^3 dx dy

        If no ``center`` is given, it assumes the origin (0, 0)

        :param center: The center to compute second moment, default (0, 0)
        :type center: tuple[float, float]
        :return: The values of Ixxx, Ixxy, Ixyy, Iyyy
        :rtype: tuple[float, float, float, float]

        Example use
        -----------

        >>> section.third_moment()
        (0., 0., 0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def geometric_center(self) -> Tuple[float]:
        """Gives the geometric center G

        G = (x_gc, y_gc)
        x_gc = (1/A) * Qy
        y_gc = (1/A) * Qx

        This center depends only on the geometry,
        not on the material

        :return: The value of geometric center G
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.geometric_center()
        (0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def bending_center(self) -> Tuple[float]:
        """Gives the bendin center B

        The bending center is the point of the
        intersection of two neutral lines, where
        the stress and strain are always zero

        :return: The value of bending center B
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section = Section(shapes, materials)
        >>> section.bending_center()
        (0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def torsion_center(self) -> Tuple[float]:
        """Gives the torsion center T

        The torsion center is the point which,
        when applied a torsion moment, the shear
        stresses at this point are zero

        :return: The value of torsion center T
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section.torsion_center()
        (0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def shear_center(self) -> Tuple[float]:
        """Gives the shear center S

        The shear center is the point which,
        when applied a transverse force, this
        force doesn't cause torsion

        :return: The value of shear center S
        :rtype: tuple[float, float]

        Example use
        -----------

        >>> section.shear_center()
        (0., 0.)

        """
        raise NotImplementedError

    @abstractmethod
    def gyradius(self) -> Tuple[float]:
        """Gives the gyradius (radii of gyration)

        R = (sqrt(Ixx/A), sqrt(Iyy/A))

        """
        raise NotImplementedError

    @abstractmethod
    def torsion_constant(self) -> float:
        """Gives the torsion constant J

        J = Ixx + Iyy - Jw

        Careful: This function solves a linear system

        :return: The value of torsion constant J
        :rtype: float

        Example use
        -----------

        >>> section.torsion_constant()
        1.

        """
        raise NotImplementedError

    @abstractmethod
    def charged_field(self) -> IField:
        """
        Gives the charged field instance to evaluate stresses

        :return: The field evaluator
        :rtype: ChargedField

        """
        raise NotImplementedError


class IBasisFunc(ABC):
    """
    Basis functions abstract class
    """

    @property
    @abstractmethod
    def ndofs(self) -> Tuple[float]:
        """
        Gives the number of the basis functions

        :getter: Returns the number of the basis functions
        :type: int
        """
        raise NotImplementedError

    @property
    @abstractmethod
    def knots(self) -> Tuple[float]:
        """
        Gives the dividing knots

        :getter: Returns the knots of the basis functions
        :type: Tuple[float]
        """
        raise NotImplementedError

    @abstractmethod
    def eval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Evaluate the (m, ) basis functions at given parameters

        :param parameters: The (n, ) parameter values to be evaluated
        :type parameters: Tuple[float]
        :return: The results in a matrix of shape (m, n)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError

    @abstractmethod
    def deval(self, parameters: Tuple[float]) -> Tuple[Tuple[float]]:
        """
        Evaluate the (m, ) derivative basis functions at given parameters

        :param parameters: The (n, ) parameter values to be evaluated
        :type parameters: Tuple[float]
        :return: The results in a matrix of shape (m, n)
        :rtype: Tuple[Tuple[float]]
        """
        raise NotImplementedError
