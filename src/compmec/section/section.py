from typing import Optional, Tuple

import numpy as np
from compmec.shape.shape import DefinedShape, IntegrateShape

from compmec.section.material import Material


class Section:
    def __init__(self, shapes: Tuple[DefinedShape], materials: Tuple[Material]):
        for shape in shapes:
            if not isinstance(shape, DefinedShape):
                raise TypeError
        for material in materials:
            if not isinstance(material, Material):
                raise TypeError
        self.__shapes = tuple(shapes)
        self.__materials = tuple(materials)
        self.__area = None
        self.__first = None
        self.__second = None
        self.__weight = True

    def __iter__(self) -> Tuple[DefinedShape, Material]:
        for shape, material in zip(self.__shapes, self.__materials):
            yield (shape, material)

    @property
    def weight(self) -> bool:
        return self.__weight

    @weight.setter
    def weight(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError
        self.__weight = value

    def area(self) -> float:
        """
        Computes the area of the section
        If weight is active, it computes
        """
        if self.__area is None:
            areas = tuple(map(IntegrateShape.area, self.__shapes))
            if self.weight:
                weights = tuple(mat.young_modulus for mat in self.__materials)
                self.__area = np.inner(areas, weights) / max(weights)
            else:
                self.__area = sum(areas)
        return self.__area

    def first_momentum(self) -> Tuple[float]:
        if self.__first is None:
            nshapes = len(self.__shapes)
            values = np.zeros((nshapes, 2), dtype="float64")
            for i, shape in enumerate(self.__shapes):
                values[i, 0] = IntegrateShape.polynomial(shape, 1, 0)
                values[i, 1] = IntegrateShape.polynomial(shape, 0, 1)
            if self.weight:
                weights = tuple(mat.young_modulus for mat in self.__materials)
                self.__first = np.einsum("ij,i->j", values, weights) / max(weights)
            else:
                self.__first = np.einsum("ij->j", values)
        return self.__first

    def second_momentum(self) -> Tuple[Tuple[float]]:
        if self.__second is None:
            nshapes = len(self.__shapes)
            values = np.zeros((nshapes, 2, 2), dtype="float64")
            for i, shape in enumerate(self.__shapes):
                values[i, 0, 0] += IntegrateShape.polynomial(shape, 2, 0)
                values[i, 0, 1] += IntegrateShape.polynomial(shape, 1, 1)
                values[i, 1, 1] += IntegrateShape.polynomial(shape, 0, 2)
            values[:, 1, 0] = values[:, 1, 0]
            if self.weight:
                weights = tuple(mat.young_modulus for mat in self.__materials)
                self.__second = np.einsum("ijk,i->jk", values, weights) / max(weights)
            else:
                self.__second = np.einsum("ijk->jk", values)
        return self.__second

    def local_second_momentum(self) -> Tuple[Tuple[float]]:
        area = self.area()
        first = self.first_momentum()
        second = self.second_momentum()
        return second - np.tensordot(first, first, axes=0) / area

    def torsion_constant(self) -> float:
        raise NotImplementedError

    def gyradius(self) -> Tuple[float]:
        area = self.area()
        second = self.second_momentum()
        return np.sqrt(np.diag(second) / area)

    def elastic_modulus(self) -> Tuple[Tuple[float]]:
        raise NotImplementedError

    def plastic_modulus(self) -> Tuple[float]:
        raise NotImplementedError

    def shear_center(self) -> Tuple[float]:
        raise NotImplementedError

    def warping_constant(self) -> float:
        raise NotImplementedError

    def monosymmetry_constants(self) -> Tuple[float]:
        raise NotImplementedError
