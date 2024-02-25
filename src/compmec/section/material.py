"""
This module contains the class 'Isotropic' to store and convert values
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Dict, Optional, Tuple


class Material(ABC):
    """
    Material abstract class, the parent of other more specific materials
    """

    instances = OrderedDict()

    @staticmethod
    def clear(names: Optional[Tuple[str]] = None):
        """
        Removes all given instances of Curve
        """
        if names is None:
            Material.instances.clear()
            return
        for name in names:
            if name in Material.instances:
                Material.instances.pop(name)

    @staticmethod
    def new_instance(tipo: str, dictionary: Dict) -> Material:
        """
        Creates a new instance of Material depending on
        the 'tipo' and the informations from the dictionary

        :param tipo: The Material-subclass to be called, in ["isotropic"]
        :type tipo: str
        :return: The created material instance
        :rtype: Material
        """
        tipo2class = {"isotropic": Isotropic}
        if tipo not in tipo2class:
            raise NotImplementedError
        return tipo2class[tipo].from_dict(dictionary)

    @staticmethod
    def __next_available_name() -> str:
        index = 1
        while True:
            name = f"custom-material-{index}"
            if name not in Material.instances:
                return name
            index += 1

    def __init__(self, name: Optional[str] = None):
        if name is None:
            name = Material.__next_available_name()
        elif name in Material.instances:
            msg = f"Cannot create material '{name}'! There's already one!"
            raise ValueError(msg)
        self.__name = name
        self.__class__.instances[name] = self

    @abstractmethod
    def to_dict(self) -> Dict:
        """
        Converts the material to a dictionary
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_dict(cls, dictionary: Dict) -> Material:
        """
        Converts the dictionary in a material instance
        """
        raise NotImplementedError

    @property
    def name(self) -> str:
        """
        Gives the material name

        :getter: Returns the material's name
        :setter: Attribuates a new name for material
        :type: str

        """
        return self.__name

    @name.setter
    def name(self, new_name: str):
        if self.name == new_name:
            return
        if new_name in self.instances:
            msg = f"Cannot set the name '{new_name}' for the material "
            msg += f"'{self.name}' cause there's already a material with "
            msg += "the same name!\n"
            msg += f"Cur: '{self.instances[new_name]}'\n"
            msg += f"New: '{self}'"
            raise ValueError(msg)
        self.instances[new_name] = self.instances.pop(self.name)


class Isotropic(Material):
    """
    Isotropic materials are materials whose properties remain
    the same when tested in different directions.
    Wikipedia link of isotropy:
    https://en.wikipedia.org/wiki/Isotropy
    """

    @staticmethod
    def __verify_young_modulus(young_modulus):
        if young_modulus <= 0:
            raise ValueError

    @staticmethod
    def __verify_poissons_ratio(poissons_ratio):
        if poissons_ratio < 0.49:
            pass
        elif poissons_ratio < 0.50:
            raise ValueError("Material is incompressible")
        else:
            raise ValueError("Cannot have poisson >= 0.50 ")

    @classmethod
    def from_dict(cls, dictionary: Dict) -> Isotropic:
        return cls(
            young_modulus=dictionary["young_modulus"],
            poissons_ratio=dictionary["poissons_ratio"],
        )

    def to_dict(self) -> Dict:
        dictionary = OrderedDict()
        dictionary["young_modulus"] = self.young_modulus
        dictionary["poissons_ratio"] = self.poissons_ratio
        return dictionary

    def __init__(
        self,
        *,
        young_modulus: float,
        poissons_ratio: float,
        name: Optional[str] = None,
    ):
        Isotropic.__verify_young_modulus(young_modulus)
        Isotropic.__verify_poissons_ratio(poissons_ratio)
        super().__init__(name)
        self.__young_modulus = young_modulus
        self.__poissons_ratio = poissons_ratio

    def __str__(self) -> str:
        values = [
            self.young_modulus,
            self.poissons_ratio,
            self.bulk_modulus,
            self.shear_modulus,
            self.lame_parameter_1,
            self.lame_parameter_2,
        ]
        # values = tuple(".3f" % v for v in values)
        msg = "Isotropic Material: "
        msg += "(E, nu, K, G, L1, L2) = (%s)"
        msg %= ", ".join(map(str, values))
        return msg

    @property
    def young_modulus(self):
        """
        Young modulus (E) is a mechanical property that measures
        the tensile or compressive stiffness of a solid material
        when the force is applied lengthwise
        https://en.wikipedia.org/wiki/Young%27s_modulus
        """
        return self.__young_modulus

    @property
    def poissons_ratio(self):
        """
        Poisson's ratio (nu) is a measure of the Poisson effect,
        the deformation (expansion or contraction) of a material
        in directions perpendicular to the specific direction of loading
        https://en.wikipedia.org/wiki/Poisson%27s_ratio
        """
        return self.__poissons_ratio

    @property
    def bulk_modulus(self):
        """
        The bulk modulus (K) of a substance is a measure
        of the resistance of a substance to compression.
        https://en.wikipedia.org/wiki/Bulk_modulus
        """
        return self.young_modulus / (3 * (1 - 2 * self.poissons_ratio))

    @property
    def shear_modulus(self):
        """
        Shear modulus (G), is a measure of the elastic shear
        stiffness of a material and is defined as the ratio
        of shear stress to the shear strain.
        https://en.wikipedia.org/wiki/Shear_modulus
        """
        return self.young_modulus / (2 * (1 + self.poissons_ratio))

    @property
    def lame_parameter_1(self):
        """
        Lamé parameters (also called the Lamé coefficients)
        are two material-dependent quantities denoted by lambda
        and mu that arise in strain-stress relationships.
        https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
        """
        shear = self.shear_modulus
        poisson = self.poissons_ratio
        return 2 * shear * poisson / (1 - 2 * poisson)

    @property
    def lame_parameter_2(self):
        """
        Lamé parameters (also called the Lamé coefficients)
        are two material-dependent quantities denoted by lambda
        and mu that arise in strain-stress relationships.
        https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
        """
        return self.shear_modulus
