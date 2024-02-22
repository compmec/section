"""
This module contains the class 'Isotropic' to store and convert values
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import OrderedDict
from typing import Dict, Optional


class Material(ABC):
    """
    Material abstract class, the parent of other more specific materials
    """

    instances = OrderedDict()

    @staticmethod
    def __next_available_name() -> str:
        index = 1
        while True:
            name = f"custom-material-{index}"
            if name not in Material.instances:
                return name

    def __new__(cls, name: Optional[str] = None) -> Material:
        if name is None:
            name = Material.__next_available_name()
        elif name in Material.instances:
            msg = f"Cannot create material '{name}'! There's already one!"
            raise ValueError(msg)
        instance = super().__new__(cls)
        instance.name = name
        Material.instances[name] = instance
        return instance

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
        try:
            cur_name = self.name
        except AttributeError:
            self.__name = new_name
            return
        if cur_name == new_name:
            return
        if new_name in self.instances:
            msg = f"Cannot set the name '{new_name}' for the material "
            msg += f"'{cur_name}' cause there's already a material with "
            msg += "the same name!\n"
            msg += f"Cur: '{self.instances[new_name]}'\n"
            msg += f"New: '{self}'"
            raise ValueError(msg)
        self.instances[new_name] = self.instances.pop(cur_name)


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

    def __init__(self, young_modulus: float, poissons_ratio: float):
        Isotropic.__verify_young_modulus(young_modulus)
        Isotropic.__verify_poissons_ratio(poissons_ratio)
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
