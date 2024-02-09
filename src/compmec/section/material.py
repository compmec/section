"""
This module contains the class 'Isotropic' to store and convert values
"""
from __future__ import annotations

from importlib import resources
from typing import Dict, Optional, Tuple, Union

from . import dataio
import jsonschema


class Material(object):
    @classmethod
    def from_json(cls, filepath: str) -> Dict[str, Material]:
        data = dataio.read_material_json(filepath)
        materials = {}
        for name, info in data.items():
            materials[name] = cls.from_dict(info)
        return materials

    @classmethod
    def from_dict(cls, info: Dict) -> Material:
        material = Isotropic(info)
        return material

    def to_dict(self) -> Dict:
        """
        Converts the current material to a dictionary
        """
        dicion = {}
        for name in dir(self):
            if name.startswith("__"):
                continue
            if name.endswith("__"):
                continue
            value = getattr(self, name)
            if not callable(value):
                dicion[name] = value
        return dicion

    def to_json(self, json_filepath: str, name: Optional[str] = "unknown"):
        material_infos = self.to_dict()
        materials_dictionary = {name: material_infos}
        saving_dictionary = {"materials": materials_dictionary}
        dataio.save_json(saving_dictionary, json_filepath)


class Isotropic(Material):
    """
    Isotropic materials are materials whose properties remain
    the same when tested in different directions.
    Wikipedia link of isotropy:
    https://en.wikipedia.org/wiki/Isotropy
    """

    def __init__(self, **kwargs):
        self._young_modulus = None
        self._poissons_ratio = None
        self._density = None
        if kwargs is not None:
            for attr, value in kwargs.items():
                setattr(self, attr, value)

    def __str__(self) -> str:
        values = [
            self.young_modulus,
            self.poissons_ratio,
            self.bulk_modulus,
            self.shear_modulus,
            self.lame_parameter_1,
            self.lame_parameter_2,
            self.density,
        ]
        # values = tuple(".3f" % v for v in values)
        msg = "Isotropic Material: "
        msg += "(E, nu, K, G, L1, L2, rho) = (%s)"
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
        return self._young_modulus

    @property
    def poissons_ratio(self):
        """
        Poisson's ratio (nu) is a measure of the Poisson effect,
        the deformation (expansion or contraction) of a material
        in directions perpendicular to the specific direction of loading
        https://en.wikipedia.org/wiki/Poisson%27s_ratio
        """
        return self._poissons_ratio

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

    @property
    def density(self):
        """
        Density (rho) is the mesure of mass/volum
        https://en.wikipedia.org/wiki/Density
        """
        return self._density

    @young_modulus.setter
    def young_modulus(self, value: float):
        if value <= 0:
            raise ValueError
        self._young_modulus = value

    @poissons_ratio.setter
    def poissons_ratio(self, value: float):
        if value < 0.49:
            pass
        elif value < 0.50:
            raise ValueError("Material is incompressible")
        else:
            raise ValueError("Cannot have poisson >= 0.50 ")
        self._poissons_ratio = value

    @density.setter
    def density(self, value: float):
        assert value > 0
        self._density = value
