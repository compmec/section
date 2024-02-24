"""
This file is responsible for INPUT and OUTPUT of data

Mainly the the two types of files are : JSON and VTK/VTU
"""

import json
from abc import ABC, abstractmethod
from collections import OrderedDict
from importlib import resources
from typing import Dict, Optional, Tuple

import jsonschema

from .curve import Curve, Nodes
from .material import Material
from .section import Section


class FileReader(ABC):
    """
    Abstract Reader class that serves as basic interface to read a file.
    The read file format is defined by child.
    """

    def __init__(self, filepath: str):
        assert isinstance(filepath, str)
        self.__filepath = filepath
        self.file = None

    @property
    def filepath(self) -> str:
        """
        Gives the json filepath

        :getter: Returns the json filepath
        :type: str

        """
        return self.__filepath

    @abstractmethod
    def is_open(self) -> bool:
        """
        Tells if the reader is open
        """
        raise NotImplementedError

    @abstractmethod
    def open(self):
        """
        Opens the file
        """
        raise NotImplementedError

    @abstractmethod
    def close(self):
        """
        Closes the file
        """
        raise NotImplementedError

    @abstractmethod
    def read_nodes(self):
        """
        Saves all the nodes from file into Nodes class
        """
        raise NotImplementedError

    @abstractmethod
    def read_curves(self):
        """
        Creates all the curves instances from file
        """
        raise NotImplementedError

    @abstractmethod
    def read_materials(self):
        """
        Creates all the materials instances from file
        """
        raise NotImplementedError

    @abstractmethod
    def read_sections(self):
        """
        Creates all the sections instances from file
        """
        raise NotImplementedError

    def read(self):
        """
        Read the file and create all instaces
        """
        self.read_nodes()
        self.read_curves()
        self.read_materials()
        self.read_sections()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class JsonIO(FileReader):
    """
    JsonIO class that serves as interface between the file
    and the data structures used in this packaged
    """


def read_json(filepath: str, schemapath: Optional[str] = None) -> Dict:
    """
    Reads a json file and returns the data inside it.
    If `schema` is given, it verifies if the `filepath`
    meets the given standard, by `jsonschema`

    If `filepath` doesn't satisfy the `schema`,
    raises the error `jsonschema.exceptions.ValidationError`

    Parameters
    -----------
    filepath: str
        json filepath to be read from
    schemapath: str, optional
        schema filepath to be checked
    return: dict
        The dictionary with all infos read from json

    """
    if not isinstance(filepath, str):
        raise TypeError
    with open(filepath, "r", encoding="ascii") as file:
        data = json.load(file)
    if schemapath:
        if not isinstance(schemapath, str):
            raise TypeError
        with open(schemapath, "r", encoding="ascii") as file:
            schema = json.load(file)
        jsonschema.validate(data, schema)
    return data


def read_section_json(filepath: str) -> Dict:
    """
    Reads a section json file and returns the data inside it.

    This file must be in accordance with the schema `section.json`

    Parameters
    ----------
    filepath: str
        section json filepath to be read from
    return: dict
        The dictionary with all infos read from json

    """
    schema_name = "schema/section.json"
    folder = resources.files("compmec.section")
    schema_path = str(folder.joinpath(schema_name))
    data = read_json(filepath, schema_path)
    return data["sections"]


def read_material_json(filepath: str) -> Dict:
    """
    Reads a material json and returns the data inside it.

    This file must be in accordance with the schema `material.json`

    Parameters
    ----------
    filepath: str
    return: dict
    """
    schema_name = "schema/material.json"
    folder = resources.files("compmec.section")
    schema_path = str(folder.joinpath(schema_name))
    data = read_json(filepath, schema_path)
    return data["materials"]


def read_curve_json(filepath: str) -> Tuple[Dict]:
    """
    Reads a curve json and returns the data inside it.

    This file must be in accordance with the schema `curve.json`

    Parameters
    ----------
    filepath: str
    return: dict
    """
    schema_name = "schema/curve.json"
    folder = resources.files("compmec.section")
    schema_path = str(folder.joinpath(schema_name))
    data = read_json(filepath, schema_path)
    curves = OrderedDict()
    for label, infos in data["curves"].items():
        curves[int(label)] = infos
    return curves


def read_nodes_json(filepath: str) -> Tuple[Dict]:
    """
    Reads a nodes json and returns the data inside it.

    This file must be in accordance with the schema `curve.json`

    Parameters
    ----------
    filepath: str
    return: dict
    """
    schema_name = "schema/curve.json"
    folder = resources.files("compmec.section")
    schema_path = str(folder.joinpath(schema_name))
    data = read_json(filepath, schema_path)
    curves = OrderedDict()
    for label, infos in data["curves"].items():
        curves[int(label)] = infos
    return data["nodes"]


def save_json(dictionary: Dict, json_filepath: str):
    """
    Saves the given dictionary in a json file.

    For now this function overwrides all the file.
    It would be nice to add the informations, keeping
    the old values (if not conflitant)

    Parameters
    ----------
    dictionary: Dict
    json_filepath: str
        The path to save the informations
    """
    json_object = json.dumps(dictionary, indent=4)
    with open(json_filepath, "w", encoding="ascii") as file:
        file.write(json_object)


def load_json(json_filepath: str):
    """
    Loads all the informations from json file,
    and create classes: Curve, Material and Section

    Parameters
    ----------
    json_filepath: str
        The path to load the informations
    """
    matrix = read_nodes_json(json_filepath)
    Nodes.insert_matrix(matrix)

    curves = read_curve_json(json_filepath)
    for label, info in curves.items():
        curve = Curve.new_instance("nurbs", info)
        curve.label = label

    materials = read_material_json(json_filepath)
    for name, info in materials.items():
        material = Material.new_instance("isotropic", info)
        material.name = name

    sections = read_section_json(json_filepath)
    for name, info in sections.items():
        section = Section.from_dict(info)
        section.name = name
