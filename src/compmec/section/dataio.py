"""
This file is responsible for INPUT and OUTPUT of data

Mainly the the two types of files are : JSON and VTK/VTU
"""

import json
from abc import ABC, abstractmethod
from importlib import resources
from typing import Dict, Optional

import jsonschema

from .curve import Curve, Nodes
from .material import Material
from .section import Section


class FileIO(ABC):
    """
    Abstract Reader class that serves as basic interface to read/save
    informations from/to given file.
    The file format is defined by child class.
    """

    def __init__(self, filepath: str):
        assert isinstance(filepath, str)
        self.__filepath = filepath

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
    def open(self, mode: str = "r"):
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


class JsonIO(FileIO):
    """
    JsonIO class that serves as interface between the file
    and the data structures used in this packaged
    """

    def __init__(self, filepath: str):
        super().__init__(filepath)
        self.__opened = False

    def is_open(self) -> bool:
        return self.__opened

    def open(self, mode: str = "r"):
        self.__opened = True

    def close(self):
        self.__opened = False

    def read_json(self, schemapath: Optional[str] = None) -> Dict:
        """
        Reads a json file and returns the data inside it.
        If `schema` is given, it verifies if the `filepath`
        meets the given standard, by `jsonschema`

        If `filepath` doesn't satisfy the `schema`,
        raises the error `jsonschema.exceptions.ValidationError`

        Parameters
        -----------
        schemapath: str, optional
            schema filepath to be checked
        return: dict
            The dictionary with all infos read from json

        """
        with open(self.filepath, "r", encoding="ascii") as file:
            data = json.load(file)
        if schemapath:
            if not isinstance(schemapath, str):
                raise TypeError
            with open(schemapath, "r", encoding="ascii") as file:
                schema = json.load(file)
            jsonschema.validate(data, schema)
        return data

    def save_json(self, dictionary: Dict):
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
        with open(self.filepath, "w", encoding="ascii") as file:
            file.write(json_object)

    def read_nodes(self):
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
        matrix = self.read_json(schema_path)["nodes"]
        Nodes.insert_matrix(matrix)

    def read_curves(self):
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
        curves = self.read_json(schema_path)["curves"]
        for label, infos in curves.items():
            label = int(label)
            curve = Curve.new_instance("nurbs", infos)
            curve.label = label

    def read_materials(self):
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
        materials = self.read_json(schema_path)["materials"]
        for name, info in materials.items():
            material = Material.new_instance("isotropic", info)
            material.name = name

    def read_sections(self):
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
        sections = self.read_json(schema_path)["sections"]
        for name, info in sections.items():
            section = Section.from_dict(info)
            section.name = name
