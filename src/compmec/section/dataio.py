"""
This file is responsible for INPUT and OUTPUT of data

Mainly the the two types of files are : JSON and VTK/VTU
"""

import json
from importlib import resources
from typing import Dict, Optional

import jsonschema


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
    with open(filepath, "r") as file:
        data = json.load(file)
    if schemapath:
        if not isinstance(schemapath, str):
            raise TypeError
        with open(schemapath, "r") as file:
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
    return data


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
    for key in data.keys():
        if key != "materials":
            data.pop(key)
    return data


def read_curve_json(filepath: str) -> Dict:
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
    for key in data.keys():
        if key not in ["nodes", "curves"]:
            data.pop(key)
    return data


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
    with open(json_filepath, "w") as file:
        file.write(json_object)
