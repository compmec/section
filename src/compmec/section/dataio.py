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

    :param filepath: json filepath to be read from
    :type filepath: str
    :param schemapath: schema filepath to be checked
    :type schemapath: str, optional
    :return: The dictionary with all infos read from json
    :rtype: dict

    """
    if not isinstance(filepath, str):
        raise TypeError
    if not isinstance(schemapath, str):
        raise TypeError
    with open(filepath, "r") as file:
        data = json.load(file)
    if schemapath:
        with open(schemapath, "r") as file:
            schema = json.load(file)
        jsonschema.validate(data, schema)
    return data


def read_section_json(filepath: str) -> Dict:
    """
    Reads a section json file and returns the data inside it.

    This file must be in accordance with the schema `section.json`

    Parameters
    -----------

    :param filepath: section json filepath to be read from
    :type filepath: str
    :return: The dictionary with all infos read from json
    :rtype: dict

    """
    schema_name = "schema/section.json"
    folder = resources.files("compmec.section")
    schema_path = str(folder.joinpath(schema_name))
    data = read_json(filepath, schema_path)
    return data
