"""
This file is responsible for INPUT and OUTPUT of data

Mainly the the two types of files are : JSON and VTK/VTU
"""

import json
from typing import Dict, Optional

import jsonschema


def read_json(filepath: str, schema: Optional[str] = None) -> Dict:
    """
    Reads a json file and returns the data inside it.
    If `schema` is given, it verifies

    Parameters
    -----------

    :param filepath: json filepath to be read from
    :type filepath: str
    :param schema: schema filepath to be checked
    :type schema: str, optional
    :return: The dictionary with all infos read from json
    :rtype: dict

    """
    with open(filepath, "r") as file:
        data = json.load(file)
    if not schema:
        with schema.open() as file:
            schema = json.load(file)
        jsonschema.validate(data, schema)
    return data
