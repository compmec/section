"""
This file contains functions and classes responsible to
deal with the boundary curves.

"""

from typing import Dict

import dataio
from compmec.shape import Point2D

from compmec import nurbs


def curves_from_json(json_filepath: str) -> Dict:
    """
    Reads a json file, returning the curves instances

    Parameters
    ----------
    json_filepath: str
        The json filename to be read
    return: dict[int, nurbs.Curve]
    """
    dictionary = dataio.read_json(json_filepath)
    return curves_from_dict(dictionary)


def curves_from_dict(dictionary: Dict) -> Dict[str, nurbs.Curve]:
    """
    Transforms the informations in the dictionary in
    nurbs.Curve instance

    """
    nodes = {}
    for label, xval, yval in dictionary["nodes"]:
        point = Point2D(xval, yval)
        nodes[label] = point
    curves = {}
    for curves_info in dictionary["curves"]:
        degree = curves_info["degree"]
        knotvector = curves_info["knotvector"]
        knotvector = nurbs.KnotVector(knotvector, degree=degree)
        ctrlpts_labels = curves_info["ctrlpoints"]
        curve = nurbs.Curve(knotvector)
        curve.ctrlpoints = [nodes[label] for label in ctrlpts_labels]
        if "weights" in curves_info:
            curve.weights = curves_info["weights"]
        label = curves_info["label"]
        curves[label] = curve
    return curves


def curves_to_dict(curves: Dict[str, nurbs.Curve]) -> Dict:
    nodes = {}
    curves_out = {}
    for name, curve in curves.items():
        node_labels = []
        for point in curve.ctrlpoints:
            idpt = id(point)
            node_labels.append(idpt)
            if idpt not in nodes:
                nodes[id(point)] = tuple(map(float, point))
        curves_out[name] = tuple(node_labels)
    nodes_out = []
    for key in sorted(nodes.keys()):
        xval, yval = nodes[key]
        nodes_out.append((key, xval, yval))
    dictionary = {}
    dictionary["nodes"] = nodes_out
    dictionary["curves"] = curves_out
    return dictionary


def curves_to_json(curves: Dict[str, nurbs.Curve], json_filepath: str):
    """
    Saves the given curves in the json file

    Parameters
    ----------
    curves: dict[str, nurbs.Curve]
        The curves to be saved, nurbs.Curve instance
    json_filepath: str
        The path to save

    """
    dictionary = curves_to_dict(curves)
    dataio.save_json(dictionary, json_filepath)
