"""
This file tests the axial properties like strain
and stress when normal force or bending moments are applied
"""
import numpy as np
import pytest
from compmec.shape import Primitive
from matplotlib import pyplot as plt

from compmec.section.material import Isotropic
from compmec.section.section import Section


@pytest.mark.order(4)
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
    ],
    scope="session",
)
@pytest.mark.dependency()
def test_begin():
    pass


class TestSinglePolygon:
    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin"])
    def test_centered_square(self):
        side = 2
        geometry = Primitive.square(side)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        strain = section.strain()
        stress = section.stress()
        points = [(0, 0)]
        sxz, syz, szz = stress.eval_interior(points)

    @pytest.mark.order(4)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_centered_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(4)
@pytest.mark.dependency(depends=["TestSinglePolygon::test_end"])
def test_end():
    pass
