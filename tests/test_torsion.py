"""
This file tests the basic geometry properties, such as
* Area
* First moment of inertia
* Second moment of inertia
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
        "tests/test_bem2d.py::test_basics",
        "tests/test_bem2d.py::test_end",
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
    @pytest.mark.skip(reason="Needs correction")
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin"])
    def test_centered_square(self):
        side = 2
        geometry = Primitive.square(side)

        npts = 11
        usample = np.linspace(0, 1, npts, endpoint=False)
        xvals = (
            list(-1 + 2 * usample) + [1] * npts + list(1 - 2 * usample) + [-1] * npts
        )
        yvals = (
            [-1] * npts + list(-1 + 2 * usample) + [1] * npts + list(1 - 2 * usample)
        )
        vertices = np.transpose([xvals, yvals])
        geometry = Primitive.polygon(vertices)

        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        torsion = section.torsion_constant()
        diff = torsion - (9 / 64) * side**4
        assert abs(diff) < 1e-3

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
