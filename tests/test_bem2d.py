"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest
from compmec.shape import Primitive

from compmec.section import Section
from compmec.section.bem2d import BEMModel
from compmec.section.material import Isotropic


@pytest.mark.order(8)
@pytest.mark.dependency(
    depends=[
        "tests/test_integral.py::test_end",
        "tests/test_material.py::test_end",
        "tests/test_curve.py::test_end",
        "tests/test_geometry.py::test_end",
        "tests/test_geomprop.py::test_end",
        "tests/test_axial.py::test_end",
        "tests/test_bending.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


class TestConstruction:
    @pytest.mark.order(8)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(8)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestConstruction::test_begin"])
    def test_centered_square(self):
        geometry = Primitive.square(2)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = Section.from_shapes(geometry, material)
        bemmodel = BEMModel(section)

        mesh = np.linspace(0, 1, 5)
        bemmodel[1] = mesh
        mesh = bemmodel[1]
        assert np.all(mesh == (0, 0.25, 0.5, 0.75, 1))

    @pytest.mark.order(8)
    @pytest.mark.dependency(
        depends=[
            "TestConstruction::test_centered_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(8)
@pytest.mark.dependency(depends=["TestConstruction::test_end"])
def test_end():
    pass
