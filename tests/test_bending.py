"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest
from compmec.shape import Primitive

from compmec.section import Section
from compmec.section.material import Isotropic


@pytest.mark.order(4)
@pytest.mark.skip(reason="Test only in the version 0.4.0")
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
        "tests/test_geomprop.py::test_end",
    ],
    scope="session",
)
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
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = Section(geometry, material)
        field = section.charged_field()
        points = [(0, 0)]
        values = field.eval(points)
        assert np.all(np.abs(values) < 1e-9)

    @pytest.mark.order(4)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_centered_square",
        ]
    )
    def test_end(self):
        pass


class TestHollowPolygon:
    @pytest.mark.order(4)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(4)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin"])
    def test_centered_square(self):
        int_side, ext_side = 1, 2
        geometry = Primitive.square(ext_side)
        geometry -= Primitive.square(int_side)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = Section(geometry, material)
        field = section.charged_field()
        points = [(0, 0)]
        values = field.eval(points)
        assert np.all(np.abs(values) < 1e-9)

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
