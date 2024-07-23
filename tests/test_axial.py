"""
File to tests cases when only axial forces are applied
"""

import numpy as np
import pytest
from shapepy import Primitive

from compmec.section import HomogeneousSection
from compmec.section.material import Isotropic


@pytest.mark.order(4)
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
        section = HomogeneousSection.from_shape(geometry, material)
        field = section.charged_field()

        points = [(0, 0)]  # origin
        points += [
            (1, 0),
            (1, 1),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (-1, -1),
            (0, -1),
            (1, -1),
        ]
        field.forces = (0, 0, 0)
        field.momentums = (0, 0, 0)
        values = field.eval(points)
        assert np.all(np.abs(values) < 1e-9)  # no charge, all zero

        field.forces = (0, 0, 4)
        field.momentums = (0, 0, 0)
        S33, S13, S23 = field.eval(points)[:3]
        assert np.all(np.abs(S13) < 1e-9)  # no shear
        assert np.all(np.abs(S23) < 1e-9)  # no shear
        assert np.all(np.abs(S33 - 1) < 1e-9)

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
        section = HomogeneousSection.from_shape(geometry, material)
        field = section.charged_field()

        points = [(0, 0)]
        points += [
            (1, 0),
            (1, 1),
            (0, 1),
            (-1, 1),
            (-1, 0),
            (-1, -1),
            (0, -1),
            (1, -1),
        ]
        field.forces = (0, 0, 0)
        field.momentums = (0, 0, 0)
        S33, S13, S23 = field.eval(points)[:3]
        assert np.all(np.abs(S33) < 1e-9)
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)

        field.forces = (0, 0, 3)
        field.momentums = (0, 0, 0)
        S33, S13, S23 = field.eval(points)[:3]
        assert np.all(np.abs(S13) < 1e-9)  # no shear
        assert np.all(np.abs(S23) < 1e-9)  # no shear
        good_normal_stress = (0, 1, 1, 1, 1, 1, 1, 1, 1)
        abs_diff = np.abs(S33 - good_normal_stress)
        assert np.all(abs_diff < 1e-9)

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
