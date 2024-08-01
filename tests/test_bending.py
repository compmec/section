"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest
from shapepy import Primitive

from compmec.section import HomogeneousSection
from compmec.section.material import Isotropic
from compmec.section.shape import shape2geometry


@pytest.mark.order(4)
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
        "tests/test_geomprop.py::test_end",
        "tests/test_axial.py::test_end",
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
        shape = Primitive.square(side)
        geometry = shape2geometry(shape)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection(geometry, material)
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
        values = field.eval(points)
        S33, S13, S23 = values[:3]
        E11, E33, E13, E23 = values[3:]
        assert np.all(np.abs(S33) < 1e-9)
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)
        assert np.all(np.abs(E11) < 1e-9)
        assert np.all(np.abs(E33) < 1e-9)
        assert np.all(np.abs(E13) < 1e-9)
        assert np.all(np.abs(E23) < 1e-9)

        field.forces = (0, 0, 0)
        field.momentums = (4 / 3, 0, 0)
        values = field.eval(points)
        S33, S13, S23 = values[:3]
        E11, E33, E13, E23 = values[3:]
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)
        good_normal_stress = [0, 0, 1, 1, 1, 0, -1, -1, -1]
        abs_diff = np.abs(S33 - good_normal_stress)
        assert np.all(abs_diff < 1e-9)

        field.forces = (0, 0, 0)
        field.momentums = (0, 4 / 3, 0)
        values = field.eval(points)
        S33, S13, S23 = values[:3]
        E11, E33, E13, E23 = values[3:]
        assert np.all(np.abs(S13) < 1e-9)  # No shear
        assert np.all(np.abs(S23) < 1e-9)  # No shear
        good_normal_stress = [0, -1, -1, 0, 1, 1, 1, 0, -1]
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
        shape = Primitive.square(ext_side)
        shape -= Primitive.square(int_side)
        geometry = shape2geometry(shape)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection(geometry, material)
        field = section.charged_field()

        points = [(0, 0)]  # origin
        points += [
            (0.5, 0),
            (0.5, 0.5),
            (0, 0.5),
            (-0.5, 0.5),
            (-0.5, 0),
            (-0.5, -0.5),
            (0, -0.5),
            (0.5, -0.5),
        ]
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
        S33, S13, S23 = values[:3]
        E11, E33, E13, E23 = values[3:]
        assert np.all(np.abs(S33) < 1e-9)
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)
        assert np.all(np.abs(E11) < 1e-9)
        assert np.all(np.abs(E33) < 1e-9)
        assert np.all(np.abs(E13) < 1e-9)
        assert np.all(np.abs(E23) < 1e-9)

        field.forces = (0, 0, 0)
        field.momentums = (1.25, 0, 0)
        values = field.eval(points)
        S33, S13, S23 = values[:3]
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)
        good_normal_stress = [0]  # origin
        good_normal_stress += [0, 0.5, 0.5, 0.5, 0, -0.5, -0.5, -0.5]
        good_normal_stress += [0, 1, 1, 1, 0, -1, -1, -1]
        abs_diff = np.abs(S33 - good_normal_stress)
        assert np.all(abs_diff < 1e-9)

        field.forces = (0, 0, 0)
        field.momentums = (0, 1.25, 0)
        values = field.eval(points)
        S33, S13, S23 = values[:3]
        assert np.all(np.abs(S13) < 1e-9)
        assert np.all(np.abs(S23) < 1e-9)
        good_normal_stress = [0]
        good_normal_stress += [-0.5, -0.5, 0, 0.5, 0.5, 0.5, 0, -0.5]
        good_normal_stress += [-1, -1, 0, 1, 1, 1, 0, -1]
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
