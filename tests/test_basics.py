"""
This file tests the basic geometry properties, such as
* Area
* First moment of inertia
* Second moment of inertia
"""
import numpy as np
import pytest
from compmec.shape import Primitive

from compmec.section.material import Isotropic
from compmec.section.section import Section


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=["tests/test_material.py::test_end"],
    scope="session",
)
@pytest.mark.dependency()
def test_begin():
    pass


class TestBuild:
    @pytest.mark.order(3)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestBuild::test_begin"])
    def test_simple(self):
        geom = Primitive.square(3)
        mat = Isotropic()
        mat.young_modulus = 210e3
        mat.poissons_ratio = 0.30
        Section([geom], [mat])

    @pytest.mark.order(3)
    @pytest.mark.dependency(depends=["TestBuild::test_simple"])
    def test_end(self):
        pass


class TestSinglePolygon:
    @pytest.mark.order(3)
    @pytest.mark.dependency(depends=["test_begin", "TestBuild::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin"])
    def test_centered_square(self):
        side = 3
        geometry = Primitive.square(side)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        assert section.area() == side**2
        first = section.first_momentum()
        assert tuple(first) == (0, 0)
        second = section.second_momentum()
        good = ((side**4 / 12, 0), (0, side**4 / 12))
        np.testing.assert_equal(second, good)

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_centered_square"])
    def test_centered_rectangle(self):
        width, height = 3, 5
        geometry = Primitive.square().scale(width, height)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        assert section.area() == width * height
        first = section.first_momentum()
        assert tuple(first) == (0, 0)
        second = section.second_momentum()
        good = ((height * width**3 / 12, 0), (0, width * height**3 / 12))
        np.testing.assert_equal(second, good)

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_centered_square"])
    def test_shifted_square(self):
        side = 3
        center = (5, -7)
        geometry = Primitive.square(side, center=center)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        assert section.area() == side**2
        first = section.first_momentum()
        assert tuple(first) == (side**2 * center[0], side**2 * center[1])
        second = section.second_momentum()
        good = side**2 * np.tensordot(center, center, axes=0)
        print(good)
        good = good + ((side**4 / 12, 0), (0, side**4 / 12))
        print(good)
        np.testing.assert_equal(second, good)

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_shifted_square",
            "TestSinglePolygon::test_centered_rectangle",
        ]
    )
    def test_shifted_rectangle(self):
        side = 3
        center = (5, -7)
        geometry = Primitive.square(side, center=center)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
        assert section.area() == side**2
        first = section.first_momentum()
        assert tuple(first) == (side**2 * center[0], side**2 * center[1])
        second = section.second_momentum()
        good = side**2 * np.tensordot(center, center, axes=0)
        good = good + ((side**4 / 12, 0), (0, side**4 / 12))
        np.testing.assert_equal(second, good)

    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_centered_square",
            "TestSinglePolygon::test_centered_rectangle",
            "TestSinglePolygon::test_shifted_square",
            "TestSinglePolygon::test_shifted_rectangle",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(3)
@pytest.mark.dependency(depends=["TestSinglePolygon::test_end"])
def test_end():
    pass
