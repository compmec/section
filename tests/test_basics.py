"""
This file tests the basic geometry properties, such as
* Area
* First moment of inertia
* Second moment of inertia
"""
import pytest
from compmec.shape import Primitive

from compmec.section.material import Isotropic
from compmec.section.section import SimpleSection


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
    def test_simple_square(self):
        geometry = Primitive.square(3)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        SimpleSection(geometry, material)

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestBuild::test_begin"])
    def test_hollow_square(self):
        geometry = Primitive.square(3) - Primitive.square(1)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        SimpleSection(geometry, material)

    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=["TestBuild::test_simple_square", "TestBuild::test_hollow_square"]
    )
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
        section = SimpleSection(geometry, material)
        area = section.area()
        Qx, Qy = section.first_moment()
        Ixx, Ixy, Iyy = section.second_moment()
        assert area == side**2
        assert Qx == 0
        assert Qy == 0
        assert Ixx == side**4 / 12
        assert Ixy == 0
        assert Iyy == side**4 / 12

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_centered_square"])
    def test_centered_rectangle(self):
        width, height = 3, 5
        geometry = Primitive.square().scale(width, height)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = SimpleSection(geometry, material)
        area = section.area()
        Qx, Qy = section.first_moment()
        Ixx, Ixy, Iyy = section.second_moment()
        assert area == width * height
        assert Qx == 0
        assert Qy == 0
        assert Ixx == width * height**3 / 12
        assert Ixy == 0
        assert Iyy == height * width**3 / 12

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
        section = SimpleSection(geometry, material)
        area = section.area()
        Qx, Qy = section.first_moment()
        Ixx, Ixy, Iyy = section.second_moment()
        assert area == side**2
        assert Qx == area * center[1]
        assert Qy == area * center[0]
        assert Ixx == side**4 / 12 + area * center[1] ** 2
        assert Ixy == area * center[0] * center[1]
        assert Iyy == side**4 / 12 + area * center[0] ** 2

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_shifted_square",
            "TestSinglePolygon::test_centered_rectangle",
        ]
    )
    def test_shifted_rectangle(self):
        width, height = 3, 5
        center = (5, -7)
        geometry = Primitive.square()
        geometry.scale(width, height)
        geometry.move(center)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = SimpleSection(geometry, material)
        area = section.area()
        Qx, Qy = section.first_moment()
        Ixx, Ixy, Iyy = section.second_moment()
        assert area == width * height
        assert Qx == area * center[1]
        assert Qy == area * center[0]
        assert Ixx == width * height**3 / 12 + area * center[1] * center[1]
        assert Ixy == area * center[0] * center[1]
        assert Iyy == height * width**3 / 12 + area * center[0] * center[0]

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
