"""
This file tests the basic geometry properties, such as
* Area
* First moment of inertia
* Second moment of inertia
"""

import math

import pytest
from shapepy import Primitive

from compmec.section.dataio import JsonIO
from compmec.section.material import Isotropic
from compmec.section.section import HomogeneousSection


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
        "tests/test_curve.py::test_end",
        "tests/test_geometry.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


class TestSinglePolygon:
    @pytest.mark.order(3)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin"])
    def test_centered_square(self):
        side = 3
        geometry = Primitive.square(side)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection.from_shape(geometry, material)
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
    @pytest.mark.dependency(
        depends=["TestSinglePolygon::test_centered_square"]
    )
    def test_centered_rectangle(self):
        width, height = 3, 5
        geometry = Primitive.square().scale(width, height)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection.from_shape(geometry, material)
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
    @pytest.mark.dependency(
        depends=["TestSinglePolygon::test_centered_square"]
    )
    def test_shifted_square(self):
        side = 3
        center = (5, -7)
        geometry = Primitive.square(side, center=center)
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection.from_shape(geometry, material)
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
        material = Isotropic(young_modulus=210e3, poissons_ratio=0.3)
        section = HomogeneousSection.from_shape(geometry, material)
        area = section.area()
        Qx, Qy = section.first_moment()
        Ixx, Ixy, Iyy = section.second_moment()
        assert area == width * height
        assert abs(Qx - area * center[1]) < 1e-9
        assert abs(Qy - area * center[0]) < 1e-9
        good_Ixx = width * height**3 / 12 + area * center[1] * center[1]
        good_Iyy = height * width**3 / 12 + area * center[0] * center[0]
        assert abs(Ixx - good_Ixx) < 1e-6
        assert Ixy == area * center[0] * center[1]
        assert abs(Iyy - good_Iyy) < 1e-6

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


class TestToFromJson:
    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=["test_begin", "TestSinglePolygon::test_end"]
    )
    def test_begin(self):
        pass

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestToFromJson::test_begin"])
    def test_read_square(self):
        json_filepath = "tests/json/steel_square.json"
        with JsonIO(json_filepath) as reader:
            reader.load_nodes()
            reader.load_curves()
            reader.load_materials()
            reader.load_sections()
        square = HomogeneousSection.instances["square"]
        area = square.area()
        Qx, Qy = square.first_moment()
        Ixx, Ixy, Iyy = square.second_moment()
        assert area == 4
        assert Qx == 0
        assert Qy == 0
        assert abs(Ixx - 4 / 3) < 1e-9
        assert Ixy == 0
        assert abs(Iyy - 4 / 3) < 1e-9

    @pytest.mark.order(3)
    @pytest.mark.skip(reason="Fix derivate of rational splines")
    @pytest.mark.timeout(20)
    @pytest.mark.dependency(depends=["TestToFromJson::test_begin"])
    def test_read_circle(self):
        json_filepath = "tests/json/steel_circle.json"
        with JsonIO(json_filepath) as reader:
            reader.load_nodes()
            reader.load_curves()
            reader.load_materials()
            reader.load_sections()
        square = HomogeneousSection.instances["square"]
        area = square.area()
        Qx, Qy = square.first_moment()
        Ixx, Ixy, Iyy = square.second_moment()
        assert abs(area - math.pi) < 2e-6
        assert abs(Qx) < 1e-6
        assert abs(Qy) < 1e-6
        assert abs(Ixx - math.pi / 4) < 1e-6
        assert abs(Ixy) < 1e-6
        assert abs(Iyy - math.pi / 4) < 1e-6

    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=[
            "TestToFromJson::test_begin",
            "TestToFromJson::test_read_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=["TestSinglePolygon::test_end", "TestToFromJson::test_end"]
)
def test_end():
    pass
