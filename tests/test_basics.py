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
        section = Section([geometry], [material])
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
        section = Section([geometry], [material])
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
        section = Section([geometry], [material])
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
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin", "TestSinglePolygon::test_centered_rectangle"])
    def test_centered_regular_polygon(self):
        radius = 1
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        for nsides in range(3, 11):
            geometry = Primitive.regular_polygon(nsides, radius=radius)
            section = Section([geometry], [material])
            area = section.area()
            Qx, Qy = section.first_moment()
            Ixx, Ixy, Iyy = section.second_moment()
            good = radius**2 * nsides * np.sin(2*np.pi/nsides)/2
            assert abs(area - good) < 1e-9
            assert abs(Qx) < 1e-9
            assert abs(Qy) < 1e-9
            inner = 4 * np.sin(2*np.pi/nsides) + np.sin(4*np.pi/nsides)
            good = radius**4 * nsides * inner /48
            assert abs(Ixx - good) < 1e-9
            assert abs(Ixy) < 1e-9
            assert abs(Iyy - good) < 1e-9

    
    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestSinglePolygon::test_begin",
                                     "TestSinglePolygon::test_shifted_rectangle",
                                     "TestSinglePolygon::test_centered_regular_polygon"])
    def test_shifted_regular_polygon(self):
        radius = 3
        center = (5, -7)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        for nsides in range(3, 11):
            geometry = Primitive.regular_polygon(nsides, radius=radius, center=center)
            section = Section([geometry], [material])
            area = section.area()
            Qx, Qy = section.first_moment()
            Ixx, Ixy, Iyy = section.second_moment()
            good = radius**2 * nsides * np.sin(2*np.pi/nsides)/2
            assert abs(area - good) < 1e-9
            assert abs(Qx - area * center[1]) < 1e-9
            assert abs(Qy - area * center[0]) < 1e-9
            inner = 4 * np.sin(2*np.pi/nsides) + np.sin(4*np.pi/nsides)
            good = radius**4 * nsides * inner /48
            assert abs(Ixx - (good + area * center[1]**2)) < 1e-9
            assert abs(Ixy - area * center[0] * center[1]) < 1e-9
            assert abs(Iyy - (good + area * center[0]**2)) < 1e-9



    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=[
            "TestSinglePolygon::test_centered_square",
            "TestSinglePolygon::test_centered_rectangle",
            "TestSinglePolygon::test_shifted_square",
            "TestSinglePolygon::test_shifted_rectangle",
            "TestSinglePolygon::test_centered_regular_polygon",
            "TestSinglePolygon::test_shifted_regular_polygon",
        ]
    )
    def test_end(self):
        pass


class TestCircle:
    """
    Shape package doesn't support exactly circles.
    Therefore, we use a tolerance

    """

    @pytest.mark.order(3)
    @pytest.mark.dependency(
        depends=[
            "test_begin",
            "TestBuild::test_end",
            "TestSinglePolygon::test_end",
        ]
    )
    def test_begin(self):
        pass

    @pytest.mark.order(3)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestCircle::test_begin"])
    def test_centered_square(self):
        side = 3
        geometry = Primitive.square(side)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
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
    @pytest.mark.dependency(depends=["TestCircle::test_centered_square"])
    def test_centered_rectangle(self):
        width, height = 3, 5
        geometry = Primitive.square().scale(width, height)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
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
    @pytest.mark.dependency(depends=["TestCircle::test_centered_square"])
    def test_shifted_square(self):
        side = 3
        center = (5, -7)
        geometry = Primitive.square(side, center=center)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        section = Section([geometry], [material])
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
            "TestCircle::test_shifted_square",
            "TestCircle::test_centered_rectangle",
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
        section = Section([geometry], [material])
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
            "TestCircle::test_centered_square",
            "TestCircle::test_centered_rectangle",
            "TestCircle::test_shifted_square",
            "TestCircle::test_shifted_rectangle",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(3)
@pytest.mark.dependency(
    depends=[
        "TestSinglePolygon::test_end",
        "TestCircle::test_end",
    ]
)
def test_end():
    pass
