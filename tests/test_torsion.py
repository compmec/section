"""
File to tests cases when only bending moments are applied
"""

import pytest

from compmec.section.basisfunc import BasisFunc
from compmec.section.bem2d import BEMModel
from compmec.section.curve import Curve
from compmec.section.geometry import ConnectedGeometry
from compmec.section.material import Isotropic
from compmec.section.section import HomogeneousSection


@pytest.mark.order(10)
@pytest.mark.dependency(
    depends=[
        "tests/test_bem2d.py::test_end",
        "tests/test_axial.py::test_end",
        "tests/test_bending.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


class TestSolvingSystem:
    @pytest.mark.order(10)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(10)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSolvingSystem::test_begin"])
    def test_full_square(self):
        vertices = [[1, 1], [-1, 1], [-1, -1], [1, -1]]
        curve = Curve.from_vertices(vertices)
        geometry = ConnectedGeometry([curve])
        material = Isotropic(young_modulus=1, poissons_ratio=0.3)
        square = HomogeneousSection(geometry, material)

        basis_function = BasisFunc.cyclic(curve.knots)
        model = BEMModel(square)
        model.basis[curve.label] = basis_function

        sources = vertices
        model.solve(sources)

    @pytest.mark.order(10)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSolvingSystem::test_begin"])
    def test_hollow_square(self):
        vertices_int = [[-1, -1], [-1, 1], [1, 1], [1, -1]]
        curve_int = Curve.from_vertices(vertices_int)
        vertices_ext = [[3, 3], [-3, 3], [-3, -3], [3, -3]]
        curve_ext = Curve.from_vertices(vertices_ext)

        geometry = ConnectedGeometry([curve_ext, curve_int])
        material = Isotropic(young_modulus=1, poissons_ratio=0.3)
        square = HomogeneousSection(geometry, material)
        model = BEMModel(square)

        basis_int = BasisFunc.cyclic(curve_int.knots)
        basis_ext = BasisFunc.cyclic(curve_ext.knots)
        model.basis[curve_int.label] = basis_int
        model.basis[curve_ext.label] = basis_ext

        sources = vertices_int + vertices_ext
        model.solve(sources)

    @pytest.mark.order(10)
    @pytest.mark.dependency(
        depends=[
            "TestSolvingSystem::test_begin",
            "TestSolvingSystem::test_full_square",
            "TestSolvingSystem::test_hollow_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(10)
@pytest.mark.dependency(depends=["TestSolvingSystem::test_end"])
def test_end():
    pass
