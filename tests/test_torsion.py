"""
File to tests cases when only bending moments are applied
"""

import pytest
from shapepy import Primitive

from compmec.section.basisfunc import SplineBasisFunction
from compmec.section.bem2d import BEMModel
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
        full_square = Primitive.square(side=2, center=(0, 0))
        steel = Isotropic(young_modulus=210, poissons_ratio=0.3)
        steel_square = HomogeneousSection.from_shape(full_square, steel)
        curve = steel_square.geometry.curves[0]

        model = BEMModel(steel_square)
        knots = (0, 1, 2, 3, 4)
        basis = SplineBasisFunction.cyclic(knots, degree=1)
        model.add_basis(curve, basis)

        sources = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        model.solve(sources)

    @pytest.mark.order(10)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestSolvingSystem::test_begin"])
    def test_hollow_square(self):
        hollow_square = Primitive.square(6) - Primitive.square(2)
        steel = Isotropic(young_modulus=1, poissons_ratio=0.3)
        section = HomogeneousSection.from_shape(hollow_square, steel)

        model = BEMModel(section)
        curve_int = section.geometry.curves[0]
        curve_ext = section.geometry.curves[1]
        basis_int = SplineBasisFunction.cyclic(curve_int.knots)
        basis_ext = SplineBasisFunction.cyclic(curve_ext.knots)
        model.add_basis(curve_int, basis_int)
        model.add_basis(curve_ext, basis_ext)

        sources = [(-3, -3), (3, -3), (3, 3), (-3, 3)]
        sources += [(-1, -1), (1, -1), (1, 1), (-1, 1)]
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


class TestAutoSolve:
    @pytest.mark.order(10)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(10)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestAutoSolve::test_begin"])
    def test_full_square(self):
        full_square = Primitive.square(side=2, center=(0, 0))
        steel = Isotropic(young_modulus=210, poissons_ratio=0.3)
        steel_square = HomogeneousSection.from_shape(full_square, steel)
        steel_square.solve()

        sources = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)]
        sources += [(2, 0), (0, 2), (-2, 0), (0, -2)]
        for source in sources:
            steel_square.warping.eval(source)

    @pytest.mark.order(10)
    @pytest.mark.dependency(
        depends=[
            "TestAutoSolve::test_begin",
            "TestAutoSolve::test_full_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(10)
@pytest.mark.dependency(depends=["TestSolvingSystem::test_end"])
def test_end():
    pass
