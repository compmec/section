"""
File to tests cases when only bending moments are applied
"""

import numpy as np
import pytest
from shapepy import Primitive

from compmec.section.basisfunc import SplineBasisFunction
from compmec.section.bem2d import BEMModel, ScalarFunction
from compmec.section.material import Isotropic
from compmec.section.section import HomogeneousSection


@pytest.mark.order(10)
@pytest.mark.skip()
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

        with BEMModel(steel_square) as model:
            curve = steel_square.geometry.curves[0]
            knots = (0, 1, 2, 3, 4)
            basis = SplineBasisFunction.cyclic(knots, degree=1)
            model.sources = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            model.solve()

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
        with BEMModel(steel_square) as model:
            model.generate_mesh()
            model.solve()

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


class TestTorsionConstant:

    @pytest.mark.order(10)
    @pytest.mark.dependency(depends=["TestAutoSolve::test_end"])
    def test_begin(self):
        pass

    @pytest.mark.order(10)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestTorsionConstant::test_begin"])
    def test_full_square(self):

        def torsion_constant_rectangle(a: float, b: float, nmax: int):
            """
            Gives the value of torsion contant for a rectangle profile
            of lenght 'a' and height 'b'
            """
            r = b / a
            temp = sum(
                np.tanh(np.pi * (2 * n + 1) / 2 * r) / (2 * n + 1) ** 5
                for n in range(nmax + 1)
            )
            tors_const = r / 3 - 64 * temp / np.pi**5
            return tors_const * a**4

        full_square = Primitive.square(side=2)
        steel = Isotropic(young_modulus=210, poissons_ratio=0.3)
        steel_square = HomogeneousSection.from_shape(full_square, steel)

        ndiv = 5
        knots = tuple(i / ndiv for i in range(4 * ndiv + 1))
        basis = SplineBasisFunction.cyclic(knots, degree=1)
        ctrlpoints = [
            4.524489888226429e-05,
            -1.761717469582290e-01,
            -6.775319276312486e-02,
            6.775319276312486e-02,
            1.761717469582290e-01,
            -4.524489888225453e-05,
            -1.761463600718403e-01,
            -6.776131453206423e-02,
            6.776131453206423e-02,
            1.761463600718403e-01,
            4.524489888224952e-05,
            -1.761717469582290e-01,
            -6.775319276312486e-02,
            6.775319276312487e-02,
            1.761717469582290e-01,
            -4.524489888225529e-05,
            -1.761463600718403e-01,
            -6.776131453206426e-02,
            6.776131453206424e-02,
            1.761463600718403e-01,
        ]
        warping_curve = ScalarFunction(basis, ctrlpoints)

        good_torconst = torsion_constant_rectangle(2, 2, 10)
        test_torconst = steel_square.torsion_constant()
        assert abs(test_torconst - good_torconst) < 1e-3

    @pytest.mark.order(10)
    @pytest.mark.dependency(
        depends=[
            "TestTorsionConstant::test_begin",
            "TestTorsionConstant::test_full_square",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(10)
@pytest.mark.dependency(
    depends=[
        "TestSolvingSystem::test_end",
        "TestAutoSolve::test_end",
        "TestTorsionConstant::test_end",
    ]
)
def test_end():
    pass
