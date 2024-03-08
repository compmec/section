"""
File to tests the basis functions
"""

import numpy as np
import pytest

from compmec.section.basisfunc import BasisFunc


@pytest.mark.order(1)
@pytest.mark.dependency(
    depends=[],
    scope="session",
)
def test_begin():
    pass


class TestCyclic:
    @pytest.mark.order(1)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(1)
    @pytest.mark.timeout(10)
    @pytest.mark.dependency(depends=["TestCyclic::test_begin"])
    def test_linear(self):
        knots = (0, 0.5, 1)
        basis = BasisFunc.cyclic(knots, degree=1)
        umesh = (0, 0.25, 0.5, 0.75, 1)
        valus = basis.eval(umesh)
        assert valus.shape == (basis.ndofs, len(umesh))
        assert np.all(valus[0] == (1, 0.5, 0, 0.5, 1))
        assert np.all(valus[1] == (0, 0.5, 1, 0.5, 0))

    @pytest.mark.order(1)
    @pytest.mark.dependency(
        depends=[
            "TestCyclic::test_linear",
        ]
    )
    def test_end(self):
        pass


@pytest.mark.order(1)
@pytest.mark.dependency(depends=["TestCyclic::test_end"])
def test_end():
    pass
