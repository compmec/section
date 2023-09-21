import pytest
from compmec.shape import Primitive

from compmec.section.material import Isotropic
from compmec.section.profile import Profile


@pytest.mark.order(2)
@pytest.mark.dependency(
    depends=[
        "tests/test_material.py::test_end",
    ],
    scope="session",
)
def test_begin():
    pass


class TestProfile:
    @pytest.mark.order(2)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestProfile::test_begin"])
    def test_construct(self):
        square = Primitive.square()
        steel = Isotropic()
        steel.young_modulus = 200e3
        steel.poissons_ratio = 0.3
        profile = Profile([square], [steel])

    @pytest.mark.order(2)
    @pytest.mark.dependency(
        depends=["TestProfile::test_begin", "TestProfile::test_construct"]
    )
    def test_end(self):
        pass


@pytest.mark.order(2)
@pytest.mark.dependency(depends=["TestProfile::test_end"])
def test_end():
    pass
