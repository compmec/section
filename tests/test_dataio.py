import pytest

from compmec.section.material import Isotropic 
from compmec.section.section import Section
from compmec.shape import Primitive

@pytest.mark.order(2)
@pytest.mark.dependency()
def test_begin():
    pass


class TestBuild:
    @pytest.mark.order(2)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestBuild::test_begin"])
    def test_simple_square(self):
        geometry = Primitive.square(3)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        Section(geometry, material)

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestBuild::test_begin"])
    def test_hollow_square(self):
        geometry = Primitive.square(3) - Primitive.square(1)
        material = Isotropic()
        material.young_modulus = 210e3
        material.poissons_ratio = 0.30
        Section(geometry, material)

    @pytest.mark.order(2)
    @pytest.mark.dependency(
        depends=["TestBuild::test_simple_square", "TestBuild::test_hollow_square"]
    )
    def test_end(self):
        pass


@pytest.mark.order(2)
@pytest.mark.dependency(depends=["TestBuild::test_end"])
def test_save_section_json():
    geometry = Primitive.square(3)
    material = Isotropic()
    material.young_modulus = 210e3
    material.poissons_ratio = 0.30
    section = Section(geometry, material)

    filename = "tests/json/steel_square.json"
    section.to_json(filename)
    


@pytest.mark.order(2)
@pytest.mark.dependency(depends=["test_save_section_json"])
def test_read_section_json():
    filename = "tests/json/steel_square.json"
    Section.from_json(filename)


@pytest.mark.order(2)
@pytest.mark.dependency(depends=["test_read_section_json"])
def test_end():
    pass
