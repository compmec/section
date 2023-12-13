import pytest

from compmec.section import dataio


@pytest.mark.order(1)
@pytest.mark.dependency()
def test_begin():
    pass


@pytest.mark.order(1)
@pytest.mark.dependency(depends=["test_begin"])
def test_read_section_json():
    filename = "tests/json/steel_square.json"
    dataio.read_section_json(filename)


@pytest.mark.order(1)
@pytest.mark.dependency(depends=["test_read_section_json"])
def test_end():
    pass
