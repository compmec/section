import pytest

from compmec.section import material


@pytest.mark.order(2)
@pytest.mark.dependency(
    depends=["tests/test_dataio.py::test_end"],
    scope="session",
)
def test_begin():
    pass


class TestIsotropic:
    @pytest.mark.order(2)
    @pytest.mark.dependency(depends=["test_begin"])
    def test_begin(self):
        pass

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestIsotropic::test_begin"])
    def test_main(self):
        Es = [100, 200, 3, 45, 900, 20, 66]
        nus = [0.3, 0.1, 0.05, 0.001, 0.489, 0.45]
        for E, nu in zip(Es, nus):
            G = E / (2 * (1 + nu))
            K = E * G / (3 * (3 * G - E))
            L = K - 2 * G / 3

            mat = material.Isotropic()
            mat.young_modulus = E
            mat.poissons_ratio = nu
            assert abs(E - mat.young_modulus) < 1e-9
            assert abs(nu - mat.poissons_ratio) < 1e-9
            assert abs(G - mat.shear_modulus) < 1e-9
            assert abs(K - mat.bulk_modulus) < 1e-9
            assert abs(L - mat.lame_parameter_1) < 1e-9
            assert abs(G - mat.lame_parameter_2) < 1e-9

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestIsotropic::test_main"])
    def test_fail_setting_young(self):
        mat = material.Isotropic()
        with pytest.raises(ValueError):
            mat.young_modulus = 0
        with pytest.raises(ValueError):
            mat.young_modulus = -100

    @pytest.mark.order(2)
    @pytest.mark.timeout(1)
    @pytest.mark.dependency(depends=["TestIsotropic::test_main"])
    def test_fail_setting_poisson(self):
        mat = material.Isotropic()
        with pytest.raises(ValueError):
            mat.poissons_ratio = 0.495
        with pytest.raises(ValueError):
            mat.poissons_ratio = 0.501

    @pytest.mark.order(2)
    @pytest.mark.dependency(depends=["TestIsotropic::test_main"])
    def test_end(self):
        pass


@pytest.mark.order(2)
@pytest.mark.dependency(depends=["TestIsotropic::test_end"])
def test_end():
    pass
