"""
Solves the poisson problem using boundary element method

nabla^2_u = f(x, y)


"""
import math
from typing import Tuple

import numpy as np


class Integration:
    """
    Nodes positions and weights to perform numerical integration
    On the interval [0, 1]
    """

    def log(npts: int) -> Tuple[Tuple[float]]:
        """Return nodes and weights to integrate:

        .. math::
            I = \\gauss_{0}^{1} f(x) \\cdot ln(x) dx
        """
        if not isinstance(npts, int) or npts <= 0:
            raise ValueError("npts must be a positive integer")
        if npts == 1:
            nodes = [0.36787_94411_71442_32159_55237_70161]
            weights = [1]
        elif npts == 2:
            nodes = [
                0.88296_86513_76530_11759_59513_85185e-1,
                0.67518_64909_09887_20103_62743_16962,
            ]
            weights = [
                0.29849_98937_05524_91470_84741_43289,
                0.70150_01062_94475_08529_15258_56711,
            ]
        elif npts == 3:
            nodes = [
                0.28811_66253_09518_31174_32844_63059e-1,
                0.30406_37296_12137_65261_08623_58639,
                0.81166_92253_44078_11686_37051_77761,
            ]
            weights = [
                0.10333_07079_64928_64676_92515_92664,
                0.45463_65259_70098_70884_06911_21425,
                0.44203_27660_64972_64439_00572_85911,
            ]
        elif npts == 4:
            nodes = [
                0.11802_59099_78449_18264_91730_11095e-1,
                0.14282_56799_77483_69513_68513_69176,
                0.48920_15226_54574_47871_90313_05699,
                0.87867_99740_69183_70280_76889_06265,
            ]
            weights = [
                0.43391_02877_84143_91101_89836_96321e-1,
                0.24045_20976_59460_67597_84500_61567,
                0.42140_34522_59775_93197_88150_24211,
                0.29475_34213_02349_00094_08365_44590,
            ]
        elif npts == 5:
            nodes = [
                0.56522_28205_08009_71359_27256_19673e-2,
                0.73430_37174_26522_73406_15889_38883e-1,
                0.28495_74044_62558_15371_45276_01926,
                0.61948_22640_84778_38140_68089_43051,
                0.91575_80830_04698_33378_46091_80928,
            ]
            weights = [
                0.21046_94579_18546_29119_00268_26421e-1,
                0.13070_55407_44446_69759_10762_54992,
                0.28970_23016_71314_15684_15903_51057,
                0.35022_03701_20398_71028_55468_04135,
                0.20832_48416_71985_80616_27839_07174,
            ]
        elif npts == 6:
            nodes = [
                0.30258_02137_54625_87097_29970_37526e-2,
                0.40978_25415_59506_15053_46596_31089e-1,
                0.17086_32955_26877_29472_51498_29786,
                0.41325_57088_44793_24766_64814_55164,
                0.70909_51467_90628_54395_00459_17084,
                0.93823_95903_77167_09135_50205_94716,
            ]
            weights = [
                0.11351_33881_72726_09440_49112_38284e-1,
                0.75241_06995_49165_22917_35628_91092e-1,
                0.18879_00416_15416_35460_95079_43772,
                0.28582_07218_27227_31198_66834_80085,
                0.28448_64278_91408_80004_51516_69844,
                0.15431_03998_93758_40100_08094_93362,
            ]
        elif npts == 7:
            nodes = [
                0.17596_52118_46577_42805_62642_84949e-2,
                0.24469_65071_25133_67427_64533_73497e-1,
                0.10674_80568_58788_95418_02597_81083,
                0.27580_76412_95917_38307_78595_12057,
                0.51785_51421_51833_71615_86689_61982,
                0.77181_54853_62384_90027_46468_69494,
                0.95284_13405_81090_55899_43065_88503,
            ]
            weights = [
                0.66326_66319_02570_51178_39049_89051e-2,
                0.45799_70797_84753_34125_57673_48120e-1,
                0.12384_02080_71318_19455_04895_64922,
                0.21210_19260_23811_93010_79148_75456,
                0.26139_06456_72007_72564_65806_06859,
                0.23163_61802_90909_38431_88155_26104,
                0.11859_86656_44451_72613_27836_41957,
            ]
        else:
            raise NotImplementedError
        nodes = np.array(nodes, dtype="float64")
        weights = np.array(weights, dtype="float64")
        return nodes, weights

    def gauss(npts: int) -> Tuple[Tuple[float]]:
        """
        Returns nodes and weights to gaussian integration
        """
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        nodes, weights = np.polynomial.legendre.leggauss(npts)
        return (1 + nodes) / 2, weights / 2

    def closed(npts: int) -> Tuple[Tuple[float]]:
        """Closed newton cotes formula"""
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        if npts == 1:
            nodes = [0.5]
            weights = [1]
        elif npts == 2:
            nodes = [0, 1]
            weights = [0.5, 0.5]
        elif npts == 3:
            nodes = [0, 0.5, 1]
            weights = [1 / 6, 2 / 3, 1 / 6]
        elif npts == 4:
            nodes = (0, 1 / 3, 2 / 3, 1 / 3)
            weights = (1 / 8, 3 / 8, 3 / 8, 1 / 8)
        elif npts == 5:
            nodes = (0, 1 / 4, 1 / 2, 3 / 4, 1)
            weights = (7 / 90, 16 / 45, 2 / 15, 16 / 45, 7 / 90)
        elif npts == 6:
            nodes = (0, 1 / 5, 2 / 5, 3 / 5, 4 / 5, 1)
            weights = (19 / 288, 25 / 96, 25 / 144, 25 / 144, 25 / 96, 19 / 288)
        elif npts == 7:
            nodes = (0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 5 / 6, 1)
            weights = (41 / 840, 9 / 35, 9 / 280, 45 / 105, 9 / 280, 9 / 35, 41 / 840)
        else:
            raise NotImplementedError
        nodes = np.array(nodes, dtype="float64")
        weights = np.array(weights, dtype="float64")
        return nodes, weights

    def chebyshev(npts: int) -> Tuple[Tuple[float]]:
        if not isinstance(npts, int) or npts < 1:
            raise ValueError(f"npts invalid: {npts}")
        nums = range(1, 2 * npts, 2)
        nums = tuple(np.float64(num) / (2 * npts) for num in nums)
        nodes = tuple(math.sin(0.5 * math.pi * num) ** 2 for num in nums)
        if npts == 1:
            weights = (1,)
        elif npts == 2:
            weights = (0.5, 0.5)
        elif npts == 3:
            weights = (2 / 9, 5 / 9, 2 / 9)
        elif npts == 4:
            root2 = np.sqrt(2)
            weights = (
                (3 - root2) / 12,
                (3 - root2) / 12,
                (3 - root2) / 12,
                (3 - root2) / 12,
            )
        elif npts == 5:
            root5 = np.sqrt(5)
            weights = (
                (13 - 3 * root5) / 75,
                (13 + 3 * root5) / 75,
                23 / 75,
                (13 + 3 * root5) / 75,
                (13 - 3 * root5) / 75,
            )
        elif npts == 6:
            root3 = np.sqrt(3)
            weights = (
                (14 - 5 * root3) / 90,
                17 / 90,
                (14 + 5 * root3) / 90,
                (14 + 5 * root3) / 90,
                17 / 90,
                (14 - 5 * root3) / 90,
            )
        else:
            raise NotImplementedError
        nodes = np.array(nodes, dtype="float64")
        weights = np.array(weights, dtype="float64")
        return nodes, weights


def IntegrateUVn(vertices: Tuple[Tuple[float]]) -> Tuple[Tuple[float]]:
    """
    Computes the integral

    2*pi * int_{Gamma} u * (dv/dn) ds

    with

    * u the objective function
    * v the green function with source at S_i
    S_i = closedcurve(source)

    It in fact returns the vector like

    B = [B_0, B_1, ..., B_{n-1}]

    (n) is the number of degree of freedom

    The integral is transformed to the equivalent

    B_{ij} = 2*pi * sum_k int_{a_k}^{b_k} phi_j (r x p')/(r^2) dt

    with r(t) = p(t) - S_i
    """
    vertices = np.array(vertices, dtype="float64")
    vectors = np.roll(vertices, -1, axis=0) - vertices
    nverts = len(vertices)
    result = np.zeros((nverts, nverts), dtype="float64")

    gauss_nodes, gauss_weights = Integration.gauss(8)
    for i, Si in enumerate(vertices):
        # Source point at Si
        for j0, dVj in enumerate(vectors):
            # Integral over the segment j
            # between points V[j] and V[j+1]
            j1 = (j0 + 1) % nverts
            Vj0 = vertices[j0] - Si
            radius = tuple(Vj0 + t * dVj for t in gauss_nodes)
            funcs = tuple(np.cross(rad, dVj) / np.inner(rad, rad) for rad in radius)
            phi1 = gauss_nodes
            phi2 = 1 - gauss_nodes
            if j0 != i:
                result[i, j0] += np.einsum("k,k,k", phi2, funcs, gauss_weights)
            if j1 != i:
                result[i, j1] += np.einsum("k,k,k", phi1, funcs, gauss_weights)
        vertices += Si
    return result


def IntegrateUnV(vertices: Tuple[Tuple[float]]) -> Tuple[float]:
    """
    Computes the integral

    2*pi * int_{Gamma} (du/dn) * v ds

    with

    * u the objective function
    * v the green function with source at S_i
    S_i = closedcurve(source)

    It in fact returns the vector like

    A = [A_0, A_1, ..., A_{n-1}]

    (n) is the number of degree of freedom

    The integral is transformed to the equivalent

    A_{ij} = 2*pi * sum_k int_{a_k}^{b_k} phi_j ln(r) * abs(p') dt

    with r(t) = p(t) - S_i
    """
    vertices = np.array(vertices, dtype="float64")
    nverts = len(vertices)
    vectors = np.roll(vertices, -1, axis=0) - vertices
    # normvectors = np.linalg.norm(vectors, axis=1)
    alphs = np.einsum("ij,ij->i", vertices, vectors)
    betas = np.einsum("ij,ij->i", vectors, vectors)
    result = np.zeros(nverts, dtype="float64")
    gauss_nodes, gauss_weights = Integration.gauss(5)
    for i, Si in enumerate(vertices):
        for j0, dVj in enumerate(vectors):
            alp, bet = alphs[j0], betas[j0]
            j1 = j0 % nverts
            Vj0 = vertices[j0] - Si
            if j0 == i:
                # Associated with basis function (t)
                result[i] -= alp + bet / 4
                absdV = np.linalg.norm(dVj)
                result[i] += (alp + bet / 2) * np.log(absdV)
            elif j1 == i:
                result[i] -= alp + 3 * bet / 4
                # Associated with basis function (1-t)
                absdV = np.linalg.norm(dVj)
                result[i] += (alp + bet / 2) * np.log(absdV)
            else:
                radius = np.array(
                    tuple(Vj0 + t * dVj for t in gauss_nodes), dtype="float64"
                )
                lnradius = np.log(np.linalg.norm(radius, axis=1))
                phis = alp + bet * gauss_nodes
                result[i] += np.einsum("i,i,i", phis, lnradius, gauss_weights)
    return result


def CornerAngles(vertices):
    """
    Returns the value of the angle of each vertice
    """
    vectors0 = np.roll(vertices, -1, axis=0) - vertices
    vectors1 = np.roll(vectors0, -1, axis=0)
    crosses = tuple(np.cross(v0, v1) for v0, v1 in zip(vectors0, vectors1))
    inners = tuple(np.inner(v0, v1) for v0, v1 in zip(vectors0, vectors1))
    angles = tuple(np.arctan2(cross, inner) for cross, inner in zip(crosses, inners))
    angles = np.array(angles, dtype="float64")
    return angles
