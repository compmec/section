.. _validation:


============
Introduction
============

The validation of the code happens in 2 steps:

1. **Numerical**: Checks if the output values correspond with given input, following the :ref:`Theory <theory>`
2. **Literature**: Compares the output values with the values found in textbooks and references

====================
Numerical Validation
====================

Circle
------

The first case is for a circular bar of radius :math:`R`, centered in the origin

Geometric properties
^^^^^^^^^^^^^^^^^^^^

Computing the integrals

.. math::
    I_{a, b} = \int_{\Omega} x^a \cdot y^b \ d\Omega

.. math::
    I = \begin{bmatrix}\pi R^2 & 0 & \frac{\pi R^4}{4} & 0 \\ 0 & 0 & 0 & 0 \\ \frac{\pi R^4}{4} & 0 & \frac{\pi R^6}{24} & 0 \\ 0 & 0 & 0 & 0\end{bmatrix}

Meaning the basic properties are

.. math::
	A = \pi R^2 
.. math::
    Q_x = 0
.. math::
    Q_y = 0
.. math::
    I_{xx} = \dfrac{\pi R^4}{4}
.. math::
    I_{xy} = 0
.. math::
    I_{yy} = \dfrac{\pi R^4}{4}
.. math::
    I_{xxx} = 0
.. math::
    I_{xxy} = 0
.. math::
    I_{xyy} = 0
.. math::
    I_{yyy} = 0

Torsion properties
^^^^^^^^^^^^^^^^^^

The warping function is given by

.. math::
	\omega = 0
	
Leads to the torsional constant :math:`J`

.. math::
	\mathbb{J}_{\omega} = 0

.. math::
    J = \dfrac{\pi R^4}{2}

The torsion center :math:`\mathbf{T}`

.. math::
    \mathbf{T} = \left(0, \ 0\right)

Shear properties
^^^^^^^^^^^^^^^^

The shear center :math:`\mathbf{S}`

.. math::
    \mathbf{S} = \left(0, \ 0\right)

Strain and Stress
^^^^^^^^^^^^^^^^^

The strain and stresses are given by:

* Normal Force

.. math::
    \sigma_{zz}(x, y) = \dfrac{F_z}{A}

.. math::
    \boldsymbol{\sigma}(x, \ y) = \dfrac{F_z}{A} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}
	
.. math::
    \boldsymbol{\varepsilon}(x, \ y) = \dfrac{F_z}{EA}\begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}

* Bending Moments

.. math::
    \sigma_{zz}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4}

.. math::
    \mathbf{\sigma}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}

* Torsion Moment

* Shear Forces

Hollowed circle
---------------

Hollowed circle of external radius :math:`R_e` and internal radius :math:`R_i`

Geometric properties
^^^^^^^^^^^^^^^^^^^^

Computing the integrals

.. math::
    I_{a, b} = \int_{\Omega} x^a \cdot y^b \ d\Omega

.. math::
    I = \begin{bmatrix}\pi \left(R_e^2-R_i^2\right) & 0 & \frac{\pi \left(R_e^4-R_i^4\right)}{4} & 0 \\ 0 & 0 & 0 & 0 \\ \frac{\pi \left(R_e^4-R_i^4\right)}{4} & 0 & \frac{\pi \left(R_e^6-R_i^6\right)}{24} & 0 \\ 0 & 0 & 0 & 0\end{bmatrix}

Meaning the basic properties are

.. math::
	A = \pi \left(R_{e}^2 -R_{i}^2\right) 
.. math::
    Q_x = 0
.. math::
    Q_y = 0
.. math::
    I_{xx} = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right) }{4}
.. math::
    I_{xy} = 0
.. math::
    I_{yy} = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right) }{4}
.. math::
    I_{xxx} = 0
.. math::
    I_{xxy} = 0
.. math::
    I_{xyy} = 0
.. math::
    I_{yyy} = 0

The bending center :math:`\mathbf{B}`

.. math::
    \mathbf{B} = \left(0, \ 0\right)

Torsion properties
^^^^^^^^^^^^^^^^^^

The warping function is given by

.. math::
	\omega = 0
	
Leads to the torsional constant :math:`J`

.. math::
	\mathbb{J}_{\omega} = 0

.. math::
    J = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right)}{2}

The torsion center :math:`\mathbf{T}`

.. math::
    \mathbf{T} = \left(0, \ 0\right)

Shear properties
^^^^^^^^^^^^^^^^

The shear center :math:`\mathbf{S}`

.. math::
    \mathbf{S} = \left(0, \ 0\right)

Strain and Stress
^^^^^^^^^^^^^^^^^

The strain and stresses are given by:

* Normal Force

* Bending Moments

* Torsion Moment

* Shear Forces

Ellipse
-------

Ellipse of major axis :math:`a` and minor axis :math:`b`, centered in origin

Geometric properties
^^^^^^^^^^^^^^^^^^^^

.. math::
    I = \begin{bmatrix}\pi ab & 0 & \frac{\pi ab^3}{4} & 0 \\ 0 & 0 & 0 & 0 \\ \frac{\pi a^3b}{4} & 0 & \frac{\pi a^3b^3}{24} & 0 \\ 0 & 0 & 0 & 0\end{bmatrix}

Meaning the basic properties are

.. math::
	A = \pi ab
.. math::
    Q_x = 0
.. math::
    Q_y = 0
.. math::
    I_{xx} = \dfrac{\pi ab^3 }{4}
.. math::
    I_{xy} = 0
.. math::
    I_{yy} = \dfrac{\pi a^3b }{4}
.. math::
    I_{xxx} = 0
.. math::
    I_{xxy} = 0
.. math::
    I_{xyy} = 0
.. math::
    I_{yyy} = 0

The bending center :math:`\mathbf{B}`

.. math::
    \mathbf{B} = \left(0, \ 0\right)

Torsion properties
^^^^^^^^^^^^^^^^^^

The warping function

.. math::
    \omega(x, y) = xy

The torsion center :math:`\mathbf{T}`

.. math::
    \mathbf{T} = \left(0, \ 0\right)

Shear properties
^^^^^^^^^^^^^^^^

The shear center :math:`\mathbf{S}`

.. math::
    \mathbf{S} = \left(0, \ 0\right)

Strain and Stress
^^^^^^^^^^^^^^^^^

The strain and stresses are given by:

* Normal Force

.. math::
    \sigma_{zz}(x, y) = \dfrac{F_z}{A}

.. math::
    \boldsymbol{\sigma}(x, \ y) = \dfrac{F_z}{A} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}
	
.. math::
    \boldsymbol{\varepsilon}(x, \ y) = \dfrac{F_z}{EA}\begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}

* Bending Moments

.. math::
    \sigma_{zz}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4}

.. math::
    \mathbf{\sigma}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}

* Torsion Moment

* Shear Forces

Rectangle
---------

The fourth is a rectangle of base :math:`b` and height :math:`g`

Geometric properties
^^^^^^^^^^^^^^^^^^^^

.. math::
    I = \begin{bmatrix}LH & 0 & \dfrac{bh^3}{12} & 0 \\ 0 & 0 & 0 & 0 \\ \frac{b^3h}{12} & 0 & \frac{b^3h^3}{144} & 0 \\ 0 & 0 & 0 & 0\end{bmatrix}

Meaning the basic properties are

.. math::
	A = bh
.. math::
    Q_x = 0
.. math::
    Q_y = 0
.. math::
    I_{xx} = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right) }{4}
.. math::
    I_{xy} = 0
.. math::
    I_{yy} = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right) }{4}
.. math::
    I_{xxx} = 0
.. math::
    I_{xxy} = 0
.. math::
    I_{xyy} = 0
.. math::
    I_{yyy} = 0

Torsion properties
^^^^^^^^^^^^^^^^^^

The warping function

.. math::
	k_n = \dfrac{\pi\left(2n+1\right)}{2}
.. math::
	\omega = xy - \dfrac{8a^2}{\pi^3}\sum_{n=0}^{\infty} \dfrac{(-1)^n}{\left(2n+1\right)^3} \cdot \dfrac{\sin (k_n \cdot x)\sinh (k_n \cdot y)}{\cosh (k_n \cdot b)}

The torsion center :math:`\mathbf{T}`

.. math::
    \mathbf{T} = \left(0, \ 0\right)

Shear properties
^^^^^^^^^^^^^^^^

The shear center :math:`\mathbf{S}`

.. math::
    \mathbf{S} = \left(0, \ 0\right)

Strain and Stress
^^^^^^^^^^^^^^^^^

The strain and stresses are given by:

* Normal Force

.. math::
    \sigma_{zz}(x, y) = \dfrac{F_z}{A}

.. math::
    \boldsymbol{\sigma}(x, \ y) = \dfrac{F_z}{A} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}
	
.. math::
    \boldsymbol{\varepsilon}(x, \ y) = \dfrac{F_z}{EA}\begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}

* Bending Moments

.. math::
    \sigma_{zz}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4}

.. math::
    \mathbf{\sigma}(x, \ y) = \dfrac{4\left(M_{y} \cdot x + M_{x} \cdot y\right)}{\pi R^4} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{bmatrix}

* Torsion Moment

* Shear Forces

