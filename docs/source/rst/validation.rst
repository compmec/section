.. _validation:


============
Introduction
============

The validation of the code happens in 2 steps:

1. Numerical:
    Checks if the output values correspond with given input, following the :ref:`Theory <theory>`
2. Literature:
    Compares the output values with the values found in textbooks and references


====================
Numerical Validation
====================


The first case is for a circular bar of radius :math:`R`, centered in the origin


The basic properties are

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
	
The warping function is given by

.. math::
	\omega = 0
	
Leads to the torsional constant :math:`J`

.. math::
	\mathbb{J}_{\omega} = 0

.. math::
    J = \dfrac{\pi R^4}{2}
	
The shear...


The strain stresses are given by:

* Normal:

.. math::
    \boldsymbol{\sigma} = \dfrac{F_z}{A} \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1\end{bmatrix}
	
.. math::
    \boldsymbol{\varepsilon} = \dfrac{F_z}{EA}\begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}
	

* Bending Moments:

.. math::
    \sigma_{zz} = \dfrac{F_z}{A}



The second case is a Hollowed Circle of external radius :math:`R_e` and internal radius :math:`R_i`



The basic properties are

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
	
The warping function is given by

.. math::
	\omega = 0
	
Leads to the torsional constant :math:`J`

.. math::
	\mathbb{J}_{\omega} = 0

.. math::
    J = \dfrac{\pi \left(R_{e}^4 -R_{i}^4\right) }{2}
	
The shear...





The third case is an ellipse

The fourth is a rectangle of base :math:`a` and height :math:`b`

.. math::
	A = ab
	
.. math::
	k_n = \dfrac{\pi\left(2n+1\right)}{2}
.. math::
	\omega = xy - \dfrac{8a^2}{\pi^3}\sum_{n=0}^{\infty} \dfrac{(-1)^n}{\left(2n+1\right)^3} \cdot \dfrac{\sin (k_n \cdot x)\sinh (k_n \cdot y)}{\cosh (k_n \cdot b)}


