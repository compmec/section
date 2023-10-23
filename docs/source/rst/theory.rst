.. _theory:


============
Introduction
============

Notations

We use the notation :math:`\vec{u}` as a bidimensional vector, while :math:`\mathbf{u}` represents a :math:`n`-dimensional vector.

Operations between vectors
--------------------------

The inner product between two vectors is defined as :math:`\langle \vec{u}, \  \vec{v}\rangle`

.. math::
    \langle \vec{u}, \  \vec{v}\rangle = u_x \cdot v_x + u_y \cdot v_y

Between two :math:`n`-dimensional vectors

.. math::
    \langle \mathbf{u}, \  \mathbf{v}\rangle = \sum_{i} u_i \cdot v_i

The cross product between two bidimensional vectors :math:`\vec{u} \times \vec{v}` gives an scalar:

.. math::
    \vec{u} \times \vec{v} = u_x \cdot v_y - u_y \cdot v_x


===============
Basic geometric
===============



The geometric properties

.. math::
    A = \int_{\Omega} \ dx \ dy
.. math::
    S_x = \int_{\Omega} x \ dx \ dy
.. math::
    S_y = \int_{\Omega} y \ dx \ dy
.. math::
    I_{xx} = \int_{\Omega} x^2 \ dx \ dy
.. math::
    I_{xy} = \int_{\Omega} x \cdot y \ dx \ dy
.. math::
    I_{yy} = \int_{\Omega} y^2 \ dx \ dy


=======
Torsion
=======

From Saint-venant theory, the warping function :math:`\phi(x, \ y)` is used to compute stress when torsion moments are applied.

.. math::
    \nabla^2 \phi = 0

.. math::
    \left\langle \nabla \phi, \ \vec{n}\right\rangle = \vec{n} \times \vec{p}

With :math:`\vec{p} = (x, \ y)`.

This function is computed by boundary element method (see bellow).
Once :math:`\phi` is found, the stress can be computed 





=======================
Boundary element method
=======================

Introduction
------------

The boundary element method is used to find numerically the Laplace equation

.. math::
    \nabla^2 u = 0

This method transforms the PDE into a boundary version

.. math::
    \xi\left(\vec{s}\right) \cdot u\left(\vec{s}\right) = \int_{\Gamma} u \cdot \dfrac{\partial v}{\partial n} \ d\Gamma - \int_{\Gamma} \dfrac{\partial u}{\partial n}  \cdot v \ d\Gamma

Which :math:`\vec{s}` is the source point of the Green function :math:`v` and :math:`\xi(\vec{s})` is the winding number.

.. math::
    v = \dfrac{1}{2\pi} \ln r = \dfrac{1}{2\pi} \ln \|\vec{p} - \vec{s}\|

Since the gives all the values of :math:`\dfrac{\partial u}{\partial n}`, the principia is find all the values of :math:`u` at the boundary. With :math:`u` and :math:`\dfrac{\partial u}{\partial n}` known at the boundary, it's possible to compute :math:`u(x, y)` at any point inside.

Parametrizing the curve :math:`\Gamma` by :math:`\vec{p}(t)`, fixing the source point :math:`\vec{s}_i = \vec{p}(t_i)` at the boundary, and setting :math:`u` as a linear combination of :math:`n` basis functions :math:`\varphi` and weights :math:`\mathbf{U}`

.. math::
    u(t) = \sum_{j=0}^{n-1} \varphi_j(t) \cdot U_j

.. math::
    \sum_{j=0}^{n-1} E_{ij} \cdot U_{j} = \sum_{j=0}^{n-1} H_{ij} \cdot U_{j} - G_{ij}

With the auxiliar values which depends only on the geometry, the source point and the basis functions

.. math::
    E_{ij} = 2\pi \cdot \xi\left(\vec{s}_i\right) \cdot \varphi_j\left(t_i\right)

.. math::
    H_{ij} = 2\pi \int_{\Gamma} \varphi_j \cdot \dfrac{\partial v_i}{\partial n} \ d\Gamma

.. math::
    G_{i} = 2\pi \int_{\Gamma} \dfrac{\partial u}{\partial n} \cdot v_i \ d\Gamma

Applying for :math:`n` different source points :math:`\vec{s}_i`, we get the matrices :math:`E`, :math:`H` and :math:`\mathbf{G}` such

.. math::
    \left(H-E\right) \cdot \mathbf{U} = \mathbf{G}

Finding the values of :math:`\mathbf{U}` means solving that system.

Matrix H
^^^^^^^^

We use

.. math::
    \dfrac{\partial v}{\partial n} ds = \dfrac{1}{2\pi} \cdot \dfrac{\vec{r} \times \vec{p}'}{\left\langle\vec{r}, \ \vec{r}\right\rangle}

to write

.. math::
    H_{ij} = \int_{t_{min}}^{t_{max}} \varphi_{j}(t) \cdot \dfrac{\vec{r} \times \vec{p}'}{\left\langle\vec{r}, \ \vec{r}\right\rangle} \ dt

Vector :math:`G`
^^^^^^^^^^^^^^^^

This depends on the Neumann's boundary condition.

* For warping function

    .. math::
        \dfrac{\partial u}{\partial n} = \left\langle \nabla u, \ \vec{n}\right\rangle = \vec{n} \times \vec{p}

    .. math::
        G_i = \int_{t_{min}}^{t_{max}} \left\langle p, \ p'\right\rangle \cdot \ln r

Winding number
^^^^^^^^^^^^^^^

This number is the mesure if a given point is inside the domain :math:`\Omega`.

.. math::
    \xi\left(\vec{s}\right) = \begin{cases}0 \ \ \ \ \ \ \ \ \text{if} \ \vec{s} \notin \Omega \\ \dfrac{\alpha}{2\pi} \ \ \ \ \text{if} \ \vec{s} \in \partial \Omega \\   1 \ \ \ \ \ \ \ \ \text{if} \ \vec{s} \in \Omega \end{cases}

Now, suppose that :math:`\vec{s}` is on the boundary. Then exists a value :math:`\tau` such :math:`\vec{p}(\tau) = \vec{s}` and the angle :math:`\alpha` is computed by

.. math::
    \vec{v}_0 = \lim_{\delta \to 0^{+}} \vec{p}'\left(\tau - \delta\right)

.. math::
    \vec{v}_1 = \lim_{\delta \to 0^{+}} \vec{p}'\left(\tau + \delta\right)

.. math::
    \alpha = \arg\left(\langle\vec{v_0}, \ \vec{v_1} \rangle + i \cdot \left(\vec{v_0} \times \vec{v_1}\right)\right)

.. note::
    In python code, it's in fact used ``alpha = arctan2(cross(v0, v1), inner(v0, v1))``

If :math:`\vec{p}\left(\tau\right)` is not a corner (it's a smooth region and first derivative of :math:`\vec{p}` is continuous), then the winding number is :math:`\dfrac{1}{2}`.


Computing matrices
^^^^^^^^^^^^^^^^^^

The matrices highly depend on the basis functions :math:`\varphi`, which were not yet defined.

In the expressions of :math:`H_{ij}` and :math:`G_{i}`, there are the terms :math:`\dfrac{1}{r}` and :math:`\ln r`, which makes these integrals singular.
The main idea to compute them is decompose the integral in intervals, use standard techniques for intervals which :math:`r\ne 0` (called outside integration), and make special methods for elements when :math:`r=0` inside an interval.





Isoparametric linear
--------------------

We restrict the geometry to polygons, present the piecewise linear basis functions and we also suppose the source points lays always in the vertices.

Let's say the parametric space :math:`t` is divided by the knots :math:`t_0`, :math:`t_1`, :math:`\cdots`, :math:`t_{n-1}`, :math:`t_n`

The function :math:`\varphi_{j}(t)` is represented by

.. math::
    \varphi_{j}(t) = \begin{cases}\frac{t - t_{j-1}}{t_{j}-t_{j-1}} \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{if} \ t_{j-1} \le t < t_{j} \\ \left(\frac{t_{j+1} - t}{t_{j+1}-t_{j}}\right) \ \ \ \ \ \text{if} \ t_{j} \le t < t_{j+1} \\ 0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{else}  \end{cases}

Call :math:`\vec{P}_{j} = (x_j, \ y_j)` the vertex :math:`j`, then

.. math::
    \vec{p}(t) = \sum_{j=0}^{n} \varphi_{j}(t) \cdot \vec{P}_{j}

For the interval :math:`\left(t_{j}, \ t_{j+1}\right)`

.. math::
    \vec{p}(t) = \dfrac{t-t_{j}}{t_{j+1}-t_{j}} \cdot \vec{P}_{j} + \dfrac{t-t_{j}}{t_{j+1}-t_{j}} \cdot \vec{P}_{j+1}

.. math::
    \vec{p}'(t) = \vec{P}_{j+1} - \vec{P}_j




=========
Integrals
=========

Polynomial integrals
--------------------

To compute area, momentums and inertias, it's needed to compute the integral

.. math::
    I_{a,b} = \int_{\Omega} x^a \cdot y^b \ dx \ dy

Which :math:`\Omega` is the defined region with closed boundary :math:`\Gamma`.

By using Green's thereom, we transform the integral

.. math::
    \int_{\Omega} \left(\dfrac{\partial Q}{\partial x} - \dfrac{\partial P}{\partial y}\right) \ dx \ dy = \int_{\Gamma} P \ dx + Q \ dy

Without loss of generality, let :math:`\alpha \in \mathbb{R}` and take

.. math::
    \dfrac{\partial Q}{\partial x} = \alpha \cdot x^a \cdot y^b \Longrightarrow Q = \dfrac{\alpha}{a+1} \cdot x^{a+1} \cdot y^b

.. math::
    \dfrac{\partial P}{\partial y} = \left(\alpha-1\right) \cdot x^a \cdot y^b \Longrightarrow P = \dfrac{\alpha - 1}{b+1} \cdot x^{a} \cdot y^{b+1}

Then

.. math::
    I_{a, b} = \dfrac{\alpha - 1}{b+1} \int_{\Gamma} x^{a} \cdot y^{b+1} \ dx + \dfrac{\alpha}{a+1} \int_{\Gamma} x^{a+1} \cdot y^b \ dy

Regular integrals
------------------

The numerical integral are computated by using quadrature schemas, rewriting

.. math::
    \int_{0}^{1} f(x) \ dx = \sum_{i=0}^{n-1} w_i \cdot f(x_i)

With specific position nodes :math:`x_i` and weights :math:`w_i`. 

Here we present some possible quadratures

* Closed Newton Cotes: Equally spaced points in interval. Degree at most :math:`p-1` with :math:`p` evaluation points

* Chebyshev: `Chebyshev nodes <https://en.wikipedia.org/wiki/Chebyshev_nodes>`_ in interval. Degree at most :math:`p-1` with :math:`p` evaluation points

* `Gauss-Legendre Quadrature <https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature>`_: 

* `Gauss-Legendre Quadrature <https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature>`_

* Lobatto Quadrature: Can be used to adaptative quadrature

* `Clenshaw–Curtis Quadrature <https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature>`_

.. _singular_integrals:

Singular integrals
------------------

There are two types of singular integrals to compute:

.. math::
    \int_{0}^{1} f(x) \cdot \ln x \ dx

.. math::
    \int_{-1}^{1} f(x) \cdot \dfrac{1}{x} \ dx

Logarithm singularity
^^^^^^^^^^^^^^^^^^^^^

We are interested in computing the integral

.. math::
    I = \int_{0}^{1} f(x) \ \cdot \ln x \ dx

If the function :math:`f(x)` is described by using series

.. math::
    f(x) = \sum_{i=0}^{\infty} a_i \cdot x^{i}

Then the integral is 

.. math::
    I = - \sum_{i=0}^{\infty} \dfrac{a_i}{\left(1+i\right)^2}

Which is well defined as long as :math:`f(x)` is a polynomial.

A logarithm quadrature was created by `Stroud and Sladek <https://www.sciencedirect.com/science/article/abs/pii/S0045782597002399>`_ with given values in table bellow

.. math::
    \int_{0}^{1} f(x)\ln x \ dx = \sum_{k=1}^{p} w_{k} \cdot f(\eta_{k})

.. list-table:: Nodes and Weights for Logarithm Quadrature 
   :widths: 20 40 40
   :header-rows: 1
   :align: center

   * - :math:`p`
     - :math:`\eta`
     - :math:`w`
   * - 2
     - 0.112008806166976
     - 0.718539319030384
   * - 
     - 0.602276908118738
     - 0.281460680969615
   * - 
     - 
     - 
   * - 3
     - 0.0638907930873254
     - 0.513404552232363
   * - 
     - 0.368997063715618
     - 0.391980041201487
   * - 
     - 0.766880303938941
     - 0.0946154065661491

    
Odd singularity
^^^^^^^^^^^^^^^

We are interested in computing the integral

.. math::
    \int_{-1}^{1} \dfrac{1}{x} \cdot f(x) \ dx

The given integral is computed as the Cauchy Principal Value

.. math::
    PV\int_{-1}^{1} \dfrac{f(x)}{x} \ dx = \lim_{\varepsilon \to 0^{+}} \int_{-1}^{-\varepsilon} \dfrac{f(x)}{x} \ dx + \int_{\varepsilon}^{1} \dfrac{f(x)}{x} \ dx 

This integral is well defined if :math:`f(x)` is a polynomial:

.. math::
    PV\int_{-1}^{1} \dfrac{1}{x} \ dx = 0
.. math::
    PV\int_{-1}^{1} \dfrac{x}{x} \ dx = 2
.. math::
    PV\int_{-1}^{1} \dfrac{x^2}{x} \ dx = 0

Expanding :math:`f(x)` by its coefficients, therefore

.. math::
    PV \int_{-1}^{1} \dfrac{1}{x} \cdot f(x) \ dx = \sum_{i=1}^{\infty} a_{i} \cdot \dfrac{1 + \left(-1\right)^{i+1}}{i} = \sum_{j=0}^{\infty} \dfrac{2}{2j+1} \cdot a_{2j+1}

It's possible to create a quadrature for it:

TO DO