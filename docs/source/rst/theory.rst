.. _theory:


============
Introduction
============

Notations

We use the notation :math:`\mathbf{u}` as a bidimensional vector, while :math:`\mathbf{u}` represents a :math:`n`-dimensional vector.

Operations between vectors
--------------------------

The inner product between two vectors is defined as :math:`\langle \mathbf{u}, \  \mathbf{v}\rangle`

.. math::
    \langle \mathbf{u}, \  \mathbf{v}\rangle = u_x \cdot v_x + u_y \cdot v_y

Between two :math:`n`-dimensional vectors

.. math::
    \langle \mathbf{u}, \  \mathbf{v}\rangle = \sum_{i} u_i \cdot v_i

The cross product between two bidimensional vectors :math:`\mathbf{u} \times \mathbf{v}` gives an scalar:

.. math::
    \mathbf{u} \times \mathbf{v} = u_x \cdot v_y - u_y \cdot v_x


===============
Basic geometric
===============



The most basic geometric properties are

* Area
.. math::
    A = \int_{\Omega} \ dx \ dy

* First moment of area

.. math::
    Q_x = \int_{\Omega} x \ dx \ dy
.. math::
    Q_y = \int_{\Omega} y \ dx \ dy

* Second moment of area

.. math::
    I_{xx} = \int_{\Omega} x^2 \ dx \ dy
.. math::
    I_{xy} = \int_{\Omega} x \cdot y \ dx \ dy
.. math::
    I_{yy} = \int_{\Omega} y^2 \ dx \ dy



=======
Torsion
=======

From Saint-venant theory, the warping function :math:`\omega(x, \ y)` is used to compute stresses when torsion is applied.

.. math::
    \nabla^2 \omega = 0

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \mathbf{n} \times \mathbf{p}

With :math:`\mathbf{p} = (x, \ y)` begin a point on the boundary. The boundary condition can be rewriten as

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \dfrac{\langle \mathbf{p}', \ \mathbf{p} \rangle}{\|\mathbf{p}'\|} 

This function is found by :ref:`boundary_element_method`.

Once :math:`\phi` is found, the torsion constant can be computed

.. math::
    J = I_{xx} + I_{yy} - \int_{\Omega} y \dfrac{\partial \omega}{\partial x} - x \dfrac{\partial \omega}{\partial y} \ dx \ dy

Since we use boundary methods, we rewrite as

.. math::
    J = I_{xx} + I_{yy} - \int_{\Gamma} \phi \cdot \dfrac{\langle \mathbf{p}, \ \mathbf{p}'\rangle}{\|\mathbf{p}'\|} \ d\Gamma


.. _shear_properties:

================
Shear properties
================

Main shear properties are:

* Shear center

.. note::
    For now, assume :math:`I_{xx} = I_{\bar{xx}}` and so on

Functions :math:`\Psi` and :math:`\Phi` are used:

.. math::
    \nabla^2 \Psi = 2\left(- I_{xx} \cdot x + I_{xy} \cdot y \right)

.. math::
    \nabla^2 \Phi = 2\left(I_{xy} \cdot x - I_{yy} \cdot y\right)

And boundary conditions

.. math::
    \left\langle\nabla \Psi, \ \mathbf{n}\right\rangle = \left\langle\mathbf{h}_{x}, \mathbf{n}\right\rangle
.. math::
    \left\langle \nabla \Phi, \ \mathbf{n}\right\rangle = \left\langle\mathbf{h}_{y}, \mathbf{n}\right\rangle
.. math::
    \mathbf{h}_{x} = \dfrac{\nu}{2}\left(I_{xx}\begin{bmatrix}1 & 0 & -1 \\ 0 & 1 & 0\end{bmatrix}+ I_{xy}\begin{bmatrix}0 & -1 & 0 \\ 1 & 0 & -1\end{bmatrix}\right)\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}
.. math::
    \mathbf{h}_{y} = \dfrac{\nu}{2}\left(I_{xy}\begin{bmatrix}-1 & 0 & 1 \\ 0 & -1 & 0\end{bmatrix}+ I_{yy}\begin{bmatrix}0 & 1 & 0 \\ -1 & 0 & 1\end{bmatrix}\right)\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}

Both equations are in fact Poisson equations.
We solve an equavalent problem by using the solution for homogeneous problem.
Set :math:`\Psi^{*}` and :math:`\Phi^{*}` as

.. math::
    \Psi^{*} = \dfrac{1}{4}\left(x^2+y^2\right)\left(-I_{xx} \cdot x + I_{xy} \cdot y\right)

.. math::
    \Phi^{*} = \dfrac{1}{4}\left(x^2+y^2\right)\left(I_{xy} \cdot x - I_{yy} \cdot y\right)

Note that :math:`\Psi^{\star} = \Psi - \Psi^{*}` and :math:`\Phi^{\star} = \Phi - \Phi^{*}` satisfy the Laplace equation with the boundary conditions

.. math::
    \nabla^2 \Psi^{\star} = 0

.. math::
    \nabla^2 \Phi^{\star} = 0

.. math::
    \left\langle \Psi^{\star}, \ \mathbf{n}\right\rangle = \left\langle \mathbf{g}_x, \mathbf{n}\right\rangle
.. math::
    \left\langle \Phi^{\star}, \ \mathbf{n}\right\rangle =\left\langle \mathbf{g}_{y}, \mathbf{n}\right\rangle

.. math::
    \mathbf{g}_x = \left(-I_{xx} \begin{bmatrix}\frac{3-2\nu}{4} & 0 & \frac{1+2\nu}{4} \\ 0 & \frac{1-2\nu}{4} & 0\end{bmatrix} + I_{xy}\begin{bmatrix}0 & \frac{1-2\nu}{4} & 0 \\ \frac{1+2\nu}{4} & 0 & \frac{3-2\nu}{4}\end{bmatrix}\right)\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}
.. math::
    \mathbf{g}_y = \left(I_{xy}\begin{bmatrix}\frac{3-2\nu}{4} & 0 & \frac{1+2\nu}{4} \\ 0 & \frac{1-2\nu}{4} & 0\end{bmatrix} - I_{yy}\begin{bmatrix}0 & \frac{1-2\nu}{4} & 0 \\ \frac{1+2\nu}{4} & 0 & \frac{3-2\nu}{4}\end{bmatrix}\right)\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}

======
Others
======

There are some other metrics like

.. math::
    Q_{\omega} = \int_{\Omega} \omega \ dx \ dy = \dfrac{1}{2}\int_{\Gamma} w \cdot \mathbf{p} \times \mathbf{p}' \ d\Gamma - \dfrac{1}{4}\int_{\Gamma}\langle \mathbf{p}, \ \mathbf{p}\rangle \cdot \dfrac{\langle \mathbf{p}, \ \mathbf{p}'\rangle}{\|\mathbf{p}'\|}

.. math::
    I_{x\omega} = \int_{\Omega} x \cdot \omega \ dx \ dy

.. math::
    I_{y\omega} = \int_{\Omega} y \cdot \omega \ dx \ dy
.. math::
    I_{\omega\omega} = \int_{\Omega} \omega^2 \ dx \ dy


.. _stress:

=================
Stress and Strain
=================

The stress in a beam is given by

.. math::
    \mathbf{\sigma} = \begin{bmatrix}\sigma_{xx} & \tau_{xy} & \tau_{xz} \\ \tau_{xy} & 0 & 0 \\ \tau_{xz} & 0 & 0\end{bmatrix}

While the strain is given by


.. math::
    \mathbf{\varepsilon} = \begin{bmatrix}\varepsilon_{xx} & \varepsilon_{xy} & \varepsilon_{xz} \\ \varepsilon_{xy} & 0 & 0 \\ \varepsilon_{xz} & 0 & 0\end{bmatrix}



.. _boundary_element_method:

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
    \xi\left(\mathbf{s}\right) \cdot u\left(\mathbf{s}\right) = \int_{\Gamma} u \cdot \dfrac{\partial v}{\partial n} \ d\Gamma - \int_{\Gamma} \dfrac{\partial u}{\partial n}  \cdot v \ d\Gamma

Which :math:`\mathbf{s}` is the source point of the Green function :math:`v` and :math:`\xi(\mathbf{s})` is the winding number.

.. math::
    v = \dfrac{1}{2\pi} \ln r = \dfrac{1}{2\pi} \ln \|\mathbf{p} - \mathbf{s}\|

Since all laplace's equations so far have only Neumann's boundary condtions, it known all the values of :math:`\dfrac{\partial u}{\partial n}`, the principia is find all the values of :math:`u` at the boundary.
Once :math:`u` and :math:`\dfrac{\partial u}{\partial n}` known at the boundary, it's possible to compute :math:`u(x, y)` at any point inside.

Parametrizing the curve :math:`\Gamma` by :math:`\mathbf{p}(t)`, fixing the source point :math:`\mathbf{s}_i = \mathbf{p}(t_i)` at the boundary, and setting :math:`u` as a linear combination of :math:`n` basis functions :math:`\varphi` and weights :math:`\mathbf{U}`

.. math::
    u(t) = \sum_{j=0}^{n-1} \varphi_j(t) \cdot U_j

.. math::
    \sum_{j=0}^{n-1} E_{ij} \cdot U_{j} = \sum_{j=0}^{n-1} H_{ij} \cdot U_{j} - G_{ij}

With the auxiliar values which depends only on the geometry, the source point and the basis functions

.. math::
    E_{ij} = 2\pi \cdot \xi\left(\mathbf{s}_i\right) \cdot \varphi_j\left(t_i\right)

.. math::
    H_{ij} = 2\pi \int_{\Gamma} \varphi_j \cdot \dfrac{\partial v_i}{\partial n} \ d\Gamma

.. math::
    G_{i} = 2\pi \int_{\Gamma} \dfrac{\partial u}{\partial n} \cdot v_i \ d\Gamma

Applying for :math:`n` different source points :math:`\mathbf{s}_i`, we get the matrices :math:`E`, :math:`H` and :math:`\mathbf{G}` such

.. math::
    \left(H-E\right) \cdot \mathbf{U} = \mathbf{G}

Finding the values of :math:`\mathbf{U}` means solving that system.

Matrix H
^^^^^^^^

We use

.. math::
    \dfrac{\partial v}{\partial n} ds = \dfrac{1}{2\pi} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle}

to write

.. math::
    H_{ij} = \int_{t_{min}}^{t_{max}} \varphi_{j}(t) \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle} \ dt

Vector :math:`G`
^^^^^^^^^^^^^^^^

This depends on the Neumann's boundary condition.

* For warping function

    .. math::
        \dfrac{\partial u}{\partial n} = \mathbf{n} \times \mathbf{p}

    .. math::
        G_i = \int_{t_{min}}^{t_{max}} \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \cdot \ln r \ dt

* For shear properties

    .. math::
        \dfrac{\partial u}{\partial n} = \left\langle \mathbf{h}, \ \mathbf{n}\right\rangle
    
    .. math::
        G_i = \int_{t_{min}}^{t_{max}} \begin{bmatrix}y' & -x'\end{bmatrix}\begin{bmatrix}\square & \square & \square \\ \square & \square & \square \end{bmatrix}\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix} \cdot \ln r \ dt

Winding number
^^^^^^^^^^^^^^^

This number is the mesure for a given point with respect to its position to the domain :math:`\Omega`.

.. math::
    \xi\left(\mathbf{s}\right) = \begin{cases}0 \ \ \ \ \ \ \ \ \text{if} \ \mathbf{s} \notin \Omega \\ \dfrac{\alpha}{2\pi} \ \ \ \ \text{if} \ \mathbf{s} \in \partial \Omega \\   1 \ \ \ \ \ \ \ \ \text{if} \ \mathbf{s} \in \Omega \end{cases}

Now, suppose that :math:`\mathbf{s}` is on the boundary. Then exists a value :math:`\tau` such :math:`\mathbf{p}(\tau) = \mathbf{s}` and the angle :math:`\alpha` is computed by

.. math::
    \mathbf{v}_0 = \lim_{\delta \to 0^{+}} \mathbf{p}'\left(\tau - \delta\right)

.. math::
    \mathbf{v}_1 = \lim_{\delta \to 0^{+}} \mathbf{p}'\left(\tau + \delta\right)

.. math::
    \alpha = \arg\left(\langle\mathbf{v_0}, \ \mathbf{v_1} \rangle + i \cdot \left(\mathbf{v_0} \times \mathbf{v_1}\right)\right)

.. note::
    In python code, it's in fact used ``alpha = arctan2(cross(v0, v1), inner(v0, v1))``

For smooth regions, the first derivative of :math:`\mathbf{p}` is continuous and therefore then :math:`\alpha = \pi`.


Computing matrices
^^^^^^^^^^^^^^^^^^

The matrices highly depend on the basis functions :math:`\varphi`.

To compute the coefficients :math:`H_{ij}` and :math:`G_{i}`, it's used numerical integration, like Gaussian-Quadrature.
Unfortunatelly, when :math:`r\approx 0` the integrants are singular and special techniques are required.

The main idea to compute them is decompose the integral in intervals, use standard techniques for intervals which :math:`r\ne 0` (called outside integration), and make special methods for elements when :math:`r=0` inside an interval (called inside integration).

* **Outside integration**: uses :ref:`regular_integrals` for elements which :math:`r\ne 0` for every point in the element

* **Inside integration**: uses :ref:`singular_integrals` for elements which :math:`r=0`







Isoparametric linear
--------------------

We restrict the geometry to polygons, present the piecewise linear basis functions and we also suppose the source points lays always in the vertices.

Let's say the parametric space :math:`t` is divided by the knots :math:`t_0`, :math:`t_1`, :math:`\cdots`, :math:`t_{n-1}`, :math:`t_n`

The function :math:`\varphi_{j}(t)` is represented by

.. math::
    \varphi_{j}(t) = \begin{cases}\frac{t - t_{j-1}}{t_{j}-t_{j-1}} \ \ \ \ \ \text{if} \ t_{j-1} \le t < t_{j} \\ \frac{t_{j+1} - t}{t_{j+1}-t_{j}} \ \ \ \ \ \text{if} \ t_{j} \le t < t_{j+1} \\ 0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{else}  \end{cases}

Call :math:`\mathbf{P}_{j} = (x_j, \ y_j)` the vertex :math:`j`, then

.. math::
    \mathbf{p}(t) = \sum_{j=0}^{n} \varphi_{j}(t) \cdot \mathbf{P}_{j}

For the interval :math:`\mathbb{I}_{j} = \left(t_{j}, \ t_{j+1}\right)`

.. math::
    \mathbf{p}(t) = \dfrac{t-t_{j}}{t_{j+1}-t_{j}} \cdot \mathbf{P}_{j} + \dfrac{t-t_{j}}{t_{j+1}-t_{j}} \cdot \mathbf{P}_{j+1}

.. math::
    \mathbf{p}'(t) = \mathbf{P}_{j+1} - \mathbf{P}_j

And then

.. math::
    \mathbf{r}(t) = \mathbf{p}(t) - \mathbf{P}_i

Matrix :math:`H`
^^^^^^^^^^^^^^^^

Since 

.. math::
    H_{ij} = \int_{t_{min}}^{t_{max}} \varphi_{j} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} dt = \sum_{k=0}^{n-1} H_{ijk} 

.. math::
    H_{ijk} = \int_{t_{k}}^{t_{k+1}} \varphi_{j} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} dt

* For outside integral, when :math:`k\ne i-1` or :math:`k \ne i`, then $H_{ij}$ is computed by :ref:`regular_integrals`.

* For the integrals over the intervals :math:`\left[t_{i-1}, \ t_{i}\right]` and :math:`\left[t_{i}, \ t_{i+1}\right]`, then we divide it into parts:


.. math::
    H_{ij} = \sum_{k=0}^{n-1} \int_{t_{k}}^{t_{k+1}} \varphi_{j} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} \ dt = \sum_{k=0}^{n-1} \int_{0}^{1} \varphi_{j} \cdot \dfrac{\left(\mathbf{p}_j\times\mathbf{p}_i\right)}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} \ dt



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

.. _regular_integrals:

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