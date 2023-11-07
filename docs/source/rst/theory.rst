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


.. _boundary_element_method:

=======================
Boundary element method
=======================

Introduction
------------

The boundary element method (BEM for short) is used to find :math:`u` numerically

.. math:: 
    :label: eq_laplace

    \nabla^2 u = 1

The BEM transforms :eq:`eq_laplace` into a boundary version :eq:`eq_bem`

.. math::
    :label: eq_bem

    \alpha\left(\mathbf{s}\right) \cdot u\left(\mathbf{s}\right) = \int_{\Gamma} u \cdot \dfrac{\partial v}{\partial n} \ d\Gamma - \int_{\Gamma} \dfrac{\partial u}{\partial n}  \cdot v \ d\Gamma

Which :math:`\mathbf{s}` is the source point of the Green function :math:`v` and :math:`\alpha(\mathbf{s})` is the angle at the point :math:`\mathbf{s}`.

.. math::
    :label: eq_source

    v(\mathbf{p}, \ \mathbf{s}) = \ln r = \ln \|\mathbf{r}\| = \ln \|\mathbf{p} - \mathbf{s}\|

Since all the PDEs used in this package have only Neumann's boundary conditions, the values of :math:`\dfrac{\partial u}{\partial n}` are known and the objective is finding all the values of :math:`u` at the boundary.
Once :math:`u` and :math:`\dfrac{\partial u}{\partial n}` are known at the boundary, it's possible to compute :math:`u(x, y)` at any point inside by using :eq:`eq_bem`, as well the other properties


Discretize solution
-------------------

Parametrize the curve :math:`\Gamma` by :math:`\mathbf{p}(t)`, fix the source point :math:`\mathbf{s}_i = \mathbf{p}(t_i)` at the boundary, and set :math:`u` as a linear combination of :math:`n` basis functions :math:`\varphi` and weights :math:`\mathbf{U}`

.. math::
    :label: eq_curve_param

    \mathbf{p}(t) = \sum_{j=0}^{m-1} \phi_{j}(t) \cdot P_{j} = \langle \mathbf{\phi}(t), \ \mathbf{P}\rangle

.. math::
    :label: eq_discret_func

    u(t) = \sum_{j=0}^{n-1} \varphi_j(t) \cdot U_j = \langle \mathbf{\varphi}(t), \ \mathbf{U}\rangle

Expanding :eq:`eq_bem` and using :eq:`eq_discret_func`, :eq:`eq_matrix_formula` is obtained

.. math::
    :label: eq_matrix_formula

    \sum_{j=0}^{n-1} A_{ij} \cdot U_{j} = \sum_{j=0}^{n-1} M_{ij} \cdot U_{j} - F_{i}

With the auxiliar values which depends only on the geometry, the source point and the basis functions

.. math::
    A_{ij} = \alpha\left(\mathbf{s}_i\right) \cdot \varphi_j\left(t_i\right)

.. math::
    M_{ij} = \int_{\Gamma} \varphi_j \cdot \dfrac{\partial v_i}{\partial n} \ d\Gamma

.. math::
    F_{i} = \int_{\Gamma} \dfrac{\partial u}{\partial n} \cdot v_i \ d\Gamma

Applying for :math:`n` different source points :math:`\mathbf{s}_i` at boundary, we get the matrices :math:`A`, :math:`M` and :math:`\mathbf{F}` such

.. math::
    :label: eq_linear_system

    \left(M-A\right) \cdot \mathbf{U} = K \cdot \mathbf{U} = \mathbf{F}

Finding the values of :math:`\mathbf{U}` means solving the linear system :eq:`eq_linear_system`

Matrix M
^^^^^^^^

We use

.. math::
    \dfrac{\partial v}{\partial n} ds = \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle}

to write

.. math::
    M_{ij} = \int_{t_{min}}^{t_{max}} \varphi_{j}(t) \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle} \ dt

Vector :math:`F`
^^^^^^^^^^^^^^^^

This vector depends on the Neumann's boundary condition.

* For warping function

    .. math::
        \dfrac{\partial u}{\partial n} = \mathbf{n} \times \mathbf{p} = \dfrac{\langle \mathbf{p}, \ \mathbf{p}'\rangle}{\|\mathbf{p}'\|}

    .. math::
        F_i = \int_{t_{min}}^{t_{max}} \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \cdot \ln \|\mathbf{r}_i\| \ dt

* For shear properties

    .. math::
        \dfrac{\partial u}{\partial n} = \left\langle \mathbf{h}, \ \mathbf{n}\right\rangle = \dfrac{\mathbf{h} \times \mathbf{p}'}{\|\mathbf{p}'\|}
    
    .. math::
        F_i = \int_{t_{min}}^{t_{max}} \begin{bmatrix}y' & -x'\end{bmatrix}\left[\bar{F}\right]\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix} \cdot \ln r \ dt

    Which :math:`\left[\bar{F}\right]` is the :math:`(2 \times 3)` matrix shown in :ref:`shear_properties`.

Angle at boundary
^^^^^^^^^^^^^^^^^

The angle :math:`\alpha` is the mesure for a given point with respect to its position to the domain :math:`\Omega`.

.. math::
    \alpha\left(\mathbf{s}\right) = \begin{cases}\in \left(0, \ 2\pi\right) \ \ \ \ \text{if} \ \mathbf{s} \in \partial \Omega \\ 0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{if} \ \mathbf{s} \notin \text{closed}\left(\Omega\right) \\   2\pi \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \text{if} \ \mathbf{s} \in \text{open}\left(\Omega\right) \end{cases}

When :math:`\mathbf{s} \in \partial \Omega`, there is a value :math:`\tau` such :math:`\mathbf{p}(\tau) = \mathbf{s}` and the angle :math:`\alpha` is computed by

.. math::
    \mathbf{v}_0 = -\lim_{\delta \to 0^{+}} \mathbf{p}'\left(\tau - \delta\right)

.. math::
    \mathbf{v}_1 = \lim_{\delta \to 0^{+}} \mathbf{p}'\left(\tau + \delta\right)

.. math::
    \alpha = \arg\left(\langle\mathbf{v_0}, \ \mathbf{v_1} \rangle + i \cdot \left(\mathbf{v_0} \times \mathbf{v_1}\right)\right)

For smooth regions, the first derivative of :math:`\mathbf{p}` is continuous and therefore then :math:`\alpha = \pi`.

.. note::
    In python code, it's in fact used ``alpha = arctan2(cross(v0, v1), inner(v0, v1))``

Computing matrices
^^^^^^^^^^^^^^^^^^

The matrices highly depend on the geometry and the basis functions :math:`\varphi`.

To compute the coefficients :math:`M_{ij}` and :math:`F_{i}`, it's used numerical integration, like Gaussian-Quadrature.
Unfortunatelly, when :math:`r = 0` at some point, the integrants are singular and special techniques are used.

The main idea to compute them is decompose the integral in intervals and use

* **Outside integration**: uses :ref:`regular_integrals` for elements which :math:`r\ne 0` for every point in the element

* **Inside integration**: uses :ref:`singular_integrals` for elements which :math:`r=0`

For polygonal domains, when :math:`\mathbf{p}(t)` is linear piecewise, the **Inside integration** is not required cause it can be done analiticaly. But for higher degrees, it's indeed necessary

.. _constraint_solution:

Constraint solution
^^^^^^^^^^^^^^^^^^^

All the PDEs have Neumann's boundary conditions, if :math:`u^{\star}` is a solution, then :math:`\left(u^{\star} + \text{const}\right)` also is a solution.
Although both functions give the same properties cause it envolves only the derivatives of :math:`u`, we restrict the solution by setting an exact value at a certain point.

* For warping function, we set :math:`u(\mathbf{c}_t) = 0`, which :math:`\mathbf{c}_t` is the twisting center

These values are obtained by minimizing the total deformation energy

The lagrange multiplier is used and the following system is solved instead:

.. math::

    \begin{bmatrix}K & \mathbf{C}^T \\ \mathbf{C} & 0\end{bmatrix} \begin{bmatrix}\mathbf{U} \\ \lambda \end{bmatrix} = \begin{bmatrix}\mathbf{F} \\ 0\end{bmatrix}

Which vector :math:`\mathbf{C}` is obtained from the constraint and depends on the problem: torsion or shear

Polygonal domain
----------------

For polygonal domains, when the basis functions :math:`\phi(t)` are piecewise linear, some computations becomes easier. Let's say the parametric space :math:`t` is divided by the knots :math:`t_0`, :math:`t_1`, :math:`\cdots`, :math:`t_{m-1}`, :math:`t_m`

For an arbitrary interval :math:`\left[t_k, \ t_{k+1}\right]`, :math:`\mathbf{p}(t)` is described as

.. math::
    \mathbf{p}(t) = \mathbf{P}_{k} + z \cdot \mathbf{V}_k
    
.. math::
    \mathbf{V}_k = \mathbf{P}_{k+1} - \mathbf{P}_{k}

.. math::
    z = \dfrac{t - t_{k}}{t_{k+1} - t_{k}} \in \left[0, \ 1\right]

Since the source point :math:`\mathbf{s}_i = \mathbf{p}(t_i)`,

* If :math:`t_i \in \left[t_{k}, \ t_{k+1}\right]` then

    .. math::
        \mathbf{r}(t) = \left(z-z_i\right) \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_{k}\right)

    .. math::
        z_i = \dfrac{t_i - t_{k}}{t_{k+1} - t_{k}}\in \left[0, \ 1\right]

* Else

    .. math::
        \mathbf{r}(t) = \left(\mathbf{P}_{k}-\mathbf{s}_i\right) + z \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_{k}\right)

Matrix :math:`M`
^^^^^^^^^^^^^^^^

.. math::
    M_{ij} = \sum_{k=0}^{m-1} \int_{t_{k}}^{t_{k+1}} \varphi_{j} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} \ dt

When :math:`t_i \in \left[t_k, \ t_{k+1}\right]`

.. math::
    \mathbf{p(t)} = \mathbf{P}_k + z \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_k\right) 
.. math::
    \mathbf{r(t)} = \left(z-z_i\right) \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_k\right)
.. math::
    \mathbf{r} \times \mathbf{p}' = 0 

Therefore, we can ignore the integration over the interval :math:`\left[t_k, \ t_{k+1}\right]`


Matrix :math:`F` for warping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For warping function, the expression :math:`F_i` is written as

.. math::
    \dfrac{\partial u}{\partial n} = \dfrac{\left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle}{\|\mathbf{p}'\|}

.. math::
    F_{i} = \sum_{k=0}^{m-1} F_{ik}
    
.. math::
    F_{ik} = \int_{0}^{1} \left(\alpha_k + z \cdot \beta_k \right) \ln\|\mathbf{r}\| \ dz

.. math::
    \alpha_k = \left\langle \mathbf{P}_k, \ \mathbf{V}_k\right\rangle

.. math::
    \beta_k = \left\langle \mathbf{V}_k, \ \mathbf{V}_k\right\rangle
    
When :math:`t_i \notin \left[t_k, \ t_{k+1}\right]`, the integral is computed by Gaussian Quadrature.

For :math:`t_i \notin \left[t_k, \ t_{k+1}\right]`

.. math::
    \mathbf{r} = \left(z-z_i\right) \cdot \mathbf{V}_k
.. math::
    F_{ik} = & \int_{0}^{1} \left(\alpha_k + z \beta_k \right) \ln\|\left(z-z_i\right) \cdot \mathbf{V}_k\| \ dz \\
           = & \left(\alpha_{k} + \dfrac{1}{2}\beta_{k}\right) \cdot \ln \|\mathbf{V}_k\| \\
             & + \alpha_{k} \int_{0}^{1} \ln |z-z_i| dz \\
             & + \beta_k \int_{0}^{1} z \cdot \ln |z-z_i| \ dz 

These two integrals are easy to compute, the expressions are complicated (`here <https://www.wolframalpha.com/input?i=int_%7B0%7D%5E%7B1%7D+ln%28abs%28x-x_0%29%29+dx%3B+0+%3C%3D+x_0+%3C%3D+1>`_ and `here <https://www.wolframalpha.com/input?i=int_%7B0%7D%5E%7B1%7D+x*ln%28abs%28x-x_0%29%29+dx%3B+0+%3C%3D+x_0+%3C%3D+1>`_) and depends on the value of :math:`z_i`. Bellow you find a table with some values

.. list-table:: Values of logarithm integrals
   :widths: 20 40 40
   :header-rows: 1
   :align: center

   * - :math:`z_i`
     - :math:`\int_0^1 \ln|z-z_i| dz`
     - :math:`\int_0^1 z\ln|z-z_i| dz`
   * - :math:`0`
     - :math:`-1`
     - :math:`\frac{-1}{4}`
   * - :math:`\frac{1}{2}`
     - :math:`-(1+\ln 2)`
     - :math:`\frac{-1}{2}\left(1+\ln 2\right)`
   * - :math:`1`
     - :math:`-1`
     - :math:`\frac{-3}{4}`

Therefore, the integral over interval which :math:`t_i` lies on is made by using analitic values, and singular integrals are not computed.



==================================
Cross-section geometric properties
==================================


Cross-section area
------------------

.. math::
    A = \int_{\Omega} \ dx \ dy


First moment of area
--------------------

.. math::
    Q_x = \int_{\Omega} x \ dx \ dy
.. math::
    Q_y = \int_{\Omega} y \ dx \ dy

Geometric centroid
------------------

.. math::
    x_{gc} = \dfrac{Q_x}{A}
.. math::
    y_{gc} = \dfrac{Q_y}{A}
.. math::
    \mathbf{c}_g = \left(x_{gc}, \ y_{gc}\right)

Second moment of area
---------------------

The global second moment of inertia are

.. math::
    I_{xx} = \int_{\Omega} x^2 \ dx \ dy
.. math::
    I_{xy} = \int_{\Omega} xy \ dx \ dy
.. math::
    I_{yy} = \int_{\Omega} y^2 \ dx \ dy
.. math::
    \mathbb{I} = \begin{bmatrix}I_{xx} & I_{xy} \\ I_{xy} & I_{yy}\end{bmatrix}

The local second moment of inertia are computed with respect to the geometric center

.. math::
    I_{\overline{xx}} = \int_{\Omega} (x-x_{gc})^2 \ dx \ dy = I_{xx} - \dfrac{Q_{x}^2}{A}
.. math::
    I_{\overline{xy}} = \int_{\Omega} (x-x_{gc})(y-y_{gc}) \ dx \ dy= I_{xy} - \dfrac{Q_{x}Q_{y}}{A}
.. math::
    I_{\overline{yy}} = \int_{\Omega} (y-y_{gc})^2 \ dx \ dy= I_{xy} - \dfrac{Q_{y}^2}{A}
.. math::
    \overline{\mathbb{I}} = \begin{bmatrix}I_{\overline{xx}} & I_{\overline{xy}} \\ I_{\overline{xy}} & I_{\overline{yy}}\end{bmatrix} = \mathbb{I} - A \cdot \mathbf{c}_g \otimes \mathbf{c}_g

    


Radius of Gyration
------------------

.. math::
    r_{x} = \sqrt{\dfrac{I_{xx}}{A}}
.. math::
    r_{y} = \sqrt{\dfrac{I_{yy}}{A}}


Principal Axis Properties
-------------------------

Let 

.. math::
    \overline{\mathbb{I}} = \begin{bmatrix}I_{\overline{xx}} & I_{\overline{xy}} \\ I_{\overline{xy}} & I_{\overline{yy}}\end{bmatrix}

The principals moment of inertia are the eigenvalues of :math:`\overline{\mathbb{I}}`.
But for a 2D matrix, :math:`I_{11}` and :math:`I_{22}` are easily calculated

.. math::
    \Delta = \sqrt{\left(\dfrac{I_{\overline{xx}}-I_{\overline{yy}}}{2}\right)^2+I_{\overline{xy}}^2}
.. math::
    I_{11} = \dfrac{I_{\overline{xx}}+I_{\overline{yy}}}{2} + \Delta
.. math::
    I_{22} = \dfrac{I_{\overline{xx}}+I_{\overline{yy}}}{2} - \Delta

The direction principal moment of inertia is the eigenvector related to the higher eigenvalue.
It's also computed as 

.. math::
    \phi = \arg\left(I_{\overline{xy}} + i \cdot \left(I_{\overline{xx}}-I_{11}\right)\right) = \text{arctan}\left(\dfrac{I_{\overline{xx}}-I_{11}}{I_{\overline{xy}}}\right)



===============================
Saint-Venant Torsion Properties
===============================

Warping Function
----------------

From Saint-venant theory, the warping function :math:`\omega(x, \ y)` is used to compute stresses when torsion is applied.

.. math::
    \nabla^2 \omega = 0

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \mathbf{n} \times \mathbf{p}

With :math:`\mathbf{p} = (x, \ y)` begin a point on the boundary. The boundary condition can be rewriten as

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \dfrac{\langle \mathbf{p}', \ \mathbf{p} \rangle}{\|\mathbf{p}'\|} 

This function is found by :ref:`boundary_element_method`.

Once :math:`\omega` is found, the torsion properties can be computed

Torsion center
---------------

As described in :ref:`constraint_solution`, we solve a Neumann's problem.
If :math:`\omega^{\star}` is a solution, then :math:`\omega^{*} = \omega^{\star} + y_0 x - x_0 y + c_0` is also a solution.

Choosing arbitrarily the values of :math:`x_0`, :math:`y_0` and :math:`c_0` doesn't change the torsion properties or the stresses due to torsion, it can be understood as a *rigid body rotation in the plane of cross-section and a displacement parallel to the axis of the bar* (from BOOK BEM).

The quantities :math:`x_0`, :math:`y_0` and :math:`c_0` can be obtained by minimizing the strain energy produced by axial normal warping stresses, which are ignored by Saint-Venant's theory.
Doing so, leads to the linear system

.. math::
    \left(\int_{\Omega} \begin{bmatrix}1 \\ x \\ y \end{bmatrix}\begin{bmatrix}1 & -y & x \end{bmatrix} \ d\Omega\right) \begin{bmatrix}c_0 \\ x_0 \\ y_0\end{bmatrix} = \int_{\Omega} \omega\begin{bmatrix}1 \\ x \\ y\end{bmatrix} \ d\Omega

The matrix on the left side is already computed by the values :math:`A`, :math:`Q_x`, :math:`Q_y`, :math:`I_{xx}`, :math:`I_{xy}`, :math:`I_{yy}`, while the values on the right side are

.. math::
    Q_{\omega} = \int_{\Omega} \omega \ dx \ dy
.. math::
    I_{x\omega} = \int_{\Omega} x \omega \ dx \ dy
.. math::
    I_{y\omega} = \int_{\Omega} y \omega \ dx \ dy

These integrals are transformed to the boundary equivalent.


.. dropdown:: Boundary reformulation of :math:`Q_{\omega}`, :math:`I_{x\omega}` and :math:`I_{y\omega}` 

    Let :math:`u` be a function such

    .. math::
        \nabla^2 u = g(x, y)

    Select :math:`u` respectivelly as
    
    .. math::
        u_{1} = \frac{1}{4}(x^2+y^2)
    
    .. math::
        u_{x} = \frac{x^3}{6}
     
    .. math::
        u_{y} = \frac{y^3}{6}
        
    and use Green's second identity

    .. math::
        \int_{\Omega} \omega \cdot g \ d\Omega & = \int_{\Omega} \omega \nabla^2 u - u \nabla^2 \omega \ d\Omega \\ & = \oint_{\Gamma} \omega \dfrac{\partial u}{\partial n} \ d\Gamma  - u \dfrac{\partial \omega}{\partial n} \ d\Gamma \\ & = \oint_{\Gamma} \omega \dfrac{\partial u}{\partial n} \ d\Gamma - \oint_{\Gamma} u \cdot \langle \mathbf{p}, \ \mathbf{p}'\rangle \ dt

    Transforming to

    .. math::
        Q_{\omega} = \dfrac{1}{2}\int_{t_{min}}^{t_{max}} \omega \cdot \mathbf{p} \times \mathbf{p}' \ dt - \dfrac{1}{4}\int_{t_{min}}^{t_{max}} \langle \mathbf{p}, \ \mathbf{p} \rangle \cdot \langle \mathbf{p}, \ \mathbf{p}' \rangle \ dt

    .. math::
        I_{x\omega} = \dfrac{1}{2} \oint_{\Gamma} \omega \cdot x^2 \ dy - \dfrac{1}{6}\int_{t_{min}}^{t_{max}} x^3 \cdot \langle \mathbf{p}, \ \mathbf{p}' \rangle  \ dt

    .. math::
        I_{y\omega} = \dfrac{-1}{2} \int_{t_{min}}^{t_{max}} \omega \cdot y^2 \ dx - \dfrac{1}{6}\int_{t_{min}}^{t_{max}} y^3 \cdot \langle \mathbf{p}, \ \mathbf{p}' \rangle  \ dt

.. note::
    Not implemented this part yet. The vector :math:`\mathbf{C}` for lagrange multiplier is a vector of ones.

Torsion constant
----------------

The torsion constant can be computed

.. math::
    J = I_{xx} + I_{yy} - \mathbb{J}_{\omega}

With

.. math::
    \mathbb{J}_{\omega} = \int_{\Omega} y \dfrac{\partial \omega}{\partial x} - x \dfrac{\partial \omega}{\partial y} \ dx \ dy

We transform this integral into a boundary one

.. math::
    \mathbb{J}_{\omega} = \int_{t_{min}}^{t_{max}} \omega \cdot \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \ dt

Since :math:`\omega = \langle \varphi, \ \mathbf{W}\rangle`, then

.. math::
    \mathbb{J}_{\omega} = \langle \Lambda, \mathbf{W} \rangle

.. math::
    \Lambda_j = \int_{t_{min}}^{t_{max}} \varphi_j \cdot \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \ dt

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
    \mathbf{g}_x = \left[F_x\right]\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}
.. math::
    \mathbf{g}_y = \left[F_y\right]\begin{bmatrix}x^2 \\ 2xy \\ y^2\end{bmatrix}

.. math::
    \left[F_x\right] = -I_{xx} \begin{bmatrix}\frac{3-2\nu}{4} & 0 & \frac{1+2\nu}{4} \\ 0 & \frac{1-2\nu}{4} & 0\end{bmatrix} + I_{xy}\begin{bmatrix}0 & \frac{1-2\nu}{4} & 0 \\ \frac{1+2\nu}{4} & 0 & \frac{3-2\nu}{4}\end{bmatrix}

.. math::
    \left[F_x\right] = I_{xy}\begin{bmatrix}\frac{3-2\nu}{4} & 0 & \frac{1+2\nu}{4} \\ 0 & \frac{1-2\nu}{4} & 0\end{bmatrix} - I_{yy}\begin{bmatrix}0 & \frac{1-2\nu}{4} & 0 \\ \frac{1+2\nu}{4} & 0 & \frac{3-2\nu}{4}\end{bmatrix}

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
    \mathbf{\sigma} = \begin{bmatrix}0 & 0 & \tau_{xz} \\ 0 & 0 & \tau_{yz} \\ \tau_{xz} & \tau_{yz} & \sigma_{zz}\end{bmatrix}

While the strain is given by

.. math::
    \mathbf{\varepsilon} = \begin{bmatrix}\varepsilon_{xx} & 0 & \varepsilon_{xz} \\ 0 & \varepsilon_{yy} & \varepsilon_{yz} \\ \varepsilon_{xz} & \varepsilon_{yz} & \varepsilon_{zz} \end{bmatrix}

And the elasticity law is 

.. math::
    \mathbf{\sigma} = \lambda \cdot \text{trace}\left(\mathbf{\varepsilon}\right) \cdot \mathbf{I} + 2\mu \cdot \mathbf{\varepsilon}
.. math::
    \mathbf{\varepsilon} = \dfrac{1}{2\mu} \cdot \mathbf{\sigma} - \dfrac{\lambda}{2\mu\left(3\lambda + 2\mu\right)} \cdot \text{trace}\left(\mathbf{\sigma}\right) \cdot \mathbf{I}

The literature shows that the displacement field :math:`\mathbf{u}` is given by

.. math::
    \mathbf{u} = 

.. _integrals:

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