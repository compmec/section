.. _theory:

============
Introduction
============

The main source of theories, hypothesis, equations are from both books of Pikley and BEM.

* Pikley contains the base theory for section
* BEM's book 
* 

While Pikley's book uses Finite Element Method (FEM for short), their equations must be adapted to boundary element method.

These transformations happens mainly from Green's Theorem and its identities, which are developed in this text. Since the main sources are from these two books, the repeted citations are ommited.

Notation and utilities
----------------------

1. The sections are defined in the :math:`xy` plane, while the perpendicular direction is :math:`z`.

2. The bold symbol :math:`\mathbf{u}` represents a bidimensional vector, sometimes represented as a column-vector

.. math::
    \mathbf{u} = \left(u_{x}, \ u_{y}\right) = \begin{bmatrix}u_{x} \\ u_{y}\end{bmatrix}

3. Inner product between two bidimensional vectors results in a scalar

.. math::
    \langle \mathbf{u}, \ \mathbf{v} \rangle = u_x v_x + u_y v_y = \mathbf{u}^{T} \cdot \mathbf{v}

4. Cross product between two bidimensional vectors results in a scalar

.. math::
    \mathbf{u} \times \mathbf{v} = u_{x} v_{y} - u_{y}v_{x}

5. The tensorial product results in a matrix

.. math::
    \mathbf{u} \otimes \mathbf{v} = \mathbf{u} \cdot \mathbf{v}^{T} = \begin{bmatrix}u_{x} \\ u_{y}\end{bmatrix} \begin{bmatrix}v_{x} & v_{y} \end{bmatrix} = \begin{bmatrix}u_{x}v_{x} & u_{x}v_{y} \\ u_{y}v_{x} & u_{y}v_{y}\end{bmatrix}

5. The modulus of a vector is denoted by :math:`\|\mathbf{u}\|`

.. math::
    \|\mathbf{u}\| = \sqrt{\langle \mathbf{u}, \ \mathbf{u}\rangle} = \sqrt{u_x^2 +u_y^2}

6. A point in the plane is :math:`\mathbf{p}`

.. math::
    \mathbf{p} = \left(x, \ y\right)

7. Parametrization of a curve is made by :math:`\mathbf{p}(t)`, with variable :math:`t`, which sometimes is ommited

.. math::
    \mathbf{p}(t) = \left(x(t), \ y(t)\right)

8. The derivative of :math:`\mathbf{p}(t)`

.. math::
    \dfrac{d}{dt} \left(x(t), \ y(t)\right) = \dfrac{d}{dt} \mathbf{p}(t) = \mathbf{p}'(t) = \left(x'(t), \ y'(t)\right)

9. The tangent vector is :math:`\mathbf{t}`, the normal vector is :math:`\mathbf{n}` and they have modulus 1 by definition

.. math::
    \mathbf{t} = \dfrac{\mathbf{p'}}{\|\mathbf{p}'\|} = \dfrac{\left(x', \ y'\right)}{\|\mathbf{p}'\|}

.. math::
    \mathbf{n} = \dfrac{\left(y', \ -x'\right)}{\|\mathbf{p}'\|}



11. The gradient of an scalar is a bidimensional vector

.. math::
    \nabla u = \left(\dfrac{\partial u}{\partial x}, \ \dfrac{\partial u}{\partial y} \right) = \begin{bmatrix}\dfrac{\partial u}{\partial x} \\ \dfrac{\partial u}{\partial y} \end{bmatrix}

12. The laplacian of an scalar is a scalar

.. math::
    \nabla^2 u = \dfrac{\partial^2 u}{\partial x^2}+\dfrac{\partial^2 u}{\partial y^2}

13. The bold symbol :math:`\mathbf{u}` also represents :math:`n`-dimensional vector

.. math::
    \mathbf{u} = \left(u_{0}, \ u_{1}, \ \cdots, \ u_{n-1}\right)


14. The special symbol :math:`\mathbb{M}` represents a matrix.

15. Matrix-vector multiplication is denoted as 

.. math::
    \mathbb{M} \cdot \mathbf{u} = \mathbf{F} \Longleftrightarrow \sum_{j} M_{ij} \cdot u_{j} = F_{i} \ \ \ \ \ \ \ \forall i


.. _boundary_element_method:

=======================
Boundary Element Method
=======================

Introduction
------------

The Boundary Element Method (BEM for short) is used to find :math:`u` numerically

.. math:: 
    :label: eq_laplace

    \nabla^2 u = 0

The BEM transforms :eq:`eq_laplace` into a boundary version :eq:`eq_bem`

.. math::
    :label: eq_bem

    \alpha\left(\mathbf{s}\right) \cdot u\left(\mathbf{s}\right) = \int_{\Gamma} u \cdot \dfrac{\partial v}{\partial n} \ d\Gamma - \int_{\Gamma} \dfrac{\partial u}{\partial n}  \cdot v \ d\Gamma

Which :math:`\mathbf{s}` is the source point of the Green function :math:`v` and :math:`\alpha(\mathbf{s})` is the angle at the point :math:`\mathbf{s}`.

.. math::
    :label: eq_source

    v(\mathbf{p}, \ \mathbf{s}) = \ln r = \ln \|\mathbf{r}\| = \ln \|\mathbf{p} - \mathbf{s}\|

Since all the PDEs used in this package have only Neumann's boundary conditions, the values of :math:`\dfrac{\partial u}{\partial n}` are known and the objective is finding all the values of :math:`u` at the boundary.

Once :math:`u` and :math:`\dfrac{\partial u}{\partial n}` are known at the boundary, it's possible to compute :math:`u(x, y)` and its derivatives at any point inside by using :eq:`eq_bem`.


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

Applying for :math:`n` different source points :math:`\mathbf{s}_i` at boundary, we get the matrices :math:`\mathbb{A}`, :math:`\mathbb{M}` and :math:`\mathbf{F}` such

.. math::
    :label: eq_linear_system

    \left(\mathbb{M}-\mathbb{A}\right) \cdot \mathbf{U} = \mathbb{K} \cdot \mathbf{U} = \mathbf{F}

Finding the values of :math:`\mathbf{U}` means solving the linear system :eq:`eq_linear_system`


Matrix :math:`\mathbb{A}`
^^^^^^^^^^^^^^^^^^^^^^^^^

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

Matrix :math:`\mathbb{M}`
^^^^^^^^^^^^^^^^^^^^^^^^^

We use

.. math::
    \dfrac{\partial v}{\partial n} ds = \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle}

to write

.. math::
    M_{ij} = \int_{t_{min}}^{t_{max}} \varphi_{j}(t) \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle\mathbf{r}, \ \mathbf{r}\right\rangle} \ dt

Vector :math:`\mathbf{F}` for warping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the warping function

.. math::
    \dfrac{\partial u}{\partial n} = \mathbf{n} \times \mathbf{p} = \dfrac{\langle \mathbf{p}, \ \mathbf{p}'\rangle}{\|\mathbf{p}'\|}

.. math::
    F_i = \int_{t_{min}}^{t_{max}} \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \cdot \ln \|\mathbf{r}_i\| \ dt


Vector :math:`\mathbf{F}` for shear
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vector :math:`\mathbf{F}` for shear are in fact 2 vectors.

We compute the value of :math:`\mathbb{X}`, which is a :math:`(n \times 6)` matrix

.. math::
    \mathbb{X}_{i} = \int_{t_{min}}^{t_{max}} \ln r \cdot \begin{bmatrix}x^2 \cdot x' \\ 2xy \cdot x' \\ y^2 \cdot x' \\ x^2 \cdot y' \\ 2xy \cdot y' \\ y^2 \cdot y' \end{bmatrix}

With this matrix, we compute the vector :math:`\mathbf{F}` and it's better explained in :ref:`shear_center`.


Evaluating matrices
^^^^^^^^^^^^^^^^^^^

The matrices highly depend on the geometry and the basis functions :math:`\varphi`.

To compute the coefficients :math:`M_{ij}` and :math:`F_{i}`, it's used numerical integration, like Gaussian-Quadrature.
Unfortunatelly, when :math:`r = 0` at some point, the integrants are singular and special techniques are used.

The main idea to compute them is decompose the integral in intervals and use

* **Outside integration**: uses :ref:`regular_integrals` for elements which :math:`r\ne 0` for all points

* **Inside integration**: uses :ref:`singular_integrals` for elements which :math:`r=0` at any point

For polygonal domains the **Inside integration** is not required cause it can be done analiticaly. But for higher degrees, it's indeed necessary

.. _constraint_solution:

Constraint solution
^^^^^^^^^^^^^^^^^^^

Although the matrix :math:`\mathbb{K}=\mathbb{M}-\mathbb{A}` is not singular, all the PDEs have Neumann's boundary conditions and has no unique solution.
If :math:`u^{\star}` is found as solution, then :math:`\left(u^{\star} + \text{const}\right)` also is a solution.

Although both functions give the same properties cause it envolves only the derivatives of :math:`u`, we restrict the solution by solving the system with Lagrange Multiplier.

.. math::
    \begin{bmatrix}K & \mathbf{C} \\ \mathbf{C}^T & 0\end{bmatrix} \begin{bmatrix}\mathbf{U} \\ \lambda \end{bmatrix} = \begin{bmatrix}\mathbf{F} \\ 0\end{bmatrix}

Which vector :math:`\mathbf{C}` is a vector of ones.

The determination exacly of the constant depends on the problem and are better treated in :ref:`torsion_center` and :ref:`shear_center`.


.. _bem_polygonal_domain:

Polygonal domain
----------------

For polygonal domains, when the basis functions :math:`\phi(t)` are piecewise linear, some computations becomes easier. Let's say the parametric space :math:`t` is divided by the knots :math:`t_0`, :math:`t_1`, :math:`\cdots`, :math:`t_{m-1}`, :math:`t_m`, which correspond to the vertices

For an arbitrary interval :math:`\left[t_k, \ t_{k+1}\right]`, :math:`\mathbf{p}(t)` is described as

.. math::
    \mathbf{p}(t) = \mathbf{P}_{k} + \tau \cdot \mathbf{V}_k
    
.. math::
    \mathbf{V}_k = \mathbf{P}_{k+1} - \mathbf{P}_{k}

.. math::
    \tau = \dfrac{t - t_{k}}{t_{k+1} - t_{k}} \in \left[0, \ 1\right]

Since the source point :math:`\mathbf{s}_i = \mathbf{p}(t_i)`,

* If :math:`t_i \in \left[t_{k}, \ t_{k+1}\right]` then

    .. math::
        \mathbf{r}(t) = \left(\tau-\tau_i\right) \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_{k}\right)

    .. math::
        \tau_i = \dfrac{t_i - t_{k}}{t_{k+1} - t_{k}}\in \left[0, \ 1\right]

* Else

    .. math::
        \mathbf{r}(t) = \left(\mathbf{P}_{k}-\mathbf{s}_i\right) + \tau \cdot \left(\mathbf{P}_{k+1} - \mathbf{P}_{k}\right)


Matrix :math:`\mathbb{A}`
^^^^^^^^^^^^^^^^^^^^^^^^^

If the source point :math:`\mathbf{s}_i` lies in the middle of the segment

.. math::
    \alpha(\mathbf{s}_i) = \pi

If the source point :math:`s_i` lies in the vertex :math:`P_{k}` then

.. math::
    \mathbf{v}_0 = \mathbf{P}_{k-1}-\mathbf{P}_{k}
.. math::
    \mathbf{v}_1 = \mathbf{P}_{k+1}-\mathbf{P}_{k}
.. math::
    \alpha = \arg\left(\langle\mathbf{v}_0, \ \mathbf{v}_1 \rangle + i \cdot \left(\mathbf{v}_0 \times \mathbf{v}_1\right)\right)


Matrix :math:`\mathbb{M}`
^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
    M_{ij} = \sum_{k=0}^{m-1} \int_{t_{k}}^{t_{k+1}} \varphi_{j} \cdot \dfrac{\mathbf{r} \times \mathbf{p}'}{\left\langle \mathbf{r}, \mathbf{r}\right\rangle} \ dt

* If :math:`t_i \notin \left[t_k, \ t_{k+1}\right]`, then the evaluation is made by :ref:`regular_integrals`

* If :math:`t_i \in \left[t_k, \ t_{k+1}\right]`

    .. math::
        \mathbf{V}_k = \mathbf{P}_{k+1} - \mathbf{P}_k
    .. math::
        \mathbf{p(t)} = \mathbf{P}_k + \tau \cdot \mathbf{V}_{k} 
    .. math::
        \mathbf{r(t)} = \left(\tau-\tau_i\right) \cdot \mathbf{V}_{k} 
    .. math::
        \mathbf{r} \times \mathbf{p}' = 0 

    Therefore, we can ignore the integration over the interval :math:`\left[t_k, \ t_{k+1}\right]`


Vector :math:`\mathbf{F}` for warping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For warping function, the expression :math:`F_i` is written as

.. math::
    \dfrac{\partial u}{\partial n} = \dfrac{\left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle}{\|\mathbf{p}'\|}
    
.. math::
    F_{i} = \sum_{k=0}^{m-1} \int_{0}^{1} \left(\alpha_k + \tau \cdot \beta_k \right) \ln\|\mathbf{r}\| \ d\tau

With :math:`\mathbf{P}_k` begin the :math:`k`-vertex and

.. math::
    \mathbf{V}_k = \mathbf{P}_{k+1} - \mathbf{P}_k
.. math::
    \alpha_k = \left\langle \mathbf{P}_k, \ \mathbf{V}_k\right\rangle
.. math::
    \beta_k = \left\langle \mathbf{V}_k, \ \mathbf{V}_k\right\rangle
    
* If  :math:`t_i \notin \left[t_k, \ t_{k+1}\right]`, :ref:`regular_integrals` are used

* If :math:`t_i \in \left[t_k, \ t_{k+1}\right]`, then
    .. math::
        \tau_i = \dfrac{t_i-t_k}{t_{k+1}-t_{k}} \in \left[0, \ 1\right]
    .. math::
        \mathbf{V}_k = \mathbf{P}_{k+1} - \mathbf{P}_k
    .. math::
        \mathbf{p(t)} = \mathbf{P}_k + \tau \cdot \mathbf{V}_{k} 
    .. math::
        \mathbf{r(t)} = \left(\tau-\tau_i\right) \cdot \mathbf{V}_{k}
    .. math::
        F_{ik} = & \int_{0}^{1} \left(\alpha_k + \tau \beta_k \right) \ln\|\left(\tau-\tau_i\right) \cdot \mathbf{V}_k\| \ d\tau \\
            = & \left(\alpha_{k} + \dfrac{1}{2}\beta_{k}\right) \cdot \dfrac{1}{2}\ln \beta_k \\
                & + \alpha_{k} \int_{0}^{1} \ln |\tau-\tau_i| dz \\
                & + \beta_k \int_{0}^{1} \tau \cdot \ln |\tau-\tau_i| \ dz 

    These two log integrals are computed analiticaly, the expressions are complicated (`here <https://www.wolframalpha.com/input?i=int_%7B0%7D%5E%7B1%7D+ln%28abs%28x-x_0%29%29+dx%3B+0+%3C%3D+x_0+%3C%3D+1>`_ and `here <https://www.wolframalpha.com/input?i=int_%7B0%7D%5E%7B1%7D+x*ln%28abs%28x-x_0%29%29+dx%3B+0+%3C%3D+x_0+%3C%3D+1>`_) and depends on the value of :math:`\tau_i`. Bellow you find a table with some values

    .. list-table:: Values of logarithm integrals
        :widths: 20 40 40
        :header-rows: 1
        :align: center

        * - :math:`\tau_i`
          - :math:`\int_0^1 \ln|\tau-\tau_i| dz`
          - :math:`\int_0^1 \tau\ln|\tau-\tau_i| dz`
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


Vector :math:`\mathbf{F}` for shear
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The evaluation of this integral is made by computing :math:`\mathbb{X}_i`

.. math::
    \mathbb{X}_{i} = \int_{t_{min}}^{t_{max}} \ln r \cdot \begin{bmatrix}x^2 \cdot x' \\ 2xy \cdot x' \\ y^2 \cdot x' \\ x^2 \cdot y' \\ 2xy \cdot y' \\ y^2 \cdot y' \end{bmatrix} \ dt


* For :math:`t_i \notin \left[t_k, \ t_{k+1}\right]`, uses :ref:`regular_integrals` to compute

* For :math:`t_i \in \left[t_k, \ t_{k+1}\right]` then

    .. math::
        \tau_i = \dfrac{t_i-t_k}{t_{k+1}-t_{k}}
    .. math::
        \mathbf{V}_k = \mathbf{P}_{k+1}-\mathbf{P}_{k}
    .. math::
        \mathbf{p}(t) = \mathbf{P}_{k}+\tau \cdot \mathbf{V}_{k}
    .. math::
        \mathbf{r}(t) = (\tau - \tau_i) \cdot \mathbf{V}_{k}
    .. math::
        \ln \|\mathbf{r}\| = \dfrac{1}{2}\ln \beta_k + \ln |\tau - \tau_i|

    Breaking into components:

    .. math::
        x(t) = x_{k} + \tau \Delta x_{k}
    .. math::
        y(t) = y_{k} + \tau \Delta y_{k}

    and let 

    

    The integrals become

    .. math::
        \mathbb{X}_{ik} = \dfrac{1}{2}\ln \beta_k \int_{0}^{1} \begin{bmatrix}\Delta x_{k} \cdot x^2 \\ \Delta x_{k} \cdot 2xy \\ \Delta x_{k} \cdot y^2 \\ \Delta y_{k} \cdot x^2 \\ \Delta y_{k} \cdot 2xy \\ \Delta y_{k} \cdot y^2\end{bmatrix} \ d\tau + \int_{0}^{1} \ln |\tau - \tau_i| \begin{bmatrix}\Delta x_{k} \cdot x^2 \\ \Delta x_{k} \cdot 2xy \\ \Delta x_{k} \cdot y^2 \\ \Delta y_{k} \cdot x^2 \\ \Delta y_{k} \cdot 2xy \\ \Delta y_{k} \cdot y^2\end{bmatrix} \ d\tau
    
    The left part is

    .. math::
        \mathbb{X}_{ik0} = \int_{0}^{1} \begin{bmatrix}x^2 \\ 2xy \\ y^2 \end{bmatrix} \ d\tau = \begin{bmatrix}x_{k}^2+x_kx_{k+1}+x_{k+1}^{2} \\ 2x_{k}y_{k} + x_{k}y_{k+1}+x_{k+1}y_{k}+2x_{k+1}y_{k+1} \\ y_{k}^2+y_ky_{k+1}+y_{k+1}^{2} \end{bmatrix}

    The right part is used logarithm integration.
    





==================================
Cross-section geometric properties
==================================


.. _cross_sectional_area:

Cross-section area
------------------

.. math::
    A = \int_{\Omega} \ dx \ dy


.. _first_moment_area:

First moment of area
--------------------

.. math::
    Q_y = \int_{\Omega} x \ dx \ dy
.. math::
    Q_x = \int_{\Omega} y \ dx \ dy



.. _geometric_center:

Geometric center
----------------

.. math::
    x_{gc} = \dfrac{Q_y}{A}
.. math::
    y_{gc} = \dfrac{Q_x}{A}

We denote the geometric centroid by :math:`\boldsymbol{G}`

.. math::
    \boldsymbol{G} = \left(x_{gc}, \ y_{gc}\right)


.. _second_moment_area:

Global Second Moment of Area
-----------------------------

The global second moment of inertia are

.. math::
    I_{yy} = \int_{\Omega} x^2 \ dx \ dy
.. math::
    I_{xy} = \int_{\Omega} xy \ dx \ dy
.. math::
    I_{xx} = \int_{\Omega} y^2 \ dx \ dy



Local Second Moment of Area
-----------------------------

The local second moment of inertia are computed with respect to the geometric center

.. math::
    I_{\overline{yy}} = \int_{\Omega} (x-x_{gc})^2 \ dx \ dy = I_{yy} - \dfrac{Q_{y}^2}{A}
.. math::
    I_{\overline{xy}} = \int_{\Omega} (x-x_{gc})(y-y_{gc}) \ dx \ dy= I_{xy} - \dfrac{Q_{x}Q_{y}}{A}
.. math::
    I_{\overline{xx}} = \int_{\Omega} (y-y_{gc})^2 \ dx \ dy= I_{xx} - \dfrac{Q_{y}^2}{A}

    


.. _radius_gyration:

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


.. _bending_center:

Bending Center
--------------

The bending center is the intersection of the two neutral lines.
It's computed by



===============================
Torsion Properties
===============================

.. _warping_function:

Warping Function
----------------

From Saint-venant theory, the warping function :math:`\omega(x, \ y)` is fundamental to compute torsion properties.
From :math:`\omega`, it's possible to find the :ref:`torsion_constant`, :ref:`torsion_center` and stresses/strains due to :ref:`torsion_moment`.

.. math::
    \nabla^2 \omega = 0

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \mathbf{n} \times \mathbf{p}

With :math:`\mathbf{p} = (x, \ y)` begin a point on the boundary.

The boundary condition can be rewriten as

.. math::
    \left\langle \nabla \omega, \ \mathbf{n}\right\rangle = \dfrac{\langle \mathbf{p}', \ \mathbf{p} \rangle}{\|\mathbf{p}'\|} 

We solve this PDE by using the :ref:`boundary_element_method`, by solving the linear system

.. math::
    \begin{bmatrix}K & \mathbf{C} \\ \mathbf{C}^{T} & 0\end{bmatrix}\begin{bmatrix}\mathbf{W} \\ \lambda\end{bmatrix} = \begin{bmatrix}\mathbf{F} \\ 0\end{bmatrix}

With :math:`\omega = \langle \varphi, \ \mathbf{W}\rangle` at the boundary

.. _torsion_constant:

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
    \mathbb{J}_{\omega} = \left\langle \mathbf{W}, \ \int_{t_{min}}^{t_{max}} \varphi_j \cdot \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \ dt \right\rangle


.. _torsion_center:

Torsion center
---------------

As described in :ref:`constraint_solution`, we solve a Neumann's problem.
If :math:`\omega^{\star}` is a solution, then :math:`\omega^{*} = \omega^{\star} + c_0` is also a solution.

This constant is arbitrary and don't change the torsion properties or the stresses due to torsion moment. Let :math:`\omega^{*}` be the solution of  
Choosing arbitrarily the values of :math:`x_0`, :math:`y_0` and :math:`c_0` doesn't change the torsion properties or the stresses due to torsion, it can be understood as a *rigid body rotation in the plane of cross-section and a displacement parallel to the axis of the bar* (from BOOK BEM).

The quantities :math:`x_0`, :math:`y_0` and :math:`c_0` can be obtained by minimizing the strain energy produced by axial normal warping stresses, which are ignored by Saint-Venant's theory.
Doing so, leads to the linear system

.. math::
    \left(\int_{\Omega} \begin{bmatrix}1 & x & y \\ x & x^2 & xy \\ y & xy & y^2 \end{bmatrix} \ d\Omega\right) \begin{bmatrix}c_0 \\ y_0 \\ -x_0\end{bmatrix} = \int_{\Omega} \omega\begin{bmatrix}1 \\ x \\ y\end{bmatrix} \ d\Omega

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
        g_{1}(x, \ y) = 1 \Longrightarrow u_{1} = \frac{1}{4}(x^2+y^2)
    
    .. math::
        g_{x}(x, \ y) = x \Longrightarrow u_{x} = \frac{x^3}{6}
     
    .. math::
        g_{y}(x, \ y) = y \Longrightarrow u_{y} = \frac{y^3}{6}
        
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



.. _shear_properties:

================
Shear properties
================


Introduction
------------

Functions :math:`\Psi` and :math:`\Phi` are used:

.. math::
    \begin{bmatrix} \nabla^2 \Psi \\ \nabla^2 \Phi \end{bmatrix} = 
    2\begin{bmatrix} -I_{\overline{xx}} & I_{\overline{xy}} \\ I_{\overline{xy}} & -I_{\overline{yy}} \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix}


And boundary conditions

.. math::
    \begin{bmatrix}\nabla \Psi \\ \nabla \Phi\end{bmatrix} \cdot \mathbf{n} = \mathbb{H} \cdot \mathbf{n}
.. math::
    \mathbb{H} = \dfrac{\nu}{2}\left((x^2-y^2)\cdot\begin{bmatrix}I_{xx} & I_{xy} \\ -I_{xy} & -I_{yy}\end{bmatrix} + 2xy \cdot \begin{bmatrix}-I_{xy} & I_{xx} \\ I_{yy} & -I_{xy}\end{bmatrix}\right)

Both equations are in fact Poisson equations. We find them by using the :ref:`boundary_element_method`, as described bellow

.. dropdown:: BEM formulation for Poisson's equation

    To find :math:`\Psi` and :math:`\Phi`, we solve an equivalent problem by transforming the Poisson PDE into a Laplace, which is easier to solve by :ref:`boundary_element_method`.
    
    Take :math:`\Psi = \Psi^{\star} + \Psi_{0}` and :math:`\Phi = \Phi^{\star} + \Phi_{0}`, the following Laplace PDE is obtained

    .. math::
        \begin{bmatrix} \nabla^2 \Psi \\ \nabla^2 \Phi \end{bmatrix} = \begin{bmatrix}0 \\ 0 \end{bmatrix}
    
    .. math::
        \begin{bmatrix}\nabla \Psi^{\star} \\ \nabla \Phi^{\star} \end{bmatrix} \cdot \mathbf{n} = \mathbb{H} \cdot \mathbf{n} - \begin{bmatrix}\nabla \Psi_{0} \\ \nabla \Phi_{0}\end{bmatrix} \cdot \mathbf{n}

    With

    .. math::
        \begin{bmatrix}\Psi_0 \\ \Phi_0 \end{bmatrix} = \dfrac{1}{4}\left(x^2+y^2\right)\begin{bmatrix} -I_{xx} & I_{xy} \\ I_{xy} & -I_{yy} \end{bmatrix} \begin{bmatrix} x \\ y \end{bmatrix}
    
    The *stiffness* matrix :math:`\mathbb{M}` is the same for the :ref:`warping_function` and the *force* vector :math:`\mathbf{F}` is computed bellow:

    .. math::
        \mathbb{H} \cdot \mathbf{n} = \dfrac{\nu}{2}\begin{bmatrix}-I_{xx} & -I_{xy} \\ I_{xy} & I_{yy}\end{bmatrix}\begin{bmatrix}2xy & x^2-y^2 \\ x^2-y^2 & 2xy\end{bmatrix} \cdot \begin{bmatrix}x' \\ y'\end{bmatrix}\dfrac{1}{\|\mathbf{p}'\|}

    
    .. math::
        \begin{bmatrix}\nabla \Psi_{0} \\ \nabla \Phi_{0}\end{bmatrix} \cdot \mathbf{n} = \dfrac{1}{4}\begin{bmatrix}-I_{xx} & I_{xy} \\ I_{xy} & -I_{yy}\end{bmatrix}\begin{bmatrix}2xy & 3x^2+y^2 \\ x^2+3y^2 & 2xy\end{bmatrix}\cdot \begin{bmatrix}-x' \\ y'\end{bmatrix}\dfrac{1}{\|\mathbf{p}'\|}

    Let

    .. math::
        \mathbb{X}_i = \int \ln r_{i} \ \begin{bmatrix}x^2 \cdot x' & x^2  \cdot y' \\ 2xy \cdot x' & 2xy \cdot y' \\ y^2 \cdot x' & y^2 \cdot  y' \end{bmatrix}\ dt

    

    .. math::
        F_i = \int \ln r \cdot 

.. _shear_center:

Shear center
------------

The shear center is defined by :math:`\boldsymbol{S} = \left(x_{sc}, \ y_{sc}\right)`

.. math::
    \boldsymbol{S} = \dfrac{\nu}{2\Delta}\begin{bmatrix}I_{yy} & I_{xy} \\ I_{xy} & I_{xx}\end{bmatrix}\begin{bmatrix}I_{yyy}+I_{xxy} \\ I_{xyy}+I_{xxx} \end{bmatrix} - \dfrac{1}{\Delta}\int \begin{bmatrix}\Psi \\ \phi\end{bmatrix} \left\langle \mathbf{p}, \ \mathbf{p}'\right\rangle \ dt

Which

.. math::
    \begin{bmatrix}I_{yyy} \\ I_{xyy} \\ I_{xxy} \\ I_{xxx} \end{bmatrix} = \int \begin{bmatrix}x^3 \\ x^2y \\ xy^2 \\ y^3 \end{bmatrix} \ dx \ dy

.. math::
    \Delta = 2(1+\nu)(I_{xx}I_{yy}-I_{xy})

.. _stress_and_strain:

=================
Stress and Strain
=================


Introduction
------------

The stress :math:`\boldsymbol{\sigma}` and strain :math:`\boldsymbol{\varepsilon}` are one of the fundamental quantities to evaluate. They arrive from 4 different phenomenums:

* :ref:`axial_force` (1 quantity: :math:`\mathrm{F}_{z}`)
* :ref:`bending_moments` (2 quantities: :math:`\mathrm{M}_{x}` and :math:`\mathrm{M}_{y}`) 
* :ref:`torsion_moment` (1 quantity: :math:`\mathrm{M}_{z}`)
* :ref:`shear_forces` (2 quantities: :math:`\mathrm{F}_{x}` and :math:`\mathrm{F}_{y}`) 

Here we develop expressions to compute stress and strain for any point :math:`\mathbf{p}` inside the section.
The stress and strain tensor in a beam are given by

.. math::
    \boldsymbol{\sigma} = \begin{bmatrix}0 & 0 & \sigma_{xz} \\ 0 & 0 & \sigma_{yz} \\ \sigma_{xz} & \sigma_{yz} & \sigma_{zz}\end{bmatrix} \ \ \ \ \ \ \ \ \ \boldsymbol{\varepsilon} = \begin{bmatrix}\varepsilon_{xx} & 0 & \varepsilon_{xz} \\ 0 & \varepsilon_{yy} & \varepsilon_{yz} \\ \varepsilon_{xz} & \varepsilon_{yz} & \varepsilon_{zz} \end{bmatrix}

The elasticity law relates both tensors 

.. math::
    \boldsymbol{\sigma} = \lambda \cdot \text{trace}\left(\boldsymbol{\varepsilon}\right) \cdot \mathbf{I} + 2\mu \cdot \boldsymbol{\varepsilon}
    
.. math::
    \boldsymbol{\varepsilon} & = \dfrac{1}{2\mu} \cdot \boldsymbol{\sigma} - \dfrac{\lambda}{2\mu\left(3\lambda +2\mu\right)} \cdot \text{trace}\left(\boldsymbol{\sigma}\right) \cdot \mathbf{I} \\ & = \dfrac{1+\nu}{E} \cdot \boldsymbol{\sigma} - \dfrac{\nu}{E} \cdot \text{trace}\left(\boldsymbol{\sigma}\right) \cdot \mathbf{I}

With :math:`\lambda` and :math:`\mu` being `Lamé Parameters <https://en.wikipedia.org/wiki/Lam%C3%A9_parameters>`_, :math:`E` beging Young Modulus and :math:`\nu` the Poisson's coefficient.

.. math::
    \lambda = \dfrac{E\nu}{(1+\nu)(1-2\nu)} \ \ \ \ \ \ \ \ \ \ \ \mu = \dfrac{E}{2(1+\nu)}

.. math::
    E = \dfrac{\mu\left(3\lambda+2\mu\right)}{\lambda+\mu} \ \ \ \ \ \ \ \ \ \ \ \nu = \dfrac{\lambda}{2(\lambda+\mu)}


.. _axial_force:

Axial Force
------------

The axial force only leads to axial stress.
Meaning, in pure axial force case, the stress tensor and strain are given by

.. math::
    \boldsymbol{\varepsilon} =  \begin{bmatrix}\varepsilon_{xx} & 0 & 0 \\ 0 & \varepsilon_{yy} & 0 \\ 0 & 0 & \varepsilon_{zz}\end{bmatrix} \ \ \ \ \ \ \ \ \ \ \ \sigma = \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & \sigma_{zz}\end{bmatrix}

The axial stress is constant when an axial force :math:`\mathrm{F}_{z}` is given by

.. math::
    \sigma_{zz} = \dfrac{\mathrm{F}_{z}}{A}

Where :math:`A` is the :ref:`cross_sectional_area`.

Hence, the strain is given by elasticity law:

.. math::
    \varepsilon_{xx} = \varepsilon_{yy} = \left(\dfrac{-\lambda}{2\mu(3\lambda+2\mu)}\right) \cdot \dfrac{\mathrm{F}_{z}}{A} = -\nu \cdot \dfrac{\mathrm{F}_{z}}{EA}
.. math::
    \varepsilon_{zz} = \dfrac{1}{3\lambda+2\mu}\left(1+\dfrac{\lambda}{\mu}\right) \cdot \dfrac{\mathrm{F}_{z}}{A} = \dfrac{\mathrm{F}_{z}}{EA}

.. math::
    \boldsymbol{\varepsilon} = \dfrac{\mathrm{F}_{z}}{EA}\begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}

.. _bending_moments:

Bending Moments
---------------

The bending moments :math:`\mathrm{M}_{x}` and :math:`\mathrm{M}_{y}` causes only axial stresses.
The tensors are 

.. math::
    \boldsymbol{\varepsilon} =  \begin{bmatrix}\varepsilon_{xx} & 0 & 0 \\ 0 & \varepsilon_{yy} & 0 \\ 0 & 0 & \varepsilon_{zz}\end{bmatrix} \ \ \ \ \ \ \ \ \ \ \ \sigma = \begin{bmatrix}0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & \sigma_{zz}\end{bmatrix}

The expression of :math:`\sigma_{zz}` depends on the position of the point :math:`\mathbf{p}` in the section. 
In the :ref:`bending_center` :math:`\boldsymbol{B} = \left(x_{bc}, \ y_{bc}\right)` the stress and the strain are zero while they increase/decrease depending on the distance to the bending center.

Let :math:`\bar{x}=x-x_{bc}` and :math:`\bar{y}=y-y_{bc}`, the function :math:`\sigma_{zz}(x, \ y)` satisfy

.. math::
    \int_{\Omega} \sigma_{zz} \cdot \begin{bmatrix}\bar{y} \\ -\bar{x}\end{bmatrix} \ d\Omega = \begin{bmatrix}M_{x} \\ M_{y}\end{bmatrix}

Add the hypothesis that :math:`\sigma_{zz}` is linear with respect to :math:`x` and :math:`y`, then 

.. math::
    \sigma_{zz}(x, \ y) & = \dfrac{1}{\det \left(\mathbb{I}_{b}\right)} \begin{bmatrix}\bar{y} & \bar{x}\end{bmatrix} \left[\mathbb{I}_{b}\right] \begin{bmatrix}M_{x} \\ M_{y}\end{bmatrix} \\
     & = -\left(\dfrac{I_{\overline{xy}}\mathrm{M}_{x} + I_{\overline{xx}}\mathrm{M}_{y}}{I_{\overline{xx}}I_{\overline{yy}}-I_{\overline{xy}}^2}\right) \cdot \bar{x} + \left(\dfrac{I_{\overline{yy}}\mathrm{M}_{x} + I_{\overline{xy}}\mathrm{M}_{y}}{I_{\overline{xx}}I_{\overline{yy}}-I_{\overline{xy}}^2}\right) \cdot \bar{y}

With constants :ref:`second_moment_area`

.. math::
    \left[\mathbb{I}_{b}\right] = \begin{bmatrix}I_{\overline{yy}} & I_{\overline{xy}} \\ I_{\overline{xy}} & I_{\overline{xx}}\end{bmatrix} = \int_{\Omega}\begin{bmatrix}\left(x-x_{bc}\right)^2 & \left(x-x_{bc}\right)\left(y-y_{bc}\right) \\ \left(x-x_{bc}\right)\left(y-y_{bc}\right) & \left(y-y_{bc}\right)^2\end{bmatrix} d\Omega

The neutral line is the set of pairs :math:`(x, \ y)` such :math:`\sigma_{zz}(x, \ y) = 0`.
That means the neutral line is the line that pass thought :math:`\boldsymbol{B}` and it's parallel to the vector :math:`\left[\mathbb{I}_{b}\right] \cdot \left(\mathrm{M}_{x}, \ \mathrm{M}_{y}\right)`

It's possible to obtain strain values from elasticity law:

.. math::
    \varepsilon_{xx} = \varepsilon_{yy} = \left(\dfrac{-\lambda}{2\mu(3\lambda+2\mu)}\right) \cdot \sigma_{zz} = -\nu \cdot \dfrac{\sigma_{zz}}{E}
.. math::
    \varepsilon_{zz} = \dfrac{1}{3\lambda+2\mu}\left(1+\dfrac{\lambda}{\mu}\right) \cdot \sigma_{zz} = \dfrac{\sigma_{zz}}{E}

.. math::
    \boldsymbol{\varepsilon} = \dfrac{\sigma_{zz}}{E} \cdot \begin{bmatrix}-\nu & 0 & 0 \\ 0 & -\nu & 0 \\ 0 & 0 & 1\end{bmatrix}



.. _torsion_moment:

Torsion Moment
--------------

The torsion moment :math:`\mathrm{M}_{z}` causes only shear stresses.
The tensors are 

.. math::
    \boldsymbol{\varepsilon} = \begin{bmatrix}0 & 0 & \varepsilon_{xz} \\ 0 & 0 & \varepsilon_{yz} \\ \varepsilon_{xz} & \varepsilon_{yz} & 0\end{bmatrix} \ \ \ \ \ \ \ \ \ \ \ \boldsymbol{\sigma} = \begin{bmatrix}0 & 0 & \sigma_{xz} \\ 0 & 0 & \sigma_{yz} \\ \sigma_{xz} & \sigma_{xz} & 0\end{bmatrix}

The :ref:`warping_function` :math:`\omega` is used to compute them

.. math::
    \sigma_{xz}(x, \ y) = \dfrac{\mathrm{M}_{z}}{J} \cdot \left(\dfrac{\partial \omega}{\partial x} - y\right)
.. math::
    \sigma_{yz}(x, \ y) = \dfrac{\mathrm{M}_{z}}{J} \cdot \left(\dfrac{\partial \omega}{\partial y} + x\right)

.. math::
    \varepsilon_{xz}(x, \ y) = \dfrac{1}{2\mu} \cdot \sigma_{xz}
.. math::
    \varepsilon_{yz}(x, \ y) = \dfrac{1}{2\mu} \cdot \sigma_{yz}

Which :math:`J` is the :ref:`torsion_constant` and :math:`\mu` is the second `Lamé Parameter <https://en.wikipedia.org/wiki/Lam%C3%A9_parameters>`_.

To compute the partial derivatives, two approaches are used:

* For a point :math:`\mathbf{p}` on the boundary

    .. math::
        \nabla \omega & = \dfrac{\partial \omega}{\partial t} \cdot \mathbf{t} + \dfrac{\partial \omega}{\partial n} \cdot \mathbf{n} \\
        & = \left\langle \mathbf{p}, \ \mathbf{t}\right\rangle \cdot \mathbf{n} + \mathbf{t} \cdot \sum_{j=0}^{n-1} \varphi_{j}'(t) \cdot W_{j}

    The derivatives by themselves don't matter, but the evaluation of :math:`\sigma_{xz}` and :math:`\sigma_{yz}`, which are rewritten as 

    .. math::
        \begin{bmatrix}\sigma_{xz} \\ \sigma_{yz}\end{bmatrix} = \dfrac{\mathrm{M}_z}{J} \cdot \left[\left\langle\mathbf{p}, \ \mathbf{n}\right\rangle + \sum_{j=0}^{n-1}\varphi_{j}'(t) \cdot W_{j}\right] \cdot \mathbf{t}
        

* For interior points, :math:`\mathbf{p} \in \text{open}\left(\Omega\right)`


.. _shear_forces:

Shear Forces
------------

The shear forces :math:`\mathrm{F}_{x}` and :math:`\mathrm{F}_{y}` causes only shear stresses. 
The tensors are

.. math::
    \boldsymbol{\varepsilon} = \begin{bmatrix}0 & 0 & \varepsilon_{xz} \\ 0 & 0 & \varepsilon_{yz} \\ \varepsilon_{xz} & \varepsilon_{yz} & 0\end{bmatrix} \ \ \ \ \ \ \ \ \ \ \ \boldsymbol{\sigma} = \begin{bmatrix}0 & 0 & \sigma_{xz} \\ 0 & 0 & \sigma_{yz} \\ \sigma_{xz} & \sigma_{xz} & 0\end{bmatrix}

Depending on the application of the shear force, it may causes torsion.

TODO

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

This integral is computed in the boundary and the expression depends on :math:`\alpha`. 

For polygonal domains, the expressions may be resumed


.. dropdown:: Integrals :math:`I_{a,b}` for polygonal domains

    Expanding the expression of :math:`I_{a,b}` we get

    .. math::
        (a+b+2)\cdot I_{a,b} & = \dfrac{\alpha}{a+1} \sum_{i=0}^{n-1}\left(\left(y_{i+1}-y_{i}\right)\sum_{j=0}^{a+1}\sum_{k=0}^{b}\dfrac{\binom{a+1}{j}\binom{b}{k}}{\binom{a+b+1}{j+k}}x_{i}^{a+1-j}x_{i+1}^{j}y_{i}^{b-k}y_{i+1}^{k}\right) \\ & + \dfrac{\alpha-1}{b+1}\sum_{i=0}^{n-1}\left(\left(x_{i+1}-x_{i}\right)\sum_{j=0}^{a}\sum_{k=0}^{b+1}\dfrac{\binom{a}{j}\binom{b+1}{k}}{\binom{a+b+1}{j+k}}x_{i}^{a-j}x_{i+1}^{j}y_{i}^{b+1-k}y_{i+1}^{k}\right)

    By setting :math:`\alpha = 1`
    
    .. math::
        I_{a,0} = \sum_{i=0}^{n-1} \dfrac{x_{i+1}^{a+2}-x_{i}^{a+2}}{x_{i+1}-x_{i}} \cdot \dfrac{y_{i+1}-y_{i}}{(a+1)(a+2)}
    
    And :math:`\alpha = 0`

    .. math::
        I_{0,b} = -\sum_{i=0}^{n-1} \dfrac{y_{i+1}^{b+2}-y_{i}^{b+2}}{y_{i+1}-y_{i}} \cdot \dfrac{x_{i+1}-x_{i}}{(b+1)(b+2)}

    For any different value, the closed formulas are too complex. I don't have much time to find a :math:`\alpha` value such :math:`I_{a,b}` becomes a simpler expression. 

    Bellow you find values for :math:`\alpha = \dfrac{1}{2}`.

    .. math::
        I_{0,0} = \dfrac{1}{2}\sum_{i=0}^{n-1} x_{i}y_{i+1}-y_{i}x_{i+1}

    .. math::
        I_{1,1} = \dfrac{1}{24} \sum_{i=0}^{n-1} \left(x_{i}y_{i+1}-y_{i}x_{i+1}\right)\left(2x_{i}y_{i}+x_{i+1}y_{i}+x_{i}y_{i+1}+2x_{i+1}y_{i+1}\right)

    .. note::
        It's possible to have :math:`x_{i+1} = x_{i}` or :math:`y_{i+1} = y_{i}` in some edge, which leads to divide by zero in :math:`I_{a,0}` and :math:`I_{0,b}`.
        
        In that case, we open the expression:

        .. math::
            \dfrac{x_{i+1}^{c+1}-x_{i}^{c+1}}{x_{i+1}-x_{i}} = \sum_{j=0}^{c} x_{i}^{c-j}x_{i+1}^{j}
        .. math::
            \dfrac{y_{i+1}^{c+1}-y_{i}^{c+1}}{y_{i+1}-y_{i}} = \sum_{j=0}^{c} y_{i}^{c-j}y_{i+1}^{j}




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

.. note::
    The current implementation allows only polygonal domains. Hence, singular integrals are evaluated analiticaly as shown in :ref:`bem_polygonal_domain`

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