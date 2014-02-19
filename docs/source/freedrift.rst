====================================================
Free Drift Solution to the Sea Ice Momentum Equation
====================================================


.. default-role:: math

:author: David Huard, Bruno Tremblay
:date: June 2008



Context
-------

The free drift solution to the momentum equation yields the velocity
of the ice when neglecting the rheology of the ice, that is, the
internal forces occurring in the ice. This solution is used in the
model as a first approximation to the full solution, a step that is
necessary to ensure that the non linear solver has an adequate guess
to being its first iteration.

Equations
---------

The free drift equation is written as:

.. math::

   0 = -\rho_i h f \vec{k} \times \vec{u}_i + \tau_a  - \tau_w - \rho_i h g \nabla H_d

where :math:`\rho_i` is the ice density, `h` the ice thickness,
`\vec{k}` the vector normal to the surface (pointing up), `\vec{u}_i`
the ice velocity vector, `\tau_a` the wind shear stress, `\tau_w` the
water shear stress and `H_d` the sea surface height. Note that the
sign of the water stress is the inverse of that of the wind (air)
stress.

.. rubric:: Coriolis Force 

Let's break this down ! The first term on the left is the Coriolis
force and its components in the `x` and `y` directions are given by:

.. math::

   -\rho_i h f \vec{k} \times \vec{u}_i =
   \left[\begin{matrix}
   \rho_i h f u_y \\
   -\rho_i h f u_x\\
   0
   \end{matrix} \right],

where `f` is the Coriolis parameter and is given by
`f=2\omega\sin\phi`, `\omega` being Earth's angular velocity and
`\phi` the latitude. Since the model is currently restricted to the
Arctic, the Coriolis parameter is taken constant.  From now on, we
will ignore the velocity in third dimension who is always 0 since we
are considering motion on a plane.

.. rubric:: Surface stress

The second and third terms are surface stresses, and are defined in
the model as quadratic laws. That is, the stress depends on the square
of the velocity.  Both the air and water stresses can be written in
the general form:

.. math::

   \tau = \rho C_d |\vec{du}|(\vec{du} \cos\theta + \vec k \times \vec{du} \sin \theta),

where `\rho` is the dragging medium density, `C_d` the drag
coefficient, `\vec{du}` the relative speed between the forcing medium
and the ice and `\theta` the turning angle. A simpler for surface
stress can be obtained if we assume that the stress is related
linearly to the forcing medium velocity : 

.. math::

   \tau = \rho C_d (\vec{du} \cos\theta + \vec k \times \vec{du} \sin \theta), 


Since wind speeds are much higher than the ice speeds, we assume for
the air stress that `\vec{u}_a -\vec{u}_i \approx \vec{u}_a`.

Specifically, for the quadratic law, the components of the air stress are:

.. math::

   \tau_a =\left[ \begin{matrix}
   \rho_{a}C_{da}\sqrt{v_a^2+u_a^2}\left(\cos\theta_a u_a - \sin\theta_a v_a \right)\\
   \rho_{a}C_{da}\sqrt{v_a^2+u_a^2}\left(\cos\theta_a v_a + \sin\theta_a u_a \right)
   \end{matrix}  \right]

and the components of the water stress are:

.. math::

   -\tau_w =
   \left[ \begin{matrix}-C_{dw}\rho_{w}\sqrt{\left(v_{i}-v_{w}\right)^2+ \left(u_{i}-u_{w}\right)^2}\left( \left(u_{i}-u_{w}\right)\cos \theta_{w} - \left(v_{i}-v_{w}\right)\sin  \theta_{w}  \right) \\
   -C_{dw}\rho_{w}\sqrt{\left(v_{i}-v_{w}\right)^2+\left(u_{i} -u_{w}\right)^2}\left( \left(v_{i}-v_{w}\right)\cos \theta_{w} + \left(u_{i}-u_{w}\right)\sin \theta_{w}   \right)
   \end{matrix} \right]


.. rubric:: Sea Surface Tilt 

The last term, the sea surface tilt, is expressed using the ocean geostrophic
currents. Geostrophic water currents are the currents that occur when we
neglect all forces except the Coriolis force and the gravity:

.. math::

   f \vec{u}_w = g \vec k \times \nabla H.

We use this relation to express the force due to sea surface tilt by the
geostrophic currents, that is,

.. math::

   \vec k \times f\vec{u}_w & = \vec k \times (g \vec k \times \nabla H) \\f \vec k \times \vec{u}_w & = -g \nabla H

so that the last term of the first equation can be written as

.. math::

   - \rho_i h g \nabla H_d =\rho_i h f \vec k \times \vec{u}_w 
   =\left [ \begin{matrix}
   -\rho_i h f v_w \\
   \rho_i h f u_w
   \end{matrix} \right].

Bringing everything together, we have:

.. math::

   0 =\left[ \begin{matrix}
   C_{dw}\rho_{w}\sqrt{v^2+u^2}\left(v\sin\theta_{w} - u\cos \theta_{w}\right)+ \rho_{i} f h v + C_{da}\rho_{a} \sqrt{v_{a}^2+u_{a}^2} \left(u_{a}\cos  \theta_{a}-v_{a}\sin  \theta_{a}\right)\\
   -C_{dw}\rho_{w}\sqrt{v^2+u^2 }\left(u\sin \theta_{w}+v\cos \theta_{w}\right) - \rho_{i} f h u + C_{da}\rho_{a} \sqrt{v_{a}^2+u_{a}^2} \left(u_{a}\sin \theta_{a}+v_{a}\cos \theta_{a}\right)
   \end{matrix} \right]

where `u = u_i-u_w` and `v = v_i-v_w`. If we call the right hand side of the
momentum equation `F(\vec{u})`, our next task is to find `\vec{u}` such
that `F(\vec{u})=0`, that is, find the roots of the equations.


Solving for `\vec u`
--------------------

To find a solution for `\vec u` the ice velocity relative to the ocean
currents, we use the Newton-Raphson (N-R) method. This method is a
generic method used to find the zero of a function, that is, the value
`x` for which `f(x)=0`. To do so, N-R uses a starting guess `x_n`, and
computes the next guess by iteration using

.. math::

   x_{n+1} = x_n - \frac{f(x_n)}{f\prime(x_n)}.

For functions of more than one variable, the idea is similar, but the derivative
is replaced by the Jacobian matrix :

.. math::

   \vec{u}_{n+1} = \vec{u}_{n} - J^{-1}(\vec{u}_n) F(\vec{u}_n).

For our particular case, we have a vector function `\vec F(u_x, u_y)` that has two
components `F_x` and `F_y`. The Jacobian can we written as follows:

.. math::

   J =
   \left[ \begin{matrix}
   J_{11}  & J_{12} \\
   J_{21} & J_{22} \end{matrix} \right]
   =
   \left[ \begin{matrix}
   \frac{\partial F_x}{\partial u} &  \frac{\partial F_x}{\partial v} \\
   \frac{\partial F_y}{\partial u} &  \frac{\partial F_y}{\partial v}
   \end{matrix}\right].

For a quadratic surface stress, the terms of the Jacobian are:

.. math::

    \partial_u F_x & =  \rho_w C_{dw} \frac{\left[ uv\sin\theta_{w} - \left( v^2 + 2 u^2\right)\cos\theta_{w}\right]}{\sqrt{v^2+u^2}} \\
    \partial_v F_x & =  \rho_w C_{dw} \frac{\left( 2 v^2 + u^2\right) \sin\theta_{w}-uv\cos \theta_{w} }{\sqrt{v^2+u^2}} + \rho_i h f\\
    \partial_u F_y & = -\rho_w C_{dw} \frac{\left( v^2 + 2u^2 \right) \sin\theta_{w}+uv\cos \theta_{w} }{\sqrt{v^2+u^2}} - \rho_i h f\\
    \partial_v F_y & = -\rho_w C_{dw} \frac{\left[ uv\sin\theta_{w} + \left( 2 v^2 + u^2\right)\cos\theta_{w}\right]}{\sqrt{v^2+u^2}}


while they are much simpler if a linear surface stress is assumed:

.. math:: 

    \partial_u F_x & = -\rho_w C_{dw} \cos\theta_{w} \\
    \partial_v F_x & = \rho_w C_{dw} \sin\theta_w + \rho_i h f\\
    \partial_u F_y & = -\rho_w C_{dw} \sin\theta_w - \rho_i h f\\
    \partial_v F_y & = -\rho_w C_{dw} \cos\theta_w

Knowing that for a 2x2 matrix `J`, the inverse is computed using:

.. math::

    J^{-1} =
    \frac{1}{\det J}\left[ \begin{matrix}
    J_{22}  & -J_{12} \\
    -J_{21} & J_{11}
    \end{matrix} \right].

where `\det J = J_{11} J_{22} - J_{21} J_{12}`. Assuming `M = J^{-1}` we have
everything we need to know to compute `\vec u`:

.. math::

    \det J & = \partial_u F_x \partial_v F_y - \partial_u F_y \partial_v F_x \\
    M_{11} & = \frac{\partial_v F_y}{\det J} \\
    M_{12} & = \frac{-\partial_y F_x}{\det J} \\
    M_{21} & = \frac{-\partial_u F_y}{\det J} \\
    M_{22} & = \frac{\partial_u F_x}{\det J}

so that

.. math::

    u_{n+1} & = u_{n} - ( M_{11} F_x(\vec u_n) + M_{12} F_y(\vec u_n))  \\
    v_{n+1} & = v_{n} - ( M_{21} F_x(\vec u_n) + M_{22} F_y(\vec u_n))

Remember that this is the solution for the relative ice velocity,not the ice velocities themselves. 

Original Numerical Implementation
---------------------------------

In the original code ``UV_solve_NR_B``, here are the steps performed to find a 
solution to the free drift equation:

  1. Set `\vec{u}=0`

  2. Compute the relative linear free drift solution on the B-grid. The linear 
     free drift solution refers to the solution to the momentum equation when
     the rheology term is absent and where the surface stresses are proportional
     to the medium's speed. That is, instead of the quadratic law for the stress, 
     we have `\tau = \rho C (\vec{u}\cos\theta + \vec k \times \vec u \sin\theta`.
     Since the function is linear, the differentation is straightforward and 
     the N-R step leads us directly to the analytical solution for `\vec u`. 
     This solution will be the first guess for the nonlinear solution. 

  3. Loop over the solution to the nonlinear free drift equation until a 
     convergence criteria is met or a maximum number of loops is attained. 
     This solution is described above, and consists in computing, for each `u` and
     `v` in the matrix, the coefficients of the inverse Jacobian matrix, the 
     solution to `\vec F(\vec u)` at the last iteration and corresponding 
     `\delta \vec u`. 

  4. Interpolate the solution on the C-grid. 

  5. Apply the boundary conditions: the velocity and its derivative along the 
     two directions are 0. 


Proposed Numerical Implementation
---------------------------------

A modular approach is to code a simple primitive procedure for each
stress component: Coriolis, linear surface stress, quadratic surface
stress, sea surface tilt. A free_drift_momentum_eq subroutine can then
be written using calls to the primitive functions. Since there are no
interactions between adjacent grid cells, the primitive procedures can
be declared ELEMENTAL[#elemental]_. The main advantage is that those
procedures can be reused in other parts of the code without
modification. A ``linear_freedrift_solve`` function would solve the
momentum equation for a linear drag law. The result would then be used
as a starting guess in an iterative ``quadratic_freedrift_solve``. The
latter subroutine would itself call a ``quadratic_freedrift_step``
which would compute a `\Delta \vec u` for one step using the momentum
equation and the computed Jacobian at the current `\vec u`. 

.. [#elemental] Elemental statements are applied to procedure
                operating on scalars to tell the compiler that 
		it has no side effect and can be safely applied
		in a loop to all the elements of an array. In other
                words, it vectorizes a scalar procedure. 

Notes
-----

In this text, `M11, M12, M21` and `M22` refer to the opposite of
variables `A11, A12, A21` and `A22` in the code. I chose this different notation
since there is mix up potential: in the paper, A is used both
for the ice covered area and a `2x2` matrix multipling the relative ice
velocity; in the code,  Aij is the negative of the inverse of the Jacobian.

