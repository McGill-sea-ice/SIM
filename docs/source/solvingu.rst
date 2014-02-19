~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Solving the Momentum Equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:author: David Huard
:date: August 5, 2008

.. default-role:: math



You should read the section on solving the free drift equation before
tackling this section on solving the complete momentum equation (ME).  


The ME is a non-linear equation. Furthermore, we do not have an analytical
formula for the Jacobian, hence, we cannot use the same technique as in 
the free drift equation (Newton's algorithm). The current strategy to solve
for `u` is to iterate over the following steps:

 1. Linearize the equation. This means that all non-linear terms in `u` are
    set to the value of `u` at the last iteration. 
 2. Once the non-linear terms are fixed, we have a linear equation in `u`. 
    We can use a linear solver to find `u` pretty efficiently. This is now
    done using GMRES. 
 3. Use the solution to the linear equation to update `u`, and return 
    to step 1. The way this is currently done is rather simple, but we are working
    on implementing a Jacobian Free Newton Method to improve this. 

The number of iteration (linearize, solve the linear equation, update) is called an 
outer loop in our jargon. To get a fully converged solution to the non-linear equation
can take hundreds of outer loop. 


Constant Terms
~~~~~~~~~~~~~~
The ME contains terms that are independent of the ice velocity `u`. These are:
 
 * The wind stress, which only depends on the wind velocity. 
 * The sea ice tilt, which only depends on the water velocity. 
 * The pressure term, `\nabla \cdot P`. 

Linear Terms
~~~~~~~~~~~~
The linear terms in the momentum equations are 

 * The Coriolis force. 
 * The linear component of the water stress.

Non-Linear Terms
~~~~~~~~~~~~~~~~

The non-linear terms in the ME are
 
 * The non-linear component of the water stress, `|\vec u - \vec{u_w}|`. 
 * The strain rate terms in the constitutive law. 

   .. math:: 

      \nabla \cdot (2 \eta \dot{\epsilon}_{ij} + [\zeta-\eta]\dot{\epsilon}_{kk}\delta_{ij})





Relation to the Code
~~~~~~~~~~~~~~~~~~~~

Although I am not 100% sure to understand everything in the code, here is rough breakdown of what 
happens before and during an outer loop: 

Before: 
 # Load wind velocity, compute wind stress and the sea ice tilt. 
 # Compute the stress related to the pressure. 
 # Put those terms (air stress, sea ice tilt, pressure) in the variables R1p and R2p.   

During: 
 # Compute the non-linear terms in the water stress (r1ppr2pp).
 # Compute the viscous coefficients (viscouscoefficients). 
 # Solve the linear equation (prep_fgmres). 
 # Update the ice velocity (r1ppr2pp).


To facilitate the implementation of other solvers and clean up the code, 
I propose to separate those chunks into well identified individual 
procedural subroutines:
 
 * Subroutine for the wind stress (done).
 * Subroutine for the sea ice tilt (done). Since this does not
   vary at all, ever, store that in a variable. 
 * Subroutine for the ice pressure (done, untested). 
 * Subroutine to gather the independent terms (todo).
 * Subroutine to compute the non linear term in the water stress. 
 * Subroutine to compute the linear term in the water stress. 
 * Subroutine to compute the viscous coefficients (viscouscoefficient, convert to procedural).
 * Subroutine to compute the left hand side (matvec, convert to procedural)
 * Subroutine to gather the right hand side (r1ppr2pp, clean up). 
 * Subroutine to compute the linear and non-linear residuals (convert to procedural). 



 