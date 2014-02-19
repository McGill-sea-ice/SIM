~~~~~~~~~~~~~~~~~
Advection Schemes
~~~~~~~~~~~~~~~~~


This section describes the upwind (upstream) advection scheme that is currently
used in the model. Next, we describe other advection schemes, their advantages and
disadvantages. 

Some of the material in this section is inspired from the `MIT GCM`_.


.. glossary::

   Courant number
      The Courant number is defined as `C=u \Delta t/\Delta x`, the
      velocity times the time step divided by the grid resolution. It
      is used to characterize advection conditions.  
      If `C=1`, it
      means that a quantity is advected over an entire grid cell in a
      single time step, which can corrupt some advection schemes. It
      is hence important to check the typical values of the Courant
      number to make sure the advection scheme is appropriate.


.. _`MIT GCM`: http://mitgcm.org/sealion/online_documents/node71.html


Advection schemes are computational methods used to solve the differential equation

.. math::

   \partial_t A + u \partial_x A = 0.

Look at Numerical Models of Oceans and Oceanic Processes for a good reference. 


Upwind Advection Scheme
~~~~~~~~~~~~~~~~~~~~~~~

This is the simplest advection scheme. It solves the advection equation using backward space difference
and forward time differences. This is what is currently used in the model to advect ice thickness and ice
covered area. The main problem with this scheme is that it is diffusive. That is, if we write the 
Taylor series of the differential equation, we see that the solving the finite difference approximation
amounts to solving an advection-diffusion equation. In other words, the scheme adds a diffusive term
that smooths the advected quantity. 