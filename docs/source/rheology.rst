========
Rheology
========



.. default-role:: math

.. glossary::


    Constitutive law
	A constitutive law is an equation relating the strain of some material to
	the imposed stress.

    Deviatoric stress tensor `\sigma\prime_{ij}`
	The difference between the stress tensor and the isotropic pressure `p`. Note
	that this mechanical pressure cannot always be identified with the thermodynamic
	pressure. 

    Isotropic pressure `p`
	The average normal pressure

	.. math:: 

           p = -\frac{1}{2}(\sigma_{kk}).

    Principal stresses 
	The principal stresses are the eigenvalues `\sigma_1`, `\sigma_2` of the 
	stress tensor with the convention that `\sigma_1 > \sigma_2`. This can be 
	understood as the diagonal elements of the stress tensor when the coordinate
	system is rotated such that the off-diagonal elements are 0:

	.. math::
	   
            \sigma^{\prime_{ij}} = \left[ \begin{matrix}
	    \sigma_{1} & 0 \\
	    0 & \sigma_{2}
	    \end{matrix}\right]

	.. figure:: principal_stress.png
	    :scale: 50
	    :align: center

	    Principal stresses.


    Rheology
	Rheology is the study of flow and matter under the influence of a stress. 


    Shear strain rate 
	The rate at which two sides of an elemental square close together:

	.. math::
	   
           \dot{\gamma}_{ij} = \left( \frac{\partial u_i}{\partial x_j} +  \frac{\partial u_j}{\partial x_i} \right)



    Shear stress `\tau_{ij}`
	The pure shearing stress. 


    Stress
	A stress is a force per unit area on a surface. The stress field is the 
	distribution of internal stresses that balance a set of external forces. 

    Stress tensor `\sigma_{ij}`
	A matrix describing the stresses on an infinitesimal surface (or cube in 
	3D). The first index `i` refers to the direction of the normal to the place 
	the stress acts on. The second index `j` refers to the direction of the
	stress. For example, `\sigma_{xx}` refers to the stress on the y-plane in
	the x direction, while \sigma_{xy}`refers to the stress on the y-plane in 
	the y direction. A normal stress refers to a stress perpendicular to the 
	surface (`\sigma_{xx}, ], \sigma_{yy}`), while a shear stress refers to a 
	stress parallel to the surface (`\sigma_{xy}`). By convention, tensile 
	stresses are positive while conpressive ones are negative. In a 2-D 
	coordinate system, stresses are described using a `2x2` matrix `\sigma_{ij}` 

	.. math::
	  
            \sigma_{ij} = \left[ \begin{matrix}
	    \sigma_{11} & \sigma_{12} \\
	    \sigma_{21} & \sigma_{22}
	    \end{matrix}\right]

	To satisfy the equilibium of moments (no rotation), the stress tensor needs
	to be symmetric, ie, `\sigma_{ij} = \sigma_{ji}`. 

	.. figure:: stress.png
	   :scale: 5
	   :align: center

	   State of stress. 

    Stress invariants
	The stress invariantes are defined as 

	.. math::
	  
            \dot{\epsilon}_I & \equiv \dot{\epsilon}_1 + \dot{\epsilon}_2 = \text{divergence} \\
	    \dot{\epsilon}_{II} & \equiv \dot{\epsilon}_1 - \dot{\epsilon}_2 = \text{maximum  shear rate} 

	If we assume that `\dot{\epsilon}_I` and `\dot{\epsilon}_{II}` are the 
	real and complex components of a complex variable, and `\theta` the angle
	between the complex vector and the real axis, then `\theta=0, \pi/4, \pi/2,
	3\pi/4` and `\pi` correspond to pure divergence, uniaxial extension, pure
	shear, uniaxial contraction and pure convergence respectively. If you find 
	the last statement puzzling, please stop and think about it before going on.   


    Strain `\epsilon`
	The deformation of a material `\frac{dl}{l}`. In 3D you can imagine the 
	strain as the displacement of the corner of a cube under stress. 

	.. figure:: types_of_motion.png
	    :align: center
	    :width: 600 


    Strain rate `\dot{\epsilon}`
	The rate of deformation over time `\frac{1}{l}\frac{dl}{dt}=\frac{v}{l}`.
	In 2-D, this is described by the strain rate tensor:

	.. math::
	  
           \dot{\epsilon}_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)

	The diagonal elements of the strain rate tensor are the normal strain rates. The 
	off-diagonal elements are one half the shear strain rate components (`\dot{\gamma}_{ij}`). 

	Note that some authors denote the strain rate simply by `\epsilon`. 

    Velocity gradient tensor
	The velocity gradient tensor is written as  `\partial u_i / \partial x_j`. The 
	symmetrical part of the tensor is the strain rate tensor, and the antisymmetrical
	part is the rotation tensor, which is related to the vorticity (`\omega`).   

	.. math::
	 
          \frac{\partial u_i}{\partial x_j}   = \dot{\epsilon}_{ij} + \Omega_{ij}

	where `\Omega = \frac{1}{2} \nabla \times \vec{u} = \frac{1}{2} \omega`. 


    Yield stress
	The minimum stress that must be applied to initiate significant flow and 
	a significant drop in viscosity. 

    Newtonian fluid
	A fluid whose viscosity is independent of the shear conditions. The 
	resistance between two layers of fluid is proportional to the difference in
	speed of these fluids. The constitutive law for a compressible  Newtonian fluid is 

	.. math::
	
           \sigma_{ij} = -2 \eta \dot{\epsilon}_{ij} + \left( \frac{2\eta}{3} - \zeta \right ) \dot{\epsilon}_{kk}\delta_{ij}

    Elastic material
	A material is elastic if the deformation follows the applied stress. 
	A perfect elastic material has a stress proportional to the strain. 

    Viscous material
	A material is viscous if the deformation increases linearly at a constant
	stress. In other words, you need to keep pushing on it to make it move. 

    Yield curve
	Curve describing the stress at which a material starts to yield.    

    Bulk viscosity `\zeta`
	Viscosity associated with volume expansion, also called volume viscosity or
	Lamé's constant.

    Shear viscosity `\eta`
	Viscosity when the applied stress is a shear stress. 


Elliptical Yield Curve
----------------------

.. math::

    \zeta & = P/2\Delta \\
    \eta & = \zeta/e^2 \\
    \Delta & = [(\dot{\epsilon}_{11}^2 + \dot{\epsilon}_{22}^2)(1 + 1/e^2) + 4e^{-2}\dot{\epsilon}_{12}^2 + 2 \dot{\epsilon}_{11} \dot{\epsilon}_{22}(1 - 1/e^2)]^{1/2}


Viscous Plastic Constitutive Law
--------------------------------

.. math::
 
   2 \eta\dot{\epsilon}_{ij} + [\zeta - \eta]\dot{\epsilon}_{kk}\delta_{ij} - P\delta_{ij}/2


Internal Stress Force
---------------------
The force components due to internal stresses in the ice are calculated from `F_i = \partial \sigma_{ij}/\partial x_j` (summation is implied on the repeated indices.)




