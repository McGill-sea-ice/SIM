!************************************************************************                                  
!     Subroutine ViscousCoefficient: computes the non-linear viscous                                         
!     coefficient at the centers (zetaC, etaC) and at the nodes (etaB) of the                                
!     grid. All the strain rate calculations are second order accurate.                                      
!     The ice strength (Pp) is collocated with zetaC, etaC.
!                                                                                                           
!     Calculation of viscous coefficients eta and zeta at the grid node (B)
!     and at the grid center (C).    
!                        
!     The user can specify the method to calculate the viscous coeff by 
!     setting the int. variable visc_method to 1, 2, 3 or 4 :
!
!     1) Tremblay and Mysak 1997
!        -ep12c=dudy+dvdx, ep11c and ep22c are calc at grid center. 
!        -zetaC (and etaC) are calc at grid center.
!        -etaB = ( etaC + etaC + etaC + etaC ) / 4 
!
!     2) Lemieux et al. 2012
!        -ep12c, ep11c and ep22c are calc at grid center. 
!        -zetaC (and etaC) are calc at grid center.
!        -ep12b, ep11b and ep22b are calc at grid node. 
!        -etaB is calc at grid node.
!     
!     3) Bouillon et al. 2013 (modified)
!        -ep11c and ep22c are calc at grid center. 
!        -ep12b are calc at the grid nodes.
!        -ep12c=( ep12b + ep12b + ep12b + ep12b ) / 4
!        -zetaC (and etaC) = function of ep11c, ep22c and ep12c.
!        -etaB = ( etaC + etaC + etaC + etaC ) / 4 
!
!     4) Bouillon et al. 2013
!        -ep11c and ep22c are calc at grid center. 
!        -ep12b are calc at the grid nodes.
!        -ep12c**2=( ep12b**2 + ep12b**2 + ep12b**2 + ep12b**2 ) / 4
!        -zetaC (and etaC) = function of ep11c, ep22c and ep12c**2.
!        -etaB = ( etaC + etaC + etaC + etaC ) / 4 
!
!     note1 : all these approaches use the replacement closure of 
!             Kreyscher et al 2000 (eq.7).                                   
!                                          
!     note2 : methods 1 and 2 give a VP solution everywhere. Methods 3 
!             leads to a VP solution everywhere except close to a coast
!             with a corner. Methods 4 leads to a VP solution everywhere 
!             except close to a coast with a corner and in the ocean 
!             where there are strong deformations.
!
!     note3 : in terms of robustness of the JFNK solver, the methods are 
!             lister in order of robustness (method 4 is the more robust, 
!             method 1 leads to more failures of the solver). 
!
!     JF Lemieux, 30 June 2015                                                                            
!                                                                                                            
!************************************************************************     


subroutine ViscousCoefficient(utp,vtp)
 
  implicit none

  include 'parameter.h'
  include 'CB_options.h'
  
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
  if (Rheology .eq. 3) then
     call MEBcoeff
  elseif (visc_method .eq. 1) then
     call ViscousCoeff_method1(utp,vtp)
  elseif (visc_method .eq. 2) then
     call ViscousCoeff_method2(utp,vtp)
  elseif (visc_method .eq. 3 .or. visc_method .eq. 4) then
     call ViscousCoeff_method3_and_4(utp,vtp)
  else
     print *, 'Wrong setting of visc_method'
     stop
  endif

  call etaB_at_open_boundaries
  
  return
end subroutine ViscousCoefficient

!************************************************************************
!     Subroutine ViscousCoeff_method1: computes the non-linear viscous
!     coefficient at the centers (etaC) and the nodes (as the average of
!     the etaC).
!
!     We also apply the replacement closure of Kreyscher et al 2000 (eq.7).
!
!************************************************************************
    
subroutine ViscousCoeff_method1(utp,vtp)
  use ellipse
  implicit none

  include 'parameter.h'
  include 'CB_DynVariables.h'
  include 'CB_const.h'
  include 'CB_mask.h'
  include 'CB_options.h'

  integer i, j, rheo, peri

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

  denomin = 2d-09 ! Hibler, 1979
  rheo = Rheology ! define local variable to speed up the code
  peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions

!------------------------------------------------------------------------
!     free slip boundary condition:
!       d(v_tangential)/d(normal) = 0 & v_normal = 0, at close boundary
!------------------------------------------------------------------------


  if ( BndyCond .eq. 'freeslip' ) then

     print *, 'WARNING: freeslip option not fully tested yet'

!------------------------------------------------------------------------
!     strain-rates calculation
!     no-slip boundary condition :  v_normal = 0 and v_tangential = 0 
!------------------------------------------------------------------------

  elseif ( BndyCond .eq. 'noslip' ) then


     do i = 0, nx+1
        do j = 0, ny+1

           etaC(i,j)  = 0d0
           zetaC(i,j) = 0d0 
           etaB(i,j)  = 0d0
 
        enddo
     enddo


     dudx       = 0d0
     dvdy       = 0d0
     dudy       = 0d0
     dvdx       = 0d0


     
     if (peri .ne. 0) call periodicBC(utp,vtp)
     

     do i = 1, nx
        do j = 1, ny

           if ( maskC(i,j) .eq. 1 ) then
                  
              dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
              dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
              if     ( maskC(i+1,j) + maskC(i-1,j) .eq. 2 ) then
                     
                 dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                      ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                      ( 4d0 * Deltax )
                     
              elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. 1 ) then
                     
                 dvdx = ( 1d0 * ( vtp(i+1,j) + vtp(i+1,j+1) ) +  &
                      3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) /  &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. -1 ) then
                     
                 dvdx = ( -1d0 * ( vtp(i-1,j) + vtp(i-1,j+1) ) - &
                      3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i+1,j) + maskC(i-1,j) .eq. 0 ) then
                     
                 print *, 'WARNING: irregular grid cell case1', i, j
                     
              endif

               
              if     ( maskC(i,j+1) + maskC(i,j-1) .eq. 2 ) then
                     
                 dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                      ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                      ( 4d0 * Deltax )
                     
              elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                     
                 dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                      3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                     
                 dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                      3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i,j+1) + maskC(i,j-1) .eq. 0 ) then
                     
                 print *, 'WARNING: irregular grid cell case2',i,j
                     
              endif
                  
!------------------------------------------------------------------------
!     Shear and bulk viscosity calculation at the grid center
!------------------------------------------------------------------------

              if ( rheo .eq. 1 ) then ! ellipse, jfl p.892
                     

                 deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                      + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                      + ell_2 * ( dvdx + dudy ) ** 2 )

                 if ( regularization .eq. 'tanh' ) then

                    deno = max( deno, 1d-20 )
                     
                    zetaC(i,j) =  ( (Pp(i,j)+Pt(i,j))/denomin ) &
                         *( tanh(denomin*(1/deno)))
                         
                    P(i,j) = ((( Pp(i,j)-Pt(i,j))/denomin ) &
                         *tanh(denomin*(1/deno))) *deno 

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin

                    zetaC(i,j) = (Pp(i,j) + Pt(i,j))/denoT 
                    
                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)

                    zetaC(i,j) = (Pp(i,j) + Pt(i,j)) / denoT

                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                        
                 endif

                 etaC(i,j)  = zetaC(i,j) * ell_2
                     
              elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                 stop
                   
              endif
                  
           endif
               
        enddo
     enddo

     do i = 1, nx+1

        if (maskC(i,0) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
            
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1 .and. Periodic_x .eq. 0) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1 .and. Periodic_x .eq. 0) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo
     
     if (peri .ne. 0) call periodicBC(etaC,zetaC)
     if (peri .ne. 0) call periodicBC(etaB,P)  !etaB = zero and is a dummy here
!------------------------------------------------------------------------
!     Shear and bulk viscosity calculation at the grid node
!------------------------------------------------------------------------
         
     do i = 1, nx+1
        do j = 1, ny+1

           deno = max ( maskC(i,j) + maskC(i-1,j) + maskC(i,j-1) +   &
                maskC(i-1,j-1), 1) 
               
           etaB(i,j)  = ( etaC(i,j) + etaC(i-1,j) + etaC(i,j-1) +    &
                etaC(i-1,j-1) ) / deno
           
        enddo
     enddo
         
  endif
  

      
  return
end subroutine ViscousCoeff_method1


!************************************************************************
!     Subroutine ViscousCoeff_method2: computes the non-linear viscous
!     coefficient at the centers (etaC) and at the nodes (etaB) of the 
!     grid. All the strain rate calculations are second order accurate. 
!     The ice strength (p) is collocated with etaC. A first order
!     expansion is done To get p at the nodes.
!
!     We also apply the replacement closure of Kreysher et al 2000(eq. 7). 
!
!************************************************************************


subroutine ViscousCoeff_method2(utp,vtp)
  use ellipse
  implicit none
  
  include 'parameter.h'
  include 'CB_DynVariables.h'
  include 'CB_const.h'
  include 'CB_mask.h'
  include 'CB_options.h'
  
  integer i, j, rheo, summaskC, peri

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin, pnode
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
  
  denomin = 2d-09 ! Hibler, 1979

  rheo = Rheology ! define local variable to speed up the code
  peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions

!------------------------------------------------------------------------
!     free slip boundary condition:
!       d(v_tangential)/d(normal) = 0 & v_normal = 0, at close boundary
!------------------------------------------------------------------------


  if ( BndyCond .eq. 'freeslip' ) then

     print *, 'WARNING: freeslip option not fully tested yet'

!------------------------------------------------------------------------
!     strain-rates calculation
!     no-slip boundary condition :  v_normal = 0 and v_tangential = 0 
!------------------------------------------------------------------------

  elseif ( BndyCond .eq. 'noslip' ) then


     do i = 0, nx+1
        do j = 0, ny+1

           etaC(i,j)  = 0d0
           zetaC(i,j) = 0d0 
           etaB(i,j)  = 0d0
           
        enddo
     enddo


     dudx       = 0d0
     dvdy       = 0d0
     dudy       = 0d0
     dvdx       = 0d0
     pnode      = 0d0


     if (peri .ne. 0) call periodicBC(utp,vtp)    
     
!------------------------------------------------------------------------
!     Shear and bulk viscosity calculation at the grid center
!------------------------------------------------------------------------   

     do i = 1, nx
        do j = 1, ny               
               
           if ( maskC(i,j) .eq. 1 ) then
                  
              dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
              dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
              if     ( maskC(i+1,j) + maskC(i-1,j) .eq. 2 ) then
                     
                 dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                      ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                      ( 4d0 * Deltax )
                     
              elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. 1 ) then
                     
                 dvdx = ( 1d0 * ( vtp(i+1,j) + vtp(i+1,j+1) ) +  &
                      3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) /  &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. -1 ) then
                     
                 dvdx = ( -1d0 * ( vtp(i-1,j) + vtp(i-1,j+1) ) - &
                      3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i+1,j) + maskC(i-1,j) .eq. 0 ) then
                     
                 print *, 'WARNING: irregular grid cell case1', i, j
                     
              endif

               
              if     ( maskC(i,j+1) + maskC(i,j-1) .eq. 2 ) then
                     
                 dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                      ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                      ( 4d0 * Deltax )
                 
              elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                     
                 dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                      3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                     
                 dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                      3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                      ( 6d0 * Deltax )
                     
              elseif ( maskC(i,j+1) + maskC(i,j-1) .eq. 0 ) then
                     
                 print *, 'WARNING: irregular grid cell case2',i,j
                 
              endif
                  

              if ( rheo .eq. 1 ) then ! ellipse, jfl p.892


                 deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                      + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                      + ell_2 * ( dvdx + dudy ) ** 2 )

 
                 if ( regularization .eq. 'tanh' ) then

                    deno = max( deno, 1d-20 )
                     
                    zetaC(i,j) =  ( (Pp(i,j)+Pt(i,j))/denomin ) &
                         *( tanh(denomin*(1/deno)))
                         
                    P(i,j) = ((( Pp(i,j)-Pt(i,j))/denomin ) &
                         *tanh(denomin*(1/deno))) *deno 

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin

                    zetaC(i,j) = (Pp(i,j)+Pt(i,j))/denoT 
                    
                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)

                    zetaC(i,j) = (Pp(i,j)+Pt(i,j))/denoT 

                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno 

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                        
                 endif

                 etaC(i,j)  = zetaC(i,j) * ell_2
                     
              elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                 stop
                   
              endif
                  
           endif
               
        enddo
     enddo

!------------------------------------------------------------------------
!     Set zetaC and etaC to 0.0 at the open boundaries (see p.32-33 EC-2)
!------------------------------------------------------------------------

     do i = 1, nx+1

        if (maskC(i,0) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
        
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1 .and. Periodic_x .eq. 0) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1 .and. Periodic_x .eq. 0) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo

     if (peri .ne. 0) call periodicBC(etaC,zetaC)
     if (peri .ne. 0) call periodicBC(etaB,P)  !etaB = zero and is a dummy here     

!------------------------------------------------------------------------
!     Shear viscosity calculation at the grid node (see p.2-118 PDF notebook)
!------------------------------------------------------------------------

         
         ! for sig12B (defined at the node)
     do j = 1, ny+1 
        do i = 1, nx+1

           summaskC = maskC(i-1,j) + maskC(i,j) + & 
                maskC(i,j-1) + maskC(i-1,j-1)

           if (summaskC .ge. 2) then

              if (summaskC .eq. 4) then
! oo
! oo normal
                 dudy = ( utp(i,j) - utp(i,j-1) ) / Deltax !case 1 
                 dvdx = ( vtp(i,j) - vtp(i-1,j) ) / Deltax
		     
                 dudx = ( (utp(i+1,j) + utp(i+1,j-1)) * maskB(i+1,j) - &
                      (utp(i-1,j) + utp(i-1,j-1)) * maskB(i-1,j) ) / (4d0*Deltax)
                 
                 dvdy = ( (vtp(i-1,j+1) + vtp(i,j+1)) * maskB(i,j+1) - &
                      (vtp(i-1,j-1) + vtp(i,j-1)) * maskB(i,j-1) ) / (4d0*Deltax)
                     
                 pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i,j-1) + Pp(i-1,j-1) + & 
                           Pt(i-1,j) + Pt(i,j) + Pt(i,j-1) + Pt(i-1,j-1) ) / 4d0 
                     
              elseif (summaskC .eq. 3) then 

                 if (maskC(i-1,j) .eq. 0) then !case 2
! xo
! oo    
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax

                    dudx = ( utp(i+1,j) + utp(i+1,j-1) - &
                         ( utp(i+2,j) + utp(i+2,j-1) ) * maskB(i+2,j)/4d0 ) / Deltax
		       
                    dvdy = ( -vtp(i-1,j-1) - vtp(i,j-1) + &
                         ( vtp(i-1,j-2) + vtp(i,j-2) ) * maskB(i,j-2)/4d0 ) / Deltax
                        
                    pnode = ( Pp(i,j) + Pp(i,j-1) + Pp(i-1,j-1) + & ! check ca 
                              Pt(i,j) + Pt(i,j-1) + Pt(i-1,j-1) ) / 3d0 ! check ca 

                 elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax
		       
                    dvdy = ( -vtp(i-1,j-1) - vtp(i,j-1) + &
                         ( vtp(i-1,j-2) + vtp(i,j-2) ) * maskB(i,j-2)/4d0 ) / Deltax

                    pnode = ( Pp(i-1,j) + Pp(i,j-1) + Pp(i-1,j-1) + & 
                              Pt(i-1,j) + Pt(i,j-1) + Pt(i-1,j-1) ) / 3d0

                 elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax

                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax

                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax
    
                    pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i-1,j-1) + &
                              Pt(i-1,j) + Pt(i,j) + Pt(i-1,j-1) ) / 3d0

                 elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                    dudx = ( utp(i+1,j) + utp(i+1,j-1) - &
                         ( utp(i+2,j) + utp(i+2,j-1) ) * maskB(i+2,j)/4d0 ) / Deltax

                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax

                    pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i,j-1) + &
                              Pt(i-1,j) + Pt(i,j) + Pt(i,j-1) ) / 3d0
                 else

                    print *, 'wowowo1'
                    stop

                 endif !summaskC .eq. 3
                 
              elseif (summaskC .eq. 2) then !case 7
                     
                 if (maskC(i-1,j) .eq. 0 .and. &
                      maskC(i-1,j-1) .eq. 0) then
! xo
! xo
                    dudy = 0d0
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                    dudx = ( utp(i+1,j) + utp(i+1,j-1) - &
                         ( utp(i+2,j) + utp(i+2,j-1) ) * maskB(i+2,j)/4d0 ) / Deltax
		      
                    dvdy = 0d0
                    pnode = ( Pp(i,j) + Pp(i,j-1) + &
                              Pt(i,j) + Pt(i,j-1) )/2d0

                 elseif(maskC(i,j) .eq. 0 .and. & !case 6
                      maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           
                    dudy = 0d0
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax

                    dvdy = 0d0
                    pnode = ( Pp(i-1,j) + Pp(i-1,j-1) + &
                              Pt(i-1,j) + Pt(i-1,j-1) )/2d0

                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                      maskC(i,j) .eq. 0) then
! xx
! oo
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = 0d0
                    dudx = 0d0
                    dvdy = ( -vtp(i-1,j-1) - vtp(i,j-1) + &
                         ( vtp(i-1,j-2) + vtp(i,j-2) ) * maskB(i,j-2)/4d0 ) / Deltax
                    
                    pnode = ( Pp(i-1,j-1) + Pp(i,j-1) + &
			      Pt(i-1,j-1) + Pt(i,j-1) )/2d0

                 elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                      maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = 0d0
                    dudx = 0d0
                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax
                        
                    pnode = ( Pp(i-1,j) + Pp(i,j) + Pt(i-1,j) + Pt(i,j)) / 2d0

                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 15
                      maskC(i,j-1) .eq. 0) then
! xo                                                                    
! ox                    ! don't do anything (could be improved)                                        
                        
			
                 elseif(maskC(i,j) .eq. 0 .and. & !case 16
                      maskC(i-1,j-1) .eq. 0) then
! ox                                                                      
! xo                    ! don't do anything (could be improved)                                                     
			
                 else

                    print *, 'wowowo2'
                    stop

                 endif  !summaskC .eq. 2

              else
                 print *, 'wowowo3'
                 stop

              endif !if (summaskC .eq. 4) then...
                  

!------------------------------------------------------------------------
!     Shear viscosity calculation at the grid node
!------------------------------------------------------------------------

              if ( rheo .eq. 1 ) then ! ellipse, jfl p.892

                 deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                      + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                      + ell_2 * ( dvdx + dudy ) ** 2 )

                 if ( regularization .eq. 'tanh' ) then

                    deno = max( deno, 1d-20 )
                    etaB(i,j) = ell_2 * ( pnode/denomin ) &
                         *(tanh( denomin*(1/deno)) )

                 elseif ( regularization .eq. 'Kreysher' ) then

                    denoT = deno + denomin
                    etaB(i,j) = ell_2 * ( pnode/denoT )

                 elseif ( regularization .eq. 'capping' ) then
                    
                    denoT = max(deno,denomin)
                    etaB(i,j) = ell_2 * ( pnode/denoT )

                 else

                    print *, 'WRONG REGULARIZATION'
                    stop

                 endif

              elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                 stop
                   
              endif
              
           endif !summaskC .ge. 2 
           
        enddo
     enddo
         
  endif
      
  return
end subroutine ViscousCoeff_method2

!************************************************************************                                  
!     Subroutine ViscousCoeff_method3_and_4: computes the non-linear viscous                                         
!     coefficient at the centers (zetaC, etaC) and at the nodes (etaB) of the                                
!     grid. All the strain rate calculations are second order accurate.                                      
!     The ice strength (Pp) is collocated with zetaC, etaC.
!                                                                                                           
!     We use the method described in Bouillon et al. 2013 (the EVP revisited):    
!                                                                                                            
!     1) ep12 are calculated at the grid nodes (ep12b).                                                      
!     2) ep12c**2 = ( ep12b**2 + ep12b**2 + ep12b**2 + ep12b**2 / 4).                                        
!     3) zetaC = Pp / 2f(ep11,ep22,ep12c**2)                                                                 
!     4) etaB = ( etaC + etaC + etaC + etaC ) / 4                                                            
!                                                                                                            
!     note1 : etaC = e^-2 * zetaC                                                                             
!
!     note2: method 4 is Bouillon's method. method 3 is a slightly modified
!            approach.
!                                                                                                            
!     We also apply the replacement closure of Kreyscher et al 2000 (eq.7).                                   
!                                                                                                            
!     JF Lemieux, 15 Febrary 2013                                                                            
!                                                                                                            
!************************************************************************     


subroutine ViscousCoeff_method3_and_4(utp,vtp)
  use ellipse
  implicit none

  include 'parameter.h'
  include 'CB_DynVariables.h'
  include 'CB_const.h'
  include 'CB_mask.h'
  include 'CB_options.h'

  integer i, j, rheo, summaskC, peri

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
  double precision ep12(0:nx+2,0:ny+2), meanep12sq, meanep12
  
  denomin = 2d-09 ! Hibler, 1979

  rheo = Rheology ! define local variable to speed up the code
  peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions

!------------------------------------------------------------------------
!     free slip boundary condition:
!       d(v_tangential)/d(normal) = 0 & v_normal = 0, at close boundary
!------------------------------------------------------------------------

  if ( BndyCond .eq. 'freeslip' ) then

     print *, 'WARNING: freeslip option not fully tested yet'

!------------------------------------------------------------------------
!     strain-rates calculation
!     no-slip boundary condition :  v_normal = 0 and v_tangential = 0 
!------------------------------------------------------------------------

  elseif ( BndyCond .eq. 'noslip' ) then


     do i = 0, nx+1
        do j = 0, ny+1

           etaC(i,j)  = 0d0
           zetaC(i,j) = 0d0 
           etaB(i,j)  = 0d0
           ep12(i,j)  = 0d0
           
        enddo
     enddo

     dudx       = 0d0
     dvdy       = 0d0
     dudy       = 0d0
     dvdx       = 0d0
     deno       = 0d0


     if (peri .ne. 0) call periodicBC(utp,vtp)
     

!------------------------------------------------------------------------
!     ep12 calculation at the grid node (see p.2-118 PDF notebook)
!------------------------------------------------------------------------

         ! for sig12B (defined at the node)
     do j = 1, ny+1 
        do i = 1, nx+1

           summaskC = maskC(i-1,j) + maskC(i,j) + & 
                maskC(i,j-1) + maskC(i-1,j-1)

           if (summaskC .ge. 2) then
                  
              if (summaskC .eq. 4) then
! oo
! oo normal
                 dudy = ( utp(i,j) - utp(i,j-1) ) / Deltax !case 1 
                 dvdx = ( vtp(i,j) - vtp(i-1,j) ) / Deltax
                     
              elseif (summaskC .eq. 3) then 

                 if (maskC(i-1,j) .eq. 0) then !case 2
! xo
! oo    
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                    
                 elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                 elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax

                 elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                 else

                    print *, 'wowowo1'
                    stop

                 endif !summaskC .eq. 3

              elseif (summaskC .eq. 2) then !case 7

                 if (maskC(i-1,j) .eq. 0 .and. &
                      maskC(i-1,j-1) .eq. 0) then
! xo
! xo
                    dudy = 0d0
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                 elseif(maskC(i,j) .eq. 0 .and. & !case 6
                      maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           
                    dudy = 0d0
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                      maskC(i,j) .eq. 0) then
! xx
! oo
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = 0d0

                 elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                      maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = 0d0

                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 15
                      maskC(i,j-1) .eq. 0) then
! xo                                                                    
! ox
                    
                 elseif(maskC(i,j) .eq. 0 .and. & !case 16
                      maskC(i-1,j-1) .eq. 0) then
! ox                      	                    
! xo                  ! don't do anything (could be improved)                                                     
                 else
                    
                    print *, 'wowowo2'
                    stop

                 endif  !summaskC .eq. 2

              else
                 print *, 'wowowo3'
                 stop

              endif !if (summaskC .eq. 4) then...

              ep12(i,j) = dudy + dvdx
                  
           endif !summaskC .ge. 2 

        enddo
     enddo


!------------------------------------------------------------------------
!     Shear and bulk viscosity calculation at the grid center
!------------------------------------------------------------------------   

     do i = 1, nx
        do j = 1, ny

               
           if ( maskC(i,j) .eq. 1 ) then
                  
              dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
              dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
              if ( rheo .eq. 1 ) then ! ellipse, jfl p.892

                 if (visc_method .eq. 3) then
                     
                    meanep12 = ( ep12(i,j) + ep12(i,j+1)+ & 
                         ep12(i+1,j+1) + ep12(i+1,j) )/4d0

                    deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                         + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                         + ell_2 * (meanep12)**2 )
                    
                 elseif (visc_method .eq. 4) then

                    meanep12sq = ( (ep12(i,j))**2 + (ep12(i,j+1))**2+ & 
                         (ep12(i+1,j+1))**2 + (ep12(i+1,j))**2 )/4d0

                    deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                         + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                         + ell_2 * meanep12sq )

                 endif

                 if ( regularization .eq. 'tanh' ) then

                    deno = max( deno, 1d-20 )
                     
                    zetaC(i,j) =  ( (Pp(i,j)+Pt(i,j))/denomin ) &
                         *( tanh(denomin*(1/deno)))
                         
                    P(i,j) = ((( Pp(i,j)-Pt(i,j))/denomin ) &
                         *tanh(denomin*(1/deno))) *deno 

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin

                    zetaC(i,j) = (Pp(i,j)+Pt(i,j))/denoT 
                    
                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)

                    zetaC(i,j) = (Pp(i,j)+Pt(i,j))/denoT 

                    P(i,j) = ((Pp(i,j) - Pt(i,j)) / denoT) * deno

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                        
                 endif

                 etaC(i,j)  = zetaC(i,j) * ell_2
                     
              elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                 stop
                   
              endif
              
           endif
               
        enddo
     enddo

!------------------------------------------------------------------------
!     Set zetaC and etaC to 0.0 at the open boundaries (see p.32-33 EC-2)
!------------------------------------------------------------------------

     do i = 1, nx+1

        if (maskC(i,0) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1 .and. Periodic_y .eq. 0) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
            
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1 .and. Periodic_x .eq. 0) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1 .and. Periodic_x .eq. 0) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo
     
     if (peri .ne. 0) call periodicBC(etaC,zetaC)
     if (peri .ne. 0) call periodicBC(etaB,P)  !etaB = zero and is a dummy here          

!------------------------------------------------------------------------
!     Shear viscosity calculation at the grid node (see p.2-118 PDF notebook)
!------------------------------------------------------------------------

     ! for sig12B (defined at the node)
     do j = 1, ny+1 
        do i = 1, nx+1
           
           summaskC = maskC(i-1,j) + maskC(i,j) + & 
                maskC(i,j-1) + maskC(i-1,j-1)

           if (summaskC .ge. 2) then

              if (summaskC .eq. 4) then
! oo
! oo normal
                     
                 etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 4d0
                     
                     
              elseif (summaskC .eq. 3) then 
                 
                 if (maskC(i-1,j) .eq. 0) then !case 2
! xo
! oo    
                        
                    etaB(i,j) = ( etaC(i,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 3d0
                    
                 elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           

                    etaB(i,j) = ( etaC(i-1,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 3d0

                 elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
    
                    etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i-1,j-1) ) / 3d0

                 elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      

                    etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i,j-1) ) / 3d0
		      
                 else

                    print *, 'wowowo1'
                    stop

                 endif !summaskC .eq. 3

              elseif (summaskC .eq. 2) then !case 7

                 if (maskC(i-1,j) .eq. 0 .and. &
                      maskC(i-1,j-1) .eq. 0) then
! xo
! xo

                    etaB(i,j) = ( etaC(i,j) + etaC(i,j-1) ) / 2d0

                 elseif(maskC(i,j) .eq. 0 .and. & !case 6
                      maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           

                    etaB(i,j) = ( etaC(i-1,j) + etaC(i-1,j-1) ) / 2d0

                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                      maskC(i,j) .eq. 0) then
! xx
! oo

                    etaB(i,j) = ( etaC(i-1,j-1) + etaC(i,j-1) ) / 2d0

                 elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                      maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                        
                    etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) ) / 2d0  
  
                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 15
                      maskC(i,j-1) .eq. 0) then
! xo                                                                    
! ox                    ! don't do anything (could be improved)                                        
                        
			
                 elseif(maskC(i,j) .eq. 0 .and. & !case 16
                      maskC(i-1,j-1) .eq. 0) then
! ox                                                                      
! xo                    ! don't do anything (could be improved)                                                     
			
                 else

                    print *, 'wowowo2'
                    stop

                 endif  !summaskC .eq. 2

              else
                 print *, 'wowowo3'
                 stop

              endif !if (summaskC .eq. 4) then...
                  
           endif !summaskC .ge. 2 

        enddo
     enddo
         
  endif
      
  return
end subroutine ViscousCoeff_method3_and_4


subroutine etaB_at_open_boundaries
  use ellipse
  implicit none

  include 'parameter.h'
  include 'CB_DynVariables.h'
  include 'CB_mask.h'
  include 'CB_options.h'

  integer i, j, summaskC, peri

!------------------------------------------------------------------------                                                                                                      
!     Set etaB to 0.0 at the open boundaries (see p.32-33 EC-2)                                                                                                                
!------------------------------------------------------------------------  

  peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions                                                                                                    
  
  if (peri .ne. 2) then
    etaB(1,1)=0d0
    etaB(1,ny+1)=0d0
    etaB(nx+1,1)=0d0
    etaB(nx+1,ny+1)=0d0
  endif  

  if (Periodic_x .eq. 0) then
      
      do j = 2, ny

	summaskC = maskC(0,j) + maskC(1,j) + &
	      maskC(1,j-1) + maskC(0,j-1)

	if (summaskC .ge. 2) then

	    if (summaskC .ge. 3) then
	      etaB(1,j) = 0d0
	    elseif (summaskC .eq. 2) then
	      if (maskC(0,j) .eq. 0 .and. maskC(1,j) .eq. 0) then
		  etaB(1,j) = 0d0
	      elseif (maskC(0,j-1) .eq. 0 .and. maskC(1,j-1) .eq. 0) then
		  etaB(1,j) = 0d0
	      elseif (maskC(0,j) .eq. 0 .and. maskC(0,j-1) .eq. 0) then
		  ! don't do anything...it can be non zero                                                                                                                   
	      else
		  print *, 'wowowo' ! ocean on the left is not possible                                                                                                         
	      endif

	    endif

	endif

      enddo

      do j = 2, ny
	    
	summaskC = maskC(nx,j) + maskC(nx+1,j) + &
	      maskC(nx+1,j-1) + maskC(nx,j-1)

	if (summaskC .ge. 2) then

	    if (summaskC .ge. 3) then
	      etaB(nx+1,j) = 0d0
	    elseif (summaskC .eq. 2) then
	      if (maskC(nx,j) .eq. 0 .and. maskC(nx+1,j) .eq. 0) then
		  etaB(nx+1,j) = 0d0
	      elseif (maskC(nx,j-1) .eq. 0 .and. maskC(nx+1,j-1) .eq. 0) then
		  etaB(nx+1,j) = 0d0
	      elseif (maskC(nx+1,j) .eq. 0 .and. maskC(nx+1,j-1) .eq. 0) then
		  ! don't do anything...it can be non zero                                                                                                                   
	      else
		  print *, 'wowowo' ! ocean on the right is not possible                                                                                                        
	      endif

	    endif

	endif

      enddo
      
  endif 

  if (Periodic_y .eq. 0) then  
  
      do i = 2, nx

	summaskC = maskC(i-1,1) + maskC(i,1) + &
	      maskC(i-1,0) + maskC(i,0)

	if (summaskC .ge. 2) then

	    if (summaskC .ge. 3) then
	      etaB(i,1) = 0d0
	    elseif (summaskC .eq. 2) then
	      if (maskC(i-1,1) .eq. 0 .and. maskC(i-1,0) .eq. 0) then
		  etaB(i,1) = 0d0
	      elseif (maskC(i,1) .eq. 0 .and. maskC(i,0) .eq. 0) then
		  etaB(i,1) = 0d0
	      elseif (maskC(i-1,0) .eq. 0 .and. maskC(i,0) .eq. 0) then
		  ! don't do anything...it can be non zero                                                                                                                   
	      else
		  print *, 'wowowo' ! ocean below is not possible                                                                                                               
	      endif

	    endif

	endif

      enddo
      
      do i = 2, nx

	summaskC = maskC(i-1,ny+1) + maskC(i,ny+1) + &
	      maskC(i-1,ny) + maskC(i,ny)
	
	if (summaskC .ge. 2) then

	    if (summaskC .ge. 3) then
	      etaB(i,ny+1) = 0d0
	    elseif (summaskC .eq. 2) then
	      if (maskC(i-1,ny+1) .eq. 0 .and. maskC(i-1,ny) .eq. 0) then
		  etaB(i,ny+1) = 0d0
	      elseif (maskC(i,ny+1) .eq. 0 .and. maskC(i,ny) .eq. 0) then
		  etaB(i,ny+1) = 0d0
	      elseif (maskC(i-1,ny+1) .eq. 0 .and. maskC(i,ny+1) .eq. 0) then
		  ! don't do anything...it can be non zero                                                                                                                      
	      else
		  print *, 'wowowo' ! ocean above is not possible                                                                                                               
	      endif
	      
	    endif

	endif

      enddo
      
  endif   
      
  return
end subroutine etaB_at_open_boundaries







subroutine MEBcoeff

!************************************************************************                                  
!     Subroutine MEBcoeff: computes the MEB coefficients in the same format as                                         
!     those of the VP model, at the centers (zetaC, etaC) and at the nodes (etaB)                               
!     of the grid. These are defined as (See Plante et al. 2019, The Cryosphere):
!     
!     zetaC = Gamma * E * (C1+C2) * Deltat
!     etaC =  Gamma * E * Deltat * C3
!                        
!     where : 
!     
!        Gamma = [ 1 + Deltat / lambda   ] ^(-1)    is a viscous dissipation factor
!        lambda = [ lambda0 * (1-d)^(1-alpha) ] / h * e^(-c(1-A))   
!             : is the viscous relaxation time as function of d, h and A
!
!        with 

!        E = Y * (1-d) * h * e^(-c(1-A))      is the Elastic Stiffness, 
!                                      as a function of Y (young modulus), d, h, A
!        C1 = 1 / (1 - poisson^2)
!        C2 = poisson / (1 - poisson^2)
!        C3 = (1-poisson) / (1 - poisson^2)
!          : are three elastic constants (a version of the Lame coefficients)
!
!    In this code, we include the Elastic stiffness in the other coefficients,
!    to decrease the number of variables. So we set:
!
!        Lame1 = Y * C3
!        Lame2 = Y * ( C1 + C2 )
!        hAfunc =  h * e^(-c(1-A)) 
!
!
!
!     Mathieu Plante, October 9 2019                                                                            
!                                                                                                            
!************************************************************************     

     use elastic
     use ellipse !we need this for the coefficient C
     implicit none

     include 'parameter.h'
     include 'CB_DynVariables.h'
     include 'CB_mask.h'
     include 'CB_options.h'
     include 'CB_const.h'

     integer i, j, peri

     double precision Lame1, Lame2,m1, m2, m3, m4
     double precision hAfunc(0:nx+2,0:ny+2), hAfuncB(0:nx+2,0:ny+2)
     peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions      
      
!---------------------------------------------------------  
!     set constants for for linear elastic Rheology (Matt, 2015)
!---------------------------------------------------------

      Lame1     =  Young/(2d0*(1d0 + Poisson)) ! Shear Modulus of sea ice
      Lame2     =  2d0*Poisson*Lame1/(1d0-Poisson) ! elastic const of sea ice
      Tdam      = max(Deltax /sqrt(Lame2/rhoice), Deltat) !Damage relaxation time scale

!------------------------------------------------------------------------
!     free slip boundary condition:
!       d(v_tangential)/d(normal) = 0 & v_normal = 0, at close boundary
!------------------------------------------------------------------------
 

      if ( BndyCond .eq. 'freeslip' ) then

               print *, 'WARNING: freeslip option not fully tested yet'

      elseif ( BndyCond .eq. 'noslip' ) then
      
         do i = 0, nx+1
            do j = 0, ny+1

               etaC(i,j)  = Lame1*Deltat
               zetaC(i,j) = (Lame1 + Lame2)*Deltat
               etaB(i,j)  = Lame1*Deltat
               hAfunc(i,j) = 1d0
               hAfuncB(i,j) = 1d0
            enddo
         enddo

         ! Apply periodic boundary conditions 
         if (peri .ne. 0) call periodicBC(dam,dfactor)  
         if (peri .ne. 0) call periodicBC(damB,dfactorB)  
         ! ---
               
!------------------------------------------------------------------------
!     Calculate the coefficient in the grid center
!------------------------------------------------------------------------

         do i = 1, nx
           do j = 1, ny

             if (IMEX .eq. 0) dfactor(i,j) = 1d0 !If IMEX=1, the following coefficients are computed 
                                                 !using the implicit damage : d = d*dfactor

             if (maskC(i,j) .eq. 1d0) then

                hAfunc(i,j)  = h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
                
                GammaMEB(i,j) = 1d0 / (1d0 + (Deltat * hAfunc(i,j)) / &
                           (lambda0 * ((dam(i,j))*dfactor(i,j))**(alpha-1d0))) 

             else
             
                hAfunc(i,j) = Lame1*Deltat                
                GammaMEB(i,j) = 1d0
                
             endif

             etaC(i,j)  = Lame1*hAfunc(i,j)*(dam(i,j))*dfactor(i,j)* &   !Elastic stiffness
                                                 Deltat * GammaMEB(i,j)
             zetaC(i,j) = hAfunc(i,j)*(Lame1 + Lame2)*Deltat * &
                                 dam(i,j)* dfactor(i,j) * GammaMEB(i,j)
             P(i,j) = 0d0
             
             CoheC(i,j) = Cohe
             sigtC(i,j) = sigt
             sigcC(i,j) = sigc

            enddo
         enddo
         
         if (peri .ne. 0) call periodicBC(etaC,zetaC)
         if (peri .ne. 0) call periodicBC(GammaMEB,P)
         
!------------------------------------------------------------------------
!     Calculate the coefficient in the grid nodes
!------------------------------------------------------------------------         
         
         do i = 1, nx+1
            do j = 1, ny+1 

               if (IMEX .eq. 0) dfactorB(i,j) = 1d0

               m1 = maskC(i,j) 
               m2 = maskC(i-1,j)
               m3 = maskC(i,j-1)
               m4 = maskC(i-1,j-1)

               !open bc
               if (Periodic_x .eq. 0) then
                   if (i .eq. 1) then
                      m2 = 0d0
                      m4 = 0d0
                   elseif (i .eq. nx+1) then
                      m1 = 0d0
                      m3 = 0d0  
                   endif
               endif                   
                   
               if (Periodic_y .eq. 0) then
                   if (j .eq. 1) then
                      m3 = 0d0
                      m4 = 0d0
                   elseif (j .eq. ny+1) then
                      m1 = 0d0
                      m2 = 0d0  
                   endif                  
               endif    

               if (m1+m2+m3+m4 .ne. 0d0) then

                  hAfuncB(i,j) = ((h(i,j)*m1 + h(i-1,j)*m2 + h(i,j-1)*m3 &
                                  + h(i-1,j-1)*m4 )/ (m1+m2+m3+m4))  &
                             * dexp(-C * ( 1d0 - ((A(i,j)*m1 + A(i-1,j)*m2 &
                              + A(i,j-1)*m3 + A(i-1,j-1)*m4 )/ (m1+m2+m3+m4) ) ) )
               endif

               GammaMEB_B(i,j) = 1d0 / ( 1d0 + (Deltat * hAfuncB(i,j)) / &
                              (lambda0 * (damB(i,j)*dfactorB(i,j))**(alpha-1d0) )  )
               etaB(i,j)  = hAfuncB(i,j)*Lame1*Deltat*&
                             (damB(i,j)*dfactorB(i,j)) *GammaMEB_B(i,j)
            enddo
         enddo
         
!------------------------------------------------------------------------ 

         if (peri .ne. 0) call periodicBC(etaB,GammaMEB_B)    
         
      endif

  return
end subroutine MEBcoeff
