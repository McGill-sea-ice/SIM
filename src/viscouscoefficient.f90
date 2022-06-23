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

  if (visc_method .eq. 1) then
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

  integer i, j, rheo

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

  denomin = 2d-09 ! Hibler, 1979

  rheo = Rheology ! define local variable to speed up the code

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
                    zetaC(i,j) =  ( Pp(i,j)/denomin ) &
                         *( tanh(denomin*(1/deno)))

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin
                    zetaC(i,j) = Pp(i,j)/denoT 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)
                    zetaC(i,j) = Pp(i,j)/denoT

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                        
                 endif

                 P(i,j) = zetaC(i,j)*deno ! replacement pressure 
                 etaC(i,j)  = zetaC(i,j) * ell_2
                     
              elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                 stop
                   
              endif
                  
           endif
               
        enddo
     enddo

     do i = 1, nx+1

        if (maskC(i,0) .eq. 1) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
            
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo

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
  
  integer i, j, rheo, summaskC

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin, pnode
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
  
  denomin = 2d-09 ! Hibler, 1979

  rheo = Rheology ! define local variable to speed up the code

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
                    zetaC(i,j) =  ( Pp(i,j)/denomin ) &
                         *( tanh(denomin*(1/deno)))

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin
                    zetaC(i,j) = Pp(i,j)/denoT 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)
                    zetaC(i,j) = Pp(i,j)/denoT

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                    
                 endif

                 P(i,j) = zetaC(i,j)*deno ! replacement pressure 
                 etaC(i,j)  = zetaC(i,j) * ell_2

                 if (Damage) then
                    damSS(i,j) = 1 - tanh(denomin*(1/deno)) ! damage steady state solution
                    
                 endif
                 
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

        if (maskC(i,0) .eq. 1) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
        
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo

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
                     
                 pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i,j-1) + Pp(i-1,j-1) ) / 4d0 
                     
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
                        
                    pnode = ( Pp(i,j) + Pp(i,j-1) + Pp(i-1,j-1) ) / 3d0 ! check ca 

                 elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax
		       
                    dvdy = ( -vtp(i-1,j-1) - vtp(i,j-1) + &
                         ( vtp(i-1,j-2) + vtp(i,j-2) ) * maskB(i,j-2)/4d0 ) / Deltax

                    pnode = ( Pp(i-1,j) + Pp(i,j-1) + Pp(i-1,j-1) ) / 3d0

                 elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax

                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax

                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax
    
                    pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i-1,j-1) ) / 3d0

                 elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                    dudx = ( utp(i+1,j) + utp(i+1,j-1) - &
                         ( utp(i+2,j) + utp(i+2,j-1) ) * maskB(i+2,j)/4d0 ) / Deltax

                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax

                    pnode = ( Pp(i-1,j) + Pp(i,j) + Pp(i,j-1) ) / 3d0

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
                    pnode = ( Pp(i,j) + Pp(i,j-1) )/2d0

                 elseif(maskC(i,j) .eq. 0 .and. & !case 6
                      maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           
                    dudy = 0d0
                    dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                    dudx = ( -utp(i-1,j) - utp(i-1,j-1) + &
                         ( utp(i-2,j) + utp(i-2,j-1) ) * maskB(i-2,j)/4d0 ) / Deltax

                    dvdy = 0d0
                    pnode = ( Pp(i-1,j) + Pp(i-1,j-1) )/2d0

                 elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                      maskC(i,j) .eq. 0) then
! xx
! oo
                    dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                    dvdx = 0d0
                    dudx = 0d0
                    dvdy = ( -vtp(i-1,j-1) - vtp(i,j-1) + &
                         ( vtp(i-1,j-2) + vtp(i,j-2) ) * maskB(i,j-2)/4d0 ) / Deltax
                    
                    pnode = ( Pp(i-1,j-1) + Pp(i,j-1) )/2d0

                 elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                      maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                    dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                    dvdx = 0d0
                    dudx = 0d0
                    dvdy = ( vtp(i-1,j+1) + vtp(i,j+1) - &
                         ( vtp(i-1,j+2) + vtp(i,j+2) ) * maskB(i,j+2)/4d0 ) / Deltax
                        
                    pnode = ( Pp(i-1,j) + Pp(i,j) ) / 2d0

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

  integer i, j, rheo, summaskC

  double precision dudx, dvdy, dudy, dvdx, deno, denoT, denomin
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
  double precision ep12(0:nx+2,0:ny+2), meanep12sq, meanep12
  
  denomin = 2d-09 ! Hibler, 1979

  rheo = Rheology ! define local variable to speed up the code

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
                    zetaC(i,j) =  ( Pp(i,j)/denomin ) &
                         *( tanh(denomin*(1/deno)))

                 elseif ( regularization .eq. 'Kreyscher' ) then
                        
                    denoT = deno + denomin
                    zetaC(i,j) = Pp(i,j)/denoT 

                 elseif ( regularization .eq. 'capping' ) then

                    denoT = max(deno,denomin)
                    zetaC(i,j) = Pp(i,j)/denoT

                 else
                        
                    print *, 'WRONG REGULARIZATION'
                    stop
                        
                 endif

                 P(i,j) = zetaC(i,j)*deno ! replacement pressure 
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

        if (maskC(i,0) .eq. 1) then
           etaC(i,1)  = 0d0
           zetaC(i,1) = 0d0
        endif
            
        if (maskC(i,ny+1) .eq. 1) then
           etaC(i,ny)  = 0d0
           zetaC(i,ny) = 0d0
        endif
            
     enddo
         
     do j = 1, ny+1
            
        if (maskC(0,j) .eq. 1) then   
           etaC(1,j)  = 0d0
           zetaC(1,j) = 0d0
        endif
            
        if (maskC(nx+1,j) .eq. 1) then  
           etaC(nx,j)  = 0d0
           zetaC(nx,j) = 0d0
        endif

     enddo

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

  integer i, j, summaskC

!------------------------------------------------------------------------                                                                                                      
!     Set etaB to 0.0 at the open boundaries (see p.32-33 EC-2)                                                                                                                
!------------------------------------------------------------------------                                                                                                      

  etaB(1,1)=0d0
  etaB(1,ny+1)=0d0
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

  etaB(nx+1,1)=0d0
  etaB(nx+1,ny+1)=0d0
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

  return
end subroutine etaB_at_open_boundaries



    



