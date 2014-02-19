
      subroutine EVP_SOLVER(tstep)
        use ellipse
        use numerical_EVP
!------------------------------------------------------------------------    
! EVP solver based on Hunke, 2001 (see CICEDOC for details)
! There are differences with Hunke 2001. 
! See p.2-152...2-156 (postdoc notebook)
! 
! author: JF Lemieux
! date: Jan 28 2011
!------------------------------------------------------------------------ 

      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'

      integer i, j, s
      integer, intent(in) :: tstep

      double precision crap(nvar)
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision, save :: sig12B(0:nx+2,0:ny+2)
      double precision, save :: sig1(0:nx+1,0:ny+1), sig2(0:nx+1,0:ny+1) 
      !sig1=sig11+sig22, sig2=sig11-sig22

      double precision uavg, vavg
      double precision hvert, Deltate, T, left1, left2, B1, B2, D

!------------------------------------------------------------------------
!  Define EVP parameters for the solver
!------------------------------------------------------------------------      

      T = Eo*Deltat 
      Deltate = Deltat / Nsub

!------------------------------------------------------------------------
!  Define constants
!------------------------------------------------------------------------

      left1 = (1d0/Deltate) + (0.5d0/T)
      left2 = (1d0/Deltate) + (0.5d0/T)*(ell2)

!------------------------------------------------------------------------
!  Subcycling loop
!------------------------------------------------------------------------
      
      do s = 1, Nsub

!------------------------------------------------------------------------ 
!  Set u^s-1 = uice, v^s-1 = vice
!  The initial iterate is the previous time step solution
!------------------------------------------------------------------------     

         utp   = uice ! utp is used as u^s-1 while uice = u^s
         vtp   = vice ! vtp is used as v^s-1 while vice = v^s   

!------------------------------------------------------------------------
!  Calculate eta, zeta, CdwC1 and CdwC2 as a function of u^s-1 and v^s-1
!------------------------------------------------------------------------ 

         call ViscousCoefficient(utp,vtp) ! there are a lot of strain rates calc 
         call bvect (utp,vtp,crap)
!        here that are also cal below...improve this!!!!!!!!!!!!

!------------------------------------------------------------------------
!     Advance stresses in time (from s-1 to s)
!------------------------------------------------------------------------

         call calc_sigs(utp, vtp, sig1, sig2, sig12B, &
                        left1, left2, Deltate, T, tstep, s)
 
!------------------------------------------------------------------------
!     Advance u velocity in time (from s-1 to s)
!------------------------------------------------------------------------

         do j = 1, ny
            do i = 2, nx

               if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

!------------------------------------------------------------------------
!     tair + term of sea surface tilt + part of water drag                  
!------------------------------------------------------------------------  

                  B1 = bu(i,j) ! calc in bvect (watchout bu here is not the
                               ! same as for the VP model (solver 1 or 2)
                  
!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------  

                  vavg = ( vtp(i,j)   + vtp(i,j+1)  &
                         + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

                  hvert = ( h(i,j) + h(i-1,j) ) / 2d0

                  B1 = B1 + rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*du/dte (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

                  D  = ( rhoice * hvert ) / Deltate

                  B1 = B1 + ( rhoice * hvert * utp(i,j) ) / Deltate

!------------------------------------------------------------------------
!    rhoice*h*du/dt (extra term to match the implicit VP solution)
!------------------------------------------------------------------------
 
!                  D  = D + ( rhoice * hvert ) / Deltat

!                  B1 = B1 + ( rhoice * hvert * un1(i,j) ) / Deltat

!------------------------------------------------------------------------
!     water drag (remaining terms)
!------------------------------------------------------------------------

                  D  = D + CdwC1(i,j) * costheta_w 

                  B1 = B1 + CdwC1(i,j) * sintheta_w * vavg 

!------------------------------------------------------------------------
!     dsig11 / dx where sig11 = (sig1 + sig2) / 2
!------------------------------------------------------------------------

                  B1 = B1 + ( ( sig1(i,j) + sig2(i,j) ) - &
                              ( sig1(i-1,j) + sig2(i-1,j) ) ) / (2d0*Deltax)

!------------------------------------------------------------------------
!     dsig12B / dy
!------------------------------------------------------------------------

                  B1 = B1 + ( sig12B(i,j+1) - sig12B(i,j) ) /  Deltax

!------------------------------------------------------------------------
!     update u component
!------------------------------------------------------------------------

                  uice(i,j) = B1 / D ! uice is u^s and utp is u^s-1

               endif
               
            enddo
         enddo

!------------------------------------------------------------------------
!     Apply boundary cond for u component (dudx=0 for open bc)
!     check p.2-106 PDF notebook 2
!------------------------------------------------------------------------
         
         do j = 1, ny
         
            if ( maskB(1,j) + maskB(1,j+1) .gt. 0 ) then

               uice(1,j) = ( 4d0 * uice(2,j)  - uice(3,j) ) / 3d0

            endif

            if ( maskB(nx+1,j) + maskB(nx+1,j+1) .gt. 0 ) then
            
               uice(nx+1,j) = ( 4d0 * uice(nx,j) - uice(nx-1,j))/ 3d0
            
            endif
         
         enddo

!------------------------------------------------------------------------
!     now the v-component
!------------------------------------------------------------------------

         do j = 2, ny
            do i = 1, nx

               if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then 

!------------------------------------------------------------------------
!     tair + term of sea surface tilt + part of water drag
!------------------------------------------------------------------------ 

                  B2 = bv(i,j) ! calc in bvect (watchout bv here is not the
                               ! same as for the VP model (solver 1 or 2)

!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------

                  uavg  = ( utp(i,j)   + utp(i+1,j) &
                          + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

                  hvert = ( h(i,j) + h(i,j-1) ) / 2d0

                  B2 = B2 - rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*dv/dte (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

                  D  = ( rhoice * hvert ) / Deltate

                  B2 = B2 + ( rhoice * hvert * vtp(i,j) ) / Deltate

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (extra term to match the implicit VP solution)
!------------------------------------------------------------------------ 
                  
!                  D  = D + ( rhoice * hvert ) / Deltat

!                  B2 = B2 + ( rhoice * hvert * vn1(i,j) ) / Deltat 

!------------------------------------------------------------------------
!     water drag (remaining terms)
!------------------------------------------------------------------------

                  D  = D  + CdwC2(i,j) * costheta_w 
           
                  B2 = B2 - CdwC2(i,j) * sintheta_w * uavg

!------------------------------------------------------------------------
!     dsig22 / dy where sig22 = (sig1 - sig2) / 2   
!------------------------------------------------------------------------

                  B2 = B2 + ( ( sig1(i,j) - sig2(i,j) ) - &
                              ( sig1(i,j-1) - sig2(i,j-1) ) ) / (2d0*Deltax)
               
!------------------------------------------------------------------------
!     dsig12B / dy
!------------------------------------------------------------------------

                  B2 = B2 + ( sig12B(i+1,j) - sig12B(i,j) ) /  Deltax

                  vice(i,j) = B2 / D

               endif

            enddo
         enddo

!------------------------------------------------------------------------
!     Apply boundary cond for v component (dvdy=0 for open bc)
!------------------------------------------------------------------------

         do i = 1, nx

            if ( maskB(i,1) + maskB(i+1,1) .gt. 0 ) then

               vice(i,1) = ( 4d0 * vice(i,2) - vice(i,3) ) / 3d0

            endif

            if ( maskB(i,ny+1) + maskB(i+1,ny+1) .gt. 0 ) then

               vice(i,ny+1) = ( 4d0 * vice(i,ny) - vice(i,ny-1) )/3d0

            endif

         enddo

      enddo ! end of subcycling loop
      
      return
    end subroutine EVP_SOLVER


    subroutine calc_sigs(utp, vtp, sig1, sig2, sig12B, &
                         left1, left2, Deltate, T, tstep, s)
      use ellipse
      use numerical_EVP

! subroutine to calculte the stresses for the EVP solver 
! and soon for the precond_evp...

      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynVariables.h' ! to get the pressure

      integer i, j, summaskC
      integer, intent(in) :: tstep, s

      double precision, intent(in) :: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision, intent(inout) :: sig12B(0:nx+2,0:ny+2)
      double precision, intent(inout) :: sig1(0:nx+1,0:ny+1)
      double precision, intent(inout) :: sig2(0:nx+1,0:ny+1)
      !sig1=sig11+sig22, sig2=sig11-sig22  

      double precision, intent(in) :: left1, left2, Deltate, T
                             
      double precision dudx, dvdy, dudy, dvdx 


    ! sig1 and sig2 (therefore sig11 and sig22) are define at the tracer point
 
         do i = 1, nx  ! is it ok to loop up to here?                    
            do j = 1, ny
               if ( maskC(i,j) .eq. 1 ) then

                  dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax ! tracer point 
                                                                
                  dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax ! tracer point     

                  if (tstep .eq. 1 .and. s .eq. 1) then 
                                      !get sig(s-1) if s=1 and tstep=1

                     if (init_stress .eq. 'VP') then

                        sig1(i,j) = 2d0 * zetaC(i,j) * ( dudx+dvdy ) &
                                   -2d0 * P(i,j) !recall that P is in fact P/2 in code

                        sig2(i,j) = 2d0 * etaC(i,j) * ( dudx-dvdy )

                     elseif (init_stress .eq. 'zero') then

                        sig1(i,j) = 0d0
                        sig2(i,j) = 0d0

                     endif

                  endif

                     sig1(i,j) = ( sig1(i,j) / Deltate -P(i,j)/T+ &     
                                 ( zetaC(i,j)/T ) * ( dudx+dvdy ) ) / left1

                     sig2(i,j) = ( sig2(i,j) / Deltate + &
                                 ( zetaC(i,j)/T ) * ( dudx-dvdy ) ) / left2

               endif

            enddo
         enddo


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
                        dudy = 0d0 ! ideally we should not calc sig12
                        dvdx = 0d0 ! for case 15      

                     elseif(maskC(i,j) .eq. 0 .and. & !case 16
                            maskC(i-1,j-1) .eq. 0) then
! ox                                                                      
! xo                                                                         
                        dudy = 0d0 ! ideally we should not calc sig12
                        dvdx = 0d0 ! for case 16
                     
                     else

                     print *, 'wowowo2'
                     stop

                     endif  !summaskC .eq. 2

                  else
                     print *, 'wowowo3'
                     stop

                  endif !if (summaskC .eq. 4) then...
                                                                           
! we don't have access to zetaB..this is why we use etaB instead
 
                  if (tstep .eq. 1 .and. s .eq. 1) then 
                                         ! to get sig(s-1) if tstep=1 and s=1
                                          
                     if (init_stress .eq. 'VP') then

                        sig12B(i,j) = etaB(i,j) * ( dudy+dvdx )

                     elseif (init_stress .eq. 'zero') then
                        

                        sig12B(i,j) = 0d0

                     endif


                  endif

                     sig12B(i,j)= ( sig12B(i,j)/ Deltate + &
                     ( etaB(i,j)*(ell2)/(2d0*T) ) * ( dudy+dvdx ) )/left2

               endif !summaskC .ge. 2 

            enddo
         enddo

    end subroutine calc_sigs
      



