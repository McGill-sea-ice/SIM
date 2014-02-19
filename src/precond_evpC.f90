
      subroutine PRECOND_EVPC( rhs, wk2 )
        use ellipse
!--------------------------------------------------------------------------
! Preconditioner: calculates x = P^-1rhs with the EVP. Note that rhs is wk1. 
! The initial guess is set to 0. The solution x is then put in wk2.
! In this approach, eta, zeta and Cdw are held at u^k-1,v^k-1. Moreover,
! f=0 and theta_w = 0 in the precondioner.See PDF notebook p2-116...
!--------------------------------------------------------------------------

      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'

      integer i, j, s, N_sub, summaskC
      
      double precision rhs(nvar), wk2(nvar)
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision rhsu(0:nx+2,0:ny+2), rhsv(0:nx+2,0:ny+2)

      double precision sig12B(0:nx+2,0:ny+2)
      double precision sig1(0:nx+1,0:ny+1), sig2(0:nx+1,0:ny+1) 
      !2*sig1=sig11+sig22, 2*sig2=sig11-sig22

      double precision dudx,dvdy,dudy,dvdx
      double precision hvert, Deltate, T, left1, left2, B1, B2, D

      call transformer(rhsu,rhsv,rhs,0)

!------------------------------------------------------------------------
!  Define EVP parameters for the preconditioner
!------------------------------------------------------------------------      

      T = 0.36d0*Deltat
      N_sub=10
      Deltate = 10d0

!------------------------------------------------------------------------
!  Define constants
!------------------------------------------------------------------------

      left1 = (1d0/Deltate) + (0.5d0/T)
      left2 = (1d0/Deltate) + (0.5d0/T)*(ell2)

!------------------------------------------------------------------------
!  Set initial guess for preconditioner
!------------------------------------------------------------------------

            utp   = 0d0
            vtp   = 0d0

            sig12B= 0d0
            sig1  = 0d0
            sig2  = 0d0

!------------------------------------------------------------------------
!  Subcycling loop
!------------------------------------------------------------------------

      do s = 1, N_sub

!------------------------------------------------------------------------
!     Advance sig* in time (from s-1 to s). Recall that in the precond
!     we use the sig* not the regular stresses (see p.2-126)
!------------------------------------------------------------------------

!         call calc_sigs(utp, vtp, sig1, sig2, sig12B, &
!                        left1, left2, Deltate, T,s) ! this subroutine is
                                                    ! in evp_solver.f90

! I cant use directly the subroutine above because we should use etaf...
! for the precond but use eta for the solver...needs some work.

! sig1 and sig2 (therefore sig11 and sig22) are define at the tracer point

         do i = 1, nx  ! is it ok to loop up to here?
            do j = 1, ny

               if ( maskC(i,j) .eq. 1 ) then
                  
                  dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax ! tracer point
                  dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax ! tracer point

                  sig1(i,j) = ( sig1(i,j) / Deltate + &
                              ( zetaCf(i,j)/T ) * ( dudx+dvdy ) ) / left1

                  sig2(i,j) = ( sig2(i,j) / Deltate + &
                              ( zetaCf(i,j)/T ) * ( dudx-dvdy ) ) / left2

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

! we don't have access to zetaBf..this is why we use etaBf instead
                  
                  sig12B(i,j)= ( sig12B(i,j)/ Deltate + &
                  ( etaBf(i,j)*(ell2)/(2d0*T) ) * ( dudy+dvdx ) )/left2
! could be improved here

               endif !summaskC .ge. 2
                     
            enddo
         enddo

!------------------------------------------------------------------------
!     Set sig12B to 0.0 at the proper open boundaries
!------------------------------------------------------------------------
      
! I think they are already zero but just to be sure

         do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then
               sig12B(i,1)  = 0.0d0
            endif

            if (maskC(i,ny+1) .eq. 1) then             
               sig12B(i,ny+1)  = 0.0d0
            endif
 
         enddo

         do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then   
               sig12B(1,j)  = 0.0d0
            endif

            if (maskC(nx+1,j) .eq. 1) then  
               sig12B(nx+1,j)  = 0.0d0
            endif

         enddo

!------------------------------------------------------------------------
!     Advance u velocity in time
!------------------------------------------------------------------------

         do j = 1, ny
            do i = 2, nx

               if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

                  B1   = rhsu(i,j)
                  
!------------------------------------------------------------------------
!     Coriolis term (f = 0 in the precond)
!------------------------------------------------------------------------  

!            vavg = ( vice(i,j)   + vice(i,j+1)                           &
!                   + vice(i-1,j) + vice(i-1,j+1) ) / 4d0

                  hvert = ( h(i,j) + h(i-1,j) ) / 2d0

!            B1 = B1 + rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*u/Deltate
!------------------------------------------------------------------------

                  D  = ( rhoice * hvert ) / Deltate

!------------------------------------------------------------------------
!    rhoice*h*du/dt (extra term to match the implicit VP solution)
!    (the rhoice*h*u^t-1/Deltat part is already included in rhs)
!    After many subcycles, the rhoice*hdu/dte term disapears and we 
!    are left with rhoice*h*u/Deltat only.
!------------------------------------------------------------------------

!                  D  = D + ( rhoice * hvert ) / Deltat

!                  B1 = B1 + ( rhoice * hvert * utp(i,j) ) / Deltate

!------------------------------------------------------------------------
!     water drag  (theta_w is considered to be zero in the precond)
!------------------------------------------------------------------------

                  D  = D + CdwC1f(i,j) !* costheta_w 

!            B1 = B1 + CdwC1f(i,j) * sintheta_w * vavg 

!------------------------------------------------------------------------
!     dsig11 / dx where sig11 = (sig1 + sig2) / 2
!------------------------------------------------------------------------

!                  B1 = B1 + (sig11(i,j) - sig11(i-1,j))/Deltax

                  B1 = B1 + ( ( sig1(i,j) + sig2(i,j) ) - &
                              ( sig1(i-1,j) + sig2(i-1,j) ) ) / (2d0*Deltax)

!------------------------------------------------------------------------
!     dsig12B / dy
!------------------------------------------------------------------------

                  B1 = B1 + ( sig12B(i,j+1) - sig12B(i,j) ) /  Deltax

!------------------------------------------------------------------------
!     update u component
!------------------------------------------------------------------------

                  utp(i,j) = B1 / D

               endif

            enddo
         enddo

!------------------------------------------------------------------------
!     Apply boundary cond for u component (dudx=0 for open bc)
!     check p.2-106 PDF notebook 2
!------------------------------------------------------------------------
         
         do j = 1, ny
         
            if ( maskB(1,j) + maskB(1,j+1) .gt. 0 ) then

               utp(1,j) = ( rhsu(1,j) + 4d0 * utp(2,j)  - utp(3,j) ) / 3d0

            endif

            if ( maskB(nx+1,j) + maskB(nx+1,j+1) .gt. 0 ) then
            
               utp(nx+1,j) = ( rhsu(nx+1,j) + 4d0 * utp(nx,j) - &
                               utp(nx-1,j))/ 3d0
            
            endif
         
         enddo


!------------------------------------------------------------------------
!     now the v-component
!------------------------------------------------------------------------

         do j = 2, ny
            do i = 1, nx

               if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then 

                  B2 = rhsv(i,j)

!------------------------------------------------------------------------
!     Coriolis term (f = 0 in the precond)
!------------------------------------------------------------------------

!            uavg  = ( uice(i,j)   + uice(i+1,j)                          &
!                    + uice(i,j-1) + uice(i+1,j-1) ) / 4d0

                  hvert = ( h(i,j) + h(i,j-1) ) / 2d0

!            B2 = B2 - rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*u/Deltate
!------------------------------------------------------------------------

                  D  = ( rhoice * hvert ) / Deltate

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (extra term to match the implicit VP solution)
!    (the rhoice*h*v^t-1/Deltat part is already included in rhs)
!    After many subcycles, the rhoice*hdv/dte term disapears and we
!    are left with rhoice*h*v/Deltat only.
!------------------------------------------------------------------------

!                  D  = D + ( rhoice * hvert ) / Deltat                        

!                  B2 = B2 + ( rhoice * hvert * vtp(i,j) ) / Deltate

!------------------------------------------------------------------------
!     water drag (theta_w is considered to be zero in the precond)
!------------------------------------------------------------------------

                  D  = D  + CdwC2f(i,j) !* costheta_w 
           
!            B2 = B2 - CdwC2f(i,j) * sintheta_w * uavg

!------------------------------------------------------------------------
!     dsig22 / dy where sig22 = (sig1 - sig2) / 2   
!------------------------------------------------------------------------

!                  B2 = B2 + ( sig22(i,j) - sig22(i,j-1) ) /  Deltax

                  B2 = B2 + ( ( sig1(i,j) - sig2(i,j) ) - &
                              ( sig1(i,j-1) - sig2(i,j-1) ) ) / (2d0*Deltax)
               
!------------------------------------------------------------------------
!     dsig12B / dy
!------------------------------------------------------------------------

                  B2 = B2 + ( sig12B(i+1,j) - sig12B(i,j) ) /  Deltax

                  vtp(i,j) = B2 / D

               endif

            enddo
         enddo

!------------------------------------------------------------------------
!     Apply boundary cond for v component (dvdy=0 for open bc)
!------------------------------------------------------------------------

         do i = 1, nx

            if ( maskB(i,1) + maskB(i+1,1) .gt. 0 ) then

               vtp(i,1) =  ( rhsv(i,1) + 4d0 * vtp(i,2)  - vtp(i,3) ) / 3d0

            endif

            if ( maskB(i,ny+1) + maskB(i+1,ny+1) .gt. 0 ) then

               vtp(i,ny+1) = ( rhsv(i,ny+1) + 4d0 * vtp(i,ny) - &
                               vtp(i,ny-1) )/3d0

            endif

         enddo

      enddo ! end of subcycling loop

      call transformer (utp,vtp,wk2,1)

      return
    end subroutine PRECOND_EVPC







