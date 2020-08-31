!****************************************************************************
!     subroutine UVsolveSOR: 
!      -Calculates the ice velocity field from the momentum equations 
!       (without the acceleration term) using the pressure field and 
!       non-linear coefficients of friction from the previous time step. 
!      -The water drag term is also linearized around the previous time
!       velocity field. 
!      -The resulting set of equations ([A] {u_i,j} = R') is linear and is
!       solved using Simultaneous Over Relaxation (SOR) Technique. When the
!       viscous term is not present, under-relaxationis necessary. 
!
!       When the ice velocity is small, the matrix A is no longer diagonally
!       dominant and the correction on u and p are less accurate 
!       (the off-diagonal term are neglected in the calculation of u' and p'). 
!       This leads to a slow convergence between the solution of Psolve and
!       UVsolve. 
!
!      -In this routine, I use the off-diagonal term from the previous time
!       step.
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             14-05-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************


      subroutine UVsolveSOR ()
        use numerical_VP
      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      integer i, j, uflag, vflag, nmax, kit
      integer ierroru, jerroru, ierrorv, jerrorv
      
      double precision tol
      double precision error, residual, D, B1, B2
      double precision hvert, uavg, vavg

      tol    = 1d-6     !1d-06  jfl watchout
      nmax   = 20000
      kit    = 0
      uflag  = 0
      vflag  = 0

!------------------------------------------------------------------------
!     Set the over(under)-relaxation parameter:
!      -For eta = zeta = zero (freedrift and cavitating), the A matrix 
!       is not always diagonally dominant and so underrelaxation is 
!       necessary ( w =< 0.1 is slow but stable, w > 0.1 not stable )
!      -For eta, zeta /= 0 (viscous-plastic), the A matrix is diagonally 
!       dominant and overrelaxation is used ( w ~= 1.4 -- 1.68, depending 
!       on how large p is).
!      -Note: Drag/Coriolis ~= rhow Cdw Ui^2 / rhoi h f Ui ~= 20 Ui
!             So for Ui << 5cm/s, the off-diagonal term dominates in the 
!             absence of viscous forces.
!      -Note: If the wind stress <<< 1 over large part of the domain
!             w << 0.1, may be necessary. This only happens in idealized
!             run with wind stress specified by the user. 
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!     Solve for the ice velocity field, usind SOR, 
!     .       [A] {u_i} = R'
!     .       [L + D + U] {u_i} = R',   where L,U = lower/upper diagonal
!     .       D ui = R' - (L + U) u_i
!     .       ui = ui_old + w * (R' - (L+U) u_i) / D
!------------------------------------------------------------------------

      
 5    continue

      error   = 0d0   !   \  
      ierroru = 0     !   |   Initialization of variables:
      jerroru = 0     !    >   error between successive iterations &  
      ierrorv = 0     !   |   location (i,j) where the max error occurs
      jerrorv = 0     !   /

      residual = 0d0

      do j = 1, ny+1
         do i = 1, nx+1


            if ( maskB(i,j) + maskB(i,j+1) .eq. 0 ) goto 100

            D = 0d0

            B1   = bu(i,j)

!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------  

            vavg = ( vice(i,j)   + vice(i,j+1)                           &
                   + vice(i-1,j) + vice(i-1,j+1) ) / 4d0

            hvert = ( h(i,j) + h(i-1,j) ) / 2d0

            B1 = B1 + rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*du/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            D  = D + ( rhoice * hvert ) / Deltat

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (u - uw) cos(theta_w) - (v - vw) sin(theta_w) ] }
!------------------------------------------------------------------------

            D  = D + CdwC1(i,j) * costheta_w
 
            B1 = B1 + CdwC1(i,j) * sintheta_w * vavg 

!------------------------------------------------------------------------
!     -d ( eta dv/dy ) / dx    B1_1  p.899... 
!------------------------------------------------------------------------

            B1  = B1 - ( etaC(i,j)   * ( vice(i,j+1)   - vice(i,j)   )   &
                       - etaC(i-1,j) * ( vice(i-1,j+1) - vice(i-1,j) )   &
                       ) /  Deltax2

!------------------------------------------------------------------------
!     d ( zeta dv/dy ) / dx    B1_2
!------------------------------------------------------------------------


            B1  = B1 + ( zetaC(i,j)   * ( vice(i,j+1)   - vice(i,j)   )  &
                       - zetaC(i-1,j) * ( vice(i-1,j+1) - vice(i-1,j) )  &
                       ) /  Deltax2


!------------------------------------------------------------------------
!     d ( eta dv/dx ) / dy    B1_3 see p.1219-1226
!
!     There could be a problem at i=2 and i=nx when there is an open bc. 
!     But the mask is built so that it does not happen and only the last 
!     case (normal) can happen.
!------------------------------------------------------------------------


            if ( maskC(i-1,j+1) .eq. 0 .and. maskC(i,j+1) .eq. 1 ) then

!     x o
!     o o  
!     o o

               B1 = B1 + ( etaB(i,j+1) *                                 &
                         ( 3d0 * vice(i,j+1) - vice(i+1,j+1)/ 3d0 )      &
                         - etaB(i,j) * ( vice(i,j) - vice(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j+1) .eq. 1 .and. maskC(i,j+1) .eq. 0) then

!     o x
!     o o  
!     o o

               B1 = B1 + ( etaB(i,j+1) *                                 &
                         ( vice(i-2,j+1)/ 3d0 - 3d0 * vice(i-1,j+1))     &
                         - etaB(i,j) * ( vice(i,j) - vice(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j-1) .eq. 0 .and. maskC(i,j-1) .eq. 1) then

!     o o
!     o o  
!     x o

               B1 = B1 + ( etaB(i,j+1) *(vice(i,j+1) - vice(i-1,j+1))    &
                         - etaB(i,j) *                                   &
                         ( 3d0 * vice(i,j) - vice(i+1,j) / 3d0 )         &
                         ) / Deltax2

         
            elseif(maskC(i-1,j-1) .eq. 1 .and. maskC(i,j-1) .eq. 0) then
                  
!     o o
!     o o  
!     o x

               B1 = B1 + ( etaB(i,j+1) *(vice(i,j+1) - vice(i-1,j+1))    &
                         - etaB(i,j) *                                   &
                         ( vice(i-2,j) / 3d0 - 3d0 * vice(i-1,j) )       &
                         ) / Deltax2

            else 

!     o o                   x x    o o
!     o o  --normal case or o o or o o
!     o o                   o o    x x
               
               B1 = B1 + ( etaB(i,j+1) * ( vice(i,j+1) - vice(i-1,j+1))  &
                         + etaB(i,j)   * ( vice(i-1,j) - vice(i,j)     ) &
                         ) / Deltax2

            endif

!------------------------------------------------------------------------
!     d [ (eta + zeta ) du/dx ] / dx      B1_4, D1_4 
!------------------------------------------------------------------------

              D  = D  + ( etaC(i,j)  + etaC(i-1,j)                       &
                        + zetaC(i,j) + zetaC(i-1,j)                      &
                        ) / Deltax2

              B1 = B1 + ( etaC(i,j)    * uice(i+1,j)                     &
                        + etaC(i-1,j)  * uice(i-1,j)                     &
                        + zetaC(i,j)   * uice(i+1,j)                     &
                        + zetaC(i-1,j) * uice(i-1,j)                     &
                        ) / Deltax2
               
!------------------------------------------------------------------------
!     d ( eta du/dy ) / dy   B1_5, D1_5
!------------------------------------------------------------------------

!     o o
!     o o   -- land just below
!     x x

           if ( maskB(i,j) .eq. 0 ) then

              D  = D  + ( etaB(i,j+1) + 3d0 * etaB(i,j)                  &
                        ) / Deltax2

              B1 = B1 + ( etaB(i,j+1) * uice(i,j+1)                      &
                        + etaB(i,j)   * uice(i,j+1) / 3d0                &
                        ) / Deltax2


!     x x
!     o o   -- land just above
!     o o

           elseif ( maskB(i,j+1) .eq. 0 ) then

              D  = D  + ( 3d0 * etaB(i,j+1) + etaB(i,j)                  &
                          ) / Deltax2

              B1 = B1 + ( etaB(i,j+1) * uice(i,j-1) / 3d0                &
                        + etaB(i,j)   * uice(i,j-1)                      &
                        ) / Deltax2


!     o o
!     o o   -- open boundary just below
!     . .

           elseif ( maskB(i,j) .eq. 1 .and. j .eq. 1) then

              D  = D  + etaB(i,j+1) / Deltax2

              B1 = B1 + ( etaB(i,j+1) * uice(i,j+1) ) / Deltax2

!     . .
!     o o   -- open boundary just above
!     o o

           elseif ( maskB(i,j+1) .eq. 1 .and. j .eq. ny) then

              D  = D  + etaB(i,j) / Deltax2

              B1 = B1 + ( etaB(i,j) * uice(i,j-1) ) / Deltax2

!     o o
!     o o   -- normal situation
!     o o

           else

              D  = D  + ( etaB(i,j+1) + etaB(i,j) ) / Deltax2

              B1 = B1 + ( etaB(i,j+1) * uice(i,j+1)                      &
                        + etaB(i,j)   * uice(i,j-1)                      &
                        ) / Deltax2

           endif

!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------


           if ( i .eq. 1 ) then
              if (maskC(3,j) .eq. 1) then ! oooo
              uice(1,j)    = ( 4d0 * uice(2,j)  - uice(3,j)    )         &
                             / 3d0
              else
              uice(1,j) = uice(2,j)
              endif

           elseif ( i .eq. nx+1 ) then
              if (maskC(nx-2,j) .eq. 1) then ! oooo                      
              uice(nx+1,j) = ( 4d0 * uice(nx,j) - uice(nx-1,j) )         &
                             / 3d0
              else
              uice(nx+1,j) = uice(nx,j)
              endif
              
           elseif ( j .eq. ny+1 ) then

              uice(i,j) = 0.0d0
              
           else 
                 
              residual = B1/D - uice(i,j)
              uice(i,j)    = uice(i,j) + wsor * residual

           endif

!------------------------------------------------------------------------
!     Maximum error between succesive iterations
!------------------------------------------------------------------------

            if ( abs( residual ) .gt. error ) then
               error = abs( residual )
               ierroru = i
               jerroru = j
               uflag   = 1
               vflag   = 0
            endif

 100        continue



!------------------------------------------------------------------------
!     Calculate the v-component
!------------------------------------------------------------------------


            if ( maskB(i,j) + maskB(i+1,j) .eq. 0 ) goto 200


            D  = 0d0

            B2 = bv(i,j)

!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------

            uavg  = ( uice(i,j)   + uice(i+1,j)                          &
                 + uice(i,j-1) + uice(i+1,j-1) ) / 4d0

            hvert = ( h(i,j) + h(i,j-1) ) / 2d0

            B2 = B2 - rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            D  = D + ( rhoice * hvert ) / Deltat

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (v - vw) * cos(theta_w) + (u - uw) * sin(theta_w) ] }
!------------------------------------------------------------------------

           D  = D  + CdwC2(i,j) * costheta_w 

           B2 = B2 - CdwC2(i,j) * sintheta_w * uavg

!------------------------------------------------------------------------
!     d ( eta du/dy ) / dx     B2_3 see p.1219-1226
!
!     There could be a problem at j=2 and j=ny when there is an open bc. 
!     But the mask is built so that it does not happen and only the last 
!     case (normal) can happen.
!------------------------------------------------------------------------

           if ( maskC(i-1,j-1) .eq. 0 .and. maskC(i-1,j) .eq. 1 ) then 

!     o o o 
!     x o o               

              B2  = B2  + ( etaB(i+1,j) * ( uice(i+1,j) - uice(i+1,j-1)) &
                          - etaB(i,j)   *                                &
                          ( 3d0 * uice(i,j) - uice(i,j+1) / 3d0 )        &
                          ) / Deltax2

           elseif (maskC(i-1,j-1) .eq. 1 .and. maskC(i-1,j) .eq. 0) then 

!     x o o 
!     o o o 
                 
              B2  = B2  + ( etaB(i+1,j) * ( uice(i+1,j) - uice(i+1,j-1)) &
                          - etaB(i,j)   *                                &
                          ( uice(i,j-2) / 3d0 - 3d0 * uice(i,j-1) )      &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 0 .and. maskC(i+1,j) .eq. 1) then 

!     o o o 
!     o o x  
              
              B2  = B2  + ( etaB(i+1,j) *                                &
                          ( 3d0 * uice(i+1,j) - uice(i+1,j+1) / 3d0 )    &
                          - etaB(i,j) * ( uice(i,j) - uice(i,j-1) )      &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 1 .and. maskC(i+1,j) .eq. 0) then 

!     o o x 
!     o o o
              
              B2  = B2  + ( etaB(i+1,j) *                                &
                          ( uice(i+1,j-2) / 3d0 - 3d0 * uice(i+1,j-1) )  &
                          - etaB(i,j) * ( uice(i,j) - uice(i,j-1) )      &
                          ) / Deltax2

           else

!     o o o -- normal case or x o o or o o x
!     o o o                   x o o    o o x

              B2  = B2  + ( etaB(i+1,j) * ( uice(i+1,j) - uice(i+1,j-1)) &
                          - etaB(i,j)   * ( uice(i,j)   - uice(i,j-1) ) &
                          ) / Deltax2

           endif
           
!------------------------------------------------------------------------
!     -d ( eta du/dx ) / dy    B2_1
!------------------------------------------------------------------------
           
           B2  = B2  - ( etaC(i,j)   * ( uice(i+1,j)   - uice(i,j)   )  &
                       - etaC(i,j-1) * ( uice(i+1,j-1) - uice(i,j-1) )  &
                       ) / Deltax2

!------------------------------------------------------------------------
!     d (zeta du/dx) / dy      B2_2
!------------------------------------------------------------------------

           B2  = B2  + ( zetaC(i,j)   * ( uice(i+1,j)   - uice(i,j)   ) &
                       - zetaC(i,j-1) * ( uice(i+1,j-1) - uice(i,j-1) ) &
                       ) / Deltax2

!------------------------------------------------------------------------
!     d [ (eta + zeta) dv/dy ] / dy   D2_4, B2_4
!------------------------------------------------------------------------

               D  = D  + ( etaC(i,j)  + etaC(i,j-1)              &
                         + zetaC(i,j) + zetaC(i,j-1)             &
                         ) / Deltax2

               B2 = B2 + ( etaC(i,j)    * vice(i,j+1)            &
                         + etaC(i,j-1)  * vice(i,j-1)            &
                         + zetaC(i,j)   * vice(i,j+1)            &
                         + zetaC(i,j-1) * vice(i,j-1)            &
                         ) / Deltax2
   
!------------------------------------------------------------------------
!     d ( eta dv/dx ) / dx   D2_5, B2_5
!------------------------------------------------------------------------

!     x o o -- land to the left

            if ( maskB(i,j) .eq. 0 ) then

               D  = D  + ( etaB(i+1,j) + 3d0 * etaB(i,j) )       &
                         / Deltax2

               B2 = B2 + ( etaB(i+1,j) * vice(i+1,j)             &
                         + etaB(i,j)   * vice(i+1,j) / 3d0       &
                         ) / Deltax2

!     o o x -- land to the right

           elseif ( maskB(i+1,j) .eq. 0 ) then

               D  = D  + ( 3d0 * etaB(i+1,j) + etaB(i,j) )       &
                         / Deltax2

               B2 = B2 + ( etaB(i+1,j) * vice(i-1,j) / 3d0       &
                         + etaB(i,j)   * vice(i-1,j)             &
                         ) / Deltax2

!     . o o -- open boundary to the left

           elseif ( maskB(i,j) .eq. 1 .and. i .eq. 1) then

               D  = D  + etaB(i+1,j) / Deltax2

               B2 = B2 + ( etaB(i+1,j) * vice(i+1,j) ) / Deltax2

!     o o . -- open boundary to the right

           elseif ( maskB(i+1,j) .eq. 1 .and. i .eq. nx) then

               D  = D  + etaB(i,j) / Deltax2

               B2 = B2 + ( etaB(i,j) * vice(i-1,j) ) / Deltax2

!     o o o -- normal situation

           else

               D  = D  + ( etaB(i+1,j) + etaB(i,j) ) / Deltax2

               B2 = B2 + ( etaB(i+1,j) * vice(i+1,j)             &
                         + etaB(i,j)   * vice(i-1,j)             &
                         ) / Deltax2
           endif


!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------


           if ( j .eq. 1 ) then
              if (maskC(i,3) .eq. 1) then
              vice(i,1)    = ( 4d0 * vice(i,2)  - vice(i,3)    ) &
                             / 3d0
              else
              vice(i,1) = vice(i,2)
              endif

           elseif ( j .eq. ny+1 ) then
              if (maskC(i,ny-2) .eq. 1) then
              vice(i,ny+1) = ( 4d0 * vice(i,ny) - vice(i,ny-1) ) &
                             / 3d0
              else
              vice(i,ny+1) = vice(i,ny)
              endif

           elseif ( i .eq. nx+1 ) then

              vice(i,j) = 0.0d0
              
           else 

              residual  = B2/D - vice(i,j)
              vice(i,j) = vice(i,j) + wsor * residual

           endif

!------------------------------------------------------------------------
!     Maximum error between succesive iterations
!------------------------------------------------------------------------

            if ( abs( residual ) .gt. error ) then
               error = abs( residual )
               ierrorv = i
               jerrorv = j
               uflag   = 0
               vflag   = 1
            endif


 200        continue

         enddo
      enddo

      kit = kit + 1

      if (kit .eq. nmax) then
         print *, 'WARNING: max # of iteration exceeded in UVsolveSOR'

         if ( uflag .eq. 1)                      &
            print *, 'max error on u = ', error, &
                   'at (i,j) =', ierroru, jerroru
         if ( vflag .eq. 1)                      &
            print *, 'max error on v = ', error, &
                   'at (i,j) =', ierrorv, jerrorv

         if ( kit .eq. nmax ) goto 10 
      endif 

!      print *, 'error', error

      if ( error .gt. tol ) goto 5  


 10   continue
      
      print *, '# of iteration in UVsolveSOR= ', kit

      return
      end







