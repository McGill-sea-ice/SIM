
      subroutine precondJacobi (rhs, wk2 )
        use numerical_VP
!--------------------------------------------------------------------------
! Preconditioner: calculates M x = rhs where rhs is wk1. x is the initial
! guess. x is set to 0 here. The solution x is then put in wk2.
! GMRES solves the mom equation fully implicitly. However, to increase the
! diag dominance, f=0 and theta_w = 0 in the precondioner.
!--------------------------------------------------------------------------


      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j, k
      
      double precision rhs(nvar), wk2(nvar)
      double precision residual, D, B1, B2
      double precision hvert    !, uavg, vavg
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision uold(0:nx+2,0:ny+2), vold(0:nx+2,0:ny+2)
      double precision rhsu(0:nx+2,0:ny+2), rhsv(0:nx+2,0:ny+2)

      call transformer(rhsu,rhsv,rhs,0)

!------------------------------------------------------------------------
!  Set initial guess for preconditioner
!------------------------------------------------------------------------
      
            utp = 0.0d0
            vtp = 0.0d0

!------------------------------------------------------------------------
!  apply Jacobi iterations
!------------------------------------------------------------------------

      do k = 1, kjac

         uold = utp
         vold = vtp
         residual = 0.0d0

      do j = 1, ny+1
         do i = 1, nx+1


            if ( maskB(i,j) + maskB(i,j+1) .eq. 0 ) goto 100

            B1   = rhsu(i,j)

!------------------------------------------------------------------------
!     Coriolis term (f = 0 in the precond)
!------------------------------------------------------------------------  

!            vavg = ( vtp(i,j)   + vtp(i,j+1)                           &
!                   + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

            hvert = ( h(i,j) + h(i-1,j) ) / 2d0

!            B1 = B1 + rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*du/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            if ( BDF .eq. 0 ) then
               D  = ( rhoice * hvert ) / Deltat
            elseif ( BDF .eq. 1 ) then
               D = ( 3d0 * rhoice * hvert ) / ( 2d0 * Deltat )
            endif

!------------------------------------------------------------------------
!     water drag  (theta_w is considered to be zero in the precond)
!------------------------------------------------------------------------

            D  = D + CdwC1f(i,j) !* costheta_w 

!            B1 = B1 + CdwC1f(i,j) * sintheta_w * vavg 
 
!------------------------------------------------------------------------
!     -d ( eta dv/dy ) / dx    B1_1  p.899... 
!------------------------------------------------------------------------

            B1  = B1 - ( etaCf(i,j)   * ( vold(i,j+1)   - vold(i,j)   )   &
                       - etaCf(i-1,j) * ( vold(i-1,j+1) - vold(i-1,j) )   &
                       ) /  Deltax2

!------------------------------------------------------------------------
!     d ( zeta dv/dy ) / dx    B1_2
!------------------------------------------------------------------------


            B1  = B1 + ( zetaCf(i,j)   * ( vold(i,j+1)   - vold(i,j)   )  &
                       - zetaCf(i-1,j) * ( vold(i-1,j+1) - vold(i-1,j) )  &
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

               B1 = B1 + ( etaBf(i,j+1) *                                 &
                         ( 3d0 * vold(i,j+1) - vold(i+1,j+1)/ 3d0 )      &
                         - etaBf(i,j) * ( vold(i,j) - vold(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j+1) .eq. 1 .and. maskC(i,j+1) .eq. 0) then

!     o x
!     o o  
!     o o

               B1 = B1 + ( etaBf(i,j+1) *                                 &
                         ( vold(i-2,j+1)/ 3d0 - 3d0 * vold(i-1,j+1))     &
                         - etaBf(i,j) * ( vold(i,j) - vold(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j-1) .eq. 0 .and. maskC(i,j-1) .eq. 1) then

!     o o
!     o o  
!     x o

               B1 = B1 + ( etaBf(i,j+1) *(vold(i,j+1) - vold(i-1,j+1))    &
                         - etaBf(i,j) *                                   &
                         ( 3d0 * vold(i,j) - vold(i+1,j) / 3d0 )         &
                         ) / Deltax2

         
            elseif(maskC(i-1,j-1) .eq. 1 .and. maskC(i,j-1) .eq. 0) then
                  
!     o o
!     o o  
!     o x

               B1 = B1 + ( etaBf(i,j+1) *(vold(i,j+1) - vold(i-1,j+1))    &
                         - etaBf(i,j) *                                   &
                         ( vold(i-2,j) / 3d0 - 3d0 * vold(i-1,j) )       &
                         ) / Deltax2

            else 

!     o o                   x x    o o
!     o o  --normal case or o o or o o
!     o o                   o o    x x
               
               B1 = B1 + ( etaBf(i,j+1) * ( vold(i,j+1) - vold(i-1,j+1))  &
                         + etaBf(i,j)   * ( vold(i-1,j) - vold(i,j)     ) &
                         ) / Deltax2

            endif

!------------------------------------------------------------------------
!     d [ (eta + zeta ) du/dx ] / dx      B1_4, D1_4 
!------------------------------------------------------------------------

              D  = D  + ( etaCf(i,j)  + etaCf(i-1,j)                       &
                        + zetaCf(i,j) + zetaCf(i-1,j)                      &
                        ) / Deltax2

              B1 = B1 + ( etaCf(i,j)    * uold(i+1,j)                     &
                        + etaCf(i-1,j)  * uold(i-1,j)                     &
                        + zetaCf(i,j)   * uold(i+1,j)                     &
                        + zetaCf(i-1,j) * uold(i-1,j)                     &
                        ) / Deltax2
               
!------------------------------------------------------------------------
!     d ( eta du/dy ) / dy   B1_5, D1_5
!------------------------------------------------------------------------

!     o o
!     o o   -- land just below
!     x x

           if ( maskB(i,j) .eq. 0 ) then

              D  = D  + ( etaBf(i,j+1) + 3d0 * etaBf(i,j)                  &
                        ) / Deltax2

              B1 = B1 + ( etaBf(i,j+1) * uold(i,j+1)                      &
                        + etaBf(i,j)   * uold(i,j+1) / 3d0                &
                        ) / Deltax2


!     x x
!     o o   -- land just above
!     o o

           elseif ( maskB(i,j+1) .eq. 0 ) then

              D  = D  + ( 3d0 * etaBf(i,j+1) + etaBf(i,j)                  &
                          ) / Deltax2

              B1 = B1 + ( etaBf(i,j+1) * uold(i,j-1) / 3d0                &
                        + etaBf(i,j)   * uold(i,j-1)                      &
                        ) / Deltax2


!     o o
!     o o   -- open boundary just below
!     . .

           elseif ( maskB(i,j) .eq. 1 .and. j .eq. 1) then

              D  = D  + etaBf(i,j+1) / Deltax2

              B1 = B1 + ( etaBf(i,j+1) * uold(i,j+1) ) / Deltax2

!     . .
!     o o   -- open boundary just above
!     o o

           elseif ( maskB(i,j+1) .eq. 1 .and. j .eq. ny) then

              D  = D  + etaBf(i,j) / Deltax2

              B1 = B1 + ( etaBf(i,j) * uold(i,j-1) ) / Deltax2

!     o o
!     o o   -- normal situation
!     o o

           else

              D  = D  + ( etaBf(i,j+1) + etaBf(i,j) ) / Deltax2

              B1 = B1 + ( etaBf(i,j+1) * uold(i,j+1)                      &
                        + etaBf(i,j)   * uold(i,j-1)                      &
                        ) / Deltax2

           endif

!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------


           if ( i .eq. 1 ) then

              D = 1.0d0
              if (maskC(3,j) .eq. 1) then ! oooo
              B1 = rhsu(i,j) + ( 4d0 * uold(2,j)  - uold(3,j)    )       &
                               / 3d0
              else ! ooox                                              
              B1 = rhsu(i,j) + uold(2,j)
              endif


              residual = B1/D - utp(i,j)
              utp(i,j)    = utp(i,j) + wjac*residual

           elseif ( i .eq. nx+1 ) then

              D = 1.0d0
              if (maskC(nx-2,j) .eq. 1) then ! oooo   
              B1 = rhsu(i,j) + ( 4d0 * uold(nx,j) - uold(nx-1,j) )       &
                               / 3d0
              else ! xooo                                      
              B1 = rhsu(i,j) + uold(nx,j)
              endif

              residual = B1/D - utp(i,j)
              utp(i,j)    = utp(i,j) + wjac*residual

           elseif ( j .eq. ny+1 ) then

              utp(i,j) = 0.0d0

           else 
                 
              residual = B1/D - utp(i,j)
              utp(i,j)    = utp(i,j) + wjac*residual

           endif


 100       continue


!------------------------------------------------------------------------
!     Calculate the v-component
!------------------------------------------------------------------------


            if ( maskB(i,j) + maskB(i+1,j) .eq. 0 ) goto 200

            B2 = rhsv(i,j)

!------------------------------------------------------------------------
!     Coriolis term (f = 0 in the precond)
!------------------------------------------------------------------------

!            uavg  = ( utp(i,j)   + utp(i+1,j)                          &
!                    + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

            hvert = ( h(i,j) + h(i,j-1) ) / 2d0

!            B2 = B2 - rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            if ( BDF .eq. 0 ) then
               D  = ( rhoice * hvert ) / Deltat
            elseif ( BDF .eq. 1 ) then
               D = ( 3d0 * rhoice * hvert ) / ( 2d0 * Deltat )
            endif

!------------------------------------------------------------------------
!     water drag (theta_w is considered to be zero in the precond)
!------------------------------------------------------------------------

            D  = D  + CdwC2f(i,j) !* costheta_w 
           
!            B2 = B2 - CdwC2f(i,j) * sintheta_w * uavg

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

              B2  = B2  + ( etaBf(i+1,j) * ( uold(i+1,j) - uold(i+1,j-1)) &
                          - etaBf(i,j)   *                                &
                          ( 3d0 * uold(i,j) - uold(i,j+1) / 3d0 )        &
                          ) / Deltax2

           elseif (maskC(i-1,j-1) .eq. 1 .and. maskC(i-1,j) .eq. 0) then 

!     x o o 
!     o o o 
                 
              B2  = B2  + ( etaBf(i+1,j) * ( uold(i+1,j) - uold(i+1,j-1)) &
                          - etaBf(i,j)   *                                &
                          ( uold(i,j-2) / 3d0 - 3d0 * uold(i,j-1) )      &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 0 .and. maskC(i+1,j) .eq. 1) then 

!     o o o 
!     o o x  
              
              B2  = B2  + ( etaBf(i+1,j) *                                &
                          ( 3d0 * uold(i+1,j) - uold(i+1,j+1) / 3d0 )    &
                          - etaBf(i,j) * ( uold(i,j) - uold(i,j-1) )      &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 1 .and. maskC(i+1,j) .eq. 0) then 

!     o o x 
!     o o o
              
              B2  = B2  + ( etaBf(i+1,j) *                                &
                          ( uold(i+1,j-2) / 3d0 - 3d0 * uold(i+1,j-1) )  &
                          - etaBf(i,j) * ( uold(i,j) - uold(i,j-1) )      &
                          ) / Deltax2

           else

!     o o o -- normal case or x o o or o o x
!     o o o                   x o o    o o x

              B2  = B2  + ( etaBf(i+1,j) * ( uold(i+1,j) - uold(i+1,j-1)) &
                          - etaBf(i,j)   * ( uold(i,j)   - uold(i,j-1) ) &
                          ) / Deltax2

           endif
           
!------------------------------------------------------------------------
!     -d ( eta du/dx ) / dy    B2_1
!------------------------------------------------------------------------
           
           B2  = B2  - ( etaCf(i,j)   * ( uold(i+1,j)   - uold(i,j)   )  &
                       - etaCf(i,j-1) * ( uold(i+1,j-1) - uold(i,j-1) )  &
                       ) / Deltax2

!------------------------------------------------------------------------
!     d (zeta du/dx) / dy      B2_2
!------------------------------------------------------------------------

           B2  = B2  + ( zetaCf(i,j)   * ( uold(i+1,j)   - uold(i,j)   ) &
                       - zetaCf(i,j-1) * ( uold(i+1,j-1) - uold(i,j-1) ) &
                       ) / Deltax2

!------------------------------------------------------------------------
!     d [ (eta + zeta) dv/dy ] / dy   D2_4, B2_4
!------------------------------------------------------------------------

               D  = D  + ( etaCf(i,j)  + etaCf(i,j-1)                 &
                         + zetaCf(i,j) + zetaCf(i,j-1)                &
                         ) / Deltax2

               B2 = B2 + ( etaCf(i,j)    * vold(i,j+1)               &
                         + etaCf(i,j-1)  * vold(i,j-1)               &
                         + zetaCf(i,j)   * vold(i,j+1)               &
                         + zetaCf(i,j-1) * vold(i,j-1)               &
                         ) / Deltax2
   
!------------------------------------------------------------------------
!     d ( eta dv/dx ) / dx   D2_5, B2_5
!------------------------------------------------------------------------

!     x o o -- land to the left

            if ( maskB(i,j) .eq. 0 ) then

               D  = D  + ( etaBf(i+1,j) + 3d0 * etaBf(i,j) )          &
                         / Deltax2

               B2 = B2 + ( etaBf(i+1,j) * vold(i+1,j)                &
                         + etaBf(i,j)   * vold(i+1,j) / 3d0          &
                         ) / Deltax2

!     o o x -- land to the right

           elseif ( maskB(i+1,j) .eq. 0 ) then

               D  = D  + ( 3d0 * etaBf(i+1,j) + etaBf(i,j) )          &
                         / Deltax2

               B2 = B2 + ( etaBf(i+1,j) * vold(i-1,j) / 3d0          &
                         + etaBf(i,j)   * vold(i-1,j)                &
                         ) / Deltax2

!     . o o -- open boundary to the left

           elseif ( maskB(i,j) .eq. 1 .and. i .eq. 1) then

               D  = D  + etaBf(i+1,j) / Deltax2

               B2 = B2 + ( etaBf(i+1,j) * vold(i+1,j) ) / Deltax2

!     o o . -- open boundary to the right

           elseif ( maskB(i+1,j) .eq. 1 .and. i .eq. nx) then

               D  = D  + etaBf(i,j) / Deltax2

               B2 = B2 + ( etaBf(i,j) * vold(i-1,j) ) / Deltax2

!     o o o -- normal situation

           else

               D  = D  + ( etaBf(i+1,j) + etaBf(i,j) ) / Deltax2

               B2 = B2 + ( etaBf(i+1,j) * vold(i+1,j)                &
                         + etaBf(i,j)   * vold(i-1,j)                &
                         ) / Deltax2
           endif


!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------


           if ( j .eq. 1 ) then

              D = 1.0d0
              if (maskC(i,3) .eq. 1) then
              B2 =  rhsv(i,j) + ( 4d0 * vold(i,2)  - vold(i,3)    ) &
                                / 3d0
              else
              B2 = rhsv(i,j) + vold(i,2)
              endif

              residual  = B2/D - vtp(i,j)
              vtp(i,j) = vtp(i,j) + wjac*residual

           elseif ( j .eq. ny+1 ) then

              D = 1.0d0
              if (maskC(i,ny-2) .eq. 1) then
              B2 =  rhsv(i,j) + ( 4d0 * vold(i,ny) - vold(i,ny-1) ) &
                                / 3d0
              else
              B2 = rhsv(i,j) + vold(i,ny)
              endif

              residual  = B2/D - vtp(i,j)
              vtp(i,j) = vtp(i,j) + wjac*residual

           elseif ( i .eq. nx+1 ) then

              vtp(i,j) = 0.0d0

           else 

              residual  = B2/D - vtp(i,j)
              vtp(i,j) = vtp(i,j) + wjac*residual

           endif


 200        continue

         enddo
      enddo


      enddo

      call transformer (utp,vtp,wk2,1)

      return
    end subroutine PRECONDJACOBI







