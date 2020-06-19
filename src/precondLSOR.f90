
      subroutine precondLSOR ( rhs, wk2 )
        use numerical_VP
!--------------------------------------------------------------------------
! Preconditioner: calculates M x = rhs where rhs is wk1. x is the initial
! guess. x is set to 0 here. The solution x is then put in wk2.
! M is almost the matrix A of the system Ax=b but with f=0 and thetaw=0.
! GMRES solves the syetm fully implicitly however because these terms are 
! implicit in Jacfreevec.f90
!--------------------------------------------------------------------------

      implicit none
      
      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_options.h'

      double precision rhs(nvar), wk2(nvar)
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision rhsu(0:nx+2,0:ny+2), rhsv(0:nx+2,0:ny+2)

      integer i, j, k
      
      call transformer(rhsu,rhsv,rhs,0)

!------------------------------------------------------------------------
!  Set initial guess for preconditioner
!------------------------------------------------------------------------

      do j = 1, ny+1
         do i = 1, nx+1

            utp(i,j) = 0.0d0
            vtp(i,j) = 0.0d0

         enddo
      enddo

      do k = 1, klsor               ! nb of times the grid is swept
      
         do j = 1, ny

            call vectu_3diag (utp, vtp, rhsu, j) 

         enddo

         do i = 1, nx

            call vectv_3diag (utp, vtp, rhsv, i) 

         enddo
         
      enddo


      call transformer (utp,vtp,wk2,1)

      return
    end subroutine precondLSOR


!************************************************************************
!     creates the 3 column vectors (u component) for the LSOR method 
!     solved with the routine tridag.f
!
!************************************************************************

      subroutine vectu_3diag (utp, vtp, rhsu, j)
        use numerical_VP
      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      integer i, j

      double precision aVEC(nx+1),bVEC(nx+1),cVEC(nx+1)
      double precision rhs(nx+1), X(nx+1), Xold(nx+1)
      double precision hvert

      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision rhsu(0:nx+2,0:ny+2)

      do i = 1, nx+1

         aVEC(i) = 0.0d0
         bVEC(i) = 0.0d0
         cVEC(i) = 0.0d0
!         rhs(k)  = 0.0d0
         Xold(i) = utp(i,j)

      enddo

      do i = 1, nx+1

         if ( maskB(i,j) + maskB(i,j+1) .eq. 0 ) goto 100
            
         rhs(i) = rhsu(i,j)

!------------------------------------------------------------------------
!     Coriolis term (f is considered to be 0 in this precond...see p.1354)
!------------------------------------------------------------------------  

!         vavg = ( vtp(i,j)   + vtp(i,j+1)                              &
!                + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

         hvert = ( h(i,j) + h(i-1,j) ) / 2d0

!         rhs(i) = rhs(i) + rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*du/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

         if ( BDF .eq. 0 ) then
            bVEC(i) = ( rhoice * hvert ) / Deltat
         elseif ( BDF .eq. 1 ) then
            bVEC(i) = ( 3d0 * rhoice * hvert ) / ( 2d0 * Deltat )
         endif

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (u - uw) cos(theta_w) - (v - vw) sin(theta_w) ] }
!     - thetaw is considered to be 0 in this precond...see p.1354
!------------------------------------------------------------------------

         bVEC(i) = bVEC(i) + CdwC1f(i,j) !* costheta_w
         bVEC(i) = bVEC(i) + Cbasal1(i,j)
!         rhs(i) = rhs(i) + CdwC1f(i,j) * sintheta_w * vavg 

!------------------------------------------------------------------------
!     -d ( eta dv/dy ) / dx    B1_1  p.899... 
!------------------------------------------------------------------------

         if (stressBC .and. i .eq. 1) then ! W

            rhs(i) = rhs(i) - ( etaCf(i,j) * ( vtp(i,j+1) - vtp(i,j) )  &
                            ) /  DxhDx

         elseif (stressBC .and. i .eq. nx+1) then ! E                                         
            rhs(i) = rhs(i) + ( etaCf(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j)) &
                                          ) /  DxhDx

         else

            rhs(i) = rhs(i) - ( etaCf(i,j) * ( vtp(i,j+1) - vtp(i,j) )  &
                         - etaCf(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j) ) &
                            ) /  Deltax2

         endif

!------------------------------------------------------------------------
!     d ( zeta dv/dy ) / dx    B1_2
!------------------------------------------------------------------------

         if (stressBC .and. i .eq. 1) then ! W

            rhs(i) = rhs(i) + ( zetaCf(i,j) * ( vtp(i,j+1) - vtp(i,j) ) &
                                          ) /  DxhDx

         elseif (stressBC .and. i .eq. nx+1) then ! E 

            rhs(i) = rhs(i) - ( zetaCf(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j) )&
                                          ) /  DxhDx

         else

            rhs(i) = rhs(i) + ( zetaCf(i,j) * ( vtp(i,j+1) - vtp(i,j) ) &
                         - zetaCf(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j)) &
                                          ) /  Deltax2

         endif

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

               rhs(i) = rhs(i) + ( etaBf(i,j+1) *                         &
                         ( 3d0 * vtp(i,j+1) - vtp(i+1,j+1)/ 3d0 )      &
                         - etaBf(i,j) * ( vtp(i,j) - vtp(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j+1) .eq. 1 .and. maskC(i,j+1) .eq. 0) then

!     o x
!     o o  
!     o o

               rhs(i) = rhs(i) + ( etaBf(i,j+1) *                         &
                         ( vtp(i-2,j+1)/ 3d0 - 3d0 * vtp(i-1,j+1))     &
                         - etaBf(i,j) * ( vtp(i,j) - vtp(i-1,j))        &
                         ) / Deltax2

            elseif(maskC(i-1,j-1) .eq. 0 .and. maskC(i,j-1) .eq. 1) then

!     o o
!     o o  
!     x o

               rhs(i) = rhs(i) + ( etaBf(i,j+1) *                         &
                         (vtp(i,j+1) - vtp(i-1,j+1))                   &
                         - etaBf(i,j) *                                   &
                         ( 3d0 * vtp(i,j) - vtp(i+1,j) / 3d0 )         &
                         ) / Deltax2

         
            elseif(maskC(i-1,j-1) .eq. 1 .and. maskC(i,j-1) .eq. 0) then
                  
!     o o
!     o o  
!     o x

               rhs(i) = rhs(i) + ( etaBf(i,j+1) *                         &
                         (vtp(i,j+1) - vtp(i-1,j+1))                   &
                         - etaBf(i,j) *                                   &
                         ( vtp(i-2,j) / 3d0 - 3d0 * vtp(i-1,j) )       &
                         ) / Deltax2

            else 

!     o o                   x x    o o
!     o o  --normal case or o o or o o
!     o o                   o o    x x
               if (stressBC .and. i .eq. 1) then ! W

                  rhs(i) = rhs(i) + 0d0

               elseif (stressBC .and. i .eq. nx+1) then ! E

                  rhs(i) = rhs(i) + 0d0

               elseif (stressBC .and. j .eq. 1 .and. i .ge. 2 .and. i .le. nx) then ! S

                  rhs(i) = rhs(i) + ( etaBf(i,j+1) *                      &
                         ( vtp(i,j+1) - vtp(i-1,j+1) ) ) / Deltax2

               elseif (stressBC .and. j .eq. ny .and. i .ge. 2 .and. i .le. nx) then ! N

                  rhs(i) = rhs(i) + ( etaBf(i,j)   *                      &
                         ( vtp(i-1,j) - vtp(i,j)     ) ) / Deltax2

               else
                     
                  rhs(i) = rhs(i) + ( etaBf(i,j+1) *                      &
                         ( vtp(i,j+1) - vtp(i-1,j+1))                  &
                         + etaBf(i,j)   * ( vtp(i-1,j) - vtp(i,j)     ) &
                         ) / Deltax2

               endif

            endif

!------------------------------------------------------------------------
!     d [ (eta + zeta ) du/dx ] / dx      B1_4, D1_4 
!------------------------------------------------------------------------

            if (stressBC .and. i .eq. 1) then ! W    


               bVEC(i) = bVEC(i) + ( etaCf(i,j)  + zetaCf(i,j) ) / DxhDx


               aVEC(i) = aVEC(i) + 0d0


               cVEC(i) = cVEC(i) - ( etaCf(i,j) + zetaCf(i,j) ) / DxhDx

            elseif (stressBC .and. i .eq. nx+1) then ! E  

               bVEC(i) = bVEC(i) + ( etaCf(i-1,j) + zetaCf(i-1,j) ) / DxhDx


               aVEC(i) = aVEC(i) - ( etaCf(i-1,j) + zetaCf(i-1,j) ) / DxhDx


               cVEC(i) = cVEC(i) + 0d0

            else

               bVEC(i) = bVEC(i) + ( etaCf(i,j)  + etaCf(i-1,j)         &
                           +   zetaCf(i,j) + zetaCf(i-1,j)        &
                             ) / Deltax2


               aVEC(i) = aVEC(i) - ( etaCf(i-1,j) +                    &
                               zetaCf(i-1,j) ) / Deltax2


               cVEC(i) = cVEC(i) - ( etaCf(i,j) +                      &
                               zetaCf(i,j) ) / Deltax2           

            endif

!------------------------------------------------------------------------
!     d ( eta du/dy ) / dy   B1_5, D1_5
!------------------------------------------------------------------------

!     o o
!     o o   -- land just below
!     x x

         if ( maskB(i,j) .eq. 0 ) then

            bVEC(i) = bVEC(i) + ( etaBf(i,j+1) + 3d0 * etaBf(i,j) &
                                ) / Deltax2

            rhs(i) = rhs(i) + ( etaBf(i,j+1) * utp(i,j+1)       &
                            +   etaBf(i,j)   * utp(i,j+1) / 3d0 &
                              ) / Deltax2


!     x x
!     o o   -- land just above
!     o o

         elseif ( maskB(i,j+1) .eq. 0 ) then

            bVEC(i) = bVEC(i) + ( 3d0 * etaBf(i,j+1) + etaBf(i,j) &
                                ) / Deltax2

            rhs(i) = rhs(i) + ( etaBf(i,j+1) * utp(i,j-1) / 3d0 &
                            +   etaBf(i,j)   * utp(i,j-1) &
                              ) / Deltax2


!     o o
!     o o   -- open boundary just below
!     . .

         elseif ( maskB(i,j) .eq. 1 .and. j .eq. 1 .and. .not. stressBC) then

            bVEC(i) = bVEC(i) + etaBf(i,j+1) / Deltax2

            rhs(i) = rhs(i) + ( etaBf(i,j+1) * utp(i,j+1) &
                               ) / Deltax2

!     . .
!     o o   -- open boundary just above
!     o o

         elseif ( maskB(i,j+1) .eq. 1 .and. j .eq. ny .and. .not. stressBC) then

            bVEC(i) = bVEC(i) +  etaBf(i,j) / Deltax2

            rhs(i) = rhs(i) + ( etaBf(i,j) * utp(i,j-1)   &
                              ) / Deltax2

!     o o
!     o o   -- normal situation
!     o o

         else

            if (stressBC .and. i .eq. 1) then ! W 

               bVEC(i) = bVEC(i) + 0d0

            elseif (stressBC .and. i .eq. nx+1) then ! E 

               bVEC(i) = bVEC(i) + 0d0

            elseif (stressBC .and. j .eq. 1 .and. i .ge. 2 .and. i .le. nx) then ! S  

               bVEC(i) = bVEC(i) + ( etaBf(i,j+1) ) / Deltax2

               rhs(i) = rhs(i) + ( etaBf(i,j+1) * utp(i,j+1) ) / Deltax2

            elseif (stressBC .and. j .eq. ny .and. i .ge. 2 .and. i .le. nx) then ! N

               bVEC(i) = bVEC(i) + ( etaBf(i,j) ) / Deltax2

               rhs(i) = rhs(i) + ( etaBf(i,j)   * utp(i,j-1) ) / Deltax2

            else

               bVEC(i) = bVEC(i) + ( etaBf(i,j+1) + etaBf(i,j) &
                    ) / Deltax2

               rhs(i) = rhs(i) + ( etaBf(i,j+1) * utp(i,j+1) &
                    +   etaBf(i,j)   * utp(i,j-1) &
                    ) / Deltax2

            endif

         endif

!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------

         if (stressBC .and. clipping) then

            if (i .eq. 1 .and. j .eq. 1) then

               aVEC(i) = 0d0
               bVEC(i) = 1.0d0 ! clipping                                    
               cVEC(i) = 0d0
               rhs(i)  = 0d0
               
            endif

         elseif (.not. stressBC) then

            if ( i .eq. 1 ) then

               if (maskC(3,j) .eq. 1) then ! oooo                               

               aVEC(i) = 0.0d0
               bVEC(i) = 1.0d0
               cVEC(i) = -4d0/3d0

               rhs(i) = rhsu(i,j) - utp(3,j) / 3d0

            else ! ooox                                                         

               aVEC(i) = 0.0d0
               bVEC(i) = 1.0d0
               cVEC(i) = -1.0d0

               rhs(i) = rhsu(i,j)
               endif

            elseif ( i .eq. nx+1 ) then

               if (maskC(nx-2,j) .eq. 1) then ! oooo                            

               aVEC(i) = -4d0/3d0
               bVEC(i) = 1.0d0
               cVEC(i) = 0.0d0

               rhs(i) = rhsu(i,j) - utp(nx-1,j) / 3d0

               else ! xooo                                                      

               aVEC(i) = -1.0d0
               bVEC(i) = 1.0d0
               cVEC(i) = 0.0d0

               rhs(i) = rhsu(i,j)
               endif

            endif

         endif


 100     continue

      enddo

      call tridag (aVEC, bVEC, cVEC, rhs, X, nx+1)

      do i = 1, nx+1

         utp(i,j) = Xold(i) + wlsor*(X(i)-Xold(i))

      enddo

      return
    end subroutine vectu_3diag


!************************************************************************
!     creates the 3 column vectors (v component) for the LSOR method 
!     solved with the routine tridag.f
!
!************************************************************************

      subroutine vectv_3diag (utp, vtp, rhsv, i)
        use numerical_VP
      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      integer i, j
      
      double precision aVEC(ny+1),bVEC(ny+1),cVEC(ny+1)
      double precision rhs(ny+1), Y(ny+1), Yold(ny+1)
      double precision hvert

      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision rhsv(0:nx+2,0:ny+2)

      do j = 1, ny+1

         aVEC(j) = 0.0d0
         bVEC(j) = 0.0d0
         cVEC(j) = 0.0d0
!         rhs(k)  = 0.0d0
         Yold(j) = vtp(i,j)

      enddo

      do j = 1, ny+1


         if ( maskB(i,j) + maskB(i+1,j) .eq. 0 ) goto 200

         rhs(j) = rhsv(i,j)

!------------------------------------------------------------------------
!     Coriolis term (f is considered to be 0 in this precond...see p.1354)
!------------------------------------------------------------------------
         
!         uavg  = ( utp(i,j)   + utp(i+1,j)                            &
!                 + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

         hvert = ( h(i,j) + h(i,j-1) ) / 2d0
            
!         rhs(j) = rhs(j) - rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

         if ( BDF .eq. 0 ) then
            bVEC(j) = ( rhoice * hvert ) / Deltat
         elseif ( BDF .eq. 1 ) then
            bVEC(j) = ( 3d0 * rhoice * hvert ) / ( 2d0 * Deltat )
         endif

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (v - vw) * cos(theta_w) + (u - uw) * sin(theta_w) ] }
!     - thetaw is considered to be 0 in this precond...see p.1354
!------------------------------------------------------------------------

         bVEC(j) = bVEC(j) + CdwC2f(i,j) !* costheta_w 
         bVEC(j) = bVEC(j) + Cbasal2(i,j)
!         rhs(j) = rhs(j) - CdwC2f(i,j) * sintheta_w * uavg
           
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

              rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                       &
                         ( utp(i+1,j) - utp(i+1,j-1))                 &
                         - etaBf(i,j)   *                                &
                         ( 3d0 * utp(i,j) - utp(i,j+1) / 3d0 )        &
                          ) / Deltax2

           elseif (maskC(i-1,j-1) .eq. 1 .and. maskC(i-1,j) .eq. 0) then 

!     x o o 
!     o o o 
                 
              rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                       &
                         ( utp(i+1,j) - utp(i+1,j-1))                 &
                         - etaBf(i,j)   *                                &
                         ( utp(i,j-2) / 3d0 - 3d0 * utp(i,j-1) )      &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 0 .and. maskC(i+1,j) .eq. 1) then 

!     o o o 
!     o o x  
              
              rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                       &
                          ( 3d0 * utp(i+1,j) - utp(i+1,j+1) / 3d0 )   &
                          - etaBf(i,j) * ( utp(i,j) - utp(i,j-1) )     &
                          ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 1 .and. maskC(i+1,j) .eq. 0) then 

!     o o x 
!     o o o
              
              rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                       &
                          ( utp(i+1,j-2) / 3d0 - 3d0 * utp(i+1,j-1) ) &
                          - etaBf(i,j) * ( utp(i,j) - utp(i,j-1) )    &
                          ) / Deltax2

           else

!     o o o -- normal case or x o o or o o x
!     o o o                   x o o    o o x

              if (stressBC .and. j .eq. 1) then ! S 
                 
                 rhs(j) = rhs(j) + 0d0

              elseif (stressBC .and. j .eq. ny+1) then ! N

                 rhs(j) = rhs(j) + 0d0

              elseif (stressBC .and. i .eq. 1 .and. j .ge. 2 .and. j .le. ny) then ! W 

                 rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                      &
                         ( utp(i+1,j) - utp(i+1,j-1) )                     &
                         ) / Deltax2

              elseif (stressBC .and. i .eq. nx .and. j .ge. 2 .and. j .le. ny) then ! E

                 rhs(j)  = rhs(j) - ( etaBf(i,j)   *                       &
                         ( utp(i,j)   - utp(i,j-1) )                       &
                         ) / Deltax2

              else

                 rhs(j)  = rhs(j)  + ( etaBf(i+1,j) *                      &
                         ( utp(i+1,j) - utp(i+1,j-1))                &
                         - etaBf(i,j)   * ( utp(i,j)   - utp(i,j-1) ) &
                         ) / Deltax2

              endif

           endif
           
!------------------------------------------------------------------------
!     -d ( eta du/dx ) / dy    B2_1
!------------------------------------------------------------------------
           
           if (stressBC .and. j .eq. 1) then ! S 

              rhs(j) = rhs(j) - ( etaCf(i,j)   *                     &
                           ( utp(i+1,j)   - utp(i,j)   ) ) / DxhDx

           elseif(stressBC .and. j .eq. ny+1) then ! N                                          
              rhs(j) = rhs(j) + ( etaCf(i,j-1) *                     &
                           ( utp(i+1,j-1) - utp(i,j-1) )      &
                           ) / DxhDx

           else

              rhs(j) = rhs(j) - ( etaCf(i,j)   *                      &
                           ( utp(i+1,j)   - utp(i,j)   )      &
                         -   etaCf(i,j-1)   *                    &
                           ( utp(i+1,j-1) - utp(i,j-1) )      &
                           ) / Deltax2

           endif

!------------------------------------------------------------------------
!     d (zeta du/dx) / dy      B2_2
!------------------------------------------------------------------------

         if (stressBC .and. j .eq. 1) then ! S     

            rhs(j) = rhs(j) + ( zetaCf(i,j)   *                     &
                           ( utp(i+1,j)   - utp(i,j)   ) ) / DxhDx

         elseif(stressBC .and. j .eq. ny+1) then ! N

            rhs(j) = rhs(j) - ( zetaCf(i,j-1) *                     &
                           ( utp(i+1,j-1) - utp(i,j-1) )      &
                           ) / DxhDx
            
         else

            rhs(j) = rhs(j) + ( zetaCf(i,j)   *                     &
                           ( utp(i+1,j)   - utp(i,j)   )      &
                         -   zetaCf(i,j-1) *                     &
                           ( utp(i+1,j-1) - utp(i,j-1) )      &
                           ) / Deltax2

         endif

!------------------------------------------------------------------------
!     d [ (eta + zeta) dv/dy ] / dy   D2_4, B2_4
!------------------------------------------------------------------------

         if (stressBC .and. j .eq. 1) then ! S

            bVEC(j) = bVEC(j) + ( etaCf(i,j) + zetaCf(i,j) ) / DxhDx

            aVEC(j) = aVEC(j) + 0d0

            cVEC(j) = cVEC(j) - ( etaCf(i,j) + zetaCf(i,j) ) / DxhDx            

         elseif(stressBC .and. j .eq. ny+1) then ! N                

            bVEC(j) = bVEC(j) + ( etaCf(i,j-1) + zetaCf(i,j-1) ) / DxhDx

            aVEC(j) = aVEC(j) - ( etaCf(i,j-1) + zetaCf(i,j-1) ) / DxhDx

            cVEC(j) = cVEC(j) + 0d0


         else

            bVEC(j) = bVEC(j) + ( etaCf(i,j)  + etaCf(i,j-1)         &
                           +   zetaCf(i,j) + zetaCf(i,j-1)        &
                             ) / Deltax2

            aVEC(j) = aVEC(j) - ( etaCf(i,j-1) + zetaCf(i,j-1)       &
                             ) / Deltax2

            cVEC(j) = cVEC(j) - ( etaCf(i,j) + zetaCf(i,j)           &
                             ) / Deltax2

         endif

!------------------------------------------------------------------------
!     d ( eta dv/dx ) / dx   D2_5, B2_5
!------------------------------------------------------------------------

!     x o o -- land to the left

         if ( maskB(i,j) .eq. 0 ) then

            bVEC(j) = bVEC(j) + ( etaBf(i+1,j) + 3d0 * etaBf(i,j) &
                                ) / Deltax2

            rhs(j) = rhs(j) + ( etaBf(i+1,j) * vtp(i+1,j)       &
                            +   etaBf(i,j)   * vtp(i+1,j) / 3d0 &
                              ) / Deltax2

!     o o x -- land to the right

         elseif ( maskB(i+1,j) .eq. 0 ) then

            bVEC(j) = bVEC(j) + ( 3d0 * etaBf(i+1,j) + etaBf(i,j) &
                                ) / Deltax2

            rhs(j) = rhs(j) + ( etaBf(i+1,j) * vtp(i-1,j) / 3d0 &
                            +   etaBf(i,j)   * vtp(i-1,j) &
                              ) / Deltax2

!     . o o -- open boundary to the left

         elseif ( maskB(i,j) .eq. 1 .and. i .eq. 1 .and. .not. stressBC) then

            bVEC(j) = bVEC(j) + etaBf(i+1,j) / Deltax2

            rhs(j) = rhs(j) + ( etaBf(i+1,j) * vtp(i+1,j) &
                              ) / Deltax2

!     o o . -- open boundary to the right

         elseif ( maskB(i+1,j) .eq. 1 .and. i .eq. nx .and. .not. stressBC) then

            bVEC(j) = bVEC(j) + etaBf(i,j) / Deltax2

            rhs(j) = rhs(j) + ( etaBf(i,j) * vtp(i-1,j)   &
                              ) / Deltax2

!     o o o -- normal situation

         else

            if (stressBC .and. j .eq. 1) then ! S       

               bVEC(j) = bVEC(j) + 0d0

            elseif (stressBC .and. j .eq. ny+1) then ! N

               bVEC(j) = bVEC(j) + 0d0

            elseif (stressBC .and. i .eq. 1 .and. j .ge. 2 .and. j .le. ny) then ! W 

               bVEC(j) = bVEC(j) + ( etaBf(i+1,j) ) / Deltax2

               rhs(j) = rhs(j) + ( etaBf(i+1,j) * vtp(i+1,j) ) / Deltax2

            elseif (stressBC .and. i .eq. nx .and. j .ge. 2 .and. j .le. ny) then ! E

               bVEC(j) = bVEC(j) + ( etaBf(i,j) ) / Deltax2

               rhs(j) = rhs(j) + ( etaBf(i,j)   * vtp(i-1,j) &
                              ) / Deltax2

            else

               bVEC(j) = bVEC(j) + ( etaBf(i+1,j) + etaBf(i,j) &
                                ) / Deltax2

               rhs(j) = rhs(j) + ( etaBf(i+1,j) * vtp(i+1,j) &
                            +   etaBf(i,j)   * vtp(i-1,j) &
                              ) / Deltax2

            endif

         endif


!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------

         if (.not. stressBC) then

         if ( j .eq. 1 ) then

            if (maskC(i,3) .eq. 1) then
            aVEC(j) = 0.0d0
            bVEC(j) = 1.0d0
            cVEC(j) = -4d0/3d0

            rhs(j)  = rhsv(i,j) - vtp(i,3) / 3d0
            else
            aVEC(j) = 0.0d0
            bVEC(j) = 1.0d0
            cVEC(j) = -1d0

            rhs(j)  = rhsv(i,j)
            endif

         elseif ( j .eq. ny+1 ) then

            if (maskC(i,ny-2) .eq. 1) then
            aVEC(j) = -4d0/3d0
            bVEC(j) = 1.0d0
            cVEC(j) = 0.0d0

            rhs(j) = rhsv(i,j) - vtp(i,ny-1) / 3d0
            else
            aVEC(j) = -1d0
            bVEC(j) = 1.0d0
            cVEC(j) = 0.0d0

            rhs(j) = rhsv(i,j)
            endif

         endif
      endif


 200     continue

      enddo

      call tridag (aVEC, bVEC, cVEC, rhs, Y, ny+1)

      do j = 1, ny+1

         vtp(i,j) = Yold(j) + wlsor*(Y(j)-Yold(j))

      enddo

      return
    end subroutine vectv_3diag


!************************************************************************
!     Subroutine tridag: Solve tri-diagonal system of equations (A x = r).
!       a : lower diagonal element of A
!       b : diagonal element of A
!       c : upper diagonal element of A
!       r : right hand side of the equation
!       u : solution x to the system of equation
!       n : dimension of the matrix A (n x n)
!************************************************************************

      subroutine tridag (a, b, c, r, u, n)
      implicit double precision (a-h, o-z)

      parameter (nmax = 4000)

      dimension gam(nmax) , a(n), b(n), c(n), r(n), u(n)

      if (b(1) .eq. 0d0) then

         u(1) = 0.0d0
         bet  = 1.0d0

      else

         bet  = b(1)
         u(1) = r(1) / bet

      endif

      do j = 2, n
         gam(j) = c(j-1) / bet
        
         if (b(j) .eq. 0d0) then

            u(j) = 0.0d0
            bet  = 1.0d0

         else
            bet = b(j) - a(j) * gam(j)
            u(j) = (r(j) - a(j) * u(j-1)) / bet
         
         endif

      enddo

      do 12 j = n-1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
 12   continue
      
      return
      end


