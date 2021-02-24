!*********************************************************************                 
!     approximates Jv by F(x+ev)-F(x)/e                                    
!*********************************************************************    

      subroutine JacfreeVec (v, Jv, Fneg, epsilon)

      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_options.h'

      integer i

      double precision xtp(nvar), x(nvar),rhs(nvar)
      double precision Fpos(nvar),Fneg(nvar)
      double precision epsilon,v(nvar),Jv(nvar)
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      call transformer (uice,vice,x,1)
      xtp = x + epsilon * v
      call transformer (utp,vtp,xtp,0)

      if ( IMEX .eq. 2 ) then ! IMEX method 2                                  
         call advection ( un1, vn1, utp, vtp, hn2, An2, hn1, An1, h, A )
         call Ice_strength()
         call bvect_ind
      endif

      call ViscousCoefficient(utp,vtp)
      call bvect (utp,vtp,rhs)
      call Funk (xtp,rhs,Fpos)

      if (Jac_finite_diff .eq. 'centred') then

         print *, 'Might not work because of replacement closure'
         stop
         xtp = x - epsilon * v
         call transformer (utp,vtp,xtp,0)

         call ViscousCoefficient(utp,vtp)
         call bvect (utp,vtp,rhs)
         call Funk (xtp,rhs,Fneg)

         do i=1, nvar

            Jv(i) = ( Fpos(i)-Fneg(i) ) / ( 2 * epsilon )

         enddo

      elseif (Jac_finite_diff .eq. 'forward') then

         ! Funeg is already in Fneg                                                                                                               
         do i=1, nvar

            Jv(i) = ( Fpos(i)-Fneg(i) ) / epsilon

         enddo

      else

         print *, 'Wrong choice of finite diff'

      endif

      return
    end subroutine JacfreeVec

!****************************************************************************                                              
!     calculates the vector given by Ax - rhs                                                                   
!**************************************************************************** 

subroutine Funk (xtp,rhs,Fout)

  implicit none

  include 'parameter.h'

  double precision, intent(in) :: xtp(nvar), rhs(nvar)
  double precision, intent(out) :: Fout(nvar)

  call matvec(xtp,Fout) ! Fout here is Ax

  Fout = Fout - rhs ! substract rhs vector                                                                     

  return
end subroutine Funk


!****************************************************************************
!     calculates matrix A times a vector x and outputs the vector Ax
!****************************************************************************

      subroutine MATVEC (x,Ax)

      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j, k

      double precision, intent(in) :: x(nvar)
      double precision, intent(out) :: Ax(nvar)
      double precision hvert, uavg, vavg
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      call transformer(utp,vtp,x,0)

      do j = 1, ny
         do i = 1, nx+1

            k = i+(j-1)*(nx+1)

            Ax(k) = 0.0d0

            if ( maskB(i,j) + maskB(i,j+1) .eq. 0 ) goto 100

!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------

            vavg = ( vtp(i,j)   + vtp(i,j+1) &
                 + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

            hvert = ( h(i,j) + h(i-1,j) ) / 2d0

            Ax(k) = Ax(k) - rhof * hvert * vavg

!------------------------------------------------------------------------
!    rhoice*h*du/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            if ( BDF .eq. 0 ) then
               Ax(k) = Ax(k) + ( rhoice * hvert * utp(i,j) ) / Deltat
            elseif ( BDF .eq. 1 ) then
               Ax(k) = Ax(k) + ( 3d0*rhoice * hvert * utp(i,j) ) /(2d0* Deltat)
            endif

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (u - uw) cos(theta_w) - (v - vw) sin(theta_w) ] }
!     drag_x does not need to be within the loop, since it does not change
!------------------------------------------------------------------------

            Ax(k) = Ax(k) + &
                      CdwC1(i,j) * costheta_w * utp(i,j) - &
                      CdwC1(i,j) * sintheta_w * vavg 

            Ax(k) = Ax(k) + Cbasal1(i,j)*utp(i,j)

!------------------------------------------------------------------------
!     -d ( eta dv/dy ) / dx    B1_1  p.899... 
!------------------------------------------------------------------------

            Ax(k) = Ax(k) + &
                       ( etaC(i,j)   * ( vtp(i,j+1)   - vtp(i,j)   ) &
                       - etaC(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j) ) &
                       ) /  Deltax2

!------------------------------------------------------------------------
!     d ( zeta dv/dy ) / dx    B1_2
!------------------------------------------------------------------------
            
            Ax(k) = Ax(k) - &
                       ( zetaC(i,j)   * ( vtp(i,j+1)   - vtp(i,j)   ) &
                       - zetaC(i-1,j) * ( vtp(i-1,j+1) - vtp(i-1,j) ) &
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
               Ax(k) = Ax(k) - &
                       ( etaB(i,j+1) * &
                       ( 3d0 * vtp(i,j+1) - vtp(i+1,j+1)/ 3d0 ) &
                       - etaB(i,j) * ( vtp(i,j) - vtp(i-1,j))   &
                        ) / Deltax2

            elseif(maskC(i-1,j+1) .eq. 1 .and. maskC(i,j+1) .eq. 0) then

!     o x
!     o o  
!     o o
               Ax(k) = Ax(k) - &
                       ( etaB(i,j+1) * &
                       ( vtp(i-2,j+1)/ 3d0 - 3d0 * vtp(i-1,j+1)) &
                       - etaB(i,j) * ( vtp(i,j) - vtp(i-1,j))    &
                        ) / Deltax2

            elseif(maskC(i-1,j-1) .eq. 0 .and. maskC(i,j-1) .eq. 1) then

!     o o
!     o o  
!     x o
               Ax(k) = Ax(k) - &
                       ( etaB(i,j+1) *(vtp(i,j+1) - vtp(i-1,j+1)) &
                       - etaB(i,j) * &
                       ( 3d0 * vtp(i,j) - vtp(i+1,j) / 3d0 ) &
                        ) / Deltax2

         
            elseif(maskC(i-1,j-1) .eq. 1 .and. maskC(i,j-1) .eq. 0) then
                  
!     o o
!     o o  
!     o x

               Ax(k) = Ax(k) - &
                       ( etaB(i,j+1) *(vtp(i,j+1) - vtp(i-1,j+1)) &
                       - etaB(i,j) * &
                       ( vtp(i-2,j) / 3d0 - 3d0 * vtp(i-1,j) ) &
                        ) / Deltax2

            else 

!     o o                   x x    o o
!     o o  --normal case or o o or o o
!     o o                   o o    x x
               
               Ax(k) = Ax(k) - &
                       ( etaB(i,j+1) * ( vtp(i,j+1) - vtp(i-1,j+1)) &
                       + etaB(i,j)   * ( vtp(i-1,j) - vtp(i,j)     )&
                        ) / Deltax2

            endif

!------------------------------------------------------------------------
!     d [ (eta + zeta ) du/dx ] / dx      B1_4, D1_4 
!------------------------------------------------------------------------

              Ax(k) = Ax(k) + &
                  utp(i,j) * ( etaC(i,j)  + etaC(i-1,j) &
                              + zetaC(i,j) + zetaC(i-1,j) &
                              ) / Deltax2   - &
                  ( etaC(i,j)    * utp(i+1,j) &
                  + etaC(i-1,j)  * utp(i-1,j) &
                  + zetaC(i,j)   * utp(i+1,j) &
                  + zetaC(i-1,j) * utp(i-1,j) &
                   ) / Deltax2
               
!------------------------------------------------------------------------
!     d ( eta du/dy ) / dy   B1_5, D1_5
!------------------------------------------------------------------------

!     o o
!     o o   -- land just below
!     x x

           if ( maskB(i,j) .eq. 0 ) then

              Ax(k) = Ax(k) + &
                  utp(i,j) * ( etaB(i,j+1) + 3d0 * etaB(i,j) &
                              ) / Deltax2  - &
                  ( etaB(i,j+1) * utp(i,j+1) &
                  + etaB(i,j)   * utp(i,j+1) / 3d0 &
                   ) / Deltax2


!     x x
!     o o   -- land just above
!     o o

           elseif ( maskB(i,j+1) .eq. 0 ) then

              Ax(k) = Ax(k) + &
                  utp(i,j) * ( 3d0 * etaB(i,j+1) + etaB(i,j) &
                              ) / Deltax2  - &
                  ( etaB(i,j+1) * utp(i,j-1) / 3d0  &
                  + etaB(i,j)   * utp(i,j-1) &
                   ) / Deltax2


!     o o
!     o o   -- open boundary just below
!     . .

           elseif ( maskB(i,j) .eq. 1 .and. j .eq. 1) then

              Ax(k) = Ax(k) + &
                  utp(i,j) * etaB(i,j+1) / Deltax2 - &
                  ( etaB(i,j+1) * utp(i,j+1) ) / Deltax2

!     . .
!     o o   -- open boundary just above
!     o o

           elseif ( maskB(i,j+1) .eq. 1 .and. j .eq. ny) then

              Ax(k) = Ax(k) + &
                  utp(i,j) * etaB(i,j) / Deltax2 - &
                  ( etaB(i,j) * utp(i,j-1) ) / Deltax2

!     o o
!     o o   -- normal situation
!     o o

           else

              Ax(k) = Ax(k) + &
                  utp(i,j) * ( etaB(i,j+1) + etaB(i,j) ) / Deltax2 - &
                  ( etaB(i,j+1) * utp(i,j+1) &
                  + etaB(i,j)   * utp(i,j-1) &
                   ) / Deltax2

           endif


!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------


              if ( i .eq. 1 ) then

                 if (maskC(3,j) .eq. 1) then ! oooo (2nd order accurate)        

                    Ax(k) = utp(1,j) - ( 4d0 * utp(2,j)  - utp(3,j) ) &
                         / 3d0

                 else ! ooox (mask should be built so it does not occur)        

                    Ax(k) = utp(1,j) - utp(2,j)

                 endif


              elseif ( i .eq. nx+1 ) then
   
                 if (maskC(nx-2,j) .eq. 1) then

                    Ax(k) = utp(nx+1,j) - ( 4d0 * utp(nx,j) &
                         - utp(nx-1,j) ) / 3d0

                 else ! xooo (1st order accurate)                               
                    Ax(k) = utp(nx+1,j) - utp(nx,j)
                 endif

              endif


 100       continue

         enddo
      enddo

!------------------------------------------------------------------------
!     Calculate the v-component
!------------------------------------------------------------------------
      
      do j = 1, ny+1
         do i = 1, nx

            k = i+(j-1)*nx+(nx+1)*ny

            Ax(k) = 0.0d0

            if ( maskB(i,j) + maskB(i+1,j) .eq. 0 ) goto 200

!------------------------------------------------------------------------
!     Coriolis term
!------------------------------------------------------------------------

            uavg  = ( utp(i,j)   + utp(i+1,j) &
                    + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

            hvert = ( h(i,j) + h(i,j-1) ) / 2d0

            Ax(k) = Ax(k) + rhof * hvert * uavg

!------------------------------------------------------------------------
!    rhoice*h*dv/dt (tendency term. Advection of momentum is neglected)
!------------------------------------------------------------------------

            if ( BDF .eq. 0 ) then
               Ax(k) = Ax(k) + ( rhoice * hvert * vtp(i,j) ) / Deltat
            elseif ( BDF .eq. 1 ) then
               Ax(k) = Ax(k) + ( 3d0*rhoice * hvert * vtp(i,j) ) /(2d0* Deltat)
            endif

!------------------------------------------------------------------------
!     Drag  { -Cdw' [ (v - vw) * cos(theta_w) + (u - uw) * sin(theta_w) ] }
!------------------------------------------------------------------------

           Ax(k) = Ax(k) + &
                     CdwC2(i,j) * costheta_w * vtp(i,j) + &
                     CdwC2(i,j) * sintheta_w * uavg

           Ax(k) = Ax(k) + Cbasal2(i,j)*vtp(i,j)

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
              Ax(k) = Ax(k) - &
                      ( etaB(i+1,j) * ( utp(i+1,j) - utp(i+1,j-1)) &
                      - etaB(i,j)   * &
                      ( 3d0 * utp(i,j) - utp(i,j+1) / 3d0 ) &
                       ) / Deltax2 

           elseif (maskC(i-1,j-1) .eq. 1 .and. maskC(i-1,j) .eq. 0) then 

!     x o o 
!     o o o 
              Ax(k) = Ax(k) - &
                      ( etaB(i+1,j) * ( utp(i+1,j) - utp(i+1,j-1)) &
                      - etaB(i,j)   * &
                      ( utp(i,j-2) / 3d0 - 3d0 * utp(i,j-1) ) &
                       ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 0 .and. maskC(i+1,j) .eq. 1) then 

!     o o o 
!     o o x  
              Ax(k) = Ax(k) - &
                      ( etaB(i+1,j) * &
                      ( 3d0 * utp(i+1,j) - utp(i+1,j+1) / 3d0 ) &
                      - etaB(i,j) * ( utp(i,j) - utp(i,j-1) ) &
                       ) / Deltax2

           elseif (maskC(i+1,j-1) .eq. 1 .and. maskC(i+1,j) .eq. 0) then 

!     o o x 
!     o o o
              Ax(k) = Ax(k) - &
                      ( etaB(i+1,j) * &
                      ( utp(i+1,j-2) / 3d0 - 3d0 * utp(i+1,j-1) ) &
                      - etaB(i,j) * ( utp(i,j) - utp(i,j-1) ) &
                       ) / Deltax2

           else

!     o o o -- normal case or x o o or o o x
!     o o o                   x o o    o o x

              Ax(k) = Ax(k) - &
                      ( etaB(i+1,j) * ( utp(i+1,j) - utp(i+1,j-1)) &
                      - etaB(i,j)   * ( utp(i,j)   - utp(i,j-1) ) &
                       ) / Deltax2

           endif

!------------------------------------------------------------------------
!     -d ( eta du/dx ) / dy    B2_1
!------------------------------------------------------------------------

           Ax(k) = Ax(k) + &
               ( etaC(i,j)   * ( utp(i+1,j)   - utp(i,j)   ) &
               - etaC(i,j-1) * ( utp(i+1,j-1) - utp(i,j-1) ) &
                ) / Deltax2

!------------------------------------------------------------------------
!     d (zeta du/dx) / dy      B2_2
!------------------------------------------------------------------------

           Ax(k) = Ax(k) - &
               ( zetaC(i,j)   * ( utp(i+1,j)   - utp(i,j)   ) &
               - zetaC(i,j-1) * ( utp(i+1,j-1) - utp(i,j-1) ) &
                ) / Deltax2

!------------------------------------------------------------------------
!     d [ (eta + zeta) dv/dy ] / dy   D2_4, B2_4
!------------------------------------------------------------------------

               Ax(k) = Ax(k) + &
                   vtp(i,j) * ( etaC(i,j)  + etaC(i,j-1) &
                               + zetaC(i,j) + zetaC(i,j-1) &
                                ) / Deltax2 - &
                 ( etaC(i,j)    * vtp(i,j+1) &
                 + etaC(i,j-1)  * vtp(i,j-1) &  
                 + zetaC(i,j)   * vtp(i,j+1) &  
                 + zetaC(i,j-1) * vtp(i,j-1) &
                  ) / Deltax2
   
!------------------------------------------------------------------------
!     d ( eta dv/dx ) / dx   D2_5, B2_5
!------------------------------------------------------------------------

!     x o o -- land to the left

            if ( maskB(i,j) .eq. 0 ) then

               Ax(k) = Ax(k) + &
                   vtp(i,j) * ( etaB(i+1,j) + 3d0 * etaB(i,j) ) &
                               / Deltax2 - &
                 ( etaB(i+1,j) * vtp(i+1,j) &
                 + etaB(i,j)   * vtp(i+1,j) / 3d0  &
                  ) / Deltax2 

!     o o x -- land to the right

           elseif ( maskB(i+1,j) .eq. 0 ) then

               Ax(k) = Ax(k) + &
                   vtp(i,j) * ( 3d0 * etaB(i+1,j) + etaB(i,j) ) &
                               / Deltax2         - &
                 ( etaB(i+1,j) * vtp(i-1,j) / 3d0  &
                 + etaB(i,j)   * vtp(i-1,j)        &
                  ) / Deltax2

!     . o o -- open boundary to the left

           elseif ( maskB(i,j) .eq. 1 .and. i .eq. 1) then

              Ax(k) = Ax(k) + &
                  vtp(i,j) * etaB(i+1,j) / Deltax2 - &
                ( etaB(i+1,j) * vtp(i+1,j) ) / Deltax2 

!     o o . -- open boundary to the right

           elseif ( maskB(i+1,j) .eq. 1 .and. i .eq. nx) then

              Ax(k) = Ax(k) + &
                  vtp(i,j) * etaB(i,j) / Deltax2 - &
                ( etaB(i,j) * vtp(i-1,j) ) / Deltax2

!     o o o -- normal situation

           else

              Ax(k) = Ax(k) + &
                  vtp(i,j) * ( etaB(i+1,j) + etaB(i,j) ) / Deltax2 - &
                ( etaB(i+1,j) * vtp(i+1,j) &
                + etaB(i,j)   * vtp(i-1,j) &
                 ) / Deltax2
           endif


!------------------------------------------------------------------------
!     Calculate the ice velocity using relaxation and apply boundary cond
!     open  boundary : du/dn = 0
!     close boundary : u_n   = 0
!     .                u_t   = 0, only for viscous plastic, used in eta calc
!------------------------------------------------------------------------
 

              if ( j .eq. 1 ) then

                 if (maskC(i,3) .eq. 1) then ! (2nd order accurate)             

                    Ax(k) = vtp(i,1) - ( 4d0 * vtp(i,2)  - vtp(i,3) ) &
                         / 3d0
                 else
                    Ax(k) = vtp(i,1)-vtp(i,2) ! (1st order accurate)            
                 endif

              elseif ( j .eq. ny+1 ) then
        
                 if (maskC(i,ny-2) .eq. 1) then ! (2nd order accurate)          
                    Ax(k) = vtp(i,ny+1) - ( 4d0 * vtp(i,ny) &
                         - vtp(i,ny-1) ) / 3d0
                 else
                    Ax(k) = vtp(i,ny+1)-vtp(i,ny) ! (1st order accurate)        
                 endif

              endif

 200       continue

         enddo
      enddo

      return
    end subroutine MATVEC






