
!***********************************************************************
!     subroutine advection (upwind scheme or upwind-RK2 or semilag):
!       calculates the tracer quantities at the next time step. 
!       Tracer #1 and #2 are ice thickness and concentration respectively.
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             14-05-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!     V4.0            20-09-2012             JF Lemieux
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!     Lemieux, J.-F., Knoll, D.A., Losch, M. and Girard, C., A second-order 
!     accurate in time IMplicit–EXplicit (IMEX) integration scheme for 
!     sea ice dynamics, Journal of Computational Physics, 2014.
!
!************************************************************************

      subroutine advection (un1, vn1, un, vn, hn2, An2, hn1, An1, hout, Aout)
      
      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j, k, caseSL, iloc, jloc
      integer isw, jsw, inw, jnw, ine, jne, ise, jse !SL SouthWest=sw, nw, ne, se corners
      double precision, intent(in)    :: un1(0:nx+2,0:ny+2), vn1(0:nx+2,0:ny+2)
      double precision, intent(in)    :: un(0:nx+2,0:ny+2), vn(0:nx+2,0:ny+2)
      double precision                :: hn1(0:nx+1,0:ny+1), An1(0:nx+1,0:ny+1)
      double precision, intent(in)    :: hn2(0:nx+1,0:ny+1), An2(0:nx+1,0:ny+1)
      double precision, intent(out)   :: hout(0:nx+1,0:ny+1), Aout(0:nx+1,0:ny+1)
      double precision                :: ustar(0:nx+2,0:ny+2), vstar(0:nx+2,0:ny+2)
      double precision                :: hstar(0:nx+1,0:ny+1), Astar(0:nx+1,0:ny+1)
      double precision                :: dFx(nx,ny), dFy(nx,ny), div(nx,ny)
      double precision                :: alphamx, alphamy, uinterp, vinterp, ftp
      double precision                :: hbef, Abef, xd, yd, xdn1, ydn1, xdn2, ydn2 
      double precision                :: fsw, fnw, fne, fse
      double precision                :: fxsw, fxnw, fxne, fxse, fysw, fynw, fyne, fyse
      double precision                :: fxysw, fxynw, fxyne, fxyse
      double precision                :: upper, lower, rhsdynh, rhsdynA
      double precision                :: fluxx, fluxy
      double precision                :: fx, fy, fxy, cubic_interp         ! function
      double precision                :: calc_fluxx, calc_fluxy, apply_lim ! functions
      logical                         :: SLlimiter

!------------------------------------------------------------------------ 
!     set dhn1/dx, dAn1/dx = 0 at the outside cell when there is an open bc 
!------------------------------------------------------------------------ 

      do i = 0, nx+1
               
         if (maskC(i,0) .eq. 1) then
            
            hn1(i,0) = ( 4d0 * hn1(i,1) - hn1(i,2) )/3d0
            hn1(i,0) = max(hn1(i,0), 0d0)
            An1(i,0) = ( 4d0 * An1(i,1) - An1(i,2) )/3d0
            An1(i,0) = max(An1(i,0), 0d0)
            An1(i,0) = min(An1(i,0), 1d0)
                  
         endif

         if (maskC(i,ny+1) .eq. 1) then
            
            hn1(i,ny+1)= ( 4d0 * hn1(i,ny) - hn1(i,ny-1) ) / 3d0
            hn1(i,ny+1)= max(hn1(i,ny+1), 0d0)
            An1(i,ny+1)= ( 4d0 * An1(i,ny) - An1(i,ny-1) ) / 3d0
            An1(i,ny+1)= max(An1(i,ny+1), 0d0)
            An1(i,ny+1)= min(An1(i,ny+1), 1d0)

         endif
 
      enddo

      do j = 0, ny+1

         if (maskC(0,j) .eq. 1) then
                  
            hn1(0,j)  = ( 4d0 * hn1(1,j) - hn1(2,j) ) / 3d0
            hn1(0,j)  = max(hn1(0,j), 0d0)
            An1(0,j)  = ( 4d0 * An1(1,j) - An1(2,j) ) / 3d0
            An1(0,j)  = max(An1(0,j), 0d0)
            An1(0,j)  = min(An1(0,j), 1d0)

         endif

         if (maskC(nx+1,j) .eq. 1) then
                 
            hn1(nx+1,j) = ( 4d0 * hn1(nx,j) - hn1(nx-1,j) ) / 3d0
            hn1(nx+1,j) = max(hn1(nx+1,j), 0d0)
            An1(nx+1,j) = ( 4d0 * An1(nx,j) - An1(nx-1,j) ) / 3d0
            An1(nx+1,j) = max(An1(nx+1,j), 0d0)
            An1(nx+1,j) = min(An1(nx+1,j), 1d0)

         endif

      enddo

      if ( adv_scheme .eq. 'upwind' ) then

!------------------------------------------------------------------------
!     compute the difference of the flux for thickness 
!------------------------------------------------------------------------

         call calc_dFx (un, hn1, dFx)
         call calc_dFy (vn, hn1, dFy)

!------------------------------------------------------------------------
!     update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hout(i,j) = hn1(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  hout(i,j) = max(hout(i,j), 0d0)

               endif
                     
            enddo
         enddo

!------------------------------------------------------------------------  
!     compute the difference of the flux for concentration                                 
!------------------------------------------------------------------------                  

         call calc_dFx (un, An1, dFx)
         call calc_dFy (vn, An1, dFy)

!------------------------------------------------------------------------ 
!     update the concentration values      
!     (in a separate do-loop to conserve mass)                                             
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  Aout(i,j) = An1(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  Aout(i,j) = max(Aout(i,j), 0d0)
                  Aout(i,j) = min(Aout(i,j), 1d0)

               endif

            enddo
         enddo

      elseif ( adv_scheme .eq. 'upwindRK2' ) then
               
!------------------------------------------------------------------------
!     predictor: compute the difference of the flux for thickness 
!------------------------------------------------------------------------

         call calc_dFx (un1, hn1, dFx)
         call calc_dFy (vn1, hn1, dFy)

!------------------------------------------------------------------------
!     predictor: update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hstar(i,j) = hn1(i,j) - (DtoverDx / 2d0) * ( dFx(i,j) + dFy(i,j) )
                  hstar(i,j) = max(hstar(i,j), 0d0)

               endif
                     
            enddo
         enddo

!------------------------------------------------------------------------  
!     predictor: compute the difference of the flux for concentration   
!------------------------------------------------------------------------                  

         call calc_dFx (un1, An1, dFx)
         call calc_dFy (vn1, An1, dFy)

!------------------------------------------------------------------------ 
!     predictor: update the concentration values      
!     (in a separate do-loop to conserve mass)                                             
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  Astar(i,j) = An1(i,j) - (DtoverDx / 2d0) * ( dFx(i,j) + dFy(i,j) )
                  Astar(i,j) = max(Astar(i,j), 0d0)
                  Astar(i,j) = min(Astar(i,j), 1d0)

               endif

            enddo
         enddo

!------------------------------------------------------------------------ 
!     set dhstar/dx, dAstar/dx = 0 at the outside cell when there is an open bc 
!------------------------------------------------------------------------ 

         do i = 0, nx+1
            
            if (maskC(i,0) .eq. 1) then

               hstar(i,0) = ( 4d0 * hstar(i,1) - hstar(i,2) )/3d0
               hstar(i,0) = max(hstar(i,0), 0d0)
               Astar(i,0) = ( 4d0 * Astar(i,1) - Astar(i,2) )/3d0
               Astar(i,0) = max(Astar(i,0), 0d0)
               Astar(i,0) = min(Astar(i,0), 1d0)
                  
            endif

            if (maskC(i,ny+1) .eq. 1) then
                  
               hstar(i,ny+1)= ( 4d0 * hstar(i,ny) - hstar(i,ny-1) ) / 3d0
               hstar(i,ny+1)= max(hstar(i,ny+1), 0d0)
               Astar(i,ny+1)= ( 4d0 * Astar(i,ny) - Astar(i,ny-1) ) / 3d0
               Astar(i,ny+1)= max(Astar(i,ny+1), 0d0)
               Astar(i,ny+1)= min(Astar(i,ny+1), 1d0)

            endif
 
         enddo

         do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then
                  
               hstar(0,j)  = ( 4d0 * hstar(1,j) - hstar(2,j) ) / 3d0
               hstar(0,j)  = max(hstar(0,j), 0d0)
               Astar(0,j)  = ( 4d0 * Astar(1,j) - Astar(2,j) ) / 3d0
               Astar(0,j)  = max(Astar(0,j), 0d0)
               Astar(0,j)  = min(Astar(0,j), 1d0)

            endif

            if (maskC(nx+1,j) .eq. 1) then
                 
               hstar(nx+1,j) = ( 4d0 * hstar(nx,j) - hstar(nx-1,j) ) / 3d0
               hstar(nx+1,j) = max(hstar(nx+1,j), 0d0)
               Astar(nx+1,j) = ( 4d0 * Astar(nx,j) - Astar(nx-1,j) ) / 3d0
               Astar(nx+1,j) = max(Astar(nx+1,j), 0d0)
               Astar(nx+1,j) = min(Astar(nx+1,j), 1d0)

            endif

         enddo

!------------------------------------------------------------------------
!     corrector: compute the difference of the flux for thickness 
!------------------------------------------------------------------------
  
         ustar = ( un1 + un ) / 2d0
         vstar = ( vn1 + vn ) / 2d0

         call calc_dFx (ustar, hstar, dFx)
         call calc_dFy (vstar, hstar, dFy)

!------------------------------------------------------------------------
!     corrector: update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hout(i,j) = hn1(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  hout(i,j) = max(hout(i,j), 0d0)

               endif
                     
            enddo
         enddo

!------------------------------------------------------------------------  
!     corrector: compute the difference of the flux for concentration                                 
!------------------------------------------------------------------------                  

         call calc_dFx (ustar, Astar, dFx)
         call calc_dFy (vstar, Astar, dFy)

!------------------------------------------------------------------------ 
!     corrector: update the concentration values      
!     (in a separate do-loop to conserve mass)                                             
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  Aout(i,j) = An1(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  Aout(i,j) = max(Aout(i,j), 0d0)
                  Aout(i,j) = min(Aout(i,j), 1d0)

               endif

            enddo
         enddo

      elseif ( adv_scheme .eq. 'semilag' ) then
         
!------------------------------------------------------------------------ 
!     Semi-Lagrangian scheme for advection. 
!     This is a 3 time level scheme (h is obtained from hn1 and hn2)
!     
!     Staniforth and Côté, Monthly Weather Review 1991.
!     Pellerin et al, Monthly Weather Review 1995. 
!
!------------------------------------------------------------------------ 

! caseSL=3 : no advection (land)
! caseSL=  : IMPROVE THIS
! TO DO : restart with semilag...need hn2, An2

         SLlimiter=.true.
         
         call calc_div(un1,vn1,div) ! calc divergence at n-1 for RHS [hdiv(u)]^{n-1}

         do i = 1, nx
            do j = 1, ny
               caseSL=3
               if (maskC(i,j) .eq. 1) then
                  caseSL=1 ! SL
                  if (i .lt. 3 .or. i .gt. nx-2 .or. j .lt. 3 .or. j .gt. ny-2) then
                     caseSL=2 ! upwind
                  else
                     do jloc=j-2, j+2
                        do iloc=i-2, i+2
                           if (maskC(iloc,jloc) .eq. 0) caseSL=2 ! upwind
                        enddo
                     enddo
                  endif
               endif
               
               if (caseSL == 1) then ! SL advection

!------------------------------------------------------------------------  
! find distances alphamx and alphamy of particle from tracer(i,j) at t=n-1
! xd, yd are calculated from sw corner of the T(i,j) cell (h, A at center). 
! The four corners are the four corners of the T-cell.
!
! nw--------ne
! |         |
!uij  Tij   | 
! |         |
! sw--vij---se
!  
!------------------------------------------------------------------------
                  alphamx=0.01d0 ! initial value
                  alphamy=0.01d0 ! initial value

                  do k = 1, 5 ! implicit loop to find alphamx and alphamy

                     xd = 0.5d0 - alphamx / Deltax ! same wether alphamx is + or -
                     yd = 0.5d0 - alphamy / Deltax

!--- interpolate u at x-alphaxm, y-alphamy

                     fsw = ( un1(i,j) + un1(i,j-1) ) / 2d0 
                     fnw = ( un1(i,j+1) + un1(i,j) ) / 2d0
                     fne = ( un1(i+1,j+1) + un1(i+1,j) ) / 2d0
                     fse = ( un1(i+1,j) + un1(i+1,j-1) ) / 2d0

                     ftp= ( un1(i-1,j) + un1(i-1,j-1) ) / 2d0
                     fxsw=fx(fse, ftp, 2d0)
                     ftp= ( un1(i-1,j+1) + un1(i-1,j) ) / 2d0
                     fxnw=fx(fne, ftp, 2d0)
                     ftp= ( un1(i+2,j+1) + un1(i+2,j) ) / 2d0
                     fxne=fx(ftp, fnw, 2d0)
                     ftp= ( un1(i+2,j) + un1(i+2,j-1) ) / 2d0
                     fxse=fx(ftp, fsw, 2d0)

                     fysw=fy(un1(i,j), un1(i,j-1), 1d0)
                     fynw=fy(un1(i,j+1), un1(i,j), 1d0)
                     fyne=fy(un1(i+1,j+1), un1(i+1,j), 1d0)
                     fyse=fy(un1(i+1,j), un1(i+1,j-1), 1d0)

                     fxysw=fxy(un1(i+1,j), un1(i-1,j), un1(i+1,j-1), un1(i-1,j-1), 2d0)
                     fxynw=fxy(un1(i+1,j+1), un1(i-1,j+1), un1(i+1,j), un1(i-1,j), 2d0)
                     fxyne=fxy(un1(i+2,j+1), un1(i,j+1), un1(i+2,j), un1(i,j), 2d0)
                     fxyse=fxy(un1(i+2,j), un1(i,j), un1(i+2,j-1), un1(i,j-1), 2d0)

                     uinterp=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                           fxsw, fxnw, fxne, fxse,  &
                                           fysw, fynw, fyne, fyse,  &
                                           fxysw,fxynw,fxyne,fxyse, &
                                           xd, yd )

!--- interpolate v at x-alphaxm, y-alphamy                                                           

                     fsw = ( vn1(i,j) + vn1(i-1,j) ) / 2d0
                     fnw = ( vn1(i,j+1) + vn1(i-1,j+1) ) / 2d0
                     fne = ( vn1(i,j+1) + vn1(i+1,j+1) ) / 2d0
                     fse = ( vn1(i+1,j) + vn1(i,j) ) / 2d0

                     fxsw=fx(vn1(i,j), vn1(i-1,j), 1d0)
                     fxnw=fx(vn1(i,j+1), vn1(i-1,j+1), 1d0)
                     fxne=fx(vn1(i+1,j+1), vn1(i,j+1), 1d0)
                     fxse=fx(vn1(i+1,j), vn1(i,j), 1d0)

                     ftp= ( vn1(i,j-1) + vn1(i-1,j-1) ) / 2d0
                     fysw=fy(fnw, ftp, 2d0)
                     ftp= ( vn1(i,j+2) + vn1(i-1,j+2) ) / 2d0
                     fynw=fy(ftp, fsw, 2d0)
                     ftp= ( vn1(i+1,j+2) + vn1(i,j+2) ) / 2d0
                     fyne=fy(ftp, fse, 2d0)
                     ftp= ( vn1(i+1,j-1) + vn1(i,j-1) ) / 2d0
                     fyse=fy(fne, ftp, 2d0)

                     fxysw=fxy(vn1(i,j+1), vn1(i-1,j+1), vn1(i,j-1), vn1(i-1,j-1), 2d0)
                     fxynw=fxy(vn1(i,j+2), vn1(i-1,j+2), vn1(i,j), vn1(i-1,j), 2d0)
                     fxyne=fxy(vn1(i+1,j+2), vn1(i,j+2), vn1(i+1,j), vn1(i,j), 2d0)
                     fxyse=fxy(vn1(i+1,j+1), vn1(i,j+1), vn1(i+1,j-1), vn1(i,j-1), 2d0)

                     vinterp=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                           fxsw, fxnw, fxne, fxse,  &
                                           fysw, fynw, fyne, fyse,  &
                                           fxysw,fxynw,fxyne,fxyse, &
                                           xd, yd )

! Can I use the latest alphamx to get vinterp and then alphamy???                                 
! (kind of Gauss-Seidel vs Jacobi) ...move it after uinterp

                     alphamx=Deltat*uinterp
                     alphamy=Deltat*vinterp

                  enddo
!------------------------------------------------------------------------
! find hbef and Abef (initial position of particle at time level n-2=n2)
!------------------------------------------------------------------------

! identify coordinates of 4 corners. These corners are 4 tracer points.

! Tnw--------Tne
!  |         |
!  |         | 
!  |         | 
! Tsw--------Tse 

! xd and yd are distances in the interval [0,1]. They are calc from the sw corner 
! xdn2, ydn2 is the position of the particle at n-2 (with respect to the sw corner)
! xdn1, ydn1 is the position of the particle at n-1 (with respect to the sw corner)

                  if (alphamx .ge. 0) then  ! particle coming from the West (u .ge. 0)
                     xdn2 = 1d0 - 2d0*alphamx / Deltax
                     xdn1 = 1d0 - alphamx / Deltax
                     if (alphamy .ge. 0) then ! particle coming from the South (v .ge. 0)
                        isw=i-1
                        jsw=j-1
                        inw=i-1
                        jnw=j
                        ine=i
                        jne=j
                        ise=i
                        jse=j-1
                        ydn2 = 1d0 - 2d0*alphamy / Deltax
                        ydn1 = 1d0 - alphamy / Deltax
                     else                     ! particle coming from the North (v .lt. 0)
                        isw=i-1
                        jsw=j
                        inw=i-1
                        jnw=j+1
                        ine=i
                        jne=j+1
                        ise=i
                        jse=j
                        ydn2 = -2d0*alphamy / Deltax
                        ydn1 = -1d0*alphamy / Deltax
                     endif
                  else                      ! particle coming from the East (u .lt. 0) 
                     xdn2 = -2d0*alphamx / Deltax
                     xdn1 = -1d0*alphamx / Deltax
                     if (alphamy .ge. 0) then ! particle coming from the South (v .ge. 0)                             
                        isw=i
                        jsw=j-1
                        inw=i
                        jnw=j
                        ine=i+1
                        jne=j
                        ise=i+1
                        jse=j-1
                        ydn2 = 1d0 - 2d0*alphamy / Deltax
                        ydn1 = 1d0 - alphamy / Deltax
                     else                     ! particle coming from the North (v .lt. 0) 
                        isw=i
                        jsw=j
                        inw=i
                        jnw=j+1
                        ine=i+1
                        jne=j+1
                        ise=i+1
                        jse=j
                        ydn2 = -2d0*alphamy / Deltax
                        ydn1 = -1d0*alphamy / Deltax
                     endif
                  endif

! find hbef using cubic interpolation 
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation

                  fsw=hn2(isw,jsw)
                  fnw=hn2(inw,jnw)
                  fne=hn2(ine,jne)
                  fse=hn2(ise,jse)
                  fxsw=fx(hn2(isw+1,jsw), hn2(isw-1,jsw), 2d0)
                  fxnw=fx(hn2(inw+1,jnw), hn2(inw-1,jnw), 2d0)
                  fxne=fx(hn2(ine+1,jne), hn2(ine-1,jne), 2d0)
                  fxse=fx(hn2(ise+1,jse), hn2(ise-1,jse), 2d0)
                  fysw=fy(hn2(isw,jsw+1), hn2(isw,jsw-1), 2d0)
                  fynw=fy(hn2(inw,jnw+1), hn2(inw,jnw-1), 2d0)
                  fyne=fy(hn2(ine,jne+1), hn2(ine,jne-1), 2d0)
                  fyse=fy(hn2(ise,jse+1), hn2(ise,jse-1), 2d0)
                  fxysw=fxy(hn2(isw+1,jsw+1),hn2(isw-1,jsw+1),hn2(isw+1,jsw-1),hn2(isw-1,jsw-1),4d0)
                  fxynw=fxy(hn2(inw+1,jnw+1),hn2(inw-1,jnw+1),hn2(inw+1,jnw-1),hn2(inw-1,jnw-1),4d0)
                  fxyne=fxy(hn2(ine+1,jne+1),hn2(ine-1,jne+1),hn2(ine+1,jne-1),hn2(ine-1,jne-1),4d0)
                  fxyse=fxy(hn2(ise+1,jse+1),hn2(ise-1,jse+1),hn2(ise+1,jse-1),hn2(ise-1,jse-1),4d0)
      
                  hbef=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                     fxsw, fxnw, fxne, fxse,  &
                                     fysw, fynw, fyne, fyse,  &
                                     fxysw,fxynw,fxyne,fxyse, &
                                     xdn2, ydn2 )

                  if (SLlimiter) then
                     upper=max(fsw, fnw, fne, fse)
                     lower=min(fsw, fnw, fne, fse)
                     hbef=apply_lim(hbef, upper, lower, xdn2, ydn2, fsw,  fnw,  fne,  fse)
                  endif

! find Abef using cubic interpolation                                                                    
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation 

                  fsw=An2(isw,jsw)
                  fnw=An2(inw,jnw)
                  fne=An2(ine,jne)
                  fse=An2(ise,jse)
                  fxsw=fx(An2(isw+1,jsw), An2(isw-1,jsw), 2d0)
                  fxnw=fx(An2(inw+1,jnw), An2(inw-1,jnw), 2d0)
                  fxne=fx(An2(ine+1,jne), An2(ine-1,jne), 2d0)
                  fxse=fx(An2(ise+1,jse), An2(ise-1,jse), 2d0)
                  fysw=fy(An2(isw,jsw+1), An2(isw,jsw-1), 2d0)
                  fynw=fy(An2(inw,jnw+1), An2(inw,jnw-1), 2d0)
                  fyne=fy(An2(ine,jne+1), An2(ine,jne-1), 2d0)
                  fyse=fy(An2(ise,jse+1), An2(ise,jse-1), 2d0)
                  fxysw=fxy(An2(isw+1,jsw+1),An2(isw-1,jsw+1),An2(isw+1,jsw-1),An2(isw-1,jsw-1),4d0)
                  fxynw=fxy(An2(inw+1,jnw+1),An2(inw-1,jnw+1),An2(inw+1,jnw-1),An2(inw-1,jnw-1),4d0)
                  fxyne=fxy(An2(ine+1,jne+1),An2(ine-1,jne+1),An2(ine+1,jne-1),An2(ine-1,jne-1),4d0)
                  fxyse=fxy(An2(ise+1,jse+1),An2(ise-1,jse+1),An2(ise+1,jse-1),An2(ise-1,jse-1),4d0)

                  Abef=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                     fxsw, fxnw, fxne, fxse,  &
                                     fysw, fynw, fyne, fyse,  &
                                     fxysw,fxynw,fxyne,fxyse, &
                                     xdn2, ydn2 )

                  if (SLlimiter) then
                     upper=max(fsw, fnw, fne, fse)
                     lower=min(fsw, fnw, fne, fse)
                     Abef=apply_lim(Abef, upper, lower, xdn2, ydn2, fsw,  fnw,  fne,  fse)
                  endif

!------------------------------------------------------------------------
! find right hand side terms rhsdynh and rhsdynA (time level n-1 = n1)
! minus sign added later in final calc of hout, Aout
!------------------------------------------------------------------------ 

! find rhsdynh using cubic interpolation                                                 
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation

                  fsw=hn1(isw,jsw)*div(isw,jsw)
                  fnw=hn1(inw,jnw)*div(inw,jnw)
                  fne=hn1(ine,jne)*div(ine,jne)
                  fse=hn1(ise,jse)*div(ise,jse)
                  fxsw=fx(fse, hn1(isw-1,jsw)*div(isw-1,jsw), 2d0)
                  fxnw=fx(fne, hn1(inw-1,jnw)*div(inw-1,jnw), 2d0)
                  fxne=fx(hn1(ine+1,jne)*div(ine+1,jne), fnw, 2d0)
                  fxse=fx(hn1(ise+1,jse)*div(ise+1,jse), fsw, 2d0)
                  fysw=fy(fnw, hn1(isw,jsw-1)*div(isw,jsw-1), 2d0)
                  fynw=fy(hn1(inw,jnw+1)*div(inw,jnw+1), fsw, 2d0)
                  fyne=fy(hn1(ine,jne+1)*div(ine,jne+1), fse, 2d0)
                  fyse=fy(fne, hn1(ise,jse-1)*div(ise,jse-1), 2d0)
                  fxysw=fxy(fne, hn1(inw-1,jnw)*div(inw-1,jnw), &
                        hn1(ise,jse-1)*div(ise,jse-1),hn1(isw-1,jsw-1)*div(isw-1,jsw-1),4d0)
                  fxynw=fxy(hn1(ine,jne+1)*div(ine,jne+1), hn1(inw-1,jnw+1)*div(inw-1,jnw+1), &
                        fse,hn1(isw-1,jsw)*div(isw-1,jsw), 4d0)
                  fxyne=fxy(hn1(ine+1,jne+1)*div(ine+1,jne+1),hn1(inw,jnw+1)*div(inw,jnw+1), &
                        hn1(ise+1,jse)*div(ise+1,jse), fsw, 4d0)
                  fxyse=fxy(hn1(ine+1,jne)*div(ine+1,jne), fnw, &
                        hn1(ise+1,jse-1)*div(ise+1,jse-1),hn1(isw,jsw-1)*div(isw,jsw-1),4d0) 

                  rhsdynh=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                        fxsw, fxnw, fxne, fxse,  &
                                        fysw, fynw, fyne, fyse,  &
                                        fxysw,fxynw,fxyne,fxyse, &
                                        xdn1, ydn1 )

! find rhsdynA using cubic interpolation                                                                                 
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation                           

                  fsw=An1(isw,jsw)*div(isw,jsw)
                  fnw=An1(inw,jnw)*div(inw,jnw)
                  fne=An1(ine,jne)*div(ine,jne)
                  fse=An1(ise,jse)*div(ise,jse)
                  fxsw=fx(fse, An1(isw-1,jsw)*div(isw-1,jsw), 2d0)
                  fxnw=fx(fne, An1(inw-1,jnw)*div(inw-1,jnw), 2d0)
                  fxne=fx(An1(ine+1,jne)*div(ine+1,jne), fnw, 2d0)
                  fxse=fx(An1(ise+1,jse)*div(ise+1,jse), fsw, 2d0)
                  fysw=fy(fnw, An1(isw,jsw-1)*div(isw,jsw-1), 2d0)
                  fynw=fy(An1(inw,jnw+1)*div(inw,jnw+1), fsw, 2d0)
                  fyne=fy(An1(ine,jne+1)*div(ine,jne+1), fse, 2d0)
                  fyse=fy(fne, An1(ise,jse-1)*div(ise,jse-1), 2d0)
                  fxysw=fxy(fne, An1(inw-1,jnw)*div(inw-1,jnw), &
                        An1(ise,jse-1)*div(ise,jse-1),An1(isw-1,jsw-1)*div(isw-1,jsw-1),4d0)
                  fxynw=fxy(An1(ine,jne+1)*div(ine,jne+1), An1(inw-1,jnw+1)*div(inw-1,jnw+1), &
                        fse,An1(isw-1,jsw)*div(isw-1,jsw), 4d0)
                  fxyne=fxy(An1(ine+1,jne+1)*div(ine+1,jne+1),An1(inw,jnw+1)*div(inw,jnw+1), &
                        An1(ise+1,jse)*div(ise+1,jse), fsw, 4d0)
                  fxyse=fxy(An1(ine+1,jne)*div(ine+1,jne), fnw, &
                        An1(ise+1,jse-1)*div(ise+1,jse-1),An1(isw,jsw-1)*div(isw,jsw-1),4d0)

                  rhsdynA=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                        fxsw, fxnw, fxne, fxse,  &
                                        fysw, fynw, fyne, fyse,  &
                                        fxysw,fxynw,fxyne,fxyse, &
                                        xdn1, ydn1 )

!                  rhsdynh = 0d0
!                  rhsdynA = 0d0

!------------------------------------------------------------------------
! find output h and A at time level n  
!------------------------------------------------------------------------

                  hout(i,j) = hbef - 2d0*Deltat*rhsdynh
                  Aout(i,j) = Abef - 2d0*Deltat*rhsdynA
       
                  hout(i,j) = max(hout(i,j), 0d0)
                  Aout(i,j) = max(Aout(i,j), 0d0)
                  Aout(i,j) = min(Aout(i,j), 1d0)

               elseif (caseSL == 2) then ! upwind for special cases  
         
                  fluxx=calc_fluxx(un1(i,j), un1(i+1,j), hn1(i-1,j), hn1(i,j), hn1(i+1,j))
                  fluxy=calc_fluxy(vn1(i,j), vn1(i,j+1), hn1(i,j-1), hn1(i,j), hn1(i,j+1))

                  hout(i,j) = hn1(i,j) - DtoverDx * ( fluxx + fluxy )
                  hout(i,j) = max(hout(i,j), 0d0)

                  fluxx=calc_fluxx(un1(i,j), un1(i+1,j), An1(i-1,j), An1(i,j), An1(i+1,j))
                  fluxy=calc_fluxy(vn1(i,j), vn1(i,j+1), An1(i,j-1), An1(i,j), An1(i,j+1))

                  Aout(i,j) = An1(i,j) - DtoverDx * ( fluxx + fluxy )
                  Aout(i,j) = max(Aout(i,j), 0d0) ! COULD BE MOVED BELOW FOR BOTH CASES SL
                  Aout(i,j) = min(Aout(i,j), 1d0)

               endif

            enddo
         enddo

      endif ! choice of method
      
      return
    end subroutine advection

    subroutine calc_dFx (utp, tracertp, dFx)

      implicit none

      include 'parameter.h'
      include 'CB_mask.h'

      integer i, j
      double precision, intent(in) :: utp(0:nx+2,0:ny+2),tracertp(0:nx+1,0:ny+1)
      double precision, intent(out):: dFx(nx,ny)
      double precision :: F1, F2
      
      do i = 1, nx
         do j = 1, ny
            
            if (maskC(i,j) .eq. 1) then

               if ( utp(i,j) .ge. 0d0 ) then
                  F1 = utp(i,j)*tracertp(i-1,j)
               else
                  F1 = utp(i,j)*tracertp(i,j)
               endif

               if ( utp(i+1,j) .ge. 0d0 ) then
                  F2 = utp(i+1,j)*tracertp(i,j)
               else
                  F2 = utp(i+1,j)*tracertp(i+1,j)
               endif

               dFx(i,j)=F2-F1

            endif

         enddo
      enddo
      
    end subroutine calc_dFx

    subroutine calc_dFy (vtp, tracertp, dFy)

      implicit none

      include 'parameter.h'
      include 'CB_mask.h'

      integer i, j
      double precision, intent(in) :: vtp(0:nx+2,0:ny+2),tracertp(0:nx+1,0:ny+1)
      double precision, intent(out):: dFy(nx,ny)
      double precision :: F1, F2

      do i = 1, nx
         do j = 1, ny

            if (maskC(i,j) .eq. 1) then

               if ( vtp(i,j) .ge. 0d0 ) then
                  F1 = vtp(i,j)*tracertp(i,j-1)
               else
                  F1 = vtp(i,j)*tracertp(i,j)
               endif

               if ( vtp(i,j+1) .ge. 0d0 ) then
                  F2 = vtp(i,j+1)*tracertp(i,j)
               else
                  F2 = vtp(i,j+1)*tracertp(i,j+1)
               endif

               dFy(i,j)=F2-F1

            endif

         enddo
      enddo

    end subroutine calc_dFy

    subroutine calc_div (utp, vtp, div)

      implicit none

      include 'parameter.h'
      include 'CB_mask.h'
      include 'CB_const.h'

      integer i, j
      double precision, intent(in) :: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision, intent(out):: div(nx,ny)

      do i = 1, nx
         do j = 1, ny

            if (maskC(i,j) .eq. 1) then

               div(i,j)=(utp(i+1,j)-utp(i,j) + vtp(i,j+1)-vtp(i,j))/Deltax 

            endif

         enddo
      enddo

    end subroutine calc_div

    function fx(fright, fleft, deno) result(dfdx)
      
      double precision, intent(in) :: fright, fleft                                 
      double precision, intent(in) :: deno
      double precision             :: dfdx ! output df/dx                                          

      dfdx = ( fright - fleft ) / deno
      
    end function fx

    function fy(fup, fdown, deno) result(dfdy)

      double precision, intent(in) :: fup, fdown 
      double precision, intent(in) :: deno
      double precision             :: dfdy ! output df/dy                     

      dfdy = ( fup - fdown ) / deno

    end function fy

    function fxy (fur, ful, fdr, fdl, deno) result(dfdxdy)

      double precision, intent(in) :: fur, ful, fdr, fdl
      double precision, intent(in) :: deno
      double precision             :: dfdxdy ! output d(dfdx)/dy                          
                                                        
    ! u:up, d:down, r:right, l:left
                                     
      dfdxdy = ( fur - ful - fdr + fdl ) / deno

    end function fxy

    function cubic_interp ( fsw,  fnw,  fne,  fse,   &
                            fxsw, fxnw, fxne, fxse,  &
                            fysw, fynw, fyne, fyse,  &
                            fxysw,fxynw,fxyne,fxyse, &
                            xd, yd ) result(finterp)

    !-------------------------------------------------------------
    ! cubic interpolation in a square of size 1x1.
    ! https://en.wikipedia.org/wiki/Bicubic_interpolation
    !
    !     (0,1)------(1,1)
    !       |          |
    !       |          |
    !       |          |
    !       |          |
    !     (0,0)------(1,0)
    !
    ! (0,0) = SouthWest (sw)
    ! (0,1) = NorthWest (nw)
    ! (1,1) = NorthEast (ne)
    ! (1,0) = SouthEast (se)
    !
    ! xd and yd are distances [0,1] from thea sw corner
    !-------------------------------------------------------------

      double precision, intent(in) :: fsw, fnw, fne, fse         ! input
      double precision, intent(in) :: fxsw, fxnw, fxne, fxse     ! input
      double precision, intent(in) :: fysw, fynw, fyne, fyse     ! input
      double precision, intent(in) :: fxysw, fxynw, fxyne, fxyse ! input  
      double precision, intent(in) :: xd, yd                     ! input
      double precision             :: finterp                    ! output    
      
      double precision a00, a10, a20, a30, a01, a11, a21, a31, a02, a12, a22, a32, a03, a13, a23, a33

      a00 = fsw
      a10 = fxsw
      a20 = -3d0*fsw + 3d0*fse - 2d0*fxsw - fxse
      a30 = 2d0*fsw - 2d0*fse + fxsw + fxse
      a01 = fysw
      a11 = fxysw
      a21 = -3d0*fysw + 3d0*fyse -2d0*fxysw - fxyse
      a31 = 2d0*fysw - 2d0*fyse + fxysw + fxyse
      a02 = -3d0*fsw + 3d0*fnw -2d0*fysw - fynw
      a12 = -3d0*fxsw + 3d0*fxnw -2d0*fxysw - fxynw
      a22 = 9d0*fsw - 9d0*fse - 9d0*fnw + 9d0*fne + 6d0*fxsw + 3d0*fxse - 6d0*fxnw - 3d0*fxne + &
            6d0*fysw - 6d0*fyse + 3d0*fynw - 3d0*fyne + 4d0*fxysw + 2d0*fxyse + 2d0*fxynw + fxyne
      a32 = -6d0*fsw + 6d0*fse + 6d0*fnw - 6d0*fne - 3d0*fxsw - 3d0*fxse + 3d0*fxnw + 3d0*fxne - &
            4d0*fysw + 4d0*fyse - 2d0*fynw + 2d0*fyne - 2d0*fxysw - 2d0*fxyse - fxynw - fxyne
      a03 = 2d0*fsw - 2d0*fnw + fysw + fynw
      a13 = 2d0*fxsw - 2d0*fxnw + fxysw + fxynw
      a23 = -6d0*fsw + 6d0*fse + 6d0*fnw - 6d0*fne - 4d0*fxsw - 2d0*fxse + 4d0*fxnw + 2d0*fxne - &
            3d0*fysw + 3d0*fyse - 3d0*fynw + 3d0*fyne - 2d0*fxysw - fxyse - 2d0*fxynw - fxyne
      a33 = 4d0*fsw - 4d0*fse - 4d0*fnw + 4d0*fne + 2d0*fxsw + 2d0*fxse - 2d0*fxnw - 2d0*fxne + &
            2d0*fysw - 2d0*fyse + 2d0*fynw - 2d0*fyne + fxysw + fxyse + fxynw + fxyne

      finterp = a00 + a01*yd + a02*(yd**2) + a03*(yd**3) + &
                a10*xd + a11*xd*yd + a12*xd*(yd**2) + a13*xd*(yd**3) + &
                a20*(xd**2) + a21*(xd**2)*yd + a22*(xd**2)*(yd**2) + a23*(xd**2)*(yd**3) + &
                a30*(xd**3) + a31*(xd**3)*yd + a32*(xd**3)*(yd**2) + a33*(xd**3)*(yd**3) 

    end function cubic_interp

      function calc_fluxx (uij, uip1j, Tim1j, Tij, Tip1j) result(fluxx)

        double precision, intent(in) :: uij, uip1j, Tim1j, Tij, Tip1j ! T=tracer
        double precision :: fluxx, F1, F2
        
        if ( uij .ge. 0d0 ) then
           F1 = uij*Tim1j
        else
           F1 = uij*Tij
        endif

        if ( uip1j .ge. 0d0 ) then
           F2 = uip1j*Tij
        else
           F2 = uip1j*Tip1j
        endif

        fluxx=F2-F1

      end function calc_fluxx

      function calc_fluxy (vij, vijp1, Tijm1, Tij, Tijp1) result(fluxy)

        double precision, intent(in) :: vij, vijp1, Tijm1, Tij, Tijp1 ! T=tracer
        double precision :: fluxy, F1, F2

        if ( vij .ge. 0d0 ) then
           F1 = vij*Tijm1
        else
           F1 = vij*Tij
        endif

        if ( vijp1 .ge. 0d0 ) then
           F2 = vijp1*Tij
        else
           F2 = vijp1*Tijp1
        endif

        fluxy=F2-F1
        
      end function calc_fluxy

      function apply_lim(var, upper, lower, xd, yd, fsw, fnw, fne, fse) result(var_lim)

        double precision, intent(in) :: var, upper, lower, xd, yd, fsw, fnw, fne, fse ! input
        double precision             :: var_lim ! output                                    
        double precision             :: fs, fn ! interp value on southern and northern boundaries
        
        var_lim=var

        if (var .gt. upper .or. var .lt. lower) then
           fs = fsw + ( fse - fsw ) * xd ! linear interp x direction...deno=dx=1
           fn = fnw + ( fne - fnw ) * xd ! linear interp x direction...deno=dx=1
           var_lim = fs + ( fn - fs ) * yd ! linear interp y direction...deno=dy=1  
        endif
        
      end function apply_lim
