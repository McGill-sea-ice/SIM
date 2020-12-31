!***********************************************************************
!     subroutine advection (upwind scheme or upwind-RungeKutta2):
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
!************************************************************************

      subroutine advection ( upts, vpts, utp, vtp, hin, Ain, hout, Aout )

      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j, k
      integer isw, jsw, inw, jnw, ine, jne, ise, jse !SL SouthWest=sw, nw, ne, se corners
      double precision, intent(in)    :: upts(0:nx+2,0:ny+2), vpts(0:nx+2,0:ny+2)
      double precision, intent(in)    :: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision, intent(inout) :: hin(0:nx+1,0:ny+1), Ain(0:nx+1,0:ny+1)
      double precision, intent(out)   :: hout(0:nx+1,0:ny+1), Aout(0:nx+1,0:ny+1)
      double precision                :: ustar(0:nx+2,0:ny+2), vstar(0:nx+2,0:ny+2)
      double precision                :: hstar(0:nx+1,0:ny+1), Astar(0:nx+1,0:ny+1)
      double precision                :: dFx(nx,ny), dFy(nx,ny)
      double precision                :: alphamx, alphamy, xd, yd
      double precision                :: fsw, fnw, fne, fse
      double precision                :: fxsw, fxnw, fxne, fxse, fysw, fynw, fyne, fyse
      double precision                :: fxysw, fxynw, fxyne, fxyse

!------------------------------------------------------------------------ 
!     set dhin/dx, dAin/dx = 0 at the outside cell when there is an open bc 
!------------------------------------------------------------------------ 

      do i = 0, nx+1
               
         if (maskC(i,0) .eq. 1) then
            
            hin(i,0) = ( 4d0 * hin(i,1) - hin(i,2) )/3d0
            hin(i,0) = max(hin(i,0), 0d0)
            Ain(i,0) = ( 4d0 * Ain(i,1) - Ain(i,2) )/3d0
            Ain(i,0) = max(Ain(i,0), 0d0)
            Ain(i,0) = min(Ain(i,0), 1d0)
                  
         endif

         if (maskC(i,ny+1) .eq. 1) then
            
            hin(i,ny+1)= ( 4d0 * hin(i,ny) - hin(i,ny-1) ) / 3d0
            hin(i,ny+1)= max(hin(i,ny+1), 0d0)
            Ain(i,ny+1)= ( 4d0 * Ain(i,ny) - Ain(i,ny-1) ) / 3d0
            Ain(i,ny+1)= max(Ain(i,ny+1), 0d0)
            Ain(i,ny+1)= min(Ain(i,ny+1), 1d0)

         endif
 
      enddo

      do j = 0, ny+1

         if (maskC(0,j) .eq. 1) then
                  
            hin(0,j)  = ( 4d0 * hin(1,j) - hin(2,j) ) / 3d0
            hin(0,j)  = max(hin(0,j), 0d0)
            Ain(0,j)  = ( 4d0 * Ain(1,j) - Ain(2,j) ) / 3d0
            Ain(0,j)  = max(Ain(0,j), 0d0)
            Ain(0,j)  = min(Ain(0,j), 1d0)

         endif

         if (maskC(nx+1,j) .eq. 1) then
                 
            hin(nx+1,j) = ( 4d0 * hin(nx,j) - hin(nx-1,j) ) / 3d0
            hin(nx+1,j) = max(hin(nx+1,j), 0d0)
            Ain(nx+1,j) = ( 4d0 * Ain(nx,j) - Ain(nx-1,j) ) / 3d0
            Ain(nx+1,j) = max(Ain(nx+1,j), 0d0)
            Ain(nx+1,j) = min(Ain(nx+1,j), 1d0)

         endif

      enddo

      if ( adv_scheme .eq. 'upwind' ) then

!------------------------------------------------------------------------
!     compute the difference of the flux for thickness 
!------------------------------------------------------------------------

         call calc_dFx (utp, hin, dFx)
         call calc_dFy (vtp, hin, dFy)

!------------------------------------------------------------------------
!     update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hout(i,j) = hin(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  hout(i,j) = max(hout(i,j), 0d0)

               endif
                     
            enddo
         enddo

!------------------------------------------------------------------------  
!     compute the difference of the flux for concentration                                 
!------------------------------------------------------------------------                  

         call calc_dFx (utp, Ain, dFx)
         call calc_dFy (vtp, Ain, dFy)

!------------------------------------------------------------------------ 
!     update the concentration values      
!     (in a separate do-loop to conserve mass)                                             
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  Aout(i,j) = Ain(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
                  Aout(i,j) = max(Aout(i,j), 0d0)
                  Aout(i,j) = min(Aout(i,j), 1d0)

               endif

            enddo
         enddo

      elseif ( adv_scheme .eq. 'upwindRK2' ) then
               
!------------------------------------------------------------------------
!     predictor: compute the difference of the flux for thickness 
!------------------------------------------------------------------------

         call calc_dFx (upts, hin, dFx)
         call calc_dFy (vpts, hin, dFy)

!------------------------------------------------------------------------
!     predictor: update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hstar(i,j) = hin(i,j) - (DtoverDx / 2d0) * ( dFx(i,j) + dFy(i,j) )
                  hstar(i,j) = max(hstar(i,j), 0d0)

               endif
                     
            enddo
         enddo

!------------------------------------------------------------------------  
!     predictor: compute the difference of the flux for concentration   
!------------------------------------------------------------------------                  

         call calc_dFx (upts, Ain, dFx)
         call calc_dFy (vpts, Ain, dFy)

!------------------------------------------------------------------------ 
!     predictor: update the concentration values      
!     (in a separate do-loop to conserve mass)                                             
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  Astar(i,j) = Ain(i,j) - (DtoverDx / 2d0) * ( dFx(i,j) + dFy(i,j) )
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
  
         ustar = ( upts + utp ) / 2d0
         vstar = ( vpts + vtp ) / 2d0

         call calc_dFx (ustar, hstar, dFx)
         call calc_dFy (vstar, hstar, dFy)

!------------------------------------------------------------------------
!     corrector: update the thickness values
!     (in a separate do-loop to conserve mass)
!------------------------------------------------------------------------
            
         do i = 1, nx
            do j = 1, ny

               if (maskC(i,j) .eq. 1) then

                  hout(i,j) = hin(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
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

                  Aout(i,j) = Ain(i,j) - DtoverDx * ( dFx(i,j) + dFy(i,j) )
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
  
!------------------------------------------------------------------------  
! find velocity at x-alphamx, y-alphamy  and t=n1
!------------------------------------------------------------------------
        alphamx=0.01d0 ! initial value
        alphamy=0.01d0 ! initial value

        do k = 1, 5
           um = (un1(i+1)+un1(i))/2d0 - (un1(i+1)-un1(i))*alpham/Deltax
           if (order .gt. 1 .and. i .gt. 1 .and. i .lt. nx) then ! O2...O3 not coded yet  
              um=um + (alpham**2d0)*(un1(i+2)-un1(i+1)-un1(i)+un1(i-1))/(4d0*Deltax2)
           endif
           alpham=Deltat*um
        enddo
!------------------------------------------------------------------------
! find hbef and Abef (initial position of particle at time level n-2=n2)
!------------------------------------------------------------------------

! 1) identify coordinates of 4 corners
! xd and yd are distances in the interval [0,1]. They are calc from the sw corner 

         if (alphamx .ge. 0) then  ! particle coming from the West (u .ge. 0)
            xd = 1d0 - 2d0*alphamx / Deltax
            if (alphamy .ge. 0) then ! particle coming from the South (v .ge. 0)
               isw=i-1
               jsw=j-1
               inw=i-1
               jnw=j
               ine=i
               jne=j
               ise=i
               jse=j-1
               yd = 1d0 - 2d0*alphamy / Deltax
            else                     ! particle coming from the North (v .lt. 0)
               isw=i-1
               jsw=j
               inw=i-1
               jnw=j+1
               ine=i
               jne=j+1
               ise=i
               jse=j
               yd = -2d0*alphamy / Deltax
            endif
         else                      ! particle coming from the East (u .lt. 0) 
            xd = -2d0*alphamx / Deltax
            if (alphamy .ge. 0) then ! particle coming from the South (v .ge. 0)                             
               isw=i
               jsw=j-1
               inw=i
               jnw=j
               ine=i+1
               jne=j
               ise=i+1
               jse=j-1
               yd = 1d0 - 2d0*alphamy / Deltax
            else                     ! particle coming from the North (v .lt. 0) 
               isw=i
               jsw=j
               inw=i
               jnw=j+1
               ine=i+1
               jne=j+1
               ise=i+1
               jse=j
               yd = -2d0*alphamy / Deltax
            endif
         endif         

! 2a) find hbef using cubic interpolation 

! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation

         fsw=hn2in(isw,jsw)
         fnw=hn2in(inw,jnw)
         fne=hn2in(ine,jne)
         fse=hn2in(ise,jse)
         fxsw=fx(hn2in(isw+1,jsw), hn2in(isw-1,jsw))
         fxnw=fx(hn2in(inw+1,jnw), hn2in(inw-1,jnw))
         fxne=fx(hn2in(ine+1,jne), hn2in(ine-1,jne))
         fxse=fx(hn2in(ise+1,jse), hn2in(ise-1,jse))
         fysw=fy(hn2in(isw,jsw+1), hn2in(isw,jsw-1))
         fynw=fy(hn2in(inw,jnw+1), hn2in(inw,jnw-1))
         fyne=fy(hn2in(ine,jne+1), hn2in(ine,jne-1))
         fyse=fy(hn2in(ise,jse+1), hn2in(ise,jse-1))
         fxysw=fxy(hn2in(isw+1,jsw+1), hn2in(isw-1,jsw+1), hn2in(isw+1,jsw-1), hn2in(isw-1,jsw-1))
         fxynw=fxy(hn2in(inw+1,jnw+1), hn2in(inw-1,jnw+1), hn2in(inw+1,jnw-1), hn2in(inw-1,jnw-1))
         fxyne=fxy(hn2in(ine+1,jne+1), hn2in(ine-1,jne+1), hn2in(ine+1,jne-1), hn2in(ine-1,jne-1))
         fxyse=fxy(hn2in(ise+1,jse+1), hn2in(ise-1,jse+1), hn2in(ise+1,jse-1), hn2in(ise-1,jse-1))
      
         hbef=cubic_interp( fsw,  fnw,  fne,  fse,   &
                            fxsw, fxnw, fxne, fxse,  &
                            fysw, fynw, fyne, fyse,  &
                            fxysw,fxynw,fxyne,fxyse, &
                            xd, yd )

! 2b) find Abef using cubic interpolation                                                                    

! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation 

         fsw=An2in(isw,jsw)
         fnw=An2in(inw,jnw)
         fne=An2in(ine,jne)
         fse=An2in(ise,jse)
         fxsw=fx(An2in(isw+1,jsw), An2in(isw-1,jsw))
         fxnw=fx(An2in(inw+1,jnw), An2in(inw-1,jnw))
         fxne=fx(An2in(ine+1,jne), An2in(ine-1,jne))
         fxse=fx(An2in(ise+1,jse), An2in(ise-1,jse))
         fysw=fy(An2in(isw,jsw+1), An2in(isw,jsw-1))
         fynw=fy(An2in(inw,jnw+1), An2in(inw,jnw-1))
         fyne=fy(An2in(ine,jne+1), An2in(ine,jne-1))
         fyse=fy(An2in(ise,jse+1), An2in(ise,jse-1))
         fxysw=fxy(An2in(isw+1,jsw+1), An2in(isw-1,jsw+1), An2in(isw+1,jsw-1), An2in(isw-1,jsw-1))
         fxynw=fxy(An2in(inw+1,jnw+1), An2in(inw-1,jnw+1), An2in(inw+1,jnw-1), An2in(inw-1,jnw-1))
         fxyne=fxy(An2in(ine+1,jne+1), An2in(ine-1,jne+1), An2in(ine+1,jne-1), An2in(ine-1,jne-1))
         fxyse=fxy(An2in(ise+1,jse+1), An2in(ise-1,jse+1), An2in(ise+1,jse-1), An2in(ise-1,jse-1))

         Abef=cubic_interp( fsw,  fnw,  fne,  fse,   &
                            fxsw, fxnw, fxne, fxse,  &
                            fysw, fynw, fyne, fyse,  &
                            fxysw,fxynw,fxyne,fxyse, &
                            xd, yd )

      endif
      
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

    function fx(fip1, fim1) result(fx)
      
      double precision, intent(in) :: fip1, fim1 ! ip1=i+1, im1=i-1                                
      double precision             :: fx         ! output df/dx                                          

      fx = ( fip1 - fim1 ) / 2d0
      
    end function fx

    function fy(fjp1, fjm1) result(fy)

      double precision, intent(in) :: fjp1, fjm1 ! jp1=j+1, jm1=j-1 
      double precision             :: fy         ! output df/dx                     

      fy = ( fjp1 - fjm1 ) / 2d0

    end function fx

    function fxy (fip1jp1, fim1jp1, fip1jm1, fim1jm1) result(fxy)

      double precision, intent(in) :: fip1jp1, fim1jp1, fip1jm1, fim1jm1 ! input 
      double precision             :: fxy    ! output                           
                                                                                             
      fxy = ( fip1jp1 - fim1jp1 - fip1jm1 + fim1jm1 ) / 4d0

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
      a22 = 9d0*fsw - 9d0*fse - 9d0*fnw + 9d0*fne + 6d0*fxsw + 3d0*fxse - 6d0*fxnw - 3d0fxne + &
            6d0*fysw - 6d0*fyse + 3d0*fynw - 3d0*fyne + 4d0*fxysw + 2d0*fxyse + 2d0*fxynw + fxyne
      a32 = -6d0*fsw + 6d0*fse + 6d0*fnw - 6d0*fne - 3d0*fxsw - 3d0*fxse + 3d0*fxnw + 3d0fxne - &
            4d0*fysw + 4d0*fyse - 2d0*fynw + 2d0*fyne - 2d0*fxysw - 2d0*fxyse - fxynw - fxyne
      a03 = 2d0*fsw - 2d0*fnw + fysw + fxnw
      a13 = 2d0*fxsw - 2d0*fxnw + fxysw + fxynw
      a23 = -6d0*fsw + 6d0*fse + 6d0*fnw - 6d0*fne - 4d0*fxsw - 2d0*fxse + 4d0*fxnw + 2d0fxne - &
            3d0*fysw + 3d0*fyse - 3d0*fynw + 3d0*fyne - 2d0*fxysw - fxyse - 2d0*fxynw - fxyne
      a33 = 4d0*fsw - 4d0*fse - 4d0*fnw + 4d0*fne + 2d0*fxsw + 2d0*fxse - 2d0*fxnw - 2d0fxne + &
            2d0*fysw - 2d0*fyse + 2d0*fynw - 2d0*fyne + fxysw + fxyse + fxynw + fxyne

      finterp = a00 + a01*yd + a02*(yd**2) + a03*(yd**3) + &
                a10*xd + a11*xd*yd + a12*xd*(yd**2) + a13*xd*(yd**3) + &
                a20*(xd**2) + a21*(xd**2)*yd + a22*(xd**2)*(yd**2) + a23*(xd**2)*(yd**3) + &
                a30*(xd**3) + a31*(xd**3)*yd + a32*(xd**3)*(yd**2) + a33*(xd**3)*(yd**3) + &

    end function cubic_interp
