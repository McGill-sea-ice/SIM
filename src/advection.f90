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

      integer i, j
      double precision, intent(in)    :: upts(0:nx+2,0:ny+2), vpts(0:nx+2,0:ny+2)
      double precision, intent(in)    :: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision, intent(inout) :: hin(0:nx+1,0:ny+1), Ain(0:nx+1,0:ny+1)
      double precision, intent(out)   :: hout(0:nx+1,0:ny+1), Aout(0:nx+1,0:ny+1)
      double precision                :: ustar(0:nx+2,0:ny+2), vstar(0:nx+2,0:ny+2)
      double precision                :: hstar(0:nx+1,0:ny+1), Astar(0:nx+1,0:ny+1)
      double precision dFx(nx,ny), dFy(nx,ny)

      
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
