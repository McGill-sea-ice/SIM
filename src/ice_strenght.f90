!***********************************************************************
!     subroutine Pressure: 
!
!     calculates the ice strenght using:
!       
!       P = Pstar * h * exp(-C (1 - A)) * (1 - dam)
!
!     using damage (added by Antoine Savard in 08/2021)
!
!************************************************************************


      subroutine Ice_strength ()
        use ellipse
        implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_mask.h'
      include 'CB_options.h'


      integer i, j


      if ( Rheology .eq. 1 ) then

         if ( Damage .eq. 0 ) then
         
            do i = 1, nx
               do j = 1, ny
                  if (maskC(i,j) .eq. 1) then
                     Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
                     Pp(i,j) = Pp(i,j) / 2d0
                  endif
               enddo
            enddo

         elseif ( Damage .eq. 1 ) then    ! if you want damage or not the model simply change the parameter

            do i = 1, nx
               do j = 1, ny
                  if (maskC(i,j) .eq. 1) then
                     Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) ) * ( 1d0 - dam(i,j) )
                     Pp(i,j) = Pp(i,j) / 2d0
                  endif
               enddo
            enddo
         endif

!------- set P = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------


         do i = 0, nx+1

            if (maskC(i,0) .eq. 1) then             
               Pp(i,1)  = 0d0
            endif

            if (maskC(i,ny+1) .eq. 1) then             
               Pp(i,ny)  = 0d0
            endif
 
         enddo

         do j = 0, ny+1

            if (maskC(0,j) .eq. 1) then             
               Pp(1,j)  = 0d0
            endif

            if (maskC(nx+1,j) .eq. 1) then             
               Pp(nx,j)   = 0d0 
            endif           

         enddo

      elseif ( Rheology .eq. 2 ) then

!         do i = 1, nx
!            do j = 1, ny
!              if (maskC(i,j) .eq. 1) then
!               Pmax(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
!               Pmin(i,j) = Cohe  * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
!               p(i,j) = 0.0d0
!              endif
!            enddo
!         enddo

      endif

      return
      end subroutine Ice_strength
      




