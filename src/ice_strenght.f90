!***********************************************************************
!     subroutine Pressure: 
!
!     calculates the ice strenght using:
!       
!       P = Pstar * h * exp(-C (1 - A))
!
!************************************************************************


      subroutine Ice_strength ()
        use ellipse
        implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_mask.h'
      include 'CB_mask_boat.h'
      include 'CB_const_stressBC.h'
      include 'CB_options.h'

      integer i, j

      if ( Rheology .eq. 1 ) then

         do i = 1, nx ! tracer (C) point
            do j = 1, ny
               if (maskC(i,j) .eq. 1) then
                  if (stressBC .and. mboat(i,j) .eq. 1) then
                     Pp(i,j) = bpfactor * Pstar / 2d0 ! h and A are considered 1
                  else
                     Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
                     Pp(i,j) = Pp(i,j) / 2d0
                  endif
               endif
            enddo
         enddo

! watchout not sure the code works if stressBC=false...

         ! e and PpB are same everywhere except on sides of boat
         do j = 1, ny+1 ! node (B) point
            do i = 1, nx+1
               if (stressBC .and. msideboat(i,j) .eq. 1) then
                  PpB(i,j) = hlevel * Pstar / 2d0 ! h and A are considered 1  
                  ell_2B(i,j) = 1d0/(ecboat**2d0)
               else
                  PpB(i,j) = ( Pp(i-1,j) + Pp(i,j) + Pp(i,j-1) + Pp(i-1,j-1) ) / 4d0 !in visco
                  ell_2B(i,j) = ell_2
               endif
            enddo
         enddo

!------- set P = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------

       if (.not. stressBC) then
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
       endif

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
      




