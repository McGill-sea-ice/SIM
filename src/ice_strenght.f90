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
      include 'CB_options.h'


      integer i, j


      if ( Rheology .eq. 1 ) then

         do i = 1, nx
            do j = 1, ny
               if (maskC(i,j) .eq. 1) then
                  Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ) )
                  Pp(i,j) = Pp(i,j) / 2d0
                  ! if (i .eq. 10) then
                  !    print *, 'PPPPPPPPPp(i,j)=  ', Pp(i,j), j
                  ! endif
               endif
            enddo
         enddo

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

      elseif ( Rheology .eq. 2 ) then !FB: working here
         !print *, 'Calculating Rheo=2'
         do i = 1, nx
            do j = 1, ny
              if (maskC(i,j) .eq. 1) then
               Pp(i,j) = Pstar * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
               if (i .eq. 10) then
                !print *, 'PPPPPPPPPp(i,j)=  ', Pp(i,j)
               endif !print
                  !Pmin(i,j) = Cohe  * h(i,j) * dexp(-C * ( 1d0 - A(i,j) ))
!               p(i,j) = 0.0d0
              endif

            enddo
         enddo

!------- set P = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------


!         do i = 0, nx+1
!
!            if (maskC(i,0) .eq. 1) then             
!               Pp(i,1)  = 0d0
!            endif
!
!            if (maskC(i,ny+1) .eq. 1) then             
!               Pp(i,ny)  = 0d0
!            endif
! 
!         enddo
!
!         do j = 0, ny+1
!
!            if (maskC(0,j) .eq. 1) then             
!               Pp(1,j)  = 0d0
!            endif
!
!            if (maskC(nx+1,j) .eq. 1) then             
!               Pp(nx,j)   = 0d0 
!            endif           
!
!         enddo

      endif !Rheology

      return
    end subroutine Ice_strength
      




