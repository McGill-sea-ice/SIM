      subroutine update_tracer

      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'

      integer i, j, k
      
!------------------------------------------------------------------------
!     Store tracer 1 and 2 into h, A
!------------------------------------------------------------------------

         do k = 1, 2
            do i = 0, nx+1
               do j = 0, ny+1
               
!                  if ( k .eq. 1 ) h(i,j) = tracer(i,j,k) 
!                  if ( k .eq. 2 ) A(i,j) = tracer(i,j,k)
                  
               enddo
            enddo
         enddo
         
      return
      end
      

