!------------------------------------------------------------------------------
! transform from matrix (the domain) to vector (u comp are stacked first and
! then the v comp). Performs also the reverse operation
!------------------------------------------------------------------------------

      subroutine transformer (xu,xv,x,kase)

      implicit none

      include 'parameter.h'

      integer i, j, kase
      double precision xu(0:nx+2,0:ny+2),xv(0:nx+2,0:ny+2),x(nvar)

!     kase = 1  : matrix --> vector
!     kase = 0  : vector --> matrix

      if (kase .eq. 1) then

         do j = 1, ny
            do i = 1, nx+1

               x(i+(j-1)*(nx+1)) = xu(i,j)

            enddo
         enddo

         do j = 1, ny+1
            do i = 1, nx

               x(i+(j-1)*nx+(nx+1)*ny) = xv(i,j)

            enddo
         enddo

      elseif (kase .eq. 0) then

         do j = 1, ny
            do i = 1, nx+1

               xu(i,j) =  x(i+(j-1)*(nx+1))

            enddo
         enddo

         do j = 1, ny+1
            do i = 1, nx

               xv(i,j) = x(i+(j-1)*nx+(nx+1)*ny)

            enddo
         enddo

      endif

      return
      end
