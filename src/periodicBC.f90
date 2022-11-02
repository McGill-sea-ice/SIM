
!************************************************************************
!
!       Subroutine periodicBC: apply periodic boundary conditions on a pair
!			of variables u and v. They come out with BC applied
!
!           0   = n
!           n+1 = 1
!           n+2 = 2
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             2015                   Mathilde Jutras
!     V2.0            30-10-2016             Mathieu Plante
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca

!************************************************************************
!************************************************************************

!
      subroutine periodicBC(u,v)

        implicit none

        include 'parameter.h'
        include 'CB_options.h'

        double precision u(0:nx+2,0:ny+2), v(0:nx+2,0:ny+2)
        integer i, j

        if (Periodic_y .eq. 1) then
        
          do i = 1, nx+1 ! periodic in y
               v(i,0)    = v(i,ny)
               v(i,ny+1) = v(i,1)
               v(i,ny+2) = v(i,2)
               u(i,0)    = u(i,ny)
               u(i,ny+1) = u(i,1)
               u(i,ny+2) = u(i,2)
          enddo
        endif
         
        if (Periodic_x .eq. 1) then        
          do j = 1, ny+1 ! periodic in x
               v(0,j)    = v(nx,j)
               v(nx+1,j) = v(1,j)
               v(nx+2,j) = v(2,j)
               u(0,j)    = u(nx,j)
               u(nx+1,j) = u(1,j)
               u(nx+2,j) = u(2,j)
          enddo
        endif

        return

      end subroutine periodicBC
