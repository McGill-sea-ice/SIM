!*************************************************************************
!     subroutine atmosphere:
!       Calculates the air temperature 
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             14-05-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!
!     JFL the dimentionalized code does not calculate the atm T (not yet)
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************


      subroutine atmosphere (date)
        USE datetime, ONLY: datetime_type, datetime_str, datetime_str_6
      implicit none

      include 'parameter.h'
      include 'CB_Thermodim.h'
      include 'CB_ThermoForcing.h'
      include 'CB_DynForcing.h'
      include 'CB_ThermoVariables.h'
      include 'CB_DynVariables.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      TYPE(datetime_type), INTENT(in) :: date

      integer day, mon, i, j

      double precision esat
      external esat


!------------------------------------------------------------------------
!     Compute the coefficient of the linear total heat flux between surface
!     and atmosphere (positive: heat flow into atmosphere);  A Tatm - B.
!     land  : 6 meters thick conductive layer with 6C base temperature
!     ocean : Tocn at t-dt is used to compute ocean/atmos heat fluxes
!     ice   : "h" meter thick conductive layer with Tfreeze base temp. 
!------------------------------------------------------------------------

      day = date%day
      mon = date%month

      call shortwave (date)

      do i = 1, nx
         do j = 1, ny
            Pvap(i,j) = esat( Ta(i,j), A(i,j) )
         enddo
      enddo



      return
      end



      double precision function esat(T, A)

      implicit none

      include 'parameter.h'

      double precision A, T, A1, A2, expo

      A1   = 9.50d0 * A +  7.50d0 * ( 1d0 - A )
      A2   = 7.66d0 * A + 35.86d0 * ( 1d0 - A )

      expo = A1 * ( T - 273.15d0 ) &
                / ( T - A2    )
      
      esat = 611d0 * 10d0 ** expo

      return
      end



! Bruno, please comment on what this does. 

      subroutine restore(a, x, icase, error, ierr, jerr)

      include 'parameter.h'

      integer i, j, k

      dimension a(0:nx+1, 0:ny+1), x(ntot)
      
      k     = 0

      if (icase .eq. 1) then

         do j = 1, ny, 2
            do i = 1, nx
               k = k + 1
               if (abs(a(i,j)-x(k)) .gt. error) then
                  error = abs(a(i,j)-x(k))
                  ierr  = i
                  jerr  = j
               endif
               a(i,j) = x(k)
            enddo
         enddo

      elseif (icase .eq. 2) then
            
         do j = 2, ny, 2
            do i = 1, nx
               k = k + 1
               if (abs(a(i,j)-x(k)) .gt. error) then
                  error = abs(a(i,j)-x(k))
                  ierr  = i
                  jerr  = j
               endif
               a(i,j) = x(k)
            enddo
         enddo

      elseif (icase .eq. 3) then

         do i = 1, nx, 2
            do j = 1, ny
               k = k + 1
               if (abs(a(i,j)-x(k)) .gt. error) then
                  error = abs(a(i,j)-x(k))
                  ierr  = i
                  jerr  = j
               endif
               a(i,j) = x(k)
            enddo
         enddo

      elseif (icase .eq. 4) then

         do i = 2, nx, 2
            do j = 1, ny
               k = k + 1
               if (abs(a(i,j)-x(k)) .gt. error) then
                  error = abs(a(i,j)-x(k))
                  ierr  = i
                  jerr  = j
               endif
               a(i,j) = x(k)
            enddo
         enddo

      endif


      return
      end



