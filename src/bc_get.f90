!***********************************************************************
!     subroutine bc_get: sets the atm T and ocn T at boundary conditions
!                        sets ocn T over all the domain if MonthlyClim
!***********************************************************************

      subroutine bc_get (date)
        USE datetime, ONLY: datetime_type, datetime_str, datetime_str_6
      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_options.h'
      
      TYPE(datetime_type), INTENT(in) :: date
      integer day, mon, i, j, k1, k2

      double precision w1

      character(LEN=6) datestr


      datestr = datetime_str_6(date)
      day = date%day
      mon = date%month

      k1 = mon + day / 16 - 1
      k2 = k1 + 1

      w1 = mod (((day + 15) * 1d0) / 31d0, 1d0)

!------------------------------------------------------------------------
!     Boundary conditions on ocean if OcnTemp .eq. 'calculated'
!------------------------------------------------------------------------

      if ( OcnTemp .eq. 'calculated' ) then

         j = 0                  ! bottom edge of domain
      
         do i = 0, nx+1
            To(i,j) = (1d0 - w1) * To_clim(i,j,k1) + w1 *To_clim(i,j,k2)
         enddo

         j = ny+1               ! top edge of domain
      
         do i = 0, nx+1
            To(i,j) = (1d0 - w1) * To_clim(i,j,k1) + w1 *To_clim(i,j,k2)
         enddo

         i = 0                  ! left edge of domain
      
         do j = 0, ny+1
            To(i,j) = (1d0 - w1) * To_clim(i,j,k1) + w1 *To_clim(i,j,k2)
         enddo

         i = nx+1               ! right edge of domain

         do j = 0, ny+1
            To(i,j) = (1d0 - w1) * To_clim(i,j,k1) + w1 *To_clim(i,j,k2)
         enddo
         
!------------------------------------------------------------------------
!     Load Ocean Temp Monthly Climatology if OcnTemp .eq. 'MonthlyClim'
!------------------------------------------------------------------------

      elseif ( OcnTemp .eq. 'MonthlyClim' ) then
      
         do j = 0, ny+1
            do i = 0, nx+1
               To(i,j) = (1d0 - w1) * To_clim(i,j,k1) + &
                                 w1 * To_clim(i,j,k2)
            enddo
         enddo
   
      endif

      return
      end


