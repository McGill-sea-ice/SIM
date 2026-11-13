
!************************************************************************
!     Subroutine load_data: load forcing data (wind stress, air temperature,
!       water temperature and ocean current.
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
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************

      subroutine ocn_Tclim

      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_mask.h'
      include 'CB_const.h'

      character(len=*), parameter :: dir = 'forcing/ocnT/'
      character(LEN=80) fname1
      character(len=2) :: cmonth
      character(len=5) :: cdelta
      character(len=4) :: cnx, cny
      integer i, j,  kmo, mo

      if ( Thermodyn ) then
         write(cnx, '(I0)') int(nx)
         write(cny, '(I0)') int(ny)
!------------------------------------------------------------------------
!     load ocean temperature created by Tocn_clim_gen.m
!------------------------------------------------------------------------

      do kmo = 0, 13

         mo = kmo

         if ( kmo .eq. 0 ) then
            mo = 12
         elseif ( kmo .eq. 13 ) then
            mo = 1
         endif

         fname1 = ''
         
         if ( OcnTemp .eq. 'MonthlyClim' .or.                     &
                   OcnTemp .eq. 'calculated' ) then
            write(cmonth, '(I2.2)') mo
            if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
                & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
                write(cdelta, '(I0)') int(Deltax/1d03)
                fname1 = dir // cdelta // '/Tocn' // cmonth
            elseif ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
                & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
                & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
                write(cdelta, '(F0.2)') Deltax/1d03
                fname1 = dir // 'Tocn' // cmonth // '_dx' // trim(cdelta) &
                    & // '_nx' // trim(cnx) // '_ny' // trim(cny)
            endif

            open(unit = 30, file = fname1, status = 'unknown')
         
            do j = 0 , ny+1
               read(30,*) ( To_clim(i,j,kmo), i = 0, nx+1 )
            enddo
         
            close(30)

         endif

            
         do j = 0, ny+1
            do i = 0, nx+1

               To_clim(i,j,kmo) = ( To_clim(i,j,kmo) + 273.15d0 ) &
                                 * maskC(i,j)

!     user specified ocean temperature  (non-dim)
!     WARNING specify To not To_clim

               if ( OcnTemp .eq. 'specified' ) then
                  To_clim(i,j,kmo) = ( -1.8d0 + 273.15d0 ) * maskC(i,j)
               endif
                  
            enddo
         enddo
         
      enddo

      endif

      return
      end
      

