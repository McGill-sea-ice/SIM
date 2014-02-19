
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

      character(LEN=60) fname1
      character(len=2) :: cdelta, cmonth
      integer i, j,  kmo, mo

      if ( Thermodyn ) then
         write(cdelta, '(I2)') int(Deltax)/1000
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
            fname1  = 'forcing/ocnT/' // cdelta // '/Tocn' // cmonth

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
      

