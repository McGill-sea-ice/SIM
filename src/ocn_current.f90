 
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


      subroutine ocn_current

      implicit none

      include 'parameter.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_mask.h'
      include 'CB_const.h'

      character(len=2) :: cdelta
      character(LEN=60) fname1, fname2

      integer startyear, endyear
      integer i, j


      read ( startdate(5:6), '(i2)' ) startyear
      read ( enddate  (5:6), '(i2)' ) endyear
      write(cdelta, '(I2)') int(Deltax)/1000

!------------------------------------------------------------------------
!     load ocean current data (time independant)
!------------------------------------------------------------------------

      fname1 = ''
      fname2 = ''

      if ( Current .eq. 'YearlyMean' ) then
         fname1= 'forcing/ocncurrent/'// cdelta // '/uwater' // cdelta // '.clim'
         fname2= 'forcing/ocncurrent/'// cdelta // '/vwater' // cdelta // '.clim'
      endif
      
      if ( fname1 .ne. '') then
          print *, 'Reading ocean currents for the ' // cdelta // ' km resolution grid.' 
          print *, 'reading ', fname1
          print *, 'reading ', fname2


         open(unit = 30, file = fname1, status = 'unknown')
         open(unit = 31, file = fname2, status = 'unknown')

         do j = 1 , ny
            read(30,*) ( uwatnd(i,j), i = 1, nx+1 ) !nd = non-divergent
         enddo
         
         do j = 1, ny+1
            read(31,*) ( vwatnd(i,j), i = 1, nx )
         enddo
         
      ! group those two do-loop ??

         do j = 1, ny+1
            do i = 1, nx+1
               if ( j.eq.1 .or. j.eq.ny+1 ) then
                  uwater(i,j) = 0d0
               else
                  uwater(i,j) = ( uwatnd(i,j-1) + uwatnd(i,j) ) / 2d0
               endif
            enddo
         enddo
         
         do j = 1, ny+1
            do i = 1, nx+1
               if ( i.eq.1 .or. i.eq.nx+1 ) then
                  vwater(i,j) = 0d0
               else
                  vwater(i,j) = ( vwatnd(i-1,j) + vwatnd(i,j) ) / 2d0
               endif
            enddo
         enddo

         close(30)
         close(31)
         close(32)

      endif

      if (Current .eq. 'specified' ) then
         stop
	 print *, 'Reading ocean currents from specified values'
         do i = 1, nx+1
            do j = 1, ny+1
                       
!     user specified ocean current

               uwater(i,j)  = 0.0d0  * maskB(i,j)
              vwater(i,j)  = 0.0d0  * maskB(i,j)
               uwatnd(i,j)  = 0.0d0  * min( maskB(i,j)   &
                                     + maskB(i,j+1), 1 )
               vwatnd(i,j)  = 0.0d0  * min( maskB(i,j) + &
                                       maskB(i+1,j), 1 )
            enddo
         enddo
      endif

      return
    end subroutine ocn_current
      

