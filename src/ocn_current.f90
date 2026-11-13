 
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

      character(len=5) :: cdelta
      character(len=4) :: cnx, cny
      character(len=*), parameter :: dir = 'forcing/ocncurrent/'
      character(LEN=80) fname1, fname2

      integer startyear, endyear
      integer i, j


      read ( startdate(5:6), '(i2)' ) startyear
      read ( enddate  (5:6), '(i2)' ) endyear

!------------------------------------------------------------------------
!     load ocean current data (time independant)
!------------------------------------------------------------------------

      fname1 = ''
      fname2 = ''
    
      if ( Current .eq. 'clim' ) then
         write(cnx, '(I0)') int(nx)
         write(cny, '(I0)') int(ny)
         if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
             & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
             write(cdelta, '(I0)') int(Deltax/1d03)
             fname1 = dir // cdelta // 'uwater' // cdelta // '.clim'
             fname2 = dir // cdelta // 'vwater' // cdelta // '.clim'
         elseif ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
             & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
             & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
             write(cdelta, '(F0.2)') Deltax/1d03
             fname1 = dir // 'uwater' // '_dx' // trim(cdelta) // '_nx' &
                    & // trim(cnx) // '_ny' // trim(cny)// '.clim'
             fname2 = dir // 'vwater' // '_dx' // trim(cdelta) // '_nx' &
                    & // trim(cnx) // '_ny' // trim(cny)// '.clim'
         endif
      endif
      
      print *, 'Reading ocean currents for the ' // cdelta // ' km resolution grid.' 
      print *, 'reading ', fname1
      print *, 'reading ', fname2

      if ( fname1 .ne. '') then

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

      do i = 1, nx+1
         do j = 1, ny+1
                       
!     user specified ocean current

            if (Current .eq. 'specified' ) then
               uwater(i,j)  = 0.0d0  * maskB(i,j)
               vwater(i,j)  = 0.0d0  * maskB(i,j)
               uwatnd(i,j)  = 0.0d0  * min( maskB(i,j)   &
                                     + maskB(i,j+1), 1 )
               vwatnd(i,j)  = 0.0d0  * min( maskB(i,j) + &
                                       maskB(i+1,j), 1 )
            endif

         enddo
      enddo

      return
    end subroutine ocn_current
      

