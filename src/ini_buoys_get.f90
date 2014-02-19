!*************************************************************************
!     subroutine ini_buoys_get:
!     Reads the initial conditions (buoy number, sink date and time and 
!     position on the mask (x,y)) of the buoys and seeds them in the model
!     In par_get110.f, choose between:
!     - 'Track'
!       Reads the initial conditions from a single file and seeds the buoys 
!       at their real initial position. The buoy velocity is computed every
!       at the model buoy position until the buoy reaches the real sink date
!       (or other conditions described in LagrangianTracer).
!     - 'Daily'
!       Reads the initial conditions from a yearly file and seeds the buoys
!       every day at their corresponding real location. The buoy velocity is
!       interpolated every day at the real buoy location. 
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             11-11-08               L.-B. Tremblay
!
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************



      subroutine ini_buoys_get(now_date)
        USE datetime, ONLY : datetime_type, datetime_str_6

      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_buoys.h'
      include 'CB_Dyndim.h'


      TYPE(datetime_type) :: now_date
      character(LEN=6) date
      character(LEN=60) fname, path
      

      integer year, month, day
      integer nn

      year = now_date%year
      month= now_date%month
      day  = now_date%day

      date = datetime_str_6(now_date)
      

!--------------------------------------------------------------------

      if ( Buoys .eq. 'Track') then

!     Open buoys initial condition and post buoy files on the first day of the run (firstday == 1)

         if ( date .eq. startdate ) then

            open (100, file = 'buoys/Buoys.txt', status = 'old')
            print *, 'opening/readingBuoys.txt'
         
            
            open (97, file = 'output/BuoyTracks.txt', status = 'unknown') 
            print *, 'opening/postingBuoyTracks.txt'
            
            nbuoys = 0
         endif   

         seeddate(nbuoys) = date

         do while ( seeddate(nbuoys) == date )

            nbuoys = nbuoys + 1
            read (100,'(a6, X, i2, X, i6, X, f10.3, & 
                 & X, f10.3, X, a6, X, i2)', END = 20) &
                 seeddate(nbuoys), seedt(nbuoys), & 
                 buoynb(nbuoys), BuoyPosition(nbuoys,1), &
                 BuoyPosition(nbuoys,2), sinkdate(nbuoys), &
                 sinkt(nbuoys)
           
            BuoyActive(nbuoys)     = 1         !Seeding the buoy

            BuoyPosition(nbuoys,1) = BuoyPosition(nbuoys,1)* 1000.0 ![m]
            BuoyPosition(nbuoys,2) = BuoyPosition(nbuoys,2)* 1000.0 ![m]

 20         continue


            if ( seeddate(nbuoys) .ne. date ) then 
            !Reading one line too far: reinitializing all parameters to 0
               BACKSPACE(100)         
               BuoyActive(nbuoys)     = 0
               seedt(nbuoys)          = 0
               buoynb(nbuoys)         = 0
               BuoyPosition(nbuoys,1) = 0d0
               BuoyPosition(nbuoys,2) = 0d0
               sinkdate(nbuoys)       = ''
               sinkt(nbuoys)          = 0
            endif   

         enddo

         nbuoys = nbuoys - 1

!     Close buoys initial conditions file on the last day of the run 
             
         if ( date .eq. enddate ) then
            close(100)
            close(97)
         endif


      endif
      
!---------------------------------------------------------------------


      if ( Buoys .eq. 'Daily') then

      
!     Resetting buoy arrays to 0 every day: avoids keeping buoys info from
!     previous day if nbuoys on previous day > nbuoys on present day
         
         if ( date .eq. startdate ) then
            nbuoys = 1
         endif
     
         do nn = 1, nbuoys         
           
            BuoyActive(nn)     = 0
            BuoyVelocity(nn,1) = 0d0
            BuoyVelocity(nn,2) = 0d0
            BuoyPosition(nn,1) = 0d0
            BuoyPosition(nn,2) = 0d0
            BuoyGridCell(nn,1) = 0
            BuoyGridCell(nn,2) = 0
            BuoyTracers(nn,1)  = 0d0
            BuoyTracers(nn,2)  = 0d0
            BuoyTracers(nn,3)  = 0d0
            BuoyTracers(nn,4)  = 0d0
            buoynb             = 0
            seeddate(0)        = date
            sinkdate(nn)       = ''
            seedt(nn)          = 0
            sinkt(nn)          = 0
                        
         enddo
      
!     Open buoys initial conditions and buoy post files on the first day of the run (startdate) 
!     or on the first day of the year (01/01)      

         if ( (( day == 1) .and. (month == 1)) .or. &
              (date .eq. startdate) ) then  
         
            
            fname = ''
            path   = 'buoys/Buoys_'                
            write  ( fname, '(a12, i4, a4)' ) path, year, '.txt'
            open  ( 10, file = fname, status = 'old' )
            print *, 'opening/reading', fname
            
            open (96, file = 'output/BuoyDailyDrift.txt', &
                 status = 'unknown')
            print *, 'opening/postingBuoyDailyDrift.txt'
                     
         endif
      
!! PRINT A LITTLE WARNING FOR THE USER TO TELL THEM ABOUT THE START OF THE RUN IF LATER
!! THAN THE FIRST BUOY DATE
      

         nbuoys = 0
      
         do while ( seeddate(nbuoys) == date )
         
            nbuoys = nbuoys + 1
            read (10,'(a6, X, i2, X, i6, X, f10.3, X, f10.3)', END=10) & 
                 seeddate(nbuoys), seedt(nbuoys), buoynb(nbuoys), &
                 BuoyPosition(nbuoys,1), BuoyPosition(nbuoys,2)
            
            BuoyActive(nbuoys)     = 1    !Seeding the buoys

            BuoyPosition(nbuoys,1) = BuoyPosition(nbuoys,1) * 1000.0 ![m]
            BuoyPosition(nbuoys,2) = BuoyPosition(nbuoys,2) * 1000.0 ![m]
 

10         continue

            if ( seeddate(nbuoys) .ne. date ) then

               BACKSPACE(10)
               BuoyActive(nbuoys)     = 0
               seedt(nbuoys)          = 0
               buoynb(nbuoys)         = 0
               BuoyPosition(nbuoys,1) = 0d0
               BuoyPosition(nbuoys,2) = 0d0


            endif
         enddo


!Closing buoys initial conditions files on the last day of the run 
!or on the last day of a year (31/12)
 
         if ( date .eq. enddate ) then
            close(10)
            close(96)
         else if ( (day == 31) .and. (month == 12) ) then
            close(10)
         endif


      endif


      return
      end
