!*************************************************************************
!     subroutine post_buoys:
!       Outputs either - the daily position and velocity of the buoys
!                        along their trajectory ('Track')
!                      - the daily drift vectors interpolated at the 
!                        real buoy position ('Daily')
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


      subroutine post_buoys

      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_buoys.h'

      integer  nn

             

      if ( Buoys .eq. 'Daily') then
         
                  
         do nn = 1, nbuoys
                                    
            if (BuoyActive(nn) .eq. 1) then
               write(96,'(i6, X, a6, X, f12.3, X, f12.3, X,  &
                    i2, X, i2, X, f10.3, X, f10.3, X,  & 
                    f10.6, X, f10.6)') &
                    buoynb(nn), &
                    Buoydate(nn), &     
                    BuoyPosition(nn,1)/1000.0, & ![km]      
                    BuoyPosition(nn,2)/1000.0, & ![km] 
                    BuoyGridCell(nn,1), &  
                    BuoyGridCell(nn,2), &
                    BuoyTracers(nn,1), &
                    BuoyTracers(nn,2), &
                    BuoyTracers(nn,3), &
                    BuoyTracers(nn,4)
            endif
            
         enddo
         
      endif


      
      if ( Buoys .eq. 'Track') then
         
         
         do nn = 1, nbuoys

            if (BuoyActive(nn) .eq. 1) then
               write(97,'(i6, X, a6, X, f12.3, X, f12.3, X,  &
                    i2, X, i2, X, f10.3, X, f10.3, X,  &
                    f10.6, X, f10.6)') &
                    buoynb(nn), &
                    Buoydate(nn), &
                    BuoyPosition(nn,1)/1000.0, & ![km] 
                    BuoyPosition(nn,2)/1000.0, & ![km]
                    BuoyGridCell(nn,1), &
                    BuoyGridCell(nn,2), &
                    BuoyTracers(nn,1), &
                    BuoyTracers(nn,2), &
                    BuoyTracers(nn,3), &
                    BuoyTracers(nn,4)
            endif
            
         enddo
         
         
      endif
      

      return
    end subroutine post_buoys

