!***********************************************************************
!     subroutine LagrangianTracer: 
!       Calculate the tracer position and various fields at the tracer
!       location (eg. h, A, ...) at each time step.
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             29-07-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************


      subroutine LagrangianTracer (date)
        USE datetime, ONLY : datetime_type, datetime_delta_type, delta_init, datetime_str_6
        USE datetime, ONLY: OPERATOR(+), OPERATOR(-)
      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_const.h'
      include 'CB_buoys.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer nn, dx, ii, dy, jj

      double precision dudx, ubuoy, dvdy, vbuoy ,u1, u2, v1, v2

      character(LEN = 6) datestr

      TYPE(datetime_type) :: date 
      TYPE(datetime_delta_type) :: date_step


      date_step = delta_init(seconds = int(Deltat))
      date = date + date_step          !Increment the date to store position at next time step
      datestr = datetime_str_6(date)
      
      do nn = 1, nbuoys
         if (BuoyActive(nn) .eq. 1) then
            
      
!------------------------------------------------------------------------
!     x, y-coord and grid cell in which buoy 'nn' is positionned
!------------------------------------------------------------------------

         dx = mod (BuoyPosition(nn,1), Deltax)
         ii = int (BuoyPosition(nn,1) / Deltax) + 1

         dy = mod (BuoyPosition(nn,2), Deltax)
         jj = int (BuoyPosition(nn,2) / Deltax) + 1

!------------------------------------------------------------------------
!     Release buoys if date = sinkdate
!------------------------------------------------------------------------

      if ( ( Buoys == 'Track' ) .and. ( datestr == sinkdate(nn) ) ) &
           BuoyActive(nn) = 0

!------------------------------------------------------------------------
!     Release buoys if outside the domain
!------------------------------------------------------------------------

      if ( (ii < 1) .or. (jj < 1) .or. (ii > nx) .or. (jj > ny) ) &
           BuoyActive(nn) = 0

!------------------------------------------------------------------------
!     Release the buoys if A(ii,jj) < 0.5 
!------------------------------------------------------------------------

      if (A(ii,jj) .le. 0.5) BuoyActive(nn) = 0

!------------------------------------------------------------------------
!     Divide the grid cell in 4 quadrants
!     Interpolate ui at the buoy position accordingly 
!     Bilinear interpolation : uses u at 4 points 
!     No slip at closed boundaries
!     Newmann condition at open boundaries (du/dn = 0)
!------------------------------------------------------------------------


         dudx = ( uice(ii+1,jj) - uice(ii,jj) ) / Deltax
         u1 = uice(ii,jj) + dudx * dx 


         if ( dy - Deltax/2d0 .ge. 0d0 ) then !QUADRANTS 1 & 4

           
            if ( jj < ny ) then
               
               u2 = uice(ii,jj+1) & 
                    + ( uice(ii+1,jj+1) - uice(ii,jj+1) ) / Deltax & 
                    * dx
               
               if ( maskC(ii,jj+1) == 1 ) then 
                  ubuoy = u1 + ( ( u2 - u1 ) / Deltax ) &
                       * ( dy - Deltax/2d0 )
               else             ! closed boundary : u2 = 0
                  ubuoy = u1 + ( ( u2 - u1 ) / ( Deltax / 2d0 ) ) &
                       * ( dy - Deltax / 2d0 )
               endif   
               
            else                ! open boundary (u2 = u1, ubuoy = u1)
               
               ubuoy = u1              
               
            endif


         elseif ( dy - Deltax/2d0 .lt. 0d0 ) then  !QUADRANTS 2 & 3                               
            
            if ( jj > 1 ) then
               
               u2 = uice(ii,jj-1) & 
                    + ( uice(ii+1,jj-1) - uice(ii,jj-1) ) / Deltax &
                    * dx

 
               if ( maskC(ii,jj-1) == 1 ) then
                  ubuoy = u1 + ( ( u1 - u2 ) / Deltax ) &
                       * ( dy - Deltax/2d0 )
               else             ! closed boundary
                  ubuoy = u1 + ( ( u1 - u2 ) / ( Deltax / 2d0 ) ) &
                  * ( dy - Deltax / 2d0 )
               endif
               
            else                ! open boundary

               ubuoy = u1

            endif
            
         endif
         

!------------------------------------------------------------------------
!     Divide the grid cell in 4 quadrants
!     Interpolate vi at the buoy position accordingly
!     Bilinear interpolation : uses v at 4 points
!     No slip at closed boundaries
!     Newmann condition at open boundaries
!------------------------------------------------------------------------


         dvdy = ( vice(ii,jj+1) - vice(ii,jj) ) / Deltax
         v1 = vice(ii,jj) + dvdy * dy
         
         if ( dx - Deltax/2d0 .ge. 0d0 ) then !QUADRANT 1 & 2                                
            
            
            if ( ii < nx ) then
               
               v2 = vice(ii+1,jj) &
                    + ( vice(ii+1,jj+1) - vice(ii+1,jj) ) / Deltax * dy
               
               if ( maskC(ii+1,jj) == 1 ) then                  
                  vbuoy = v1 + ( ( v2 - v1) / Deltax ) & 
                       * ( dx - Deltax/2d0 )
               else             ! closed boundary ( v2 = 0 )
                  vbuoy = v1 + ( ( v2 - v1) / ( Deltax/2d0 ) ) &
                       * ( dx - Deltax/2d0 )                        
               endif   
               
            else                ! open boundary (v2 = v1)
               
               vbuoy = v1
               
            endif
            
            
         elseif ( dx - Deltax/2d0 .lt. 0d0 ) then !QUADRANTS 3 & 4
            
            
            if ( ii > 1 ) then
               
               v2 = vice(ii-1,jj) &
                    + ( vice(ii-1,jj+1) - vice(ii-1,jj) ) / Deltax * dy
               
               if ( maskC(ii-1,jj) == 1 ) then 
                  vbuoy = v1 + ( ( v1 - v2) / Deltax ) & 
                       * ( dx - Deltax/2d0 )                
               else             ! closed boundary

                  vbuoy = v1 + ( ( v1 - v2) / ( Deltax/2d0 ) ) &
                       * ( dx - Deltax/2d0 )
               endif
               
            else                ! open boundary (v2 = v1)
               
               vbuoy = v1 
               
            endif   
            
            
         endif
         
       
!------------------------------------------------------------------------               

         BuoyPosition(nn,1) = BuoyPosition(nn,1) + Deltat * ubuoy 
         BuoyPosition(nn,2) = BuoyPosition(nn,2) + Deltat * vbuoy
         BuoyGridCell(nn,1) = int (BuoyPosition(nn,1) / Deltax) + 1
         BuoyGridCell(nn,2) = int (BuoyPosition(nn,2) / Deltax) + 1

         if (date%hour == 0) then
            Buoydate(nn) = datestr
         endif
         
        
!------------------------------------------------------------------------
!     Additional tracers : h, A, sigma1, sigma2
!     Stores the average h, a, sigma1, sigma2 of the buoy's grid cell
!     These quantities are not interpolated at the buoy location
!------------------------------------------------------------------------

         BuoyTracers(nn,1) = h(ii,jj)
         BuoyTracers(nn,2) = A(ii,jj)
!         BuoyTracers(nn,3) = sigmaI(ii,jj)
!         BuoyTracers(nn,4) = sigmaII(ii,jj)
        
         endif
                                         
      enddo

      date = date - date_step          !Decrement the date to thepresent date


      return
    end subroutine LagrangianTracer


