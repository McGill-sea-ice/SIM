!========================================================================
!     Common block buoys: Buoys related variables
!========================================================================

      integer nbuoys  



      double precision                           &
               BuoyPosition  (nbuoy, 2),         &
               BuoyVelocity  (nbuoy, 2),         &
               BuoyTracers   (nbuoy, 4) 
      
      character(LEN=6)              &  
               seeddate (0:nbuoy),  &
               sinkdate (nbuoy),    &
               Buoydate (nbuoy)

      integer                            &
               BuoyGridCell (nbuoy, 2),  &
               BuoyActive (nbuoy),       &
               seedt (nbuoy),            &
               sinkt (nbuoy),            &
               buoynb (nbuoy)  

      common/buoys/             &
               BuoyPosition,    & ! buoy position [x1,y1; x2,y2 ...]
               BuoyGridCell,    & ! grid cell (i,j) in which the buoy is positioned
               BuoyVelocity,    & ! buoy velocity [u1,v1; u2,v2 ...]
               BuoyTracers,     & ! h, A, sigmaI, sigmaII in grid cell of buoy
               BuoyActive,      & ! state of buoy (seeded or not)
               Buoydate,        & ! time history of date of buoy activity
               nbuoys,          & ! number of buoys
               seedt,           & ! seeding time (00Z or 12Z)
               sinkt,           & ! sinking time (00Z or 12Z)
               buoynb,          & ! buoy argos nb 
               seeddate,        & ! seeding date
               sinkdate           ! sinking date
