!========================================================================
!     Common block ForcingFields: forcing field variables
!     (defined on the B-grid)
!========================================================================

      double precision                     &
                uair     (0:nx+2,0:ny+2),  &
                vair     (0:nx+2,0:ny+2),  &
                speeda   (0:nx+2,0:ny+2),  &
                uwater   (0:nx+2,0:ny+2),  &
                vwater   (0:nx+2,0:ny+2),  &
                uwatnd   (0:nx+2,0:ny+2),  &
                vwatnd   (0:nx+2,0:ny+2),  &
                speediw  (0:nx+2,0:ny+2)

      double precision                     &
                R1       (0:nx+2,0:ny+2),  &
                R2       (0:nx+2,0:ny+2),  &
                R1n      (0:nx+1,0:ny+1),  &
                R2n      (0:nx+1,0:ny+1),  &
                bu_ind   (0:nx+2,0:ny+2),  &
                bv_ind   (0:nx+2,0:ny+2),  &
                bu     (0:nx+2,0:ny+2),  &
                bv     (0:nx+2,0:ny+2)

      common/ForcingFields/     &
                uair,           & ! x-comp wind velocity 
                vair,           & ! y-comp wind velocity
                speeda,         & ! wind speed
                uwater,         & ! x-comp ocean current velocity
                vwater,         & ! y-comp ocean current velocity
                uwatnd,         & ! x-comp ocean current velocity (non-div)
                vwatnd,         & ! y-comp ocean current velocity (non-div)
                speediw           ! ice speed relative to ocean current
     
      common/ForcingFields/     &
                R1,             & ! Forcing on the C-grid
                R2,             & ! Forcing on the C-grid
                R1n,            & ! Forcing on the B-grid 
                R2n,            & ! Forcing on the B_grid
                bu_ind,         &
                bv_ind,         &
                bu,      &
                bv






