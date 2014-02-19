!========================================================================
!     Common block mask: masks for the node (B-grid) and grid center (C-grid)
!========================================================================


      integer                                &
                maskB    ( 0:nx+2, 0:ny+2 ), &
                maskC    ( 0:nx+1, 0:ny+1 )

      common/mask/          &
                maskB,      & ! defined at the node   ( 0: land, 1: ocean )
                maskC         ! defined at the center ( 0: land, 1: ocean )






