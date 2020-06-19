!========================================================================
!     Common block mask: masks for the boat 
!========================================================================


      integer                                   &
                msideboat   ( 0:nx+2, 0:ny+2 ), & 
                mboat       ( 0:nx+1, 0:ny+1 )    

      common/bmask/         &
                msideboat,  & ! defined at the node   ( 0: ice, 1: boat )
                mboat         ! defined at the center ( 0: ice, 1: boat )






