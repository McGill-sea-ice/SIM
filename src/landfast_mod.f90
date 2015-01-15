!MODULE bathymetry

!  IMPLICIT NONE
!  include 'parameter.h'

!  DOUBLE PRECISION :: bathy(0:nx+1,0:ny+1)

!END MODULE bathymetry

MODULE basal_param

  IMPLICIT NONE
  LOGICAL :: BasalStress 
  DOUBLE PRECISION :: CC, umin, k1, k2, crit

END MODULE basal_param



