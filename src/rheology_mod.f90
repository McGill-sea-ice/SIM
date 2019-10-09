MODULE ellipse
!
! C : ice strength parameter
! Pstar : Ice compression strength parameter
! ell2  : ellipticity**2 
! ell_2 : 1/ellipticity**2
  IMPLICIT NONE
  DOUBLE PRECISION :: C, Pstar, ell2, ell_2, Tens

END MODULE ellipse



MODULE triangle
!
! Cohe : Ice tensile strength
! phi : internal angle of friction
! delta : angle of dilatancy
! etamax : maximum shear viscosity
  IMPLICIT NONE
  DOUBLE PRECISION :: Cohe, phi, delta, etamax

END MODULE triangle



MODULE elastic

! phi : internal angle of friction
! sigC : tensile strength cut-off
! Poisson : Poisson ratio of sea ice
! Young : Young's ration of sea ice 
! lambda0: viscous relaxation time scale

  IMPLICIT NONE

  DOUBLE PRECISION :: Young, Poisson, lambda0
  DOUBLE PRECISION :: alpha, Theal 
  DOUBLE PRECISION :: Cohe, sigt, sigc, phi, Tdam

END MODULE elastic
