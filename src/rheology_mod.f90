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
