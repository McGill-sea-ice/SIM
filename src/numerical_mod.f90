MODULE numerical_VP

  IMPLICIT NONE
  INTEGER :: kjac, ksor, klsor
  INTEGER :: NLmax ! max nb of Newton loop for JFNK
  INTEGER :: OLmax ! max nb of Outer loop for Picard
  INTEGER :: klinesearch ! linesearch is applied for JFNK for k .ge. klinesearch
  DOUBLE PRECISION :: wjac, wsor, wlsor, gamma_nl, dropini, res_t

END MODULE Numerical_VP

MODULE numerical_EVP

  IMPLICIT NONE
  INTEGER :: Nsub
  DOUBLE PRECISION :: Eo ! T = Eo*Deltat 
  character(LEN=20) :: init_stress ! VP or zero (stress for tstep=1, s=1)

END MODULE numerical_EVP

MODULE solver_choice

  IMPLICIT NONE
  INTEGER :: solver

END MODULE solver_choice
