!========================================================================
!     Common block stressBC: define stress at the boundaries 
!                            sigma = normal stress  (sigma11 or sigma22)
!                            tau   = shear stress (sigma12) 
!========================================================================


      double precision                    &
                sigmaS         (1:nx)   , &
                sigmaN         (1:nx)   , &
                sigmaW         (1:ny)   , &
                sigmaE         (1:ny)   , &
                tauS           (1:nx+1) , &
                tauN           (1:nx+1) , &
                tauW           (1:ny+1) , &
                tauE           (1:ny+1)



      common/stressBC/     &
                sigmaS,    &  ! BC sigma on South side
                sigmaN,    &  ! BC sigma on North side
                sigmaW,    &  ! BC sigma on West side
                sigmaE,    &  ! BC sigma on East side
                tauS,      &  ! BC tau on South side
                tauN,      &  ! BC tau on North side
                tauW,      &  ! BC tau on West side
                tauE          ! BC tau on East side






