!========================================================================
!     Common block Dim: parameters in the mom equations
!========================================================================

      double precision                    &
                Cdw, Cda

      double precision                    &
                CdwC1    (0:nx+2,0:ny+2), &
                CdwC2    (0:nx+2,0:ny+2), &
                Cbasal1  (0:nx+2,0:ny+2), &
                Cbasal2  (0:nx+2,0:ny+2), &
                CdwC1f   (0:nx+2,0:ny+2), &
                CdwC2f   (0:nx+2,0:ny+2)

      common/Dim/               &
                Cdw,            &
                Cda,            &
                CdwC1,          & ! water drag coefficient (C-grid)
                CdwC2,          & ! water drag coefficient (C-grid)
                Cbasal1,        & ! basal stress coefficient (u)
                Cbasal2,        & ! basal stress coefficient (v)
                CdwC1f,         & ! water drag coefficient (C-grid)
                CdwC2f            ! water drag coefficient (C-grid)

