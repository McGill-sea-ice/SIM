!========================================================================
!     Common block Dim: parameters in the mom equations
!========================================================================

      double precision                    &
                Ktracer(10), Cdw, Cda

      double precision                    &
                CdwC1    (0:nx+2,0:ny+2), &
                CdwC2    (0:nx+2,0:ny+2), &
                Cbasal1   (0:nx+2,0:ny+2), &
                Cbasal2   (0:nx+2,0:ny+2), &
                CdwC1f   (0:nx+2,0:ny+2), &
                CdwC2f   (0:nx+2,0:ny+2)


      common/Dim/               &
                Cdw,            &
                Cda,            &
                CdwC1,          & ! water drag coefficient (C-grid)
                CdwC2,          & ! water drag coefficient (C-grid)
                Cbasal1,        &
                Cbasal2,        &
                CdwC1f,         & ! water drag coefficient (C-grid)
                CdwC2f,         & ! water drag coefficient (C-grid)
                Ktracer           ! diffusion coefficient for the tracers


