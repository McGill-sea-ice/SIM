!========================================================================
!     Common block ThermoVariables: Thermodynamic variables
!     (defined on the C-grid) 
!========================================================================


      double precision                    &
                Ta       (0:nx+1,0:ny+1), &
                Ti       (0:nx+1,0:ny+1), &
                Tl       (0:nx+1,0:ny+1), &
                To       (0:nx+1,0:ny+1), &
                Qoa      (0:nx+1,0:ny+1), &
                Qoa_f    (0:nx+1,0:ny+1), &
                Qia      (0:nx+1,0:ny+1), &
                Qsh_io   (0:nx+1,0:ny+1), &
                Qadvdiff (0:nx+1,0:ny+1), &
                Pvap     (0:nx+1,0:ny+1)



      common/ThermoVariables/   &
                Ta,             & ! air temperature
                Ti,             & ! ice surface temperature
                Tl,             & ! land surface temperature
                To,             & ! ocean mixed layer temperature
                Qoa,            & ! ocean-atmosphere heat flux
                Qoa_f,          & ! ocean-atmosphere heat flux (ice growth)
                Qia,            & ! conductive heat flux through ice
                Qsh_io,         & ! sensible heat flux (ocn/ice)
                Qadvdiff,       & ! Advection-Diffusion heat transfer
                Pvap           ! atmospheric vapour pressure





