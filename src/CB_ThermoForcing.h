!========================================================================
!     Common block ThermoForcing: Thermodynamic Forcing (air temp, ocn temp 
!       and shortwave radiation (defined on the C-grid).
!========================================================================


      double precision                             &
                To_clim     (0:nx+1,0:ny+1,0:13), &
                Qsw         (0:nx+1,0:ny+1)


      common/ThermoForcing/     &
                To_clim,        & ! monthly clim ocean ML temperature
                Qsw               ! shortave radiation from the sun





