!=======================================================================
!     Common block DamageVariables: Damage variables
!     (defined on the C-grid) 
!=======================================================================


      double precision                    &
                Sdam       (0:nx+1,0:ny+1)



      common/ThermoVariables/   &
                Sdam                ! thermo source term (dam, continuity equation)