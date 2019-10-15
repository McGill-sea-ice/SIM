!========================================================================
!     Common block DynVariables: dynamic variables (defined on the C-grid)
!========================================================================


      double precision                       &
                h        (0:nx+2,0:ny+2),    &
                A        (0:nx+2,0:ny+2),    &
		hn1      (0:nx+2,0:ny+2),    &
		An1      (0:nx+2,0:ny+2),    &
                hn2      (0:nx+2,0:ny+2),    &
                An2      (0:nx+2,0:ny+2),    & 
                hn1      (0:nx+2,0:ny+2),    &
                An1      (0:nx+2,0:ny+2),    &
                uice     (0:nx+2,0:ny+2), &
                vice     (0:nx+2,0:ny+2), &
                un1      (0:nx+2,0:ny+2), &
                vn1      (0:nx+2,0:ny+2), &
                un2      (0:nx+2,0:ny+2), &
                vn2      (0:nx+2,0:ny+2), &
                Pp       (0:nx+2,0:ny+2), &
                Pt       (0:nx+2,0:ny+2), & 
                P        (0:nx+2,0:ny+2), &
                Pf       (0:nx+2,0:ny+2), &
                etaC     (0:nx+2,0:ny+2), &
                etaB     (0:nx+2,0:ny+2), &
                zetaC    (0:nx+2,0:ny+2), &
                etaCf    (0:nx+2,0:ny+2), &
                etaBf    (0:nx+2,0:ny+2), &
                zetaCf   (0:nx+2,0:ny+2), &
                GammaMEB (0:nx+2,0:ny+2), &
                GammaMEB_B (0:nx+2,0:ny+2), &
                dam      (0:nx+2,0:ny+2), &
                damB     (0:nx+2,0:ny+2), &
                dam1     (0:nx+2,0:ny+2), &
                damB1    (0:nx+2,0:ny+2), &
                dfactor  (0:nx+2,0:ny+2), &
                dfactorB (0:nx+2,0:ny+2), &
                CoheC    (0:nx+2,0:ny+2), &
                sigcC    (0:nx+2,0:ny+2), &
                sigtC    (0:nx+2,0:ny+2), &
                sigxy    (0:nx+2,0:ny+2), &
                sigxyB   (0:nx+2,0:ny+2), &
                sigxx    (0:nx+2,0:ny+2), &
                sigyy    (0:nx+2,0:ny+2)


      common/DynVariables/      &
                h,              & ! ice thickness
                A,              & ! ice concentration
                hn1,            & ! previous time step ice thickness
                An1,            & ! previous time step ice concentration
                hn2,            & ! ice thickness at time level n-2
                An2,            & ! ice concentration at time level n-2
                uice,           & ! x-comp ice velocity 
                vice,           & ! y-comp ice velocity
                un1,            & ! previous time level solution 
                vn1,            & ! previous time level solution
                un2,            & ! previous time level solution
                vn2,            & ! previous time level solution
                Pp,             & ! ice compressive strength
                Pt,             & ! ice tensile strength
                P,              & ! replacement pressure
                Pf,             & ! replacement pressure
                etaC,           & ! coefficient of shear viscosity (C-grid)
                etaB,           & ! coefficient of shear viscosity (B-grid)
                zetaC,          & ! coefficient of bulk viscosity (C-grid)
                etaCf,          & ! coefficient of shear viscosity (C-grid)
                etaBf,          & ! coefficient of shear viscosity (B-grid)
                zetaCf,         & ! coefficient of bulk viscosity (C-grid)
                GammaMEB,       & ! visco-elastic dissipation coefficient (C-grid)
                GammaMEB_B,     & ! visco-elastic dissipation coefficient (B-grid)
                dam,            & ! damage parameter (C-grid)
                damB,           & ! damage parameter (B-grid)
                dam1,           & ! previous time level damage parameter (C-grid)
                damB1,          & ! previous time level damage parameter (B-grid)
                dfactor,        & ! damage factor (C-grid)
                dfactorB,       & ! damage factor (B-grid)
                CoheC,          & ! Cohesion at grid center
                sigCC,          & ! Compressive strength at grid center
                sigTC,          & ! Tensile strength at grid center
                sigxy,          & ! shear stress at center
                sigxyB,         & ! shear stress at node
                sigxx,          & ! normal stress at center
                sigyy             ! normal stress at center





