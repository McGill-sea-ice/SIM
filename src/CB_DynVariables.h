!========================================================================
!     Common block DynVariables: dynamic variables (defined on the C-grid)
!========================================================================


      double precision                       &
                h        (0:nx+1,0:ny+1),    &
                A        (0:nx+1,0:ny+1),    &
                hn1      (0:nx+1,0:ny+1),    &
                An1      (0:nx+1,0:ny+1),    &
                tracer   (0:nx+1,0:ny+1,2), &
                uice     (0:nx+2,0:ny+2), &
                vice     (0:nx+2,0:ny+2), &
                un1      (0:nx+2,0:ny+2), &
                vn1      (0:nx+2,0:ny+2), &
                un2      (0:nx+2,0:ny+2), &
                vn2      (0:nx+2,0:ny+2), &
                Pp       (0:nx+1,0:ny+1), &
                P        (0:nx+1,0:ny+1), &                
                Pf       (0:nx+1,0:ny+1), &
                etaC     (0:nx+1,0:ny+1), &
                etaB     (0:nx+2,0:ny+2), &
                zetaC    (0:nx+1,0:ny+1), &
                etaCf    (0:nx+1,0:ny+1), &
                etaBf    (0:nx+2,0:ny+2), &
                zetaCf   (0:nx+1,0:ny+1)

      common/DynVariables/      &
                h,              & ! ice thickness ( = tracer(i,j,1) )
                A,              & ! ice concentration ( = tracer(i,j,2) )
                hn1,            & ! previous time step ice thickness
                An1,            & ! previous time step ice concentration
                tracer,         & ! passive tracer 
                uice,           & ! x-comp ice velocity 
                vice,           & ! y-comp ice velocity
                un1,            & ! previous time level solution 
                vn1,            & ! previous time level solution
                un2,            & ! previous time level solution
                vn2,            & ! previous time level solution
                Pp,             & ! ice strength
                P,              & ! replacement pressure
                Pf,             & ! replacement pressure
                etaC,           & ! coefficient of shear viscosity (C-grid)
                etaB,           & ! coefficient of shear viscosity (B-grid)
                zetaC,          & ! coefficient of bulk viscosity (C-grid)
                etaCf,          & ! coefficient of shear viscosity (C-grid)
                etaBf,          & ! coefficient of shear viscosity (B-grid)
                zetaCf           ! coefficient of bulk viscosity (C-grid)





