!------------------------------------------------------------------------
! Defines the size and position of the domain and num parameters
!
! Multi resolution 2008 version 
!
! Resolution [km]      nx        ny
! ---------------------------------
!         10          518       438
!         20          258       218
!         40          128       108
!         80           63        53
! ---------------------------------
!------------------------------------------------------------------------

integer nx, ny, nxh, nyh, ntot, nbuoy, nvar, img, img1

      double precision beta, dx_pole, dy_pole, S0

      parameter (                             &
                nx    = 8,                  & ! x-dim of the domain
                ny    = 8,                  & ! y-dim of the domain
                nxh   = 512,                  & ! define high res grid for coarsegraining
                nyh   = 512,                  & ! define high res grid for coarsegraining
                dx_pole= 2500d3,              & ! tracer point (0,0) x distance from North pole [m] 
                dy_pole= 2250d3,              & ! tracer point (0,0) y distance from North pole [m] 
                beta  = 32.0d0,               & ! angle of the dom wr to Greenwich
                ntot  = nx * ny,              & ! total number of cells
                nvar  = (nx+1)*ny+(ny+1)*nx,  & ! total nb of u,v  var
                img   = 50,                   & ! rstr value for FGMRES
                img1  = img + 1,              & !
                nbuoy = 3000,                 & ! number of buoys
                S0 = 1340d0                   & ! Solar constant [W / m2]
                )



