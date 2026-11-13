!------------------------------------------------------------------------
! Defines the size and position of the domain and num parameters
!
! Multi resolution 2026 version 
! (backwards compatibility is ensured, old grids still work)
! jkr: Definition of dx_pole and dy_pole is now done dynamically
! in par_get.f90
!
! base 32 grids
!
! Resolution dx [km]   nx        ny
! ---------------------------------
!        1.00         5120      4320
!        2.00         2560      2160
!        4.00         1280      1080
!        8.00          640       540
!       16.00          320       270
!       32.00          160       135
! ---------------------------------
!
! base 80 grids
! old grids kept for backwards compatibility
!
! Resolution dx [km]   nx        ny
! ---------------------------------
!       10.00          512       432
!       20.00          256       216
!       40.00          128       108
!       80.00           64        54
! ---------------------------------
!
! Uniaxial loading experiments
! uniaxial=.true.
!
! Resolution dx [km]   nx        ny
! ---------------------------------
!          1          100       250
! ---------------------------------
!
! Ideal ice bridge experiments
! ideal_bridge=.true.
!
! Resolution dx [km]   nx        ny
! ---------------------------------
!          2          102       402
! ---------------------------------
!
!
!------------------------------------------------------------------------

      integer nx, ny, ntot, nbuoy, nvar, img, img1
      logical uniaxial, ideal_bridge
      double precision Deltax, beta, S0

      parameter (                            &
                nx = 160,                    & ! x-dim of the domain
                ny = 135,                    & ! y-dim of the domain
                Deltax = 32d03,              & ! grid resolution
                ntot  = nx * ny,             & ! total number of cells
                nvar  = (nx+1)*ny+(ny+1)*nx, & ! total nb of u,v  var
                uniaxial = .false.,          & ! uniaxial experiment
                ideal_bridge = .false.,      & ! ideal ice bridge experiment
                beta = 32.0d0,               & ! angle of the dom wr to Greenwich
                img = 50,                    & ! rstr value for FGMRES
                img1 = img + 1,              & !
                nbuoy = 3000,                & ! number of buoys
                S0 = 1340d0                  & ! Solar constant [W / m2]
                )
      
