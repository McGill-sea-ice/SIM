!========================================================================
!     Common block const: program constants
!========================================================================

      integer ntracer

      double precision                 &
                sinlat(0:nx+1,0:ny+1), &
                coslat(0:nx+1,0:ny+1)

      double precision relhum, Deltax, Deltax2, theta_a, theta_w
      double precision sintheta_a, costheta_a, sintheta_w, costheta_w
      double precision Deltat, DtoverDx, rhof, Tif, Tof
      double precision emisice, emisatml,emisocn
      double precision AMR, hmin, Cpwater, rhowater, rhoice, Lfusion


      common/const/             &
                sinlat,         & ! sin (lat) of every grid point
                coslat,         & ! cos (lat) of every grid point
                relhum,         & ! atmosphere relative humidity
                Deltax,         & ! grid size [m]
		Deltax2,        & ! Deltax**2
		Deltat,         & ! time step
                DtoverDx,       & ! Deltat / Deltax
                theta_a,        & ! air turning angle
                theta_w,        & ! water turning angle
                sintheta_a,     & ! sin (theta_a)
                costheta_a,     & ! cos (theta_a)
                sintheta_w,     & ! sin (theta_w)
                costheta_w,     & ! cos (theta_w)
                rhof,           & ! rhoice * Coriolis parameter
                Tif,            & ! ice freezing point temperature
                Tof               ! ocn freezing point temperature



      common/const/             &
                emisice,        & ! ice emissivity
                emisatml,       & ! atmosphere emissivity
                emisocn,        & ! ocean emissivity
                AMR,            & ! atomic mass ration
                hmin,           & ! open water equivalent thickness
                Cpwater,        & ! Specific Heat of water
                rhowater,       & ! water density
                rhoice,         & ! ice density
                Lfusion           ! Latent heat of fusion [J/kg]

      common/const/  &
                ntracer           ! number of tracer ( max = 10 )

