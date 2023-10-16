!========================================================================
!     Common block options: program options
!========================================================================

      logical                   &
                Dynamic,        &
                Thermodyn,      &
                BuoyTrack,      &
                calc_month_mean,&
                RampupWind,     &
                RampupForcing,  &
                runoff

      character(LEN=20)         &
                BndyCond,       &
                linearization,  &
                regularization, &
                ini_guess,      &
                AirTemp,        &
                OcnTemp,        &
                Wind,           &
                Current,        &
                adv_scheme,     &
                startdate,      &
                enddate,        &
                Buoys,          &
                Jac_finite_diff

      integer   Rheology, IMEX, BDF, visc_method, Periodic_x, Periodic_y

      common/options/           &
                Dynamic,        & ! sea ice dynamic model (yes or no)
                Thermodyn,      & ! sea ice thermodynamic (yes or no)
                BuoyTrack,      & ! Track buoy position (yes or no)
                calc_month_mean,& ! to calc monthly mean fields (yes or no)
                Buoys,          & ! 'Track' or 'Daily'
                runoff,         & ! River runoff switch
                AirTemp,        & ! specified, MonthlyMean
                OcnTemp,        & ! specified, MonthlyMean
                Wind,           & ! specified, MonthlyMean
                RampupWind,     & ! smooth increase of specified wind
                RampupForcing,   & ! smooth increase of specified wind forcing
                Current,        & ! specified, YearlyMean
                adv_scheme,     & ! advection scheme: upwind or upwindRK2
                BndyCond,       & ! noslip or freeslip
                Periodic_x,     & ! open or periodic condition in x
                Periodic_y,     & ! open or periodic condition in y
		Rheology,       & ! ellipse, triangle, MEB
		IMEX,           & ! 0: standard (splitting in time), 1 and 2: IMEX
                BDF,            & ! 0: standard, 1: Backward diff formula (2nd order)
                visc_method,    & ! choice of calc of visc coeff
                startdate,      & ! starting date
                enddate           ! end date

      common/options/           &
                linearization,  & ! Tremblay, Zhang
                regularization, & ! tanh, Kreyscher
                ini_guess,      & ! freedrift, previous time step
                Jac_finite_diff

