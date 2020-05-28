!========================================================================
!     Common block options: program options
! FB: flag for domain Idealized in per_get.f90
!========================================================================

      logical                   &
                Dynamic,        &
                Thermodyn,      &
                BuoyTrack,      &
                calc_month_mean,&
                runoff,         &
                Idealized_domain
      character(LEN=20)         &
                BndyCond,       &
                DragLaw,        &
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
                          
      integer   Rheology, IMEX, BDF, visc_method

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
                Current,        & ! specified, YearlyMean
                Idealized_domain, & ! FB: true when domain is Idealized, false when domain is the Arctic 
		adv_scheme,     & ! advection scheme: upwind or upwindRK2
                BndyCond,       & ! noslip or freeslip
                DragLaw,        & ! linear, square, linearH or squareH
		Rheology,       & ! FB: ellipse=1, triangle=2
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

