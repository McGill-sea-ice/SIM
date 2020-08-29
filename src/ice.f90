!************************************************************************
!
!                        PROGRAM ICE
!
!     Sea ice dynamic model: viscous plastic rheology (ellipse, Mohr-Coulomb)
!                            
!     Advection scheme     : upstream scheme
!
!     Divergence free ocean currents, geostrophic winds and atmospheric 
!     temperatures need to be prescribed. The ocean temperature is calculated.
!     
!     Appropriate references: 
!
!     Lemieux, J.-F. et. al., A second-order accurate in time IMEX 
!     integration scheme for sea ice dynamics, J. of Comp. Phys. (submitted). 
!
!     Lemieux, J.-F. et. al. (2012), A comparison of the Jacobian-free Newton-Krylov 
!     method and the EVP model for solving the sea ice momentum equation with a viscous-
!     plastic formulation: a serial algorithm study, J. of Comp. Phys., 231, 
!     5926-5944, doi:10.1016/j.jcp.2012.05.024.  
!
!     Lemieux, J.-F. et. al. (2010), Improving the numerical convergence of
!     viscous-plastic sea ice models with the Jacobian-free Newton-Krylov method,
!     J. of Comp. Phys., 229, 2840-2852, doi:10.1016/j.jcp.2009.12.011.
!          
!     Lemieux, J.-F., and B. Tremblay (2009), Numerical convergence of 
!     viscous-plastic sea ice models, J. Geophys. Res., 114, C05009, 
!     doi:10.1029/2008JC005017.
!
!     Lemieux, J.-F., B. Tremblay, S. Thomas, J. Sedlacek, and L.Mysak (2008),
!     Using the preconditioned Generalized Minimum RESidual (GMRES) method to 
!     solve the sea-ice momentum equation, J. Geophys. Res., 113, C10004,
!     doi:10.1029/2007JC004680.
!
!     Tremblay, L.-B. and Mysak, L.A., "Modelling sea ice as a granular 
!     material including the dilatancy effect", JPO, 27, 2342-2360, 1997. 
!
!     Hibler, W. "A dynamic thermodynamic sea ice model", JPO, 9,815-846, 1979.
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Authors
!
!     V01             14-05-97               B. Tremblay
!     V2.0            16-10-06               B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & B. Tremblay
!     V4.0            22-10-08               D. Huard, JF Lemieux, B. Tremblay
!
!     V2.0 is the first version of the dimentionalized code. V01 was 
!     undimentionalized. V3.0 is the 1st official code with GMRES.
!
!     V4.0 Multi-dimension grid (10,20,40,80) with new forcings (see sim_tools).
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!========================================================================
!
!     grid numbering is as follows: 'x' are extra grid cells used to
!     define the boundary conditions
!
!
!                          B-grid                       C-grid
!
!
!                     -|--------------|-         -|--------------|-
!                      |              |           |              |
!                      |              |           |              |
! variable location    |   Psi(i,j)   |        u(i,j) Psi(i,j)   |
!                      |              |           |              |
!                u(i,j)|              |           |              |
!                     -|--------------|-         -|----v(i,j)----|-
!                       v(i,j)
!
!
!                       |---|---|---|-             |---|---|---|-
!                       | x |   |   |              | x | . | . |
!                       |---|---|---|-             |---|---|---|-
! numbering             | x |   |   |              |0,1|1,1| . |
!                      0,1-1,1--|---|-             |---|---|---|-
!                       | x | x | x |              |0,0|1,0| x |
!                      0,0-1,0--|---|-             |---|---|---|-
!
!========================================================================


PROGRAM ice
    USE datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    USE datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    USE datetime, ONLY: str2dt, datetime_str_6, now, delta_str

      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_const.h'
      include 'CB_ThermoForcing.h'
      include 'CB_DynVariables.h'
      include 'CB_mask.h'
      include 'CB_buoys.h'
      include 'CB_Dyndim.h' 

      character(len=40) :: datestring
      integer ndays(12)
      INTEGER i, k, m_current, ts_m_counter, tstep
      INTEGER, parameter :: delta = 10
      INTEGER expno_r, expno, restart, readnamelist
      TYPE(datetime_type) :: post_date(1000), start_date, end_date
      TYPE(datetime_type) :: now_date, restart_date, tic, tac
      type(datetime_delta_type) :: date_step
     
      double precision, DIMENSION(:,:), ALLOCATABLE :: hmean,Amean,umean,vmean

      data ndays   /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

!------------------------------------------------------------------------
!     Read run information
!------------------------------------------------------------------------

      tic = now()

      print *, 'Read namelist?'
      read *, readnamelist
      PRINT *, readnamelist

      print *, 'Use restart file for (h,A,u,v)?'
      read  *, restart

      if (restart .eq. 1) then
         PRINT *, "Restart date :"
         READ (*,'(a19)') datestring
         restart_date = str2dt(datestring)
         READ(*, '(i2)') expno_r
         WRITE(*,*) "Restart on: ", datetime_str(restart_date), " -- exp #", expno_r
      endif

      print *, 'experiment #?'
      read  *, expno 
      PRINT *, expno

      print *, 'Starting date?'
      read (*,'(a19)') datestring
      start_date = str2dt(datestring)
      PRINT *, 'start date : ', datetime_str(start_date)

      
      print *, 'End date?'
      read (*,'(a19)') datestring
      end_date = str2dt(datestring)
      PRINT *, 'end date : ', datetime_str(end_date)

      print *, 'Posting date, time?'

      datestring=''
      DO i=1,1000
         READ (*,'(a19)') datestring
         IF (datestring == 'stop') EXIT         
         post_date(i) = str2dt(datestring)
         print *, 'postdate ', datetime_str(post_date(i))
      END DO

      startdate =  datetime_str_6(start_date)
      enddate =   datetime_str_6(end_date)
      ! e.g. start_date = 1990-01-01-00-00, startdate = 010190 (old stuff)

      call get_default           ! get default settings and parameters
      if (readnamelist .eq. 1) then
         call read_namelist      ! overwrite default based on namelist
      endif
      call verify_options        ! verify validity of options
      call get_mask_and_bathy

      ! This is a datetime delta. It can be added to a
      ! datetime object. 
      date_step = delta_init(seconds=int(Deltat))

      call ocn_current               ! load current data
      call ocn_Tclim                 ! load monthly clim ocean T
      call ini_get (restart, expno_r,restart_date)! ini conds

      if ( Wind .eq. 'specified' .or. Wind .eq. '60yrs_clim') then      
         print *, 'calling specified or 60yrs_clim wind'
         call wind_forcing (start_date)
      endif

!------------------------------------------------------------------------
!     Allocate variables to calc monthly mean fields if specified by user
!------------------------------------------------------------------------

      if ( calc_month_mean ) then

         m_current = start_date%month ! current month of mean fields
         ts_m_counter = 0 ! counts nb of time step in the current month
         
         call prep_monthly_mean_fields( end_date, start_date, date_step )

         ALLOCATE( hmean(0:nx+1,0:ny+1), Amean(0:nx+1,0:ny+1) )
         ALLOCATE( umean(0:nx+2,0:ny+2), vmean(0:nx+2,0:ny+2) )

         hmean = 0d0 ! initialize to zero
         Amean = 0d0
         umean = 0d0
         vmean = 0d0

      endif

!------------------------------------------------------------------------
!     Step the model forward from stardate to enddate and print the
!     results every postdate
!------------------------------------------------------------------------

      k = 1
     
      now_date = start_date + date_step
      tstep = 0
      do while (now_date < end_date)
         
            if (BuoyTrack) then

               if ( (now_date%year >= 1979) .and. (now_date%hour == 0) ) then
                 call ini_buoys_get(now_date)
               endif

            endif

         tstep = tstep + 1
         print *, datetime_str(now_date), tstep
    
         call stepper (now_date, tstep, expno)

         if (tstep .eq. 1) call info_file (expno)
         if (now_date .eq. post_date(k)) then
            write(*,*) k, ' -- Calling sea-ice_post'
            call sea_ice_post (now_date, expno)
            k = k + 1
         endif

         now_date = now_date + date_step

         if (calc_month_mean) then   ! calc monthly mean fields if specified
            call monthly_mean_fields ( hmean, Amean, umean, vmean, expno, &
                                       m_current, ts_m_counter, end_date, &
                                       now_date )
         endif

      enddo

      tac = now()
      WRITE(*,*) "Total execution time ", delta_str(tac-tic)
      
    END PROGRAM ice


