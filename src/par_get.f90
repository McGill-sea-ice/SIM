!************************************************************************
!     Subroutine get_default: set options, constants, physical and grid parameters 
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             14-05-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************

      subroutine get_default

        use ellipse
        use ice_albedo
        use numerical_VP
        use numerical_EVP
        use solver_choice
        use basal_param
        
        implicit none
        
      include 'parameter.h' 

      include 'CB_const.h'
      include 'CB_Dyndim.h'
      include 'CB_Thermodim.h'
      include 'CB_options.h'
      include 'CB_mask.h'
      include 'CB_ThermoVariables.h'
      include 'CB_buoys.h'
      include 'CB_bathymetry.h'

      double precision f, Clat_ia, Clat_oa 
      double precision Csens_ia, Csens_oa, Csens_oi
      double precision pi, r_earth
      double precision StefanB, Cpair, rhoair
      double precision deg2rad, rad2deg
      double precision Levap, Lsubli, Psurf
      double precision e_ratio, Cdair, Cdwater
      double precision x1, y1, r1, rs, tanteta
      double precision lat(0:nx+1,0:ny+1), long(0:nx+1,0:ny+1)
      integer i, j

!------------------------------------------------------------------------
!     set run parameters (dynamic - thermodynamic - options - domain)
!------------------------------------------------------------------------

      f          =  1.46d-04         ! Coriolis parameter [1/s] 
      Cdair      =  1.2d-03          ! air-ice drag coeffient []1.2e-03
      Cdwater    =  5.5d-03          ! water-ice drag coeffient[]5.5e-03
      theta_a    =  25d0             ! wind turning angle [degree] 
      theta_w    =  25d0             ! water turning angle [degree]
      Pstar      =  27.5d03          ! ice yield stress [N/m2] 
      C          =  20d0             ! ice concentration parameter  
      !phi        =  30d0             ! internal angle of friction
      !delta      =  10d0             ! angle of dilatancy
      !Cohe       =  0d0 !4d03        ! cohesion (tensile strght) [N/m2]
      !etamax     =  1.0d12           ! max shear viscosity
      e_ratio    = 2.0d0             ! ellipse aspect ratio

      td         =  1d0 * 24*60*60   ! damage time scale [s] 
      th         =  3d1 * 24*60*60   ! healing time scale [s]
      nexp       =  5d0              ! exponent of damage source term
 
      Clat_ia    =  1d-03            ! LH transfer coeff (ice/atm) []
      Clat_oa    =  1d-03            ! LH transfer coeff (ocn/atm) []
      Csens_ia   =  1d-03            ! SH transfer coeff (ice/atm) []
      Csens_oa   =  1d-03            ! SH transfer coeff (ocn/atm) []
      Csens_oi   =  1d-03            ! SH transfer coeff (ocn/ice) []

      hmin       =  0.5d0            ! open water equivalent ice thick

      Hocn       = 100d0             ! Ocean mixed layer depth [m]

      Kocn       =  1d11             ! ocean  diffusion coeff [m2/s]

      relhum     =  0.8d0            ! atmosphere relative humidity

      BndyCond   = 'noslip'          ! noslip
      Rheology   = 1                 ! ellipse = 1, triangle = 2
      Damage     = .true.            ! damage parameter 
      linearization = 'Zhang'        ! Tremblay, Zhang
      regularization = 'tanh'        ! tanh, Kreyscher, capping
      visc_method = 2                ! see viscousCoeff routine for details
      ini_guess  = 'previous time step' ! freedrift, previous time step
      adv_scheme = 'upwind'          ! upwind, upwindRK2, semilag 
      IMEX       = 0                 ! 0:split in time, 1:Picard, 2:JFNK
      BDF        = 0                 ! 0: back. Euler, 1: 2nd order back. diff. formula
      Dynamic    = .true.            ! ice model type
      Thermodyn  = .true.            ! ice model type
      BuoyTrack  = .false.
      Buoys      = 'Daily'           ! Buoy traj: 'Track' or 'Daily'
      Current    = 'YearlyMean'      ! YearlyMean, specified
      Wind       = '6hours'          ! 6hours, 60yrs_clim, specified
      RampupWind = .false.           ! smooth increase of specified wind
      AirTemp    = 'MonthlyMean'     ! MonthlyMean, specified (-10C)
      OcnTemp    = 'calculated'      ! MonthlyClim, specified,calculated
      calc_month_mean = .false.      ! to calc monthly mean fields
      runoff     = .false.

      if ((nx == 518) .and. (ny == 438)) then
         Deltax     =  10d03           ! grid size [m] 
      elseif  ((nx == 258) .and. (ny == 218)) then
         Deltax     =  20d03           ! grid size [m] 
      elseif  ((nx == 128) .and. (ny == 108)) then
         Deltax     =  40d03           ! grid size [m] 
      elseif  ((nx == 63) .and. (ny == 53)) then
         Deltax     =  80d03           ! grid size [m] 
      else
         write(*,*) "Wrong grid size dimenions.", nx, ny
         STOP
      endif

      Deltax2 = Deltax**2d0

!------------------------------------------------------------------------
!     Numerical parameters
!------------------------------------------------------------------------      
      
      solver = 2               ! 1 = Picard, 2 = JFNK, 3 = EVP   

      wjac  = 0.575d0
      wsor  = 0.95d0           ! relaxation parameter for SOR precond
      wlsor = 1.30d0           ! relaxation parameter for SOR precond
      kjac  = 10               !
      ksor  = 10               ! nb of ite of precond SOR
      klsor = 10               ! nb of ite of precond line SOR

      gamma_nl = 1d-03         ! nonlinear convergence criterion for JFNK 
      dropini = 1.5d0          ! res_t = L2norm_ini/dropini (L2norm_ini: beg of Newton loop)
      NLmax = 200              ! max nb of Newton loop for JFNK
      OLmax = 500              ! max nb of Outer loop for Picard
      klinesearch = 1          ! linesearch is applied for JFNK for k .ge. klinesearch
      Jac_finite_diff = 'forward' ! forward, centred (for JFNK solver)

      Nsub  = 120              ! nb of EVP subcycles
      Eo    = 0.36d0           ! T = Eo*Deltat    
      init_stress = 'VP'       ! init stress (VP or zero) EVP(tstep=1, s=1)

!------------------------------------------------------------------------
!     Program constants, physical constants and conversion factors
!------------------------------------------------------------------------

      AMR         =  0.622d0           ! Atomic mass ratio [mH20 / mair]
      pi          =  4d0 * datan(1d0)  ! pi []
      Psurf       =  101.3d03          ! atmosphere surface pressure
      r_earth     =  6370d03           ! radius of the earth [m]
      StefanB     =  5.67d-08          ! StefanBoltzmann const[W/(m2K4)]

      deg2rad     =  pi / 180d0        ! [rad] / [deg]
      rad2deg     =  180d0/pi          ! [deg] / [rad] 

!------------------------------------------------------------------------
!     Time step
!------------------------------------------------------------------------

      Deltat     =  1200d0
      DtoverDx   = Deltat / Deltax

!------------------------------------------------------------------------
!     Material properties (in alphabetical order)
!------------------------------------------------------------------------

      Cpair    = 1d03                ! Specific heat of air.
      Cpwater  = 4d03                ! Specific heat of water

      emisatml = 0.90d0 !0.96d0      ! atmosphere emisivity, downward 
      emisice  = 0.97d0              ! ice emissivity
      emisocn  = 0.96d0              ! ocean emissivity

      Kice     =  2d0                ! ice thermal conductivitiy [W/m K]

      Levap    = 2.50d06             ! Latent heat of evaporation [J/kg]
      Lfusion  = 3.34d05             ! Latent heat of fusion [J/kg]
      Lsubli   = 2.83d06             ! Latent heat of sublimation [J/kg]

      Tof      = 273.15d0 - 1.8d0    ! Freezing point of salty water [K]
      Tif      = 273.15d0            ! Freezing point of fresh water [K]

      rhoair   =  1.3d0              ! air density [kg/m3]
      rhoice   =  9d02               ! ice density [kg/m3]
      rhowater =  1026d0             ! water density [kg/m3]

!------------------------------------------------------------------------                
!     Landfast ice parameters
!------------------------------------------------------------------------      

      CC=20d0
      k1=8d0
      k2=15d0
      umin=5d-05
      crit=5d-04      
      BasalStress = .false. ! T if LF ice basal stress param is used
      
!------------------------------------------------------------------------
!     Parameters (dynamic and thermodynamic)
!------------------------------------------------------------------------

      rhof      =  rhoice * f             
      Cdw       =  rhowater * Cdwater     
      Cda       =  rhoair * Cdair         
      theta_a   =  theta_a * deg2rad      ! wind turning angle [rad]
      theta_w   =  theta_w * deg2rad      ! water turning angle [rad]

      sintheta_a = sin( theta_a )             
      costheta_a = cos( theta_a ) 
      sintheta_w = sin( theta_w )             
      costheta_w = cos( theta_w ) 

! Triangle rheology ... coming soon. 
!      phi       =  phi * deg2rad          ! angle of friction [rad]
!      delta     =  delta * deg2rad        ! angle of dilatancy [rad]
      !tandelta   = tan ( delta )
      !sinphi     = sin (phi)

      ell2       = e_ratio**2
      ell_2      = 1/(e_ratio**2)

      Kemis_i   = emisice  * StefanB
      Kemis_al  = emisatml * StefanB
      Kemis_o   = emisocn  * StefanB


      Klat_ia   = rhoair * Clat_ia * Lsubli * AMR / Psurf
      Klat_oa   = rhoair * Clat_oa * Levap  * AMR / Psurf  
 

      Kadvo     = rhowater * Cpwater             ! cts of ocean adv
      Ksens_io  = rhowater * Csens_oi * Cpwater  ! cts of sensible heat
      Ksens_ai  = rhoair   * Csens_ia * Cpair    ! cts of sensible heat
      Ksens_ao  = rhoair   * Csens_oa * Cpair    ! cts of sensible heat

!------------------------------------------------------------------------
!     latitude and longitude of mask's tracer points
!     same calculation as in mask_gen.f (see p.1017-1018)
!------------------------------------------------------------------------

       do j = 0, ny+1
         do i = 0, nx+1
            
!---------- position on the stereographic plane -------------------------------

            x1 = ( i * Deltax) - dx_pole

            y1 = ( j * Deltax) - dy_pole

            r1 = sqrt ( x1**2 + y1**2 )

!---------- angle of the cone -------------------------------------------------

            tanteta  = r1/(2d0 * r_earth)

!----------- short radius on the sphere ---------------------------------------

            rs = 2d0*r_earth * tanteta / (1+tanteta**2)

!----------- get the latitude -------------------------------------------------
               
            lat(i,j)  = acos ( rs / r_earth )

            sinlat(i,j)   = sin( lat(i,j) )
            coslat(i,j)   = cos( lat(i,j) )

            lat(i,j) = lat(i,j) * rad2deg 
            
            

            if ( x1 .ge. 0.0d0 ) then

               if (x1 .lt. 1.0d-08) x1 = 1.0d-08
              
               if ( y1 .ge. 0.0d0 ) then
     
                  long(i,j) = atan ( y1 / x1 ) * rad2deg + beta

               else
                  
                  long(i,j) = atan ( y1 / x1 ) * rad2deg + beta + &
                              360.0d0

               endif

            else

                  long(i,j) = atan ( y1 / x1 ) * rad2deg + beta + &
                              180.0d0

            endif

            if ( long(i,j) .ge. 360.0d0 ) long(i,j) =long(i,j) - 360.0d0

         enddo
      enddo

      return
    end subroutine get_default

!************************************************************************                    
!     Subroutine read_namelist: read option choices and parameter values
!                               from namelist file                        
!                                                                                            
!     Author: JF Lemieux                                           
!                                                                                            
!************************************************************************ 

subroutine read_namelist

        use ellipse
        use numerical_VP
        use numerical_EVP
        use solver_choice
        use basal_param

      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_const.h'
      include 'CB_Dyndim.h'
      
      integer :: nml_error, filenb
      double precision :: e_ratio
      double precision :: rhoair, Cdair, Cdwater, f
      character filename*32

      !---- namelist variables -------
             
      namelist /option_nml/ &
           Dynamic, Thermodyn,                                  &
           linearization, regularization, ini_guess,            &
           adv_scheme, AirTemp, OcnTemp, Wind, RampupWind,      &
           Current, Rheology, Damage, IMEX, BDF, visc_method,   &
           solver, BasalStress

      namelist /numerical_param_nml/ &
           Deltat, gamma_nl, NLmax, OLmax, Nsub

      namelist /phys_param_nml/ &
           Pstar, C, e_ratio, k1, k2, rhoair, rhoice, rhowater, &
           Cdair, Cdwater, f, td, th, nexp

      filename ='namelistSIM'
      filenb = 10

      print *, 'Reading namelist values'
        
      open (filenb, file=filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
         
      do while (nml_error > 0)
         print*,'Reading option_nml'
         read(filenb, nml=option_nml,iostat=nml_error)
         if (nml_error /= 0) exit
         print*,'Reading other_nml'
         read(filenb, nml=numerical_param_nml,iostat=nml_error)
         if (nml_error /= 0) exit
         print*,'Reading phys_param_nml'
         read(filenb, nml=phys_param_nml,iostat=nml_error)
         print *, nml_error
      enddo

      close(filenb)

      DtoverDx   = Deltat / Deltax
      ell2       = e_ratio**2
      ell_2      = 1/(e_ratio**2)
      Cdw        =  rhowater * Cdwater
      Cda        =  rhoair * Cdair
      rhof       =  rhoice * f
      td         = td * 24*60*60
      th         = th * 24*60*60

      end subroutine read_namelist

      subroutine verify_options

        use ellipse

        use numerical_VP
        use numerical_EVP
        use solver_choice
        use basal_param
        
        implicit none

        include 'parameter.h'
        include 'CB_options.h'
        include 'CB_const.h'

!-------------------------------------------------------------------------                         
!     Verify validity of some inputs                                                               
!-------------------------------------------------------------------------                         

      if (BndyCond .ne. 'noslip') then
         print *, 'Wrong BndyCond chosen by user'
         stop
      endif

      if ( Rheology .ne. 1 .and. Rheology .ne. 2) then
         print *, 'Wrong Rheology chosen by user'
         stop
      endif
      
      if ( linearization .ne. 'Zhang' .and.                            &
           linearization .ne. 'Tremblay' ) then
         print *, 'Wrong linearization chosen by user'
         stop
      endif

      if ( solver .eq. 1 .and. IMEX .eq. 2 ) then
         print *, 'IMEX 2 does not work with Picard solver'
         stop
      endif

      if ( solver .eq. 3 .and. IMEX .gt. 0) then
         print *, 'IMEX does not work with EVP solver'
         stop
      endif

      if ( adv_scheme .ne. 'upwind' .and.                              &
           adv_scheme .ne. 'upwindRK2' .and.                           &
           adv_scheme .ne. 'semilag') then
         print *, 'adv_scheme should be upwind or upwindRK2 or semilag'
         stop
      endif

      if ( ini_guess .ne. 'freedrift' .and.                            &
           ini_guess .ne. 'previous time step' ) then
         print *, 'Wrong initial guess chosen by user'
         stop
      endif

      if ( Current .ne. 'YearlyMean' .and.                             &
           Current .ne. 'specified' ) then
         print *, 'Wrong Current chosen by user'
         stop
      endif

      if ( Wind .ne. '6hours' .and. Wind .ne. 'specified' .and.        &
           Wind .ne. '60yrs_clim' ) then
         print *, 'Wrong Wind chosen by user'
         stop
      endif

      if ( AirTemp .ne. 'MonthlyMean' .and.  AirTemp .ne. 'specified') then
         print *, 'Wrong AirTemp chosen by user'
         stop
      endif

      if ( OcnTemp .ne. 'MonthlyClim' .and. OcnTemp .ne. 'specified' &
           .and. OcnTemp .ne. 'calculated') then
         print *, 'Wrong OcnTemp chosen by user'
         stop
      endif

      if (1d0*Deltat .gt. Deltax) then
         print *, 'CFL condition not respected. Reduce time step'
         stop
      endif

!-------------------------------------------------------------------------                                                  
! print info of the run for the output txt file                                                                                        
!-------------------------------------------------------------------------                                                      
 
      print *,
      print *, 'Rheology      =   ', Rheology
      print *, 'Damage        =   ', Damage
      print *, 'Pstar         =   ', Pstar
      print *, 'linearization =   ', linearization
      print *, 'regularization=   ', regularization
      print *, 'initial guess =   ', ini_guess
      print *, 'Dynamic       =   ', Dynamic
      print *, 'Thermodyn     =   ', Thermodyn
      print *, 'Current       =   ', Current
      print *, 'Wind          =   ', Wind
      print *, 'AirTemp       =   ', AirTemp
      print *, 'OcnTemp       =   ', OcnTemp
      print *, 'Healing [s]   =   ', th
      print *, 'exp of dam    =   ', nexp
      print *,
      print *, 'time step [s] =   ', Deltat
      print *, 'Basal', BasalStress

    end subroutine verify_options

    subroutine get_mask_and_bathy

      use basal_param

      implicit none

      include 'parameter.h'

      include 'CB_const.h'
      include 'CB_options.h'
      include 'CB_mask.h'
      include 'CB_bathymetry.h'

      integer :: i,j
      character(len=2) :: cdelta

!------------------------------------------------------------------------                                     
!     Grid parameter: land mask (grid center), velocity mask (node)                                           
!------------------------------------------------------------------------                                     

      write(cdelta, '(I2)') int(Deltax)/1000

      open (unit = 20, file = 'src/mask'//cdelta//'.dat', status = 'old')

      do j = 0, ny+1               ! land mask                                                                
         read (20,10) ( maskC(i,j), i = 0, nx+1 )
      enddo

      close (unit = 20)
      
10    format (1x,1000(i1)) ! different format because of the grid                                              

!-----------------------------------------                                                                    

      do j = 0, ny+2                   ! velocity mask                                                        
         do i = 0, nx+2
            maskB(i,j) = 0
         enddo
      enddo
      

      do j = 1, ny+1
         do i = 1, nx+1
            
            maskB(i,j) = ( maskC(i,j)   + maskC(i-1,j) +          &
                 maskC(i,j-1) + maskC(i-1,j-1) ) / 4
            
         enddo
      enddo

!------------------------------------------------------------------------                                                              
!     load bathymetry for seabed (basal) stress paramaterization                                                                   
!------------------------------------------------------------------------   

      if (BasalStress) then ! LF ice basal stress param is used                                                                        

         open (unit=21,file='src/bathymetry'//cdelta//'km.dat', status = 'old')

         do j = 0, ny+1               ! bathy                                                                                        
            read (21,*) ( bathy(i,j), i = 0, nx+1 )
         enddo

         close (unit = 21)

         do j=0,ny+1
            do i=0,3
               if (maskC(i,j) .eq. 1) then
                  bathy(i,j)=9999d0
               endif
            enddo
         enddo

         do j=0,ny+1
            do i=nx-2,nx+1
               if (maskC(i,j) .eq. 1) then
                  bathy(i,j)=9999d0
               endif
            enddo
         enddo


         do i=0,nx+1
            do j=0,3
               if (maskC(i,j) .eq. 1) then
                  bathy(i,j)=9999d0
               endif
            enddo
         enddo

         do j = 0, ny+1 ! bathy should be gt 5m (+) for ocean and -10 for land  
            do i = 0, nx+1

               if (maskC(i,j) .eq. 0) then
                  if (bathy(i,j) .ne. -10d0) then
                     print *, 'wrong bathy on land'
                     stop
                  endif
               else
                  if (bathy(i,j) .lt. 4.9999d0) then
                     print *, 'wrong bathy on ocean'
                     stop
                  endif
               endif

            enddo
         enddo

         if ( Current .ne. 'specified' ) then
            print *, 'Currents should be zero (specified) with basal stress param'
            stop
         endif

      endif
      
    end subroutine get_mask_and_bathy
