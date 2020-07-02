!*************************************************************************
!     subroutine stepper:
!       Calculate the ice thickness (h), concentration (A) and ice velocity
!       (u_ice, v_ice) at the next time level, for a given forcing.
!       The ocean currents and air temperatures are prescribed/calculated. The 
!       standard time integration scheme is based on a splitting in time:
!       
!       momentum equation:  rho*h^n-1(u^n-u^n-1)/Deltat = f(u^n,h^n-1,A^n-1)
!                           This equation is solved implicitly for u^n
!                           (uice). u^n-1 is the previous time level solution.
!                           h and A are explicit as they are defined at n-1.  
!
!       continuity equation:h^n = f(h^n-1, u^n)  
!                           The new value of h (and A) at time n is obtained 
!                           by advecting (+thermo) h^n-1 with u^n.
!
!       Note that this splitting in time (of the continuity and momentum 
!       equations) leads to an instability if Deltat is too large. This 
!       problem is more severe at high resolution. See the following paper 
!       for details:
!
!       Lipscomb et al., Ridging, strength and stability in high-resolution
!       sea ice models, 112, C03S91, JGR, 2007.
!
!       By setting IMEX=2 (for JFNK), these equations are solved at the same time
!       for u^n, h^n and A^n (without thermo however which is done after). Using 
!       the upwindRK2 advection scheme and setting BDF=1 (backward difference 
!       formula) leads to second-order accuracy in time.
!
!************************************************************************


      subroutine stepper (date, tstep, expno)
        USE datetime, ONLY: datetime_type, datetime_str, datetime_str_6, time_init, time_set_from_datetime
        use datetime, only: operator(==), operator(/=)
        use io, only: daily_air_temp_from_monthly_mean
        use numerical_VP
        use solver_choice
      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_DynVariables.h'
      include 'CB_const.h'
      include 'CB_buoys.h'
      include 'CB_ThermoVariables.h'
      
      TYPE(datetime_type), INTENT(in) :: date

      character(LEN=6) :: datestr
      
      integer, intent(in) :: tstep, expno
      integer ::  year    ! current year

      integer :: k, tot_its, i, j
      integer, save :: sumtot_its, nbfail

      double precision :: h2sec = 3.6d03            ! [sec] / [hour]
      double precision :: xtp(nvar), rhs(nvar), Fu(nvar)
      double precision :: uk2(0:nx+2,0:ny+2), vk2(0:nx+2,0:ny+2)
      double precision :: ul(0:nx+2,0:ny+2), vl(0:nx+2,0:ny+2)
      double precision :: res, resk_1, time1, time2, timecrap
      double precision :: maxdu, tpcalc, thmaxdu

      double precision, save :: NLtol

      year = date%year
      datestr = datetime_str_6(date)

      if (tstep .eq. 1) nbfail = 0 ! number of failures of the nonlinear solver during the run

      if ( Dynamic ) then

         if ( Wind .eq. '6hours' ) then

            call wind_forcing (date) ! get wind forcing field

         endif

         if ( BDF .eq. 1 ) then
            un2 = un1
            vn2 = vn1
         endif

         un1 = uice ! previous time step solution
         vn1 = vice ! previous time step solution
         hn1 = h
         An1 = A

!------- Set initial guess to freedrift if required (if not, PTS is used)-
      
         if (ini_guess .eq. 'freedrift') call UVsolveNR_B

!------- Calc ice strength and part of b vector if IMEX=0 -----------------

         if ( IMEX .eq. 0 ) then
            call Ice_strength()
            if (solver .le. 2) then ! Picard or JFNK
               call bvect_ind ! function of h ( not directly f(u) )
            endif
         endif

!------- Solves NL mom eqn at specific time step with solver1 or solver2
!        F(u) = A(u)u - b(u) = 0, u is the solution vector

         call cpu_time(timecrap)
         call cpu_time(time1)

         if (solver .eq. 1) then !standard solver

!------- Set initial value of uk2, vk2, ul and vl to value of ini guess --

            uk2 = uice
            vk2 = vice
            ul  = uice
            vl  = vice

!------- Initialize stuff for tolerance and FGMRES its count -------------

            sumtot_its = 0

!------- Beginning of outer loop (OL) iterations -------------------------

            do k = 1, OLmax 
               
               call uv_linearization (uk2, vk2, ul, vl) ! get ul for lin syst

               uk2  = uice ! uk2 is u^k-2, uk1 is u^k-1 = uice here
               vk2  = vice

               call transformer (uice,vice,xtp,1)

               if ( IMEX .eq. 1 ) then ! IMEX 1 (2 doesn't work with Picard) 
                  call advection ( un1, vn1, uice, vice, hn1, An1, h, A )
                  call Ice_strength()
                  call bvect_ind
               endif

               call ViscousCoefficient(uice,vice)
               call bvect (ul,vl,rhs) ! forms b(ul)
               call Funk(xtp,rhs,Fu)

               res = sqrt(DOT_PRODUCT(Fu,Fu)) ! L2norm  
               
               if (k.eq. 1) NLtol = gamma_nl * res

               if (res .lt. NLtol .or. res .lt. 1d-06) then
                  print *, 'L2norm is', k,res,'(final)'
                  print *, 'nb outer ite, FGMRES ite =',k-1, sumtot_its
                  exit
               endif

               call PrepFGMRES (k,res,tot_its,xtp,rhs)  ! FGMRES solver
!               call UVsolveSOR              ! SOR linear solver
               ! solves A(ul^k)u^k = b(ul^k), uice = u^k
               sumtot_its = sumtot_its + tot_its

               if (k .eq. OLmax) then
                  print *, 'WARNING Picard solver DID NOT CONVERGED'
                  nbfail = nbfail + 1
               endif

            enddo


!------- End of outer loop (OL) iterations ------------------------------

         elseif (solver .eq. 2) then ! JFNK solver

!------- Initialize stuff for tolerance and FGMRES its count ------------
            
            sumtot_its = 0 ! grand total of gmres ite for Newton loop

!------- Beginning of Newton loop ---------------------------------------

            do k = 1, NLmax

               call transformer (uice,vice,xtp,1)

               if ( IMEX .gt. 0 ) then ! IMEX method 1 or 2                     
                  call advection ( un1, vn1, uice, vice, hn1, An1, h, A )
                  call Ice_strength()
                  call bvect_ind
               endif

               if ( k .le. klinesearch ) then        ! if k .gt. klinesearch zeta,
                  call ViscousCoefficient(uice,vice) ! eta,Cw,Fu were already      
                  call bvect (uice,vice,rhs)         ! calc in linesearch at
                  call Funk (xtp,rhs,Fu)             ! at end of last Newton ite.
               endif

               res = sqrt(DOT_PRODUCT(Fu,Fu)) ! L2norm
               if (k.eq. 1) then
                  NLtol = gamma_nl * res
                  res_t = res / dropini !transition between fast & slow phases 
               endif
               
               if (res .lt. NLtol .or. res .lt. 1d-06) then
                  print *, 'L2norm is', k,res,'(final)'
                  print *, 'nb Newton ite, FGMRES ite =',k-1, sumtot_its
                  exit
               endif

!               if (k .eq. 499) call failure(xtp,rhs,date,k,expno) !debug
               
               call PrepFGMRES_NK (k,res,resk_1,tot_its,xtp,Fu) 
               ! solves J(u^k-1)du^k = -F(u^k-1), du^k = u^k - u^k-1
               sumtot_its = sumtot_its + tot_its

               if (k .eq. NLmax) then
                  print *, 'WARNING JFNK DID NOT CONVERGED'
                  nbfail = nbfail + 1
               endif

            enddo

!------- End of Newton loop ----------------------------------------------            
         
         elseif (solver .eq. 3) then ! EVP solver 

            call evp_solver(tstep)

         endif

         call cpu_time(time2)
         print *, 'cpu time =', time2-time1
         print *, 'Total nb of failures during the simulation =', nbfail

      endif

!------- calc max(|uice-un1|) to check if steady-state is obtained ----------      

      maxdu=0d0
      do j=1,ny
         do i=1,nx+1
            tpcalc=abs(uice(i,j)-un1(i,j))
            if ( tpcalc .gt. maxdu ) maxdu=tpcalc
         enddo
      enddo

      do j=1,ny+1
         do i=1,nx
            tpcalc=abs(vice(i,j)-vn1(i,j))
            if ( tpcalc .gt. maxdu ) maxdu=tpcalc
         enddo
      enddo

      print *, '**** max delta uice ****=', maxdu

!            if (tstep .eq. 24) then
      thmaxdu=1d-09
      if (maxdu .lt. thmaxdu) then
      !               call check_if_plastic(uice, vice)
         call stress_strain (uice, vice, date, 9, expno)                
      !               call stress_strainB (uice, vice, expno)
         call var_analysis
         stop                                                       
      endif

!------------------------------------------------------------------------
!     Integrate the continuity equations to get  h,A^t (and other tracers) 
!     from h,A^t-1. We use uice and vice (u^t and v&t) to advect the tracers.
!------------------------------------------------------------------------

      if ( Dynamic ) then

         if (IMEX .eq. 0 .and. .not. stressBC) then ! already done with IMEX 1 and 2
            call advection ( un1, vn1, uice, vice, hn1, An1, h, A )
         endif
         tracer(:,:,1) = h
         tracer(:,:,2) = A
            
      endif

!         call diffusion ( tracer )   ! not modified yet

! we should have monthly winds fr thermo forcing...modify load_forcing.f
         
      if ( Thermodyn ) then
           
         ! This doesn't deal with the climatological temperatures or the specified temperatures.
         select case (AirTemp)
         case( "MonthlyMean") 
            CALL daily_air_temp_from_monthly_mean(date, int(Deltax)/1000, Ta)
         case("Specified")
            Ta = 263.15d0
         end select
         
         call thermodynamic (date)
         
         call update_tracer
      
      endif

!------- Get some statistics for h, A, u and v --------------------------

      call var_analysis
      print *, ''
!------------------------------------------------------------------------
!    Advect buoys using uice & vice to the position at the next time step
!------------------------------------------------------------------------


      if ( BuoyTrack ) then

         call LagrangianTracer (date)
          
         if ( (datestr .ne. enddate) .and.  & 
            ( (date%hour + int(Deltat/h2sec) ) .eq. 24 ) )  & 
             
              call post_buoys

      endif

!      call diagnostics (date)
      
      return
    end subroutine stepper
      

