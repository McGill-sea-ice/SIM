!************************************************************************
!     Subroutine stress_strain: computes the stresses (invariants), the 
!     normalized stresses (invariants) and the strain rates (shear + div) 
!     for the ellipse yield curve.
!
!     Note: calc only if A(i,j) > 0.5
!     
!     land poitns have values of -999 and ocean points with A(i,j)< 0.5
!     have values of -888.
!
!     Should be put in stepper after the calculation of u and v. u and v 
!     are therefore at iteration k while the visc. coeff. are at k-1.
!
!     JF Lemieux, 4 August 2009
!
!************************************************************************


      subroutine stress_strain(utp, vtp, date, k, expno)
        use datetime, only: datetime_type
        use ellipse
      implicit none

      include 'parameter.h'
      include 'CB_Dyndim.h'
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      type(datetime_type), intent(in) :: date

      character filename*45

      integer, intent(in) :: k, expno
      integer i, j, m, year, month, day, hour, minute

      double precision dudx, dvdy, dudy, dvdx, land, lowA

      double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision shear(0:nx+1,0:ny+1), div(0:nx+1,0:ny+1)
      double precision sigI(0:nx+1,0:ny+1), sigII(0:nx+1,0:ny+1)
      double precision sigInorm(0:nx+1,0:ny+1), sigIInorm(0:nx+1,0:ny+1)! stress invariants
      double precision sig1norm(0:nx+1,0:ny+1), sig2norm(0:nx+1,0:ny+1) ! princ stresses
      double precision zetaCout(0:nx+1,0:ny+1)

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute

      land=-999d0
      lowA=-888d0

!      land=0d0
!      lowA=0d0

      div       = land
      shear     = land
      sigI      = land
      sigII     = land
      sigInorm  = land
      sigIInorm = land
      sig1norm  = land
      sig2norm  = land
      zetaCout  = land

! note: to study the numerical convergence of the stress, zeta and eta should be calculated 
!       with u^{k-1} and the deformations with u^k. This is why we use here zetaCf and etaCf
!       Calculating both with u^k leads to stress that are all VP whatever the level of 
!       convergence of the solution.

         do i = 1, nx
            do j = 1, ny

               dudx       = 0d0
               dvdy       = 0d0
               dudy       = 0d0
               dvdx       = 0d0
               
               if ( maskC(i,j) .eq. 1 ) then

                  zetaCout(i,j) = zetaCf(i,j) ! no special value for A<0.5

                  if (A(i,j) .lt. 0.5d0) then
                     div(i,j)       = lowA
                     shear(i,j)     = lowA
                     sigI(i,j)      = lowA
                     sigII(i,j)     = lowA
                     sigInorm(i,j)  = lowA
                     sigIInorm(i,j) = lowA

                  elseif (A(i,j) .ge. 0.5d0) then
                  
                  dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                  dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
                  if     ( maskC(i+1,j) + maskC(i-1,j) .eq. 2 ) then
                     
                     dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                              ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                              ( 4d0 * Deltax )
                     
                  elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. 1 ) then
                     
                     dvdx = ( 1d0 * ( vtp(i+1,j) + vtp(i+1,j+1) ) +  &
                              3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) /  &
                              ( 6d0 * Deltax )
                     
                  elseif ( maskC(i+1,j) - maskC(i-1,j) .eq. -1 ) then
                     
                     dvdx = ( -1d0 * ( vtp(i-1,j) + vtp(i-1,j+1) ) - &
                               3d0 * ( vtp(i,j)   + vtp(i,j+1) ) ) / &
                             ( 6d0 * Deltax )
                     
                  elseif ( maskC(i+1,j) + maskC(i-1,j) .eq. 0 ) then
                     
                     print *, 'WARNING: irregular grid cell case1', i, j
                     
                  endif

               
                  if     ( maskC(i,j+1) + maskC(i,j-1) .eq. 2 ) then
                     
                     dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                               ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                               ( 4d0 * Deltax )
                     
                  elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. 1 ) then
                     
                     dudy = ( 1d0 * ( utp(i,j+1) + utp(i+1,j+1) ) +  &
                               3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                               ( 6d0 * Deltax )
                     
                  elseif ( maskC(i,j+1) - maskC(i,j-1) .eq. -1 ) then
                     
                     dudy = ( -1d0 * ( utp(i,j-1) + utp(i+1,j-1) ) - &
                               3d0 * ( utp(i,j)   + utp(i+1,j) ) ) / &
                               ( 6d0 * Deltax )
                     
                  elseif ( maskC(i,j+1) + maskC(i,j-1) .eq. 0 ) then
                     
                     print *, 'WARNING: irregular grid cell case2',i,j
                     
                  endif
                  
!----- stresses and strain rates at the grid center -------------------------   


                  shear(i,j) = sqrt(( dudx - dvdy )**2d0 &  ! in day-1
                       + ( dudy + dvdx )**2d0 )*86400d0

                  div(i,j)   = ( dudx + dvdy )*86400d0      ! in day-1

! watchout p in our code is in fact p/2 in Hibler's equations

                  sigI(i,j)   = -1d0*( dudx + dvdy )*zetaCf(i,j)+Pf(i,j)

                  sigII(i,j) = sqrt(( dudx - dvdy )**2d0 &
                       + ( dudy + dvdx )**2d0 )*etaCf(i,j)

                  sigInorm(i,j)  = sigI(i,j)  / (2d0*max(Pp(i,j),1d-10))
                  
                  sigIInorm(i,j) = sigII(i,j) / (2d0*max(Pp(i,j),1d-10))

                  sig1norm(i,j) = -1d0*sigInorm(i,j) + sigIInorm(i,j)
                  sig2norm(i,j) = -1d0*sigInorm(i,j) - sigIInorm(i,j)

                  endif

               endif
               
               
            enddo
         enddo

! Take care of 1st row, 1st column, last row and last column

         j=0
         do i=0,nx+1
            if ( maskC(i,j) .eq. 1 ) then
               div(i,j)       = lowA
               shear(i,j)     = lowA
               sigI(i,j)      = lowA
               sigII(i,j)     = lowA
               sigInorm(i,j)  = lowA
               sigIInorm(i,j) = lowA
               sig1norm(i,j)  = lowA
               sig2norm(i,j)  = lowA
               zetaCout(i,j)  = 0d0
            endif
         enddo

         j=ny+1
         do i=0,nx+1
            if ( maskC(i,j) .eq. 1 ) then
               div(i,j)       = lowA
               shear(i,j)     = lowA
               sigI(i,j)      = lowA
               sigII(i,j)     = lowA
               sigInorm(i,j)  = lowA
               sigIInorm(i,j) = lowA
               sig1norm(i,j)  = lowA
               sig2norm(i,j)  = lowA
               zetaCout(i,j)  = 0d0
            endif
         enddo

         i=0
         do j=0,ny+1
            if ( maskC(i,j) .eq. 1 ) then
               div(i,j)       = lowA
               shear(i,j)     = lowA
               sigI(i,j)      = lowA
               sigII(i,j)     = lowA
               sigInorm(i,j)  = lowA
               sigIInorm(i,j) = lowA
               sig1norm(i,j)  = lowA
               sig2norm(i,j)  = lowA
               zetaCout(i,j)  = 0d0
            endif
         enddo

         i=nx+1
         do j=0,ny+1
            if ( maskC(i,j) .eq. 1 ) then
               div(i,j)       = lowA
               shear(i,j)     = lowA
               sigI(i,j)      = lowA
               sigII(i,j)     = lowA
               sigInorm(i,j)  = lowA
               sigIInorm(i,j) = lowA
               sig1norm(i,j)  = lowA
               sig2norm(i,j)  = lowA
               zetaCout(i,j)  = 0d0
            endif
         enddo

      write (filename,'("output/sigI",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (12, file = filename, status = 'unknown')

      write (filename,'("output/sigII",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (13, file = filename, status = 'unknown')

      write (filename,'("output/sigInorm",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (14, file = filename, status = 'unknown')

      write (filename,'("output/sigIInorm",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (15, file = filename, status = 'unknown')

      write (filename,'("output/sig1norm",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (16, file = filename, status = 'unknown')

      write (filename,'("output/sig2norm",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (17, file = filename, status = 'unknown')

      write (filename,'("output/div",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (18, file = filename, status = 'unknown')

      write (filename,'("output/shear",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (19, file = filename, status = 'unknown')

      write (filename,'("output/zetaC",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (20, file = filename, status = 'unknown')

      do j = 0, ny+1
         write(12,100) ( sigI(i,j), i = 0, nx+1 )
         write(13,100) ( sigII(i,j), i = 0, nx+1 )
         write(14,200) ( sigInorm(i,j), i = 0, nx+1 )
         write(15,200) ( sigIInorm(i,j), i = 0, nx+1 )
         write(16,200) ( sig1norm(i,j), i = 0, nx+1 )
         write(17,200) ( sig2norm(i,j), i = 0, nx+1 )
         write(18,200) ( div(i,j), i = 0, nx+1 )
         write(19,200) ( shear(i,j), i = 0, nx+1 )
         write(20,300) ( zetaCout(i,j), i = 0, nx+1 )
      enddo

      do m=12,20
         close(m)
      enddo

100            format (1x, 1000(f20.10, 1x))
200            format (1x, 1000(f15.10, 1x))
300            format (1x, 1000(f20.4,  1x))
      
      return
    end subroutine stress_strain




      subroutine stress_strain_MEB(utp, vtp, date, k, expno)

!************************************************************************
!     This subroutine contains the damage parameterization computations.
!     First, the stress tensor is computed from the strain rates. In the
!     MEB case, the momentum equations uses the shear stress memory component
!     defined at the nodes: so we compute it here at the
!     nodes. The shear strain rate at the grid center is then computed using
!     simple averages of neighbouring shear strain rate. 
!     
!     Then, the yield function is used to find whether we have overshooting.
!     The damage factor is computed to keep the stress state on the yield curve
!
!     Note that the stress memory is defined by the corrected stress state, and
!     is only updated if outside of the outer loop. The damage parameter is
!     also ony updated outside of the outer loop.
!
!     For details on the damage equations, see Plante et al. 2019, 
!     The Cryosphere Discussion
!
!     Mathieu Plante, October 10 2019
!
!************************************************************************


      use datetime, only: datetime_type
      use elastic
      use ellipse
      implicit none

      include 'parameter.h'
      include 'CB_Dyndim.h'
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      character filename*100

      type(datetime_type), intent(in) :: date

      integer, intent(in) :: k, expno
      integer i, j, m, year, month, day, hour, minute, second, milli, peri

      double precision dudx, dvdy, dudy, dvdx, land, lowA, deg2rad
      double precision sigI_uncor, sigII_uncor, eI, eII, a_overshoot
      double precision pi, m1, m2, m3, m4, m5, frict
      double precision Rmax, Rfactor, Ndam, Rabsmax, Rabsfactor

      double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision sigI(0:nx+2,0:ny+2), sigII(0:nx+2,0:ny+2), Fn(0:nx+2,0:ny+2)
      double precision sig1(0:nx+2,0:ny+2), sig2(0:nx+2,0:ny+2), theta(0:nx+2,0:ny+2)
      double precision Dsigxx(0:nx+2,0:ny+2), Dsigyy(0:nx+2,0:ny+2), Dsigxy(0:nx+2,0:ny+2)
      double precision sigxx_it(0:nx+2,0:ny+2), sigyy_it(0:nx+2,0:ny+2), sigxy_it(0:nx+2,0:ny+2)
      double precision DsigxyB(0:nx+2,0:ny+2), Dflag(0:nx+2,0:ny+2)
      double precision Dexx(0:nx+2,0:ny+2), Deyy(0:nx+2,0:ny+2), Dexy(0:nx+2,0:ny+2)
      double precision DexyB(0:nx+2,0:ny+2), R(0:nx+2,0:ny+2)
      double precision sigIM(0:nx+2,0:ny+2), sigIIM(0:nx+2,0:ny+2)

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute
      second = date%second
      milli = date%milli

      peri = Periodic_x + Periodic_y ! =1 if we have periodic conditions
      
      land=-999d0
      lowA=0d0
      pi        = 4d0 * datan(1d0)
      deg2rad   =  pi /180d0

      sigI      = land
      sigII     = land
      Fn        = land
      sig1      = land
      sig2      = land
      frict     = sin(phi*deg2rad)
      Dexx      = 0d0
      Deyy      = 0d0
      Dexy      = 0d0
      DexyB     = 0d0
      Dflag     = 0d0
      dfactor   = 1d0
      dfactorB   = 1d0

      R = 0d0
      Rfactor = 0d0
      Ndam = 0d0
      Rmax = 0d0
      Rabsmax = 0d0
      Rabsfactor = 0d0
      if (peri .ne. 0) call periodicBC(utp,vtp)

      call ViscousCoefficient(utp,vtp)

     
!----------------------------------------
!------------------------------------------
!Calculate the shear strain at the node
!------------------------------------------
!-----------------------------------------

    do i = 1, nx+1
      do j = 1, ny+1

                dvdx  = 0d0
                dudy  = 0d0
                dudx  = 0d0
                dvdy  = 0d0

        if ( maskB(i,j) .eq. 0) then

            if ( (maskC(i,j) + maskC(i-1,j) + &
                maskC(i-1,j-1) + maskC(i,j-1)) .lt. 2) then
            !Do nothing
        ! x o
        ! x x
        !


            elseif ( ((maskC(i-1,j) + maskC(i,j-1)) .eq. 0) .or. &
                ((maskC(i-1,j-1) + maskC(i,j)) .eq. 0) ) then
        ! x o
        ! o x
        !    Do nothing
              print *, 'WARNING: stress-strain irregular grid cell', i, j

            elseif ( (maskC(i,j) + maskC(i,j-1)) .eq. 0) then

                    if ( (maskB(i-1,j) .eq. 0) .or. ((i .eq. 2) &
                    .and. (Periodic_x .eq. 0)) ) then
        ! x o x     !first order app. should not happen
        ! x o x
                        dvdx  = -2d0*vtp(i-1,j) / Deltax
                        dudy  = 0d0

                    else
        ! o o x
        ! o o x
                        dvdx  = ( vtp(i-2,j) - 9d0*vtp(i-1,j)) / (3d0*Deltax)
                        dudy  = 0d0
                    endif

             elseif ( (maskC(i-1,j) + maskC(i-1,j-1)) .eq. 0) then

                    if ( (maskB(i+1,j) .eq. 0) .or. ((i .eq. nx) &
                    .and. (Periodic_x .eq. 0)) ) then
        ! x o x
        ! x o x
                        dvdx  =  2d0*vtp(i,j) / Deltax
                        dudy  = 0d0
                    else
        ! x o o
        ! x o o
                        dvdx  = (9d0*vtp(i,j) - vtp(i+1,j)) / (3d0*Deltax)
                        dudy  = 0d0
                    endif

             elseif ( (maskC(i,j) + maskC(i-1,j)) .eq. 0) then

                    if ( (maskB(i,j-1) .eq. 0) .or. ((j .eq. 2) &
                    .and. (Periodic_y .eq. 0)) ) then
        ! x x
        ! o o
        ! x x
                        dudy  =  -2d0*utp(i,j-1) / Deltax
                        dvdx  = 0d0
                    else
        ! x x
        ! o o
        ! o o
                        dudy  = ( utp(i,j-2) - 9d0*utp(i,j-1)) / (3d0*Deltax)
                        dvdx  = 0d0

                    endif

             elseif ( (maskC(i-1,j-1) + maskC(i,j-1)) .eq. 0) then

                if ( (maskB(i,j+1) .eq. 0) .or. ((j .eq. ny) &
                    .and. (Periodic_y .eq. 0)) ) then
        ! x x
        ! o o
        ! x x
                        dudy  = 2d0*utp(i,j) / Deltax
                        dvdx  = 0d0
                else
        ! o o
        ! o o
        ! x x
                        dudy  = (9d0*utp(i,j) - utp(i,j+1)) / (3d0*Deltax)
                        dvdx  = 0d0

                endif

             elseif ( maskC(i,j) .eq. 0) then

        ! o o x
        ! o o o
                        dvdx  = ( vtp(i-2,j) - 9d0*vtp(i-1,j)) / (3d0*Deltax)
                        dudy  = ( utp(i,j-2) - 9d0*utp(i,j-1)) / (3d0*Deltax)

             elseif ( maskC(i,j-1) .eq. 0) then

        ! o o o
        ! o o x
                    	dvdx  = ( vtp(i-2,j) - 9d0*vtp(i-1,j)) / (3d0*Deltax)
                    	dudy  = (9d0*utp(i,j) - utp(i,j+1)) / (3d0*Deltax)

             elseif ( maskC(i-1,j) .eq. 0) then

        ! x o o 
        ! o o o
                    	dvdx  = (9d0*vtp(i,j) - vtp(i+1,j)) / (3d0*Deltax)
                    	dudy  = ( utp(i,j-2) - 9d0*utp(i,j-1)) / (3d0*Deltax)

             elseif ( maskC(i-1,j-1) .eq. 0) then

        ! o o o 
        ! x o o
                    	dvdx  = (9d0*vtp(i,j) - vtp(i+1,j)) / (3d0*Deltax)
                    	dudy  = (9d0*utp(i,j) - utp(i,j+1)) / (3d0*Deltax)
             endif

        else

        ! o o o 
        ! o o o
                    dvdx  = (vtp(i,j) - vtp(i-1,j)) / Deltax
                    dudy  = (utp(i,j) - utp(i,j-1)) / Deltax

        endif ! endif maskB = 0


        !Open boundaries
        
        if ((i .eq. 1) .and. (Periodic_x .eq. 0)) then

            dvdx = 0d0

        elseif (( i .eq. nx+1) .and. (Periodic_x .eq. 0)) then

            dvdx = 0d0

        endif

        if ((j .eq. 1) .and. (Periodic_y .eq. 0)) then

            dudy = 0d0

        elseif ((j .eq. ny+1) .and. (Periodic_y .eq. 0)) then

            dudy = 0d0

        endif

        !increments of shear strain and stress at nodes
        DexyB(i,j) = ( dvdx + dudy ) / 2d0
        DsigxyB(i,j) = 2d0*etaB(i,j)*DexyB(i,j) !&
                            ! + sigxyB(i,j)/GammaMEB_B(i,j)

      enddo
    enddo
    if (peri .ne. 0) call periodicBC(DexyB,DsigxyB)



!----------------------------------------
!------------------------------------------
!Calculate the normal strain and stress at center
!------------------------------------------
!-----------------------------------------

    do i = 1, nx
      do j = 1, ny

        if ( maskC(i,j) .ge. 1 ) then

            dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
            dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax

            Dexx(i,j) = dudx
            Deyy(i,j) = dvdy
            Dexy(i,j) = ( DexyB(i,j) + DexyB(i,j+1) + &
                        DexyB(i+1,j) + DexyB(i+1,j+1)) / 4d0

            Dsigxx(i,j)  = ( dudx *(zetaC(i,j)+etaC(i,j)) &
                        + (zetaC(i,j)-etaC(i,j))*dvdy )
            Dsigyy(i,j)  = ( dudx *(zetaC(i,j)-etaC(i,j)) &
                        + (zetaC(i,j)+etaC(i,j))*dvdy )
            Dsigxy(i,j)  = (DsigxyB(i,j) + DsigxyB(i,j+1) + &
                        DsigxyB(i+1,j) + DsigxyB(i+1,j+1) ) / 4d0

        endif
      enddo
    enddo
    
    if (peri .ne. 0) call periodicBC(Dexx,Deyy)
    if (peri .ne. 0) call periodicBC(Dsigxx,Dsigyy)
    if (peri .ne. 0) call periodicBC(Dsigxy,Dexy)

!----------------------------------
! Calculate the total normal stress (prev. + increm.) and yield fct at center
!----------------------------------
    do i = 1, nx
       do j = 1, ny

          if ( maskC(i,j) .ge. 1 ) then

            sigxx_it(i,j)  = sigxx(i,j)*GammaMEB(i,j) + Dsigxx(i,j) 
            sigyy_it(i,j)  = sigyy(i,j)*GammaMEB(i,j) + Dsigyy(i,j)
            sigxy_it(i,j)  = sigxy(i,j)*GammaMEB(i,j) + Dsigxy(i,j) 

            sigI(i,j)   = (1/2d0)*( sigxx_it(i,j) + sigyy_it(i,j))
            sigII(i,j) = sqrt( ( (sigyy_it(i,j) - sigxx_it(i,j))/2d0 )**2d0 &
                                     + sigxy_it(i,j)**2d0 )
            sigIM(i,j)   = (1/2d0)*( sigxx(i,j) + sigyy(i,j))
            sigIIM(i,j) = sqrt( ( (sigyy(i,j) - sigxx(i,j))/2d0 )**2d0 &
                                     + sigxy(i,j)**2d0 )

            Fn(i,j)  = sigII(i,j)  + frict * sigI(i,j) - CoheC(i,j)


          endif
      enddo
    enddo

    if (peri .ne. 0) call periodicBC(sigxx_it,sigyy_it)
    if (peri .ne. 0) call periodicBC(sigII,sigxy_it)
    if (peri .ne. 0) call periodicBC(Fn,sigI)


!--------------------------------------------------------
!-------------COMPUTE THE DAMAGE PARAMETER at the center
!--------------------------------------------------------

    do i = 1, nx
       do j = 1, ny

        dfactor(i,j) = 1d0
        dfactorB(i,j) = 1d0

        if ( maskC(i,j) .eq. 1 .and. etaC(i,j) .gt. 1d1) then
             

             ! 1. compression capping
             if ( (sigI(i,j)-sigII(i,j)) .lt. sigcC(i,j)*(1+frict) ) then
                Dflag(i,j) = 1d0
                dfactor(i,j) = ( ( Deltat / Tdam)* &
                     ((sigcC(i,j)*(1+frict) / (sigI(i,j) - sigII(i,j)) ) - 1d0)+1d0)

 
             ! 2. tensile capping
             elseif ( (sigI(i,j)+sigII(i,j)) .gt. sigtC(i,j)) then
                Dflag(i,j) = 1d0
                dfactor(i,j) = ( ( Deltat / Tdam)* &
                        ((sigtC(i,j) / (sigI(i,j)+sigII(i,j))) - 1d0)+1d0)

             ! 3.  Mohr-Coulomb
             elseif ( Fn(i,j) .gt. 0d0 .and.  MEB_flow_rule .eq. 'none') then

               Dflag(i,j) = 1d0

               if (Dam_correction .eq. 'standard') then
		!Standard line to origin 

                   dfactor(i,j) = (( Deltat / Tdam)* &
                          ((CoheC(i,j) / (frict*sigI(i,j) + sigII(i,j))) - 1d0)+1d0)

                   R(i,j) = ((sigII(i,j)**2d0 + frict*frict*sigI(i,j)**2d0) / &
                             (sigII(i,j) + frict*sigI(i,j))**2d0)**5d-1
                   Rfactor = Rfactor + R(i,j)
                   Ndam = Ndam + 1
                   Rmax = max(Rmax,R(i,j))
                   Rabsfactor = Rabsfactor + R(i,j)*sigII(i,j)*dfactor(i,j)             
                   Rabsmax = max(Rabsmax,R(i,j)*sigII(i,j))*dfactor(i,j)

               elseif (Dam_correction .eq. 'specified') then

		 if (sigI(i,j) .lt. tan(theta_cor*deg2rad)*sigII(i,j)) then
		   !Stress correction scheme following theta_cor

                   dfactor(i,j) = (( Deltat / Tdam)* &
                          (((CoheC(i,j) + frict*(tan(theta_cor*deg2rad))*sigII(i,j) &
                            -frict*sigI(i,j)) / &
                            ((1+frict*tan(theta_cor*deg2rad))*sigII(i,j))) - 1d0)+1d0)

                   R(i,j) = (( (CoheC(i,j)- frict*sigI(i,j) )**2d0 + &
                                                 frict*frict*sigI(i,j)**2d0) / &
                       ( CoheC(i,j) + frict*tan(theta_cor*deg2rad)*sigII(i,j) &
                                                - frict*sigI(i,j))**2d0)**5d-1
                   Rfactor = Rfactor + R(i,j)
                   Ndam = Ndam + 1
                   Rmax = max(Rmax,R(i,j))
                   Rabsfactor = Rabsfactor + R(i,j)*sigII(i,j)*dfactor(i,j)             
                   Rabsmax = max(Rabsmax,R(i,j)*sigII(i,j))*dfactor(i,j)

		 else
                   !Standard line to origin

                   dfactor(i,j) = (( Deltat / Tdam)* &
                          ((CoheC(i,j) / (frict*sigI(i,j) + sigII(i,j))) - 1d0)+1d0)

                   R(i,j) = ((sigII(i,j)**2d0 + frict*frict*sigI(i,j)**2d0) / &
                             (sigII(i,j) + frict*sigI(i,j))**2d0)**5d-1
                   Rfactor = Rfactor + R(i,j)
                   Ndam = Ndam + 1
                   Rmax = max(Rmax,R(i,j))
                   Rabsfactor = Rabsfactor + R(i,j)*sigII(i,j)*dfactor(i,j)             
                   Rabsmax = max(Rabsmax,R(i,j)*sigII(i,j))*dfactor(i,j)

		 endif

               elseif (Dam_correction .eq. 'associated') then
		!not done yet
                  print *, 'Associated stress correction not coded yet'
                  stop
               else
                  print *, 'Wrong choice of stress correction'
                  stop
               endif

             ! 3.  Mohr-Coulomb with flow rule
             elseif ( Fn(i,j) .gt. 0d0 .and.  &
                              (MEB_flow_rule .eq. 'specified')) then
               Dflag(i,j) = 1d0 
               dfactor(i,j) = (( Deltat / Tdam)* &
                          ((sigIIM(i,j) / sigII(i,j)) - 1d0)+1d0)

             endif

             dfactor(i,j) = min(dfactor(i,j), 1d0)
             dfactor(i,j) = max(dfactor(i,j), 0d0)

        endif ! (if maskC = 1)

      enddo
    enddo

    if (peri .ne. 0) call periodicBC(dam,dfactor)

    do i = 1, nx+1
       do j = 1, ny+1

          m1 = maskC(i,j) 
          m2 = maskC(i-1,j)
          m3 = maskC(i,j-1)
          m4 = maskC(i-1,j-1)

          !open bc
          if (Periodic_x .eq. 0) then
             if (i .eq. 1) then
                m2 = 0d0
                m4 = 0d0
             elseif (i .eq. nx+1) then
                m1 = 0d0
                m3 = 0d0  
             endif
          endif                   
                   
          if (Periodic_y .eq. 0) then
             if (j .eq. 1) then
                m3 = 0d0
                m4 = 0d0
             elseif (j .eq. ny+1) then
                m1 = 0d0
                m2 = 0d0  
             endif                  
          endif 

          if (m1+m2+m3+m4 .ne. 0d0) then

            dfactorB(i,j) = (dam(i,j)*dfactor(i,j)*m1 + &
                             dam(i-1,j)*dfactor(i-1,j)*m2 + &
                             dam(i,j-1)*dfactor(i,j-1)*m3 + &
                             dam(i-1,j-1)*dfactor(i-1,j-1)*m4) &
                               / (damB(i,j)*(m1+m2+m3+m4))

!                    dfactorB(i,j) = (dfactor(i,j)*m1 + dfactor(i-1,j)*m2 + &
!                                         dfactor(i,j-1)*m3 + dfactor(i-1,j-1)*m4) &
!                                               / (m1+m2+m3+m4)

          else
          
            dfactorB(i,j) = 1d0
            
          endif



       enddo
    enddo
    
    if (peri .ne. 0) call periodicBC(damB,dfactorB)

    if ((k .eq. 1d0)) then
    ! this insures that we only update stress history outside IMEX

!------------------------------------------
!--------UPDATING THE DAMAGE AND HISTORY FIELDS
!------------------------------------------

        print *, "instability factor (final): ", Rmax, Rfactor, Ndam
        print *, "absolute error (final): ", Rabsmax, Rabsfactor
        Rfield = Rfield + R


        do i = 1, nx+1
          do j = 1, ny+1
            if ( maskC(i,j) .eq. 1 .and. Dflag(i,j) .eq. 1d0) then

 	      if ((sigxx_it(i,j)-sigI(i,j)) .eq. 0d0) then 
                theta(i,j) = pi/2d0
              else
                theta(i,j) = atan(sigxy_it(i,j)/ &
                            (sigxx_it(i,j)-sigI(i,j)))
              endif
                
              sigI_uncor = sigI(i,j)
              sigII_uncor = sigII(i,j)
              eI = (1/2d0)*( Dexx(i,j) + Deyy(i,j))
              eII = sqrt( ( (Deyy(i,j) - Dexx(i,j))/2d0 )**2d0 &
                                     + Dexy(i,j)**2d0 )

              if (MEB_flow_rule .eq. 'none') then

                if (Dam_correction .eq. 'standard') then !line to origin

                  sigII(i,j) = dfactor(i,j)*sigII(i,j)
                  sigI(i,j)  = dfactor(i,j)*sigI(i,j)

                elseif(Dam_correction .eq. 'specified') then !following theta_cor
 		  if (sigI(i,j) .lt. tan(theta_cor*deg2rad)*sigII(i,j)) then
                    sigI(i,j)  = sigI(i,j) - (sigII(i,j)*(1-dfactor(i,j)) &
                                             *tan(theta_cor*deg2rad))
                    sigII(i,j) = dfactor(i,j)*sigII(i,j)
		  else
                    sigII(i,j) = dfactor(i,j)*sigII(i,j)
                    sigI(i,j)  = dfactor(i,j)*sigI(i,j)
		  endif

                else
                  print *, 'wrong stress correction path chosen by user'
                  stop
                endif

              else
                  print *, 'wrong flow rule chosen by user'
                  stop
              endif




              sigxy(i,j) = dfactor(i,j)*sigxy_it(i,j)
              dam(i,j) = dam(i,j)*dfactor(i,j)

              if (sigxx_it(i,j) .ge. sigyy_it(i,j)) then
                  sigxx(i,j) = sigI(i,j) + sigII(i,j)*cos(theta(i,j))
                  sigyy(i,j) = sigI(i,j) - sigII(i,j)*cos(theta(i,j))
              else
                  sigxx(i,j) = sigI(i,j) - sigII(i,j)*cos(theta(i,j))
                  sigyy(i,j) = sigI(i,j) + sigII(i,j)*cos(theta(i,j))
              endif

!              if (dam(i,j) .gt. 9d-1 .and. A(i,j) .ge. 5d-1) then 

!                write(50,100) i*1d0, j*1d0, sigI_uncor, sigII_uncor, &
!                               sigI(i,j), sigII(i,j), eI, eII, &
!                               sigIM(i,j), sigIIM(i,j), &
!                               day*1d0, hour*1d0, minute*1d0, &
!                               second*1d0, milli*1d0 

!              endif

              sigxyB(i,j) = dfactorB(i,j)*( DsigxyB(i,j) &
                          + sigxyB(i,j)* GammaMEB_B(i,j) )
              damB(i,j) = damB(i,j)*dfactorB(i,j)
 
            else

              sigxx(i,j) = sigxx_it(i,j)
              sigyy(i,j) = sigyy_it(i,j)
              sigxy(i,j) = sigxy_it(i,j) 
              sigI(i,j) = sigI(i,j)
              sigII(i,j) = sigII(i,j)

              sigxyB(i,j) = dfactorB(i,j)*( DsigxyB(i,j) &
                          + sigxyB(i,j)* GammaMEB_B(i,j) )
              damB(i,j) = damB(i,j)*dfactorB(i,j)

            endif




          enddo
        enddo


        dfactor = 1d0
        dfactorB = 1d0

        do i = 1, nx+1
          do j = 1, ny+1
             dam(i,j) = dam(i,j) + Deltat*Theal
             damB(i,j) = damB(i,j) + Deltat*Theal
             dam(i,j) = min(dam(i,j), 1d0)
             damB(i,j) = min(damB(i,j), 1d0)
          enddo
        enddo


        if (peri .ne. 0) call periodicBC(sigxyB,sigxy)
        if (peri .ne. 0) call periodicBC(sigxx,sigyy)
        if (peri .ne. 0) call periodicBC(dam,damB)


        !------------------------------------------
        ! post the sea ice stress state (every 5 min)
        !------------------------------------------


        if (((milli .eq. 0) .and. (second .eq. 0))) then! &
!	    .and. ( hour .lt. 5 ) .and. (day .lt. 2) &
!            .and. ((minute .eq. 0) .or. (minute .eq. 10) &
!            .or. (minute .eq. 5) .or. (minute .eq. 15) &
!            .or. (minute .eq. 20) .or. (minute .eq. 30) &
!            .or. (minute .eq. 25) .or. (minute .eq. 35) &
!            .or. (minute .eq. 45) .or. (minute .eq. 55) &
!            .or. (minute .eq. 40) .or. (minute .eq. 50))) then

!          call post_MEB_stress(utp, vtp, sigI, sigII, Dexx, Deyy, Dexy, date, expno)

 !       elseif ((milli .eq. 0) .and. (second .eq. 0) &
 !           .and. ( minute .eq. 0 ) ) then

          call post_MEB_stress(utp, vtp, sigI, sigII, Dexx, Deyy, Dexy, date, expno)

        endif





    endif ! if k == 1

100            format (1x, 1000(f20.10, 1x))

      return
    end subroutine stress_strain_MEB



    subroutine post_MEB_stress(utp, vtp, sigI, sigII, exx, eyy, exy, date, expno)

    use datetime, only: datetime_type
    use ellipse
    use elastic

    implicit none

    include 'parameter.h'
    include 'CB_Dyndim.h'
    include 'CB_DynVariables.h'
    include 'CB_DynForcing.h'
    include 'CB_const.h'
    include 'CB_mask.h'
    include 'CB_options.h'

    type(datetime_type), intent(in) :: date

    character filename*100
    character outputpath*100
    integer, intent(in) :: expno
    integer i, j, m, year, month, day, hour, minute, second, milli

    double precision land, lowA

    double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

    double precision sigI(0:nx+2,0:ny+2), sigII(0:nx+2,0:ny+2)
    double precision exx(0:nx+2,0:ny+2), eyy(0:nx+2,0:ny+2), exy(0:nx+2,0:ny+2)

    year = date%year
    month = date%month
    day = date%day
    hour = date%hour
    minute = date%minute
    second = date%second
    milli = date%milli

    land=-999d0
    lowA=0d0


    write (filename,'("output/sigI",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (12, file = filename, status = 'unknown')

    write (filename,'("output/sigII",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (13, file = filename, status = 'unknown')

    write (filename,'("output/u",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (14, file = filename, status = 'unknown')

    write (filename,'("output/v",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (15, file = filename, status = 'unknown')

    write (filename,'("output/dam",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (16, file = filename, status = 'unknown')

    write (filename,'("output/exx",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (17, file = filename, status = 'unknown')

    write (filename,'("output/eyy",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (18, file = filename, status = 'unknown')

    write (filename,'("output/exy",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (19, file = filename, status = 'unknown')

    write (filename,'("output/sigxx",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (20, file = filename, status = 'unknown')

    write (filename,'("output/sigyy",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (21, file = filename, status = 'unknown')

    write (filename,'("output/sigxy",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (22, file = filename, status = 'unknown')

    write (filename,'("output/h",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (23, file = filename, status = 'unknown')

    write (filename,'("output/A",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (24, file = filename, status = 'unknown')
    write (filename,'("output/R",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
                year, month, day, hour, minute, milli, expno
    open (25, file = filename, status = 'unknown')

    do j = 0, ny+1
        write(12,100) ( sigI(i,j), i = 0, nx+1 )
        write(13,100) ( sigII(i,j), i = 0, nx+1 )
        write(14,200) ( utp(i,j), i = 0, nx+1 )
        write(15,200) ( vtp(i,j), i = 0, nx+1 )
        write(16,100) ( dam(i,j), i = 0, nx+1 )
        write(17,100) ( exx(i,j),   i = 0, nx+1 )
        write(18,100) ( eyy(i,j),   i = 0, nx+1 )
        write(19,100) ( exy(i,j),   i = 0, nx+1 )
        write(20,100) ( sigxx(i,j), i = 0, nx+1 )
        write(21,100) ( sigyy(i,j), i = 0, nx+1 )
        write(22,100) ( sigxy(i,j), i = 0, nx+1 )
        write(23,200) ( h(i,j), i = 0, nx+1 )
        write(24,200) ( A(i,j), i = 0, nx+1 )
        write(25,200) ( Rfield(i,j), i = 0, nx+1 )
    enddo

    do m=12,25
        close(m)
    enddo



100            format (1x, 1000(f20.10, 1x))
200            format (1x, 1000(f15.10, 1x))
300            format (1x, 1000(f20.4,  1x))

return
end subroutine post_MEB_stress



