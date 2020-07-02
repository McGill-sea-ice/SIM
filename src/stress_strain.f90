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
      include 'CB_const_stressBC.h'
      include 'CB_mask.h'
      include 'CB_options.h'
      include 'CB_stressBC.h'

      type(datetime_type), intent(in) :: date

      character filename*80

      integer, intent(in) :: k, expno
      integer i, j, m, year, month, day, hour, minute, ncell, int11, int22, int12

      double precision dudx, dvdy, dudy, dvdx, land, lowA, pbc, maxdevp
      double precision maxp, minp, meanp, stdp, maxs, tempval, meanstrength

      double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision shear(0:nx+1,0:ny+1), div(0:nx+1,0:ny+1)
      double precision sig11(0:nx+1,0:ny+1), sig22(0:nx+1,0:ny+1) ! at t-point
      double precision sigI(0:nx+1,0:ny+1), sigII(0:nx+1,0:ny+1) ! at t-point
      double precision sigInorm(0:nx+1,0:ny+1), sigIInorm(0:nx+1,0:ny+1)! stress invariants at t-point
      double precision sig1norm(0:nx+1,0:ny+1), sig2norm(0:nx+1,0:ny+1) ! princ stresses at t-point
      double precision zetaCout(0:nx+1,0:ny+1), Pfout(0:nx+1,0:ny+1)

      double precision dsig11dx(0:nx+1,0:ny+1) ! at u-point
      double precision dsig22dy(0:nx+1,0:ny+1) ! at v-point

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute

      land=-999d0
      lowA=-888d0
      maxp=-1000000d0
      maxs=-1000000d0
      meanp=0d0
      meanstrength=0d0
      ncell=0
      minp=100000d0

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
      Pfout     = land

      sig11     = land
      sig22     = land
      dsig11dx  = land
      dsig22dy  = land

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
                  Pfout(i,j)    = Pf(i,j)

                  if (A(i,j) .lt. 0d0) then !changed value for stressBC
                     div(i,j)       = lowA
                     shear(i,j)     = lowA
                     sigI(i,j)      = lowA
                     sigII(i,j)     = lowA
                     sigInorm(i,j)  = lowA
                     sigIInorm(i,j) = lowA

                  elseif (A(i,j) .ge. 0d0) then
                  
                  dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                  dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax

                  if (stressBC) then

                     if (i .eq. 1) then ! 1st order

                        dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) )/2d0 -  &
                             ( vtp(i,j)   + vtp(i,j+1)   )/2d0 ) / Deltax

                     elseif (i .eq. nx) then

                        dvdx = ( ( vtp(i,j) + vtp(i,j+1) )/2d0 -  &
                             ( vtp(i-1,j)   + vtp(i-1,j+1)   )/2d0 ) / Deltax

                     else

                        dvdx = ( ( vtp(i+1,j) + vtp(i+1,j+1) ) -        &
                             ( vtp(i-1,j) + vtp(i-1,j+1) ) ) /      &
                             ( 4d0 * Deltax )
                        
                     endif


                  elseif (.not. stressBC) then
                  
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

                  endif

                  if (stressBC) then

                     if (j .eq. 1) then ! 1st order

                        dudy = ( ( utp(i,j+1) + utp(i+1,j+1) )/2d0 -  &
                             ( utp(i,j)   + utp(i+1,j)   )/2d0 ) / Deltax

                     elseif (j .eq. ny) then ! 1st order

                        dudy = ( ( utp(i,j) + utp(i+1,j) )/2d0 -  &
                             ( utp(i,j-1)   + utp(i+1,j-1)   )/2d0 ) / Deltax

                     else

                        dudy = ( ( utp(i,j+1) + utp(i+1,j+1) ) -        &
                             ( utp(i,j-1) + utp(i+1,j-1) ) ) /     &
                             ( 4d0 * Deltax )

                     endif

                     
                  elseif (.not. stressBC) then

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

                  endif
                  
!----- stresses and strain rates at the grid center -------------------------   


                  shear(i,j) = sqrt(( dudx - dvdy )**2d0 &  ! in day-1
                       + ( dudy + dvdx )**2d0 )*86400d0

                  div(i,j)   = ( dudx + dvdy )*86400d0      ! in day-1

! watchout p in our code is in fact p/2 in Hibler's equations

                  sigI(i,j)   = -1d0*( dudx + dvdy )*zetaCf(i,j)+Pf(i,j)
                  sig11(i,j) = zetaCf(i,j)*( dudx + dvdy ) + etaCf(i,j)*( dudx - dvdy ) -Pf(i,j)
                  sig22(i,j) = zetaCf(i,j)*( dudx + dvdy ) + etaCf(i,j)*( dvdy - dudx ) -Pf(i,j)

                  if (sigI(i,j) .gt. maxp) maxp=sigI(i,j)
                  if (sigI(i,j) .lt. minp) minp=sigI(i,j)
                  meanp=meanp+sigI(i,j)
                  meanstrength=meanstrength+Pp(i,j)
                  ncell=ncell+1
                  
                  sigII(i,j) = sqrt(( dudx - dvdy )**2d0 &
                       + ( dudy + dvdx )**2d0 )*etaCf(i,j)
                  
                  if (sigII(i,j) .gt. maxs) maxs=sigII(i,j)

                  sigInorm(i,j)  = sigI(i,j)  / (2d0*max(Pp(i,j),1d-10))
                  
                  sigIInorm(i,j) = sigII(i,j) / (2d0*max(Pp(i,j),1d-10))

                  sig1norm(i,j) = -1d0*sigInorm(i,j) + sigIInorm(i,j)
                  sig2norm(i,j) = -1d0*sigInorm(i,j) - sigIInorm(i,j)

                  endif

               endif
               
               
            enddo
         enddo

!----- calc dsig11dx at u point -----------------------------------------------
         
         do i = 2, nx
            do j = 1, ny
               dsig11dx(i,j) = ( sig11(i,j) - sig11(i-1,j) ) / Deltax
            enddo
         enddo
         
         if (stressBC) then
            do j=1,ny
               dsig11dx(1,j) = ( sig11(1,j) - 1000d0*sig11bc ) / Deltaxh
               dsig11dx(nx+1,j) = ( 1000d0*sig11bc - sig11(nx+1,j) ) / Deltaxh
            enddo
         endif

!----- calc dsig22dy at v point -----------------------------------------------                                                                     

         do i = 1, nx
            do j = 2, ny
               dsig22dy(i,j) = ( sig22(i,j) - sig22(i,j-1) ) / Deltax
            enddo
         enddo

         if (stressBC) then
            do i=1,nx
               dsig22dy(i,1) = ( sig22(i,1) - 1000d0*sig22bc ) / Deltaxh
               dsig22dy(i,ny+1) = ( 1000d0*sig22bc - sig22(i,ny+1) ) / Deltaxh
            enddo
         endif

!----- stats for pressure ------------------------------------------------------

         meanp=meanp/(ncell*1d0)
         meanstrength=2d0*meanstrength/(ncell*1d0)
         
         stdp=0d0
         do i = 1, nx
            do j = 1, ny
               if ( maskC(i,j) .eq. 1 ) then
                  if (A(i,j) .ge. 0.5d0) then
                     stdp=stdp+(sigI(i,j)-meanp)**2d0
                  endif
               endif
            enddo
         enddo
         stdp=sqrt(stdp/((ncell-1)*1d0))

         print *, '0.995 shear strength = ', 0.995d0*Pp(2,2)/sqrt(ell2)
         print *, 'mean, max, min, var = ', meanp, maxp, minp, stdp
         print *, 'max shear = ', maxs
         print *, 'mean strength = ',  meanstrength
   
!----- verify pressure at boundaries -------------------------------------------

         pbc=( -1d0*sigmaW(ny-1) - 1d0*sigmaN(nx-1) )/2d0
         maxdevp=0d0
         do i = 1, nx
            if (abs(sigI(i,1)-pbc) .gt. maxdevp) maxdevp=abs(sigI(i,1)-pbc)
            if (abs(sigI(i,ny)-pbc) .gt. maxdevp) maxdevp=abs(sigI(i,ny)-pbc)
         enddo

         do j = 1, ny
            if (abs(sigI(1,j)-pbc) .gt. maxdevp) maxdevp=abs(sigI(1,j)-pbc)
            if (abs(sigI(nx,j)-pbc) .gt. maxdevp) maxdevp=abs(sigI(nx,j)-pbc)
         enddo

         print *, 'max p deviation (%)', maxdevp*100d0/pbc

         if (maxdevp*100d0/pbc .gt. 10d0) print *, 'WARNING p deviation!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!-------------------------------------------------------------------------------
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
               zetaCout(i,j)  = lowA
               Pfout(i,j)     = lowA
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
               zetaCout(i,j)  = lowA
               Pfout(i,j)     = lowA
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
               zetaCout(i,j)  = lowA
               Pfout(i,j)     = lowA
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
               zetaCout(i,j)  = lowA
               Pfout(i,j)     = lowA
            endif
         enddo

         tempval=abs(sig11bc)
         int11=int(tempval)
         tempval=abs(sig22bc)
         int22=int(tempval)
         tempval=abs(sig12bc)
         int12=int(tempval)

      write (filename,'("output/sigI","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (12, file = filename, status = 'unknown')

      write (filename,'("output/sigII","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (13, file = filename, status = 'unknown')

      write (filename,'("output/sigInorm","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (14, file = filename, status = 'unknown')

      write (filename,'("output/sigIInorm","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (15, file = filename, status = 'unknown')

      write (filename,'("output/sig1norm","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (16, file = filename, status = 'unknown')

      write (filename,'("output/sig2norm","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (17, file = filename, status = 'unknown')

      write (filename,'("output/div","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (18, file = filename, status = 'unknown')

      write (filename,'("output/shear","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (19, file = filename, status = 'unknown')

      write (filename,'("output/zetaC","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (20, file = filename, status = 'unknown')

      write (filename,'("output/Prep","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (21, file = filename, status = 'unknown')

      write (filename,'("output/dsig11dx","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (22, file = filename, status = 'unknown')

      write (filename,'("output/dsig22dy","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (23, file = filename, status = 'unknown')

      write (filename,'("output/sig11","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (24, file = filename, status = 'unknown')

      write (filename,'("output/sig22","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (25, file = filename, status = 'unknown')

      write (filename,'("output/pinfo","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (32, file = filename, status = 'unknown')

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
         write(21,300) ( Pfout(i,j), i = 0, nx+1 )
         write(22,100) ( dsig11dx(i,j), i = 0, nx+1 )
         write(23,100) ( dsig22dy(i,j), i = 0, nx+1 )
         write(24,100) ( sig11(i,j), i = 0, nx+1 )
         write(25,100) ( sig22(i,j), i = 0, nx+1 )      
      enddo

      do m=12,25
         close(m)
      enddo
      
      write(32,400) nx, meanp, maxp, stdp, maxs
      close(32)

100            format (1x, 1000(f20.10, 1x))
200            format (1x, 1000(f15.10, 1x))
300            format (1x, 1000(f20.4,  1x))
400            format (i4, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10)      
      return
    end subroutine stress_strain


    subroutine stress_strainB(utp, vtp, expno)
      use ellipse
      implicit none

      include 'parameter.h'
      include 'CB_Dyndim.h'
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_const_stressBC.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      character filename*80

      integer, intent(in) :: expno
      integer i, j, m, int11, int22, int12

      double precision dudy, dvdx, land, lowA, tempval

      double precision, intent(in):: utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      double precision sig12(0:nx+1,0:ny+1), eps12(0:nx+1,0:ny+1)
      double precision dsig12dx(0:nx+1,0:ny+1), dsig12dy(0:nx+1,0:ny+1)

      land=-999d0
      lowA=-888d0

      eps12     = land
      sig12     = land
      dsig12dx  = land
      dsig12dy  = land

      do i = 2, nx
         do j = 2, ny

            dudy = ( utp(i,j) - utp(i,j-1) ) / Deltax
            dvdx = ( vtp(i,j) - vtp(i-1,j) ) / Deltax

            eps12(i,j)=(dudy+dvdx)*86400d0/2d0     ! in day-1 
            sig12(i,j)=etaBf(i,j)*(dudy+dvdx)

         enddo
      enddo

!----- dsig12dy at u point ---------------

      do i = 2, nx
         do j = 2, ny-1
            dsig12dy(i,j) = ( sig12(i,j+1) - sig12(i,j) ) / Deltax
         enddo
      enddo
      if (stressBC) then
         do j = 1, ny
            dsig12dy(1,j) = 0d0 ! sig12bc - sig12bc (constant)
            dsig12dy(nx+1,j) = 0d0
         enddo
         do i = 2, nx
            dsig12dy(i,1) = ( sig12(i,2) - 1000d0*sig12bc ) / Deltax
            dsig12dy(i,ny) = ( 1000d0*sig12bc - sig12(i,ny) ) / Deltax
         enddo
      endif

!----- dsig12dx at v point ---------------                                                                                      

      do i = 2, nx-1
         do j = 2, ny
            dsig12dx(i,j) = ( sig12(i+1,j) - sig12(i,j) ) / Deltax
         enddo
      enddo
      if (stressBC) then
         do i = 1, nx
            dsig12dx(i,1) = 0d0 ! sig12bc - sig12bc (constant)
            dsig12dx(i,ny+1) = 0d0
         enddo
         do j = 2, ny
            dsig12dx(1,j) = ( sig12(2,j) - 1000d0*sig12bc ) / Deltax
            dsig12dx(nx,j) = ( 1000d0*sig12bc - sig12(nx,j) ) / Deltax
         enddo
      endif

!-------------------------------------

      tempval=abs(sig11bc)
      int11=int(tempval)
      tempval=abs(sig22bc)
      int22=int(tempval)
      tempval=abs(sig12bc)
      int12=int(tempval)

      write (filename,'("output/sig12","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (15, file = filename, status = 'unknown')

      write (filename,'("output/eps12","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (16, file = filename, status = 'unknown')

      write (filename,'("output/dsig12dy","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (17, file = filename, status = 'unknown')

      write (filename,'("output/dsig12dx","_nx",i4.4,"_sbc11_",i2.2,"_sbc22_",i2.2,"_sbc12_",i2.2,".",i2.2)') &
           nx, int11, int22, int12, expno
      open (18, file = filename, status = 'unknown')

      do j = 0, ny+1
         write(15,100) ( sig12(i,j), i = 0, nx+1 )
         write(16,200) ( eps12(i,j), i = 0, nx+1 )
         write(17,200) ( dsig12dy(i,j), i = 0, nx+1 )
         write(18,200) ( dsig12dx(i,j), i = 0, nx+1 )      
      enddo

      do m=15,18
         close(m)
      enddo

100   format (1x, 1000(f20.10, 1x))
200   format (1x, 1000(f15.10, 1x))


    end subroutine stress_strainB
