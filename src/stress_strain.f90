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

                  sigI(i,j)   = -1d0*( dudx + dvdy )*zetaCf(i,j)+P(i,j)

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



