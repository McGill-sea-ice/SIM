!***********************************************************************
!     Subroutine ini_get: set the initial conditions
!***********************************************************************
subroutine ini_get (restart, expno_r, restart_date)
    use datetime, only: datetime_type
    use ellipse
      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_DynVariables.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_mask.h'
      include 'CB_mask_boat.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_const_stressBC.h'
      include 'CB_buoys.h'
      include 'CB_Dyndim.h'

      character(LEN=32) filename, crack_orientation  ! restart file name

      integer, intent(in) :: restart, expno_r
      integer :: i, j, k, year, month, day, hour, minute, ncc, ncount
      integer :: ii, jj, nr, icrack, jcrack, nwidth, nlength
      integer :: i1, i2, i3, i4, j1, j2, j3, j4, ncw, ncl
      integer :: jbeg(nxh), jend(nxh)
      TYPE(datetime_type), intent(in) :: restart_date

      double precision :: rdnb, hfactor, hmean, hminimum, sumh, sumA, sumb
      double precision :: hh(1:nxh,1:nyh), Ah(1:nxh,1:nyh)

      double precision :: r8_normal_ab, gmean, gstd
      external r8_normal_ab

      year = restart_date%year
      month = restart_date%month
      day = restart_date%day
      hour = restart_date%hour
      minute = restart_date%minute

!------------------------------------------------------------------------
!     Set initial conditions for h,A,u,v,Ta,Ti,Tl
!------------------------------------------------------------------------

      if ( restart .eq. 0 .and. .not. stressBC) then

         do i = 0, nx+1
            do j = 0, ny+1

!
!     already zeros on edges

               h(i,j)   =  1d0 * maskC(i,j) ! initial ice thick
               A(i,j)   =  1d0 * maskC(i,j) ! initial ice conc

!     h and A set to zero on open boundaries

               if (i.eq.0 .or. i.eq.nx+1 .or. j.eq.0 .or. j.eq.ny+1) &
                         h(i,j) = 0d0
               if (i.eq.0 .or. i.eq.nx+1 .or. j.eq.0 .or. j.eq.ny+1) &
                         A(i,j) = 0d0

               tracer(i,j,1) = h(i,j)
               tracer(i,j,2) = A(i,j)

               Pp(i,j) = 0d0 
               P(i,j)  = 0d0

               etaC(i,j)= 0d0
               zetaC(i,j) = 0d0
               etaB(i,j)  = 0d0

               Ta(i,j)  =  273.15d0              ! air temp

               To(i,j)  =  Tof * maskC(i,j)      ! ocean temp
               Ti(i,j)  =  273.15d0 * maskC(i,j)    ! jfl ice surface temp
               Tl(i,j)  =  1d0                   ! land surface temp
               
            enddo
         enddo

         uice = 0d0 
         vice = 0d0
         un1  = 0d0  ! u-velocity pts 
         vn1  = 0d0  ! v-velocity pts 

      endif

!------------------------------------------------------------------------
!     Load restart files for h,A,u,v,Ta,Ti,Tl initial conditions
!------------------------------------------------------------------------

      if ( restart .eq. 1 .and. .not. stressBC) then ! load restart files

!------------------------------------------------------------------------
!     Open file and assign unit to it 
!------------------------------------------------------------------------
         
         if ( Dynamic ) then

            write (filename,'("output/h",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (16, file = filename, status = 'old')
            
            write (filename,'("output/A",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (17, file = filename, status = 'old')
            

            write (filename,'("output/u",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (18, file = filename, status = 'old')
  
            write (filename,'("output/v",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (19, file = filename, status = 'unknown')
            
            do j = 1, ny
               read (18,*) ( uice(i,j), i = 1, nx+1 )
            enddo
            
            do j = 1, ny+1
               read (19,*) ( vice(i,j), i = 1, nx )
            enddo
            
            do j = 0, ny+1

               read (16,*) ( h(i,j),    i = 0, nx+1)
               read (17,*) ( A(i,j),    i = 0, nx+1)

               do i = 0, nx+1
                  tracer(i,j,1) = h(i,j)
                  tracer(i,j,2) = A(i,j)
               enddo

            enddo

            do k = 16, 19
               close(k)
            enddo

            un1 = uice
            vn1 = vice

         endif
         
	 Cbasal1 = 0d0
	 Cbasal2 = 0d0
         
         if ( Thermodyn ) then

            write (filename,'("output/Ta",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (21, file = filename, status = 'unknown')

            write (filename,'("output/To",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r            
            open (22, file = filename, status = 'unknown')

            write (filename,'("output/Ti",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (24, file = filename, status = 'unknown')

            write (filename,'("output/Qoa",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (73, file = filename, status = 'unknown')

            write (filename,'("output/Qsh_io",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (74, file = filename, status = 'unknown')

            write (filename,'("output/Pvap",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
                 year, month, day, hour, minute, expno_r
            open (75, file = filename, status = 'unknown')
            
            do j = 1, ny
               read (24,*) ( Ti(i,j),   i = 1, nx )
            enddo
            

            do j = 0, ny+1
               read (21,*) ( Ta(i,j), i = 0, nx+1 )
               read (22,*) ( To(i,j), i = 0, nx+1 )
               read (73,*) ( Qoa(i,j), i = 0, nx+1 )
               read (74,*) ( Qsh_io(i,j), i = 0, nx+1 )
               read (75,*) ( Pvap(i,j), i = 0, nx+1 )
            enddo
            
            do k = 21, 25
               close(k)
            enddo

            do k = 73, 75
               close(k)
            enddo

         endif


      endif

      mboat=0 ! set mask to zero everywhere                                                                         
      msideboat=0

      if (stressBC) then

         if (initialh .eq. 'constant') then

            do ii = 1, nxh
               do jj = 1, nyh
                  hh(ii,jj)   =  hlevel ! initial ice thick        
                  Ah(ii,jj)   =  1d0    ! initial ice conc
               enddo
            enddo

         elseif (initialh .eq. 'icebridge') then

            open (62, file = 'input/IDCSI4_20190406_a.dat', status = 'old')
            do j = 1, nyh
               read (62,*) ( hh(i,j),   i = 1, nxh )
            enddo
            close(62)

            do ii = 1, nxh
               do jj = 1, nyh
                  Ah(ii,jj)   =  1d0    ! initial ice conc        
               enddo
            enddo


         elseif (initialh .eq. 'whitenoise') then

            hfactor=1.2d0
            hmean=1d0

            do i = 1,7948

               call random_number(rdnb)

            enddo

            do ii = 1, nxh
               do jj = 1, nyh
                  
                  call random_number(rdnb) ! 0 to 1

                  hh(ii,jj)   =  ((0.5d0-rdnb)*hfactor + hmean) ! initial ice thick  
                  Ah(ii,jj)   =  1d0                            ! initial ice conc              
                  
               enddo
            enddo

         elseif (initialh .eq. 'normalDist') then

            gmean=1d0
            gstd=0.25d0
            hminimum=10000d0
            ncount=0

            do ii = 1, nxh ! define high res field (20 m res)                                                                      
               do jj = 1, nyh
                  Ah(ii,jj)=1d0
                  hh(ii,jj)=r8_normal_ab (gmean, gstd, iseed )
                  hh(ii,jj)=max(hh(ii,jj), 0d0)
                  if (hh(ii,jj) .eq. 0d0) then
                     ncount=ncount+1
                     print *, ncount, nxh*nyh, 'thickness is zero'
                     Ah(ii,jj)=0d0
                  endif
                  if (hh(ii,jj) .lt. hminimum) hminimum=hh(ii,jj)
               enddo
            enddo
            print *, 'hmin=', hminimum

         elseif (initialh .eq. 'crack') then ! zzz

            do ii = 1, nxh
               do jj = 1, nyh
                  hh(ii,jj)   =  hlevel ! initial ice thick  
                  Ah(ii,jj)   =  1d0    ! initial ice conc
               enddo
            enddo

            crack_orientation='horizontal' ! horizontal, vertical, diagonal
            nlength=100 ! yop
            nwidth=1
            icrack=(nxh-nlength)/2
            jcrack=nyh/2 - 1

            if (crack_orientation .eq. 'horizontal') then
               do ii=1,nlength
                  do jj=1,nwidth
                     hh(icrack+ii,jcrack+jj)=0d0
                     Ah(icrack+ii,jcrack+jj)=0d0
                  enddo
               enddo

            elseif (crack_orientation .eq. 'vertical') then
               do jj=1,nlength
                  do ii=1,nwidth
                     hh(icrack+ii,jcrack+jj)=0d0
                     Ah(icrack+ii,jcrack+jj)=0d0
                  enddo
               enddo

            elseif (crack_orientation .eq. 'diagonal') then ! voir EC-5 p.60
               ncl=10
               ncw=7
               i2=71
               j2=71
               i1=i2-ncw
               j1=j2+ncw
               i3=i1+ncl
               j3=j1+ncl
               i4=i2+ncl
               j4=j2+ncl
               jbeg=1
               jend=1
               
               do ii=i1,i4
                  if (ii .le. i2) jbeg(ii)=(i1-ii)+j1
                  if (ii .gt. i2) jbeg(ii)=(ii-i2)+j2
                  if (ii .le. i3) jend(ii)=(ii-i1)+j1
                  if (ii .gt. i3) jend(ii)=(i3-ii)+j3
               enddo

               do ii=i1,i4
                  do jj=jbeg(ii),jend(ii)
                     hh(ii,jj)=0d0
                     Ah(ii,jj)=0d0
                  enddo
               enddo

            endif
         
            else

               print *, 'wrong setting of initialh'
               stop
      
         endif

!---------- coarsegraining ------------------------------------------                                                                   

            h=hlevel
            A=1d0

            nr=nxh/nx ! ratio of nxh and nx                                                                                      

            do i = 1, nx
               do j = 1, ny
                  sumh=0d0
                  sumA=0d0
                  ncount=0
                  do ii = (i-1)*nr+1, i*nr
                     do jj = (j-1)*nr+1, j*nr
                        ncount=ncount+1
                        sumh=sumh+hh(ii,jj)
                        sumA=sumA+Ah(ii,jj)
                     enddo
                  enddo
                  h(i,j)=sumh/(ncount*1d0)
                  A(i,j)=sumA/(ncount*1d0)
               enddo
            enddo

!---------- end of coarsegraining -----------------------------------


!------- DEFINE boat shape and dimensions ---------------------------
         if (boat) then

            if (nx .ne. nxh) then
               print *, 'nx should be set to nxh when boat is true'
               stop
            endif
            ! boat pixels...                                                

            mboat(71,80)=1
            do jj = 1, 6
               mboat(70,80+jj)=1 ! C-point        
               mboat(71,80+jj)=1
               mboat(72,80+jj)=1
            enddo
            mboat(71,87)=1
!            mboat=0 ! set mask to zero everywhere        
!            msideboat=0

            do j = 1, nyh+1
               do i = 1, nxh+1

                  sumb = mboat(i,j)   + mboat(i-1,j) + mboat(i,j-1) + mboat(i-1,j-1)
                  if (sumb .eq. 0) then
                     msideboat(i,j)=0
                  elseif (sumb .gt. 0 .and. sumb .lt. 4) then
                     msideboat(i,j)=1
                  elseif (sumb .eq. 4) then
                     msideboat(i,j)=0
                  endif
                  if (msideboat(i,j) .eq. 1) print *, 'side of boat', j,i
               enddo
            enddo

         endif
!--------------------------------------------------------------------  

            uice = 0d0
            vice = 0d0
            un1  = 0d0  ! u-velocity pts      
            vn1  = 0d0  ! v-velocity pts 

      endif

    end subroutine ini_get


!***********************************************************************     
!     Subroutine get_BC_stress: set B! stress when stressB! = .true. 
!***********************************************************************     
                                     
 
    subroutine get_BC_stress
      
      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_const.h'
      include 'CB_const_stressBC.h'
      include 'CB_stressBC.h'
      
      integer :: i, j

      print *, sig11bc, sig22bc, sig12bc

      do i = 1, nx
         sigmaS(i)= sig22bc*1d03 ! from N/m to kN/m
         sigmaN(i)= sig22bc*1d03
      enddo

      do j = 1, ny
         sigmaW(j)= sig11bc*1d03
         sigmaE(j)= sig11bc*1d03
      enddo

      do i = 1, nx+1
         tauS(i)=sig12bc*1d03
         tauN(i)=sig12bc*1d03
      enddo

      do j = 1, ny+1
         tauW(j)=sig12bc*1d03
         tauE(j)=sig12bc*1d03
      enddo      

    end subroutine get_BC_stress

!*****************************************************************************
! library of fortran function to get random number from Gaussian dist
! Author: John Burkardt
! https://people.sc.fsu.edu/~jburkardt/f_src/normal/normal.f90
!*****************************************************************************

function c4_normal_01 ( seed )

!*****************************************************************************80
!
!! C4_NORMAL_01 returns a unit pseudonormal C4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, complex ( kind = 4 ) C4_NORMAL_01, a unit pseudonormal value.
!
  implicit none

  complex ( kind = 4 ) c4_normal_01
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2
  real ( kind = 4 ) x_c
  real ( kind = 4 ) x_r

  v1 = r4_uniform_01 ( seed )
  v2 = r4_uniform_01 ( seed )

  x_r = sqrt ( - 2.0E+00 * log ( v1 ) ) * cos ( 2.0E+00 * r4_pi * v2 )
  x_c = sqrt ( - 2.0E+00 * log ( v1 ) ) * sin ( 2.0E+00 * r4_pi * v2 )

  c4_normal_01 = cmplx ( x_r, x_c )

  return
end
function c8_normal_01 ( seed )

!*****************************************************************************80
!
!! C8_NORMAL_01 returns a unit pseudonormal C8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, complex ( kind = 8 ) C8_NORMAL_01, a sample of the PDF.
!
!    Output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
  implicit none

  complex ( kind = 8 ) c8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) x_c
  real ( kind = 8 ) x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end
function i4_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! I4_NORMAL_AB returns a scaled pseudonormal I4.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!    The result is then rounded to the nearest integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the mean of the PDF.
!
!    Input, real ( kind = 4 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the
!    random number generator.
!
!    Output, integer ( kind = 4 ) I4_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) i4_normal_ab
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  i4_normal_ab = nint ( a + b * x )

  return
end
function i8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! I8_NORMAL_AB returns a scaled pseudonormal I8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!    The result is then rounded to the nearest integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the
!    random number generator.
!
!    Output, integer ( kind = 8 ) I8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 8 ) i8_normal_ab
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  i8_normal_ab = nint ( a + b * x )

  return
end
function r4_normal_01 ( seed )

!*****************************************************************************80
!
!! R4_NORMAL_01 returns a unit pseudonormal R4.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) R4_NORMAL_01, a sample of the standard
!    normal PDF.
!
  implicit none

  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r4_normal_01
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  r4_normal_01 = x

  return
end
function r4_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R4_NORMAL_AB returns a scaled pseudonormal R4.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the mean of the PDF.
!
!    Input, real ( kind = 4 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) R4_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r4_normal_ab
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x

  r1 = r4_uniform_01 ( seed )
  r2 = r4_uniform_01 ( seed )
  x = sqrt ( - 2.0E+00 * log ( r1 ) ) * cos ( 2.0E+00 * r4_pi * r2 )

  r4_normal_ab = a + b * x

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r4_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r4vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is an array of real ( kind = 4 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 4 ) * 4.656612875E-10

  end do

  return
end
subroutine r4vec_normal_ab ( n, a, b, seed, x )

!*****************************************************************************80
!
!! R4VEC_NORMAL_AB returns a scaled pseudonormal R4VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    An R4VEC is a vector of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input, real ( kind = 4 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 4 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 4 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) m
  real ( kind = 4 ) r(n+1)
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  real ( kind = 4 ) r4_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r4_uniform_01 ( seed )

    if ( r(1) == 0.0E+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4VEC_NORMAL_AB - Fatal error!'
      write ( *, '(a)' ) '  R4_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r4_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0E+00 * log ( r(1) ) ) * cos ( 2.0E+00 * r4_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r4vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r4vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0E+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0E+00 * r4_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0E+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0E+00 * r4_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end
function r8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, the mean of the PDF.
!
!    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 8 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_normal_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_01 ( m * n, seed, r )

  return
end
subroutine r8mat_normal_ab ( m, n, a, b, seed, r )

!*****************************************************************************80
!
!! R8MAT_NORMAL_AB returns a scaled pseudonormal R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the array.
!
!    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudonormal values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  call r8vec_normal_ab ( m * n, a, b, seed, r )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_normal_01 ( n, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) m
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  return
end
subroutine r8vec_normal_ab ( n, a, b, seed, x )

!*****************************************************************************80
!
!! R8VEC_NORMAL_AB returns a scaled pseudonormal R8VEC.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values desired.
!
!    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
!    random values.  Its dimension is N+1, but really it is only needed
!    to be the smallest even number greater than or equal to N.
!
!    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
!    of entries of X that we need to compute. 
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) m
  real ( kind = 8 ) r(n+1)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) x_hi_index
  integer ( kind = 4 ) x_lo_index
!
!  Record the range of X we need to fill in.
!
  x_lo_index = 1
  x_hi_index = n
!
!  If we need just one new value, do that here to avoid null arrays.
!
  if ( x_hi_index - x_lo_index + 1 == 1 ) then

    r(1) = r8_uniform_01 ( seed )

    if ( r(1) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC_NORMAL_AB - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop 1
    end if

    r(2) = r8_uniform_01 ( seed )

    x(x_hi_index) = &
      sqrt ( - 2.0D+00 * log ( r(1) ) ) * cos ( 2.0D+00 * r8_pi * r(2) )
!
!  If we require an even number of values, that's easy.
!
  else if ( mod ( x_hi_index - x_lo_index, 2 ) == 1 ) then

    m = ( x_hi_index - x_lo_index + 1 ) / 2

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-1:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m:2) )
!
!  If we require an odd number of values, we generate an even number,
!  and handle the last pair specially, storing one in X(N), and
!  saving the other for later.
!
  else

    x_hi_index = x_hi_index - 1

    m = ( x_hi_index - x_lo_index + 1 ) / 2 + 1

    call r8vec_uniform_01 ( 2*m, seed, r )

    x(x_lo_index:x_hi_index-1:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(x_lo_index+1:x_hi_index:2) = &
      sqrt ( - 2.0D+00 * log ( r(1:2*m-3:2) ) ) &
      * sin ( 2.0D+00 * r8_pi * r(2:2*m-2:2) )

    x(n) = sqrt ( - 2.0D+00 * log ( r(2*m-1) ) ) &
      * cos ( 2.0D+00 * r8_pi * r(2*m) )

  end if

  x(1:n) = a + b * x(1:n)

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
