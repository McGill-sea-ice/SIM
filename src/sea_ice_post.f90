!************************************************************************
!     subroutine sea_ice _post:
!       Outouts the results.
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


      subroutine sea_ice_post (date, expno)
        use datetime, only: datetime_type
      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_DynVariables.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_DynForcing.h'
      include 'CB_buoys.h'

      type(datetime_type), intent(in) :: date
      
      integer, intent(in) :: expno
      integer i, j, k, kk, year, month, day, hour, minute

      character filename*32

      double precision hactual(0:nx+1,0:ny+1)

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute
     
      write (filename,'("output/h",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
           year, month, day, hour, minute, expno
      open (10, file = filename, status = 'unknown')

      write (filename,'("output/A",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
           year, month, day, hour, minute, expno
      open (11, file = filename, status = 'unknown')
         
      write (filename,'("output/dam",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
           year, month, day, hour, minute, expno
      open (40, file = filename, status = 'unknown')

      if ( Dynamic ) then

         write (filename,'("output/u",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (12, file = filename, status = 'unknown')
         
         write (filename,'("output/v",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (13, file = filename, status = 'unknown')

      endif

      if ( Thermodyn ) then

         write (filename,'("output/Ta",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (20, file = filename, status = 'unknown')

         write (filename,'("output/To",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (21, file = filename, status = 'unknown')

         write (filename,'("output/Ti",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (22, file = filename, status = 'unknown')

         write (filename,'("output/Pvap",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (24, file = filename, status = 'unknown')

         write (filename,'("output/Qsh_io",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (31, file = filename, status = 'unknown')

         write (filename,'("output/Qoa",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (33, file = filename, status = 'unknown')

      endif

      if ( BuoyTrack) then

         write (filename,'("output/Tracer",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,".",i2.2)') &
              year, month, day, hour, minute, expno
         open (35, file = filename, status = 'unknown')

      endif

      do j = 1, ny
         do i = 1, nx

            if (A(i,j) .lt. 0.0000001) then
               hactual(i,j) = 0d0
            else
               hactual(i,j) = h(i,j) / A(i,j)
            endif

         enddo
      enddo

      do j = 0, ny+1
         write(10,10) ( h(i,j),       i = 0, nx+1 )
         write(11,10) ( A(i,j),       i = 0, nx+1 )
      enddo

      do j = 0, ny + 1
         write(40,10) ( dam(i,j), i = 0, nx+1 )
      enddo

      if ( Dynamic ) then

         do j = 1, ny
            write(12,20) ( uice(i,j), i = 1, nx+1 )
!            write(17,20) ( uwatnd(i,j), i = 1, nx+1 )
         enddo

         do j = 0, ny+1
!            write(15,30) ( etaC(i,j), i = 0, nx+1 )
!            write(38,30) ( zetaC(i,j), i = 0, nx+1 )
         enddo

         do j = 1, ny+1
            write(13,10) ( vice(i,j), i = 1, nx )
!            write(18,10) ( vwatnd(i,j), i = 1, nx )
         enddo

      endif

      if ( ThermoDyn ) then

         do j = 0, ny+1
            write(20,30) ( Ta(i,j), i = 0, nx+1 )
            write(21,30) ( To(i,j), i = 0, nx+1 )
            write(24,30) ( Pvap(i,j), i = 0, nx+1 )
!            write(25,30) ( Qlwup(i,j), i = 0, nx+1 )
!            write(26,30) ( Qlwdo(i,j), i = 0, nx+1 )
!            write(27,30) ( Qlh(i,j), i = 0, nx+1 )
!            write(28,30) ( Qsh(i,j), i = 0, nx+1 )
!            write(29,30) ( Qsw(i,j), i = 0, nx+1 )
!            write(30,30) ( Qia(i,j), i = 0, nx+1 )
            write(31,30) ( Qsh_io(i,j), i = 0, nx+1 )
!            write(32,30) ( Qoa_f(i,j), i = 0, nx+1 )
            write(33,30) ( Qoa(i,j), i = 0, nx+1 )
!            write(34,30) ( Qadvdiff(i,j), i = 0, nx+1 )
         enddo

         do j = 1, ny
            write(22,30) (Ti(i,j), i = 1, nx)
!           write(23,30) (Tl(i,j), i = 1, nx)
         enddo

      endif


      if ( BuoyTrack ) then
         do kk = 1, nbuoys
!            ntime = int((time-0d0) / record_time_buoy)
!            write(35,80) (BuoyPositions(kk,kt,1), kt = 1, ntime)
!            write(35,80) (BuoyPositions(kk,kt,2), kt = 1, ntime)
!            write(35,80) (BuoyPositions(kk,kt,3), kt = 1, ntime)
!            write(35,80) (BuoyPositions(kk,kt,4), kt = 1, ntime)
         enddo
      endif


      do k = 9, 40
         close(k)
      enddo

 10   format (1x, 1000(f20.16, 1x))
 20   format (1x, 1000(f20.16, 1x))
 30   format (1x, 1000(f12.6, 1x))
 !40   format (1x, 1000(f12.6, 1x))
 !50   format (1x, 1000(f14.6, 1x))
 !80   format (1x, 1000(f12.6, 1x))
 !100  format (1x, 1000(e12.4, 1x))
      return
    end subroutine sea_ice_post

subroutine info_file (expno)

  use ellipse
  use numerical_VP
  use numerical_EVP
  use solver_choice

  implicit none

  include 'parameter.h'
  include 'CB_const.h'
  include 'CB_options.h'

  character filename*30

  integer, intent(in) :: expno

  write (filename, '("output/info.",i2.2)') expno
  open (10, file = filename, status = 'unknown')

  write(10,*) ('PHYSICAL PARAMETERS')
  write(10,*) ('')
  write(10,*) ('C =')
  write(10,10) (C)
  write(10,*) ('Pstar =')
  write(10,10) (Pstar)
  write(10,*) ('e =') 
  write(10,10) (sqrt(ell2))

  write(10,*) ('')
  write(10,*) ('RUN SETUP')
  write(10,*) ('')
  write(10,*) ('Deltax (m) =')
  write(10,10) (Deltax)
  write(10,*) ('Deltat (s) =')
  write(10,10) (Deltat)
  write(10,*) ('Damage = ')
  write(10,*) (Damage)
  write(10,*) ('Dynamic =')
  write(10,*) (Dynamic)
  write(10,*) ('Thermodynamic =')
  write(10,*) (Thermodyn)
  write(10,*) ('Current =')
  write(10,*) (Current)
  write(10,*) ('Wind =')
  write(10,*) (Wind)
  write(10,*) ('Method to calculate the viscous coefficients')
  write(10,*) (visc_method)
  
  write(10,*) ('')
  write(10,*) ('NUMERICAL')
  write(10,*) ('')
  if (solver .eq. 1) then
     write(10,*) ('Picard solver')
     write(10,*) ('')
     write(10,*) ('max nb of Outer loops')
     write(10,*) (OLmax)
     write(10,*) ('nonlinear convergence criterion for Picard solver')
     write(10,10) (gamma_nl)
  elseif (solver .eq. 2) then
     write(10,*) ('JFNK solver')
     write(10,*) ('')
     write(10,*) ('max nb of Newton loops')
     write(10,*) (NLmax)
     write(10,*) ('nonlinear convergence criterion for JFNK')
     write(10,10) (gamma_nl)
  elseif (solver .eq. 3) then
     write(10,*) ('nb of EVP subcycles')
     write(10,*) (Nsub)
     write(10,*) ('Eo EVP tuning parameter')
     write(10,10) (Eo)
  endif

10 format (f20.10)

end subroutine info_file
