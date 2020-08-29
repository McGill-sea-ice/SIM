!***********************************************************************
!     Subroutine ini_get: set the initial conditions
!***********************************************************************
subroutine ini_get (restart, expno_r, restart_date)
    use datetime, only: datetime_type
      implicit none

      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_DynVariables.h'
      include 'CB_ThermoVariables.h'
      include 'CB_ThermoForcing.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_buoys.h'
      include 'CB_Dyndim.h'

      character(LEN=32) filename  ! restart file name

      integer, intent(in) :: restart, expno_r
      integer :: i, j, k, year, month, day, hour, minute
      TYPE(datetime_type), intent(in) :: restart_date

      year = restart_date%year
      month = restart_date%month
      day = restart_date%day
      hour = restart_date%hour
      minute = restart_date%minute

!------------------------------------------------------------------------
!     Set initial conditions for h,A,u,v,Ta,Ti,Tl
!------------------------------------------------------------------------

      if ( restart .eq. 0 ) then

         do i = 0, nx+1
            do j = 0, ny+1

!
!     already zeros on edges

               h(i,j)   =  1d0 * maskC(i,j) ! initial ice thick
               A(i,j)   =  1d0 * maskC(i,j) ! initial ice conc
!
!     h and A set to zero on open boundaries 
!	       if (i.lt.11 .or. i.gt.nx-10 .or.j.eq.ny+1) then !lt.3 gt.2
!                   h(i,j) = 0d0
!                   A(i,j) = 0d0
!	       endif 

!     FB: This is to see effects on heterogineity check Ringeisen et al. 2019 
      if (Uni_Load_Exp .eqv. .true.) then          

               if (i.lt.11 .or. i.gt.nx-10 .or. j.eq.ny+1) then !FB: lt. 3 gt.2
                  h(i,j) = 0d0  !FB: initial ice thick. Check Ringeisen et al. 2019 section 3.2.3
                  A(i,j) = 0d0  !FB: initial ice concentration
               endif
      elseif (Uni_Load_Exp .eqv. .false.) then
		  h(i,j) = 1d0 * maskC(i,j) ! initial ice thickness  
		  A(i,j) = 1d0 * maskC(i,j) ! initial ice concentration	
      endif !Uni_Load_Exp


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

      endif !restart

!------------------------------------------------------------------------
!     Load restart files for h,A,u,v,Ta,Ti,Tl initial conditions
!------------------------------------------------------------------------

      if ( restart .eq. 1 ) then               ! load restart files

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


    end subroutine ini_get

