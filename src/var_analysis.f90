
      subroutine var_analysis(date, expno)
        use datetime, only: datetime_type, datetime_str
      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      type(datetime_type), intent(in) :: date
      integer, intent(in) :: expno
      integer i, j, ncell, hour, year, month, day
      integer ihmax, jhmax, iumax, jumax, ivmax, jvmax

      character filename*30

      double precision umax,vmax,hmax, havg, VOL

      umax = 0d0
      vmax = 0d0
      hmax = 0d0
      ihmax = 0
      jhmax = 0
      iumax = 0
      jumax = 0
      ivmax = 0
      jvmax = 0

      ncell = 0
      havg =0d0

      do j=1,ny
         do i=1,nx

            if (h(i,j) .gt. hmax) then 
               hmax = h(i,j) 
               ihmax = i
               jhmax = j
            endif

            if (abs(uice(i,j)) .gt. umax) then
               umax = abs(uice(i,j) )
               iumax = i
               jumax = j
            endif

            if (abs(vice(i,j)) .gt. vmax) then
               vmax = abs(vice(i,j) )
               ivmax = i
               jvmax = j
            endif
            
            if (maskC(i,j) .eq. 1) then
               ncell = ncell + 1
               havg = havg + h(i,j)
            endif
                      
         enddo
      enddo
      
      print *, 'havg (m)   = ', havg/ncell
      print *, 'hmax (m)   = ', hmax, ihmax, jhmax
      print *, 'umax (m/s) = ', umax, iumax, jumax
      print *, 'vmax (m/s) = ', vmax, ivmax, jvmax

!-------------------------------------------------------
!     Calc and output total volume in km^3
!------------------------------------------------------- 
      hour = date%hour

      if (hour == 0) then
         VOL=havg*Deltax2*1d-09
         print *, 'Volume [km^3]=', VOL
         year = date%year
         month = date%month
         day = date%day
         write (filename, '("output/sea_ice_volume.",i2.2)') expno
         open (10, file = filename, status = 'unknown')
         write(10,10) year, month, day, VOL
      endif

10 format (i4.4, i2.2, i2.2, f20.10)

      return
    end subroutine var_analysis







