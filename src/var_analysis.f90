
      subroutine var_analysis

      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      integer i, j, ncell
      integer ihmax, jhmax, iumax, jumax, ivmax, jvmax
      
      double precision umax,vmax,hmax, havg

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

      return
    end subroutine var_analysis







