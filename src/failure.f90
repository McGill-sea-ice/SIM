!************************************************************************
!     Subroutine failure: Routine used to investigate failures of JFNK. 
!                         It computes and outputs on the residual 
!                         for the u and v components on the grid. It also
!                         outputs h, A, R1, R2 and the velocity field.
!
!     JF Lemieux, 9 September 2011
!
!************************************************************************

      subroutine failure(xtp, rhs, date, k, expno)
        use datetime, only: datetime_type
      implicit none

      include 'parameter.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_DynVariables.h'
      
      type(datetime_type), intent(in) :: date

      character filename*36

      integer, intent(in) :: k, expno
      integer i, j, l, year, month, day, hour, minute

      double precision, intent(in):: xtp(nvar), rhs(nvar)
      double precision :: Ftp(nvar)
      double precision :: resu(0:nx+2,0:ny+2), resv(0:nx+2, 0:ny+2)
      double precision :: res(0:nx+1,0:ny+1) ! at tracer point
      double precision :: Aout(0:nx+1,0:ny+1), hout(0:nx+1,0:ny+1)
      double precision :: land

      year = date%year
      month = date%month
      day = date%day
      hour = date%hour
      minute = date%minute

      land = -999d0
      res = 0d0
      resu =0d0
      resv = 0d0

      call Funk (xtp, rhs, Ftp)
      call transformer(resu,resv,Ftp,0)
      
      do j = 0, ny+1
         do i = 0, nx+1
      
            if (maskC(i,j) .eq. 1) then

               res(i,j) = ( resu(i,j)**2d0 + resu(i+1,j)**2d0 + &
                            resv(i,j)**2d0 + resv(i,j+1)**2d0 ) / 4d0

            elseif (maskC(i,j) .eq. 0) then
               res(i,j) = land
            else
               print *, 'Wo wo wo wrong mask'
            endif

         enddo
      enddo

      do j=0,ny+1
         do i=0,nx+1
         
            if (maskC(i,j) .eq. 1) then
               Aout(i,j) = A(i,j)
               hout(i,j) = h(i,j)
            elseif (maskC(i,j) .eq. 0) then
               Aout(i,j) = land
               hout(i,j) = land
            else
               print *, 'Wo wo wo wrong mask'
            endif

         enddo
      enddo

      write (filename,'("output/resu",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (12, file = filename, status = 'unknown')

      write (filename,'("output/resv",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (13, file = filename, status = 'unknown')

      write (filename,'("output/res",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (14, file = filename, status = 'unknown')

      write (filename,'("output/R1_",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (15, file = filename, status = 'unknown')

      write (filename,'("output/R2_",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (16, file = filename, status = 'unknown')
      
      write (filename,'("output/A",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (17, file = filename, status = 'unknown')

      write (filename,'("output/h",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (18, file = filename, status = 'unknown')

      write (filename,'("output/u",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (19, file = filename, status = 'unknown')

      write (filename,'("output/v",i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_k",i4.4,".",i2.2)') &
           year, month, day, hour, minute, k, expno
      open (20, file = filename, status = 'unknown')
      
      do j = 1, ny
         write(12,100) ( resu(i,j), i = 1, nx+1 )
         write(15,100) ( R1(i,j), i = 1, nx+1 )
         write(19,100) ( uice(i,j), i = 1, nx+1 )
      enddo

      do j = 1, ny+1
         write(13,100) ( resv(i,j), i = 1, nx )
         write(16,100) ( R2(i,j), i = 1, nx )
         write(20,100) ( vice(i,j), i = 1, nx )
      enddo

      do j = 0, ny+1
         write(14,100) ( res(i,j), i = 0, nx+1 )
         write(17,100) ( Aout(i,j), i = 0, nx+1 )
         write(18,100) ( hout(i,j), i = 0, nx+1 )
      enddo

      do l=12,20
         close(l)
      enddo

100            format (1x, 1000(f20.12, 1x))
      
      return
    end subroutine failure



