!****************************************************************************
! 
! The following subroutines are used to create monthly mean fields. There is 
! a WARNING posted if the initial and last months are incomplete (the monthly
! mean fields are still calculated anyway but did not use use all the time 
! steps in a month for the calculation. For the first and months to be 
! complete, the date should be XXXX-XX-01 00.
!
! Let's say there are m months calculated, the output files are of the form
! 
! hmean, Amean: 1...nx in the 'horixzontal'
!               m x (1...ny) in the 'vertical'
!
! umean:        1...nx+1 in the 'horixzontal'
!               m x (1...ny) in the 'vertical' 
!
! vmean:        1...nx in the 'horixzontal'
!               m x (1...ny+1) in the 'vertical' 
!
! Date: Sept 29, 2009
! Author: JF Lemieux
! Intitution: McGill University
!
!****************************************************************************

subroutine prep_monthly_mean_fields ( end_date, start_date )
 
  USE datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
  USE datetime, ONLY: datetime_delta_type, OPERATOR(+),OPERATOR(-)
 
  implicit none
  
  TYPE(datetime_type) :: start_date, end_date
  
  print *,
  print *, '********* monthly mean fields will be calculated ********'

  if ( start_date%day .ne. 1 .or. start_date%hour .ne. 0 .or. &
       start_date%minute .ne. 0 )  then
     print *, 'WARNING: 1st month is incomplete for monthly mean fields'

  endif

  if ( end_date%day .ne. 1 .or. end_date%hour .ne. 0 .or. &
       end_date%minute .ne. 0 ) then
     print *, 'WARNING: last month is incomplete for monthly mean fields'

  endif

  print *,

  return
end subroutine prep_monthly_mean_fields


!****************************************************************************
!     calculates monthly mean fields
!****************************************************************************

subroutine monthly_mean_fields (htp, Atp, utp, vtp, expno, m_current, &
                                ts_m_counter, end_date, now_date )
  USE datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
  USE datetime, ONLY: datetime_delta_type, OPERATOR(+),OPERATOR(-),OPERATOR(==)

  implicit none
      
  include 'parameter.h'
  include 'CB_DynVariables.h'

  integer expno, m_current, ts_m_counter, ts_m_counter_1

  double precision htp(0:nx+1,0:ny+1), Atp(0:nx+1,0:ny+1)
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

  TYPE(datetime_type) :: end_date, now_date

  ts_m_counter_1 = ts_m_counter
  ts_m_counter = ts_m_counter + 1

  htp = ( htp * ts_m_counter_1 + h )    / ts_m_counter
  Atp = ( Atp * ts_m_counter_1 + A )    / ts_m_counter
  utp = ( utp * ts_m_counter_1 + uice ) / ts_m_counter
  vtp = ( vtp * ts_m_counter_1 + vice ) / ts_m_counter
  
  if ( now_date .eq. end_date ) then

     call output_monthly_mean_fields ( htp, Atp, utp, vtp, expno)

  else

     if ( m_current .ne. now_date%month) then

        call output_monthly_mean_fields ( htp, Atp, utp, vtp, expno)
        
        htp = 0d0 ! initialize to zero because of new month
        Atp = 0d0
        utp = 0d0
        vtp = 0d0

        ts_m_counter = 0 ! reset counter because of new month
     
        if ( m_current .eq. 12 ) then 
           m_current = 1             ! year changed
        else
           m_current = m_current + 1 ! month changed
        endif
        
     endif

  endif

  return
end subroutine monthly_mean_fields


subroutine output_monthly_mean_fields ( htp,Atp, utp, vtp, expno)

  implicit none
  
  include 'parameter.h'
  include 'CB_DynVariables.h'

  integer i,j,expno
  character filename*30
  double precision htp(0:nx+1,0:ny+1), Atp(0:nx+1,0:ny+1)
  double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

  print *, 'outputing monthly mean fields'

  write (filename, '("output/hmean.",i2.2)') expno
  open (11, file = filename, status = 'unknown', position = 'append')
  write (filename, '("output/Amean.",i2.2)') expno
  open (12, file = filename, status = 'unknown', position = 'append')
  write (filename, '("output/umean.",i2.2)') expno
  open (13, file = filename, status = 'unknown', position = 'append')
  write (filename, '("output/vmean.",i2.2)') expno
  open (14, file = filename, status = 'unknown', position = 'append')

     do j = 1, ny
        write(11,10) ( htp(i,j),       i = 1, nx )
        write(12,10) ( Atp(i,j),       i = 1, nx )
        write(13,10) ( utp(i,j),       i = 1, nx+1 )
     enddo

     do j = 1, ny+1
        write(14,10) ( vtp(i,j),       i = 1, nx )
     enddo

  do i=11, 14
     close(i)
  enddo

10 format (1x, 1000(f10.6, 1x))

  return
end subroutine output_monthly_mean_fields
