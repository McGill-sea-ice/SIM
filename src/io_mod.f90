MODULE IO

!  ===================
!      IO MODULE
!  ===================
!
!  Conains all input/output operations. 
!
!  Procedures:
!    * load_ocn_current(delta, u, v)
!    * load ocn_temperature(delta, t)
!    * load_air_temperature(date, delta, t)
!    * load_geostrophic_wind(date, delta, u, v)
!
!  Utilities:
!    * read_array(unit, out)

  USE datetime, DTRP=>RP
  USE netcdf

  IMPLICIT NONE

  include 'parameter.h'

  INTEGER, PARAMETER :: RP = SELECTED_REAL_KIND(12)
  INTEGER :: STD_OUT = 6    ! Default standard output
  INTEGER :: LOG_OUT = 9    ! Log

  TYPE NETCDF_ID_TYPE
     ! A simple container for variables ID.
     ! These IDs are used to identify variables in netCDF files and the file itself.
     INTEGER :: file                 ! File ID
     TYPE(datetime_type) :: since    ! Starting date for output file
     CHARACTER(len=6)    :: units    ! Counting units for the index.
     INTEGER :: step                 ! Number of units in each time step.     
     

     ! Dimension IDs
     INTEGER :: dim_time, dim_x, dim_y

     ! Variables IDs
     INTEGER :: time        ! time (units since)
     INTEGER :: h           ! Ice thickness
     INTEGER :: A           ! Ice covered area
     INTEGER :: u,v         ! Ice velocity 
     INTEGER :: ta, to, ti  ! Air, ocean and ice temperature
     INTEGER :: pvap        ! Atmospheric vapor pressure
     INTEGER :: qsh_io      ! Sensible heat flux (ice/ocean) 
     INTEGER :: Qoa         ! Heat flux (ocean/atmosphere)

     
  END TYPE NETCDF_ID_TYPE



  INTERFACE read_array
     ! Read an array from a file. 
     ! 
     ! If the `out` array has shape (nx, ny), the first ny rows and nx columns 
     ! will be read from the file. Remaining rows and columns are ignored. 
     !
     ! :Input:
     ! unit : integer
     !   The unit file number.
     ! 
     ! :Output:
     ! out : 2D array (real or integer)
     !   The output array. 
     MODULE PROCEDURE read_array_real__, read_array_integer__
  END INTERFACE


  PRIVATE :: read_array_real__, read_array_integer__


CONTAINS


  !---------------------------------------------------------------------
  !             Type specific subroutines for read_array.              !
  !---------------------------------------------------------------------
  

   SUBROUTINE read_array_real__(unit, out)
    ! Read a real array in unit.

    INTEGER, INTENT(in) :: unit
    REAL(KIND=RP), INTENT(out) :: out(:,:)
    INTEGER :: nx, ny, i,j
    
    nx = SIZE(out,1)
    ny = SIZE(out,2)

    DO j=1, ny
       READ(unit, *) (out(i,j), i = 1, nx)
    END DO
  END SUBROUTINE read_array_real__


  SUBROUTINE read_array_integer__(unit, out)
    ! Read a integer array in unit.

    INTEGER, INTENT(in) :: unit
    INTEGER, INTENT(out) :: out(:,:)
    INTEGER :: nx, ny, i,j
    
    nx = SIZE(out,1)
    ny = SIZE(out,2)

    DO j=1, ny
       READ(unit, *) (out(i,j), i = 1, nx)
    END DO
  END SUBROUTINE read_array_integer__

  SUBROUTINE read_mask(unit, out)
    ! Read a mask file. 

    INTEGER, INTENT(in) :: unit
    INTEGER, INTENT(out) :: out(:,:)
    INTEGER :: nx, ny, i, j

    nx = SIZE(out, 1)
    ny = SIZE(out, 2)
    WRITE(*,*) 'ny:', ny, 'nx: ', nx
    DO j=1, ny
       READ(unit, '(i1)') (out(i,j), i=1,nx)
    END DO
  END SUBROUTINE read_mask





  !---------------------------------------------------------------------
  !                          OCEAN CURRENTS                            !
  !---------------------------------------------------------------------

  SUBROUTINE load_ocn_current(delta, u, v)

    ! Load the climatological ocean currents.
    !
    ! The values are specified on a C-grid (staggered). 
    !
    ! :Input:
    ! delta : {10,20,40,80}
    !   The model grid resolution. 
    ! 
    ! :Output:
    ! u : 2D real array (nx+1, ny)
    !   Current velocity in the x direction [m/s].
    ! v : 2D real array (nx, ny+1)
    !   Current velocity in the y direction [m/s].
    

    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: delta
    REAL(KIND=RP), DIMENSION(:,:), INTENT(out) :: u,v

    CHARACTER(len=*), PARAMETER :: dir = 'forcing/ocncurrent/' 
    CHARACTER(len=80) :: name_u, name_v
    CHARACTER(len=5) :: cdelta
    CHARACTER(len=4) :: cnx, cny
    INTEGER :: nx,ny

    nx = SIZE(v,1)
    ny = SIZE(u,2)

    write(cnx, '(I0)') int(nx)
    write(cny, '(I0)') int(ny)

    ! Define the names of the file.
    if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
        & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
        write(cdelta, '(I0)') int(Deltax/1d03)
        name_u = dir // cdelta // 'uwater' // cdelta // '.clim'
        name_v = dir // cdelta // 'vwater' // cdelta // '.clim'
    elseif ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
        & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
        & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
        write(cdelta, '(F0.2)') Deltax/1d03
        name_u = dir // 'uwater_dx' // trim(cdelta) // '_nx' &
               & // trim(cnx) // '_ny' // trim(cny)// '.clim'
        name_v = dir // 'vwater_dx' // trim(cdelta) // '_nx' &
               & // trim(cnx) // '_ny' // trim(cny)// '.clim'
    endif

    ! Old names
    ! name_u= 'forcing/current/uwater.clim'
    ! name_v= 'forcing/current/vwater.clim'


    OPEN(unit=30, file=TRIM(name_u), status='old')
    OPEN(unit=31, file=TRIM(name_v), status='old')

    ! Read in the arrays.
    CALL read_array(30, u)
    CALL read_array(31, v)

    CLOSE(30)
    CLOSE(31)

  END SUBROUTINE load_ocn_current


  !---------------------------------------------------------------------
  !                         OCEAN TEMPERATURES                         !
  !---------------------------------------------------------------------


  SUBROUTINE load_ocn_temperature(delta, temp)

    ! Load the climatological ocean temperature.
    !
    ! The values are specified on a C-grid (tracer point). 
    !
    ! :Input:
    ! delta : {10,20,40,80}
    !   The model grid resolution [km]. 
    ! 
    ! :Output:
    ! temp : 3D real array (nx+2, ny+2, 14)
    !   Ocean temperature [C] on the model grid from december to january. 
    !   temp(:,:,1) -> december
    !   temp(:,:,2) -> january
    !   temp(:,:,13) -> december
    !   temp(:,:,14) -> january
    
    IMPLICIT NONE

    INTEGER, INTENT(in) :: delta
    REAL(KIND=RP), DIMENSION(:,:,:), INTENT(out) :: temp

    CHARACTER(len=*), PARAMETER :: dir = 'forcing/ocnT/'
    CHARACTER(len=80) :: name
    CHARACTER(len=2) :: cmonth
    CHARACTER(len=5) :: cdelta
    CHARACTER(len=4) :: cnx, cny
    INTEGER :: nx,ny,month

    nx = SIZE(temp,1)
    ny = SIZE(temp,2)
    IF (SIZE(temp,3) /= 14) STOP

    write(cnx, '(I0)') int(nx)
    write(cny, '(I0)') int(ny)

    ! Read in the climatological temperature for each month.
    ! In the files, land values are indicated with -4.
    ! Replace those -4 by -9999 so errors can be catched. 
    DO month=1,12
      
       WRITE(cmonth, '(I2.2)') month
       if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
           & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
           write(cdelta, '(I0)') int(Deltax/1d03)
           name = dir // cdelta // '/Tocn' // cmonth
           OPEN(unit=32, file=TRIM(name), status='old')
           CALL read_array(32, temp(:,:,month+1))
           WHERE (temp(:,:,month+1) <= -4) temp(:,:,month+1) = -9999.
           CLOSE(32)
       elseif  ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
           & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
           & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
           write(cdelta, '(F0.2)') Deltax/1d03
           name = dir // 'Tocn' // cmonth // '_dx' // trim(cdelta) &
               & // '_nx' // trim(cnx) // '_ny' // trim(cny)
           OPEN(unit=32, file=TRIM(name), status='old')
           CALL read_array(32, temp(:,:,month+1))
           WHERE (temp(:,:,month+1) <= -4) temp(:,:,month+1) = -9999.
           CLOSE(32)
       endif
    END DO


    ! Fill the end points. 
    temp(:,:,1) = temp(:,:,13)
    temp(:,:,14)= temp(:,:,2)
    
    ! Convert to Kelvins
    WHERE(temp >= -9998) temp = temp + 273.15_RP


  END SUBROUTINE load_ocn_temperature


  !---------------------------------------------------------------------
  !                              MASKS                                 !
  !---------------------------------------------------------------------

  SUBROUTINE load_mask(delta, mask)
    ! Load the mask (0:land, 1:ocean) at the given
    ! model resolution. 
    !
    ! :Input:
    ! delta : {10,20,40,80}
    !   The model grid resolution [km].
    !
    ! :Output:
    ! mask : 2D integer array
    !  Array of 0s (land) and 1s (ocean). 

    IMPLICIT NONE

    INTEGER, INTENT(in) :: delta
    INTEGER, DIMENSION(:,:), INTENT(out) :: mask
    INTEGER :: nx, ny
    CHARACTER(len=5) :: cdelta
    CHARACTER(len=4) :: cnx, cny
    CHARACTER(len=80) :: name
    
    nx = SIZE(mask,1)
    ny = SIZE(mask,2)
    
    write(cnx, '(I0)') int(nx)
    write(cny, '(I0)') int(ny)

    if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
        & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
        write(cdelta, '(I0)') int(Deltax/1d03)
        name = 'src/mask' // cdelta // '.dat'
    elseif  ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
        & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
        & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
        write(cdelta, '(F0.2)') Deltax/1d03
        name = 'src/mask_dx' // trim(cdelta) &
                & // '_nx' // trim(cnx) // '_ny' // trim(cny) // '.dat'
    endif

    OPEN(unit=33, file=TRIM(name), status='old')

    CALL read_array(33, mask)
   
  END SUBROUTINE load_mask


  !---------------------------------------------------------------------
  !                                WINDS                               !
  !---------------------------------------------------------------------
 
    SUBROUTINE load_geostrophic_wind(dt, delta, u, v)
      ! Return geostrophic winds at the given date and time.
      !
      ! :Input:
      ! dt : datetime_type
      !   The date and time at which the winds are desired. Winds are available
      !   from 1979-01-01:00 to 2006-12-31:18
      ! delta : integer
      !   The model grid resolution [km].
      !
      ! :Inout:
      ! u : (0:nx+2, 0:ny+2) array
      !   x-component of the geostrophic wind.
      ! v : (0:nx+2, 0:ny+2) array
      !   y-component of the geostrophic wind.
      ! 
      ! Notes
      ! -----
      ! `uair` and `vair` are the wind forcing arrays used in the model and 
      ! are of dimensions(0:nx+2, 0:ny+2).
      ! The geostrophic winds are computed from the sea level pressure from ERA5
      ! and interpolated using cubic interpolation from xarray. 
      ! 
      ! :author: Huard
      ! :date: Sept. 2008
      
      USE DATETIME, ONLY: datetime_type, datetime_init, datetime2num
      USE NETCDF, ONLY: nf90_open, nf90_nowrite, nf90_inq_varid
      USE NETCDF, ONLY: nf90_close, nf90_get_var

      IMPLICIT NONE

      TYPE(datetime_type), INTENT(in) :: dt
      INTEGER, INTENT(in) :: delta
      REAL(KIND=RP), DIMENSION(:,:), INTENT(inout) :: u, v

      CHARACTER(LEN=*), PARAMETER :: dir = "forcing/wind/"
      CHARACTER(len=80) :: file_name
      CHARACTER(len=5) :: cdelta
      CHARACTER(len=4) :: cnx, cny

      INTEGER :: ncid     ! This file ID.
      INTEGER :: u_id, v_id, time_id ! Variable IDs.
 
      INTEGER, PARAMETER :: ndims=3  ! Number of dimensions of the variables (y, x, time)
      INTEGER, PARAMETER :: steps=6  ! Time resolution of winds in the netCDF file.
      INTEGER, DIMENSION(ndims) :: start, count  ! Lower and upper bounds for indices.
      INTEGER :: time_index, nx, ny
      REAL(KIND=RP) :: time0, n_units

      TYPE(datetime_type):: since

      ! Check that the wind arrays have the right shape. 
      CALL nxny(delta, nx, ny)
      IF ( nx /= SIZE(u,1)-3) STOP
      IF ( ny /= SIZE(u,2)-3) STOP
      write(cnx, '(I0)') int(nx)
      write(cny, '(I0)') int(ny)

      if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
          & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
          write(cdelta, '(I2)') delta
          file_name = dir // cdelta // '/geowinds.nc'
      elseif ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
          & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
          & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
          write(cdelta, '(F0.2)') Deltax/1d03
          file_name = dir // 'wind_dx' // trim(cdelta) // '_nx' &
               & // trim(cnx) // '_ny' // trim(cny)// '.nc'
      endif

      ! Open the file. 
      CALL check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
      
      ! Compute the time index corresponding to the given date.
      ! Start date of file. 
      since = datetime_init(1900,1,1,0) 
      CALL check( nf90_inq_varid(ncid, "time", time_id) )

      CALL check( nf90_get_var(ncid, time_id, time0) ) 
      n_units = datetime2num(dt, 'hours', since)
      time_index = INT((n_units-time0)/steps)+1

      start = (/ 1, 1,time_index /)
      count = (/ nx+1, ny+1, 1 /)

      ! Get the variable id of uwnd and vwnd
      CALL check( nf90_inq_varid(ncid, "uwnd", u_id) )
      CALL check( nf90_inq_varid(ncid, "vwnd", v_id) )


      ! Read the data.
      u = 0.; v = 0. 
      CALL check(nf90_get_var(ncid, u_id, u(2:nx+2, 2:ny+2), start=start, count=count))
      CALL check(nf90_get_var(ncid, v_id, v(2:nx+2, 2:ny+2), start=start, count=count))

      CALL check( nf90_close(ncid) )

    END SUBROUTINE load_geostrophic_wind
      






  SUBROUTINE nxny(delta, nx, ny)
    ! Return nx and ny (interior grid points) given the model resolution.

    INTEGER, INTENT(in) :: delta
    INTEGER, INTENT(out) :: nx, ny

    ! Total number of cells, with boundaries.
    SELECT CASE (delta)
       CASE (10)
          nx = 512
          ny = 432
       CASE (20)
          nx = 256
          ny = 216
       CASE (40)
          nx = 128
          ny = 108
       CASE (80)
          nx = 64
          ny = 54
       CASE(1)
          nx = 5120
          ny = 4320
       CASE(2)
          nx = 2560
          ny = 2160
       CASE(4)
          nx = 1280
          ny = 1080
       CASE(8)
          nx = 640
          ny = 540
       CASE(16)
          nx = 320
          ny = 270
       CASE(32)
          nx = 160
          ny = 135
     END SELECT
       
   END SUBROUTINE nxny




    SUBROUTINE daily_air_temp_from_monthly_mean(dt,delta, Ta)

      ! Given a datetime, return the daily air temperature
      ! interpolated from monthly means.
      ! 
      ! 
      ! The routine ``load_monthly_mean_air_temperature``
      ! returns mean monthly air temperature a the middle of each
      ! month. In order to interpolate the temperature at any given
      ! moment, we need to find the months that occur before and after
      ! the given date. So to interpolate the temperature on February 12, 
      ! we need to load the temperatures for January (15) and February (14).
      ! 
      ! :Input:
      ! dt : datetime_type
      !  Date and time at which the air temperature is to be
      !  interpolated.
      !
      ! :Output:
      ! Ta : 2D array
      !   Air temperature [K].
      !
      ! :Warning:
      ! This subroutine assumes that it will be called 
      ! successively with increasing dt.

      USE datetime
      USE utils, ONLY : interpolate

      TYPE(datetime_type), INTENT(in) :: dt
      INTEGER, INTENT(in) :: delta
      REAL(KIND=RP), DIMENSION(:,:), INTENT(out) :: Ta

      TYPE(datetime_type), SAVE :: last_dt2
      INTEGER, SAVE :: last_delta=-10
      TYPE(datetime_type):: dt1, dt2, mid1, mid2
      TYPE(datetime_delta_type) :: delta_t

      REAL(KIND=RP), DIMENSION(:,:), SAVE, ALLOCATABLE :: Ta1, Ta2
      INTEGER :: ndays   ! Number of days in current month
      INTEGER :: dh, nx, ny
      REAL(KIND=RP) :: w

      
      ! Checking that Ta has the correct shape.
      CALL nxny(delta, nx, ny)
      IF (SIZE(Ta, 1) /= nx+2) STOP
      IF (SIZE(Ta, 2) /= ny+2) STOP


      IF (delta /= last_delta) THEN
         ALLOCATE(Ta1(nx+2,ny+2), Ta2(nx+2,ny+2))
         last_dt2 = datetime_init(1,1,1)
      END IF


      ndays = days_in_month(dt%year, dt%month)
      dh = ndays * 12  ! Half the number of hours in the month
      delta_t = delta_init(hours=dh)

      ! Setting dt1 and dt2 at the beginning of the months used to interpolate
      ! Ta.

      dt1 = dt - delta_t
      dt1%day = 1
      dt1%hour = 0

      dt2 = dt + delta_t
      dt2%day = 1
      dt2%hour = 0 

      ! Computing the middle of those months
      mid1 = dt1 + delta_init(hours=12*days_in_month(dt1%year, dt1%month))
      mid2 = dt2 + delta_init(hours=12*days_in_month(dt2%year, dt2%month))

      IF (dt2 == last_dt2) THEN
         ! Do nothing. 
      ELSEIF (dt1 == last_dt2)  THEN
         ! Advance one month.
         Ta1 = Ta2
         CALL load_monthly_mean_air_temp(dt2, delta, Ta2)
      ELSE
         ! Fetch everything.
         CALL load_monthly_mean_air_temp(dt1, delta, Ta1)
         CALL load_monthly_mean_air_temp(dt2, delta, Ta2)
      END IF

      w = hours(dt-mid1) / hours(mid2 - mid1)

      Ta  = interpolate(Ta1, Ta2, w)

      last_dt2 = dt2
      last_delta = delta

    END SUBROUTINE daily_air_temp_from_monthly_mean


    SUBROUTINE load_monthly_mean_air_temp(dt,delta, Ta)
      ! Put the monthly mean air temperature in Ta.
      !
      ! :Input:
      ! dt : datetime_type
      !   The date and time at which the winds are desired. Available
      !   air temperature cover 1978-12-1 to 2008-2-1. All dates in a given
      !   month will return the same Ta.
      ! delta : {10,20,40,80}
      !   Model grid resolution [km].
      !
      ! :Output:
      ! Ta : 2D real array
      !   Monthly mean temperature from NCEP1 interpolated on the model C-grid
      !   using bicubic Akima interpolation.
      !
      ! :author: D.Huard
      ! :date: June 16, 2008. Adapted on Sept 25. 2008 to new model grid.
      ! 


      USE DATETIME, ONLY: datetime_type, datetime_init, datetime2num, hours, OPERATOR(-)
      USE NETCDF, ONLY: nf90_open, nf90_nowrite, nf90_inq_varid
      USE NETCDF, ONLY: nf90_close, nf90_get_var
      USE utils, ONLY: find

      IMPLICIT NONE

      TYPE(datetime_type), INTENT(in) :: dt
      INTEGER, INTENT(in) :: delta
      REAL(KIND=RP), DIMENSION(:,:), INTENT(inout) :: Ta

      CHARACTER(len=80) :: file_name
!      CHARACTER(LEN=*), PARAMETER :: file_name = "air.mon.mean.nc"
      CHARACTER(len=*), PARAMETER :: dir = "forcing/airT/"
      CHARACTER(len=5) :: cdelta
      CHARACTER(len=4) :: cnx, cny

      TYPE(datetime_type) :: dtstart   ! The beginning of the month of dt.
      INTEGER :: ncid     ! This file ID.
      INTEGER :: air_id   ! Variable IDs for airT
      INTEGER :: time_id, time_dim_id  
      INTEGER, PARAMETER :: ndims=3  ! Number of dimensions of the variables (x, y, time)
      INTEGER, DIMENSION(ndims) :: start, count  ! Lower and upper bounds for indices.
      INTEGER :: time_index, time_dim, time
      INTEGER, DIMENSION(:), ALLOCATABLE :: times
      
      INTEGER :: nx, ny

      TYPE(datetime_type):: since

      since = datetime_init(1900,1,1,0)  ! Start date of file. 

      CALL nxny(delta, nx, ny)
      IF (SIZE(Ta, 1) /= nx+2) STOP
      IF (SIZE(Ta, 2) /= ny+2) STOP

      ! Beginning of the month
      dtstart = datetime_init(dt%year, dt%month, 1, 0)

      ! Initialize the output variables:
      Ta = 0.0

      ! Open the file. 
      ! NF90_NOWRITE tells netCDF we want READ-ONLY access to the file.
      write(cnx, '(I0)') int(nx)
      write(cny, '(I0)') int(ny)

      if  ( ( Deltax == 80d03 ) .or. ( Deltax == 40d03 ) .or. &
          & ( Deltax == 20d03 ) .or. ( Deltax == 10d03 ) ) then
          write(cdelta, '(I2)') delta
          file_name = dir // cdelta // '/air.mon.mean.nc'
      elseif  ( ( Deltax == 32d03 ) .or. ( Deltax == 16d03 ) .or. &
          & ( Deltax == 8d03 )  .or. ( Deltax == 4d03 ) .or. &
          & ( Deltax == 2d03 )  .or. ( Deltax == 1d03 ) ) then
          write(cdelta, '(F0.2)') Deltax/1d03
          file_name = dir // 'air.mon.mean_dx' // trim(cdelta) &
                & // '_nx' // trim(cnx) // '_ny' // trim(cny) // '.nc'
      endif

      CALL check( nf90_open(file_name, NF90_NOWRITE, ncid) )

      ! Get the variable ids
      CALL check( nf90_inq_varid(ncid, "air", air_id) )
      CALL check( nf90_inq_varid(ncid, "time", time_id) )
      
      ! Get the time dimension id
      CALL check( nf90_inq_dimid(ncid, 'time', time_dim_id) )
      CALL check( nf90_inquire_dimension(ncid, time_dim_id, len=time_dim) )
      
      ! Get the time variable
      ALLOCATE(times(time_dim))
      CALL check(nf90_get_var(ncid, time_id, times) )
      
      ! Compute the time index corresponding to the given date
      time = INT(hours(dtstart - since))
      time_index = find(times, time)

      start = (/ 1,1,time_index /)
      count = (/ nx+1, ny+1, 1 /)
      
      CALL check(nf90_get_var(ncid, air_id, Ta, start=start, count=count))

      CALL check( nf90_close(ncid) )

      DEALLOCATE(times)

    END SUBROUTINE load_monthly_mean_air_temp


    SUBROUTINE check(status)
      INTEGER, INTENT ( in) :: status

      
      IF(status /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(status))
         STOP 2
      END IF
    END SUBROUTINE check


  END MODULE IO
