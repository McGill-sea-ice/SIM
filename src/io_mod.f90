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

    CHARACTER(len=*), PARAMETER :: dir = 'forcing10/current/' 
    CHARACTER(len=80) :: name_u, name_v
    CHARACTER(len=2) :: cdelta
    INTEGER :: nx,ny

    nx = SIZE(v,1)
    ny = SIZE(u,2)

    WRITE(cdelta, '(I2)') delta


    ! Define the names of the file.
    name_u = dir // cdelta // '/uwater' // cdelta // '.clim'
    name_v = dir // cdelta // '/vwater' // cdelta // '.clim'

    ! Old names
    name_u= 'forcing/current/uwater.clim'
    name_v= 'forcing/current/vwater.clim'



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
    CHARACTER(len=2) :: cdelta, cmonth
    INTEGER :: nx,ny,month

    nx = SIZE(temp,1)
    ny = SIZE(temp,2)
    IF (SIZE(temp,3) /= 14) STOP

    WRITE(cdelta, '(I2)') delta


    ! Read in the climatological temperature for each month.
    ! In the files, land values are indicated with -4.
    ! Replace those -4 by -9999 so errors can be catched. 
    DO month=1,12
      
       WRITE(cmonth, '(I2.2)') month
       name = dir // cdelta // '/Tocn' // cmonth

       name = 'forcing/Tocn/monthlyclim/Tocn' // cmonth

       OPEN(unit=32, file=TRIM(name), status='old')
       CALL read_array(32, temp(:,:,month+1))
       WHERE (temp(:,:,month+1) <= -4) temp(:,:,month+1) = -9999.
       CLOSE(32)
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
    CHARACTER(len=2) :: cdelta
    CHARACTER(len=80) :: name
    
    nx = SIZE(mask,1)
    ny = SIZE(mask,2)
    
    WRITE(cdelta, '(I2)') delta
    name = 'newforcing/mask/mask' // cdelta // '_.dat'
    name = 'src/maskC.dat'

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
      ! The data in the netCDF file is of dimension (nx+1, ny+1), so we need 
      ! to set the borders of the arrays at 0. 
      ! The geostrophic winds are computed from the sea level pressure from NCEP1
      ! and interpolated using Akima bicubic interpolation (see ncep2model.py in windblend/). 
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
      CHARACTER(len=2) :: cdelta


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


      ! Open the file. 
      WRITE(cdelta, "(I2)") delta
      file_name = dir // cdelta // '/geowinds.nc'
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
          nx = 520
          ny = 440
       CASE (20)
          nx = 260
          ny = 220
       CASE (40)
          nx = 130
          ny = 110
       CASE (80)
          nx = 65
          ny = 55
     END SELECT
       
     ! nx, ny is the number of interior grid points.
     nx = nx-2
     ny = ny-2
       
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

      CHARACTER(LEN=*), PARAMETER :: file_name = "air.mon.mean.nc"
      CHARACTER(len=*), PARAMETER :: dir = "forcing/airT/"
      CHARACTER(len=2) :: cdelta

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
      WRITE(cdelta, '(I2)') delta
      CALL check( nf90_open(dir//cdelta//"/"//FILE_NAME, NF90_NOWRITE, ncid) )

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
      count = (/ nx+2, ny+2, 1 /)

      CALL check(nf90_get_var(ncid, air_id, Ta, start=start, count=count))

      CALL check( nf90_close(ncid) )

      DEALLOCATE(times)

    END SUBROUTINE load_monthly_mean_air_temp
!!$      
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$    SUBROUTINE nc_load_NARR_wind(dt, u, v)
!!$      ! Return NARR-NCEP2 merged winds u and v at the date and time dt.
!!$      !
!!$      ! Input
!!$      ! -----
!!$      ! dt : datetime_type
!!$      !   The date and time at which the winds are desired. Winds are available
!!$      !   from 1979-01-01:00 to 2007-12-31:18
!!$      !
!!$      ! Inout
!!$      ! -----
!!$      ! u : (0:522, 0:442) array
!!$      !   x-component of the 10m wind.
!!$      ! v : (0:522, 0:442) array
!!$      !   y-component of the 10m wind.
!!$      ! 
!!$      ! Notes
!!$      ! -----
!!$      ! `uair` and `vair` are the wind forcing arrays used in the model and 
!!$      ! are of dimensions(0:nx+2, 0:ny+2), with nx=520 and ny=440.
!!$      ! The data in the netCDF file is of dimension (522, 442), so we need 
!!$      ! to set the borders of the arrays at 0. 
!!$      ! 
!!$      ! The NARR and NCEP2 winds are interpolated on the model grid then merged.
!!$      ! 
!!$      ! The Fortran ordering of dimensions is the reversed of the Python ordering. So while in 
!!$      ! Python the dimensions of the variables are (time, lat, lon), in fortran they become
!!$      ! (lon, lat, time). Since I didn't know about this when I created the dataset, we have
!!$      ! to transpose the arrays. 
!!$      !
!!$      ! :author: D. Huard
!!$      ! :date: May 2008
!!$      
!!$      USE DATETIME, ONLY: datetime_type, datetime_init, datetime2num
!!$      USE NETCDF, ONLY: nf90_open, nf90_nowrite, nf90_inq_varid
!!$      USE NETCDF, ONLY: nf90_close, nf90_get_var
!!$
!!$      IMPLICIT NONE
!!$
!!$      TYPE(datetime_type), intent(in) :: dt
!!$      REAL(KIND=RP), DIMENSION(0:522,0:442), intent(inout) :: u, v
!!$      REAL(KIND=RP), DIMENSION(1:441, 1:521) :: utmp, vtmp
!!$      CHARACTER(LEN=*), PARAMETER :: file_name = "forcing/NARR_NCEP2_uvwind.nc"
!!$
!!$      
!!$      INTEGER :: ncid     ! This file ID.
!!$      INTEGER :: nc_u_id, nc_v_id, time_id   ! Variable IDs for uwnd, vwnd and time.
!!$ 
!!$      INTEGER, PARAMETER :: ndims=3  ! Number of dimensions of the variables (y, x, time)
!!$      INTEGER, PARAMETER :: steps=6  ! Number of hours between each data. 
!!$      INTEGER, DIMENSION(ndims) :: start, count  ! Lower and upper bounds for indices.
!!$      INTEGER :: time_index
!!$      REAL(KIND=RP) :: time0, n_units, timei
!!$      TYPE(datetime_type):: since
!!$      
!!$      since = datetime_init(1900,1,1,0)  ! Start date of time unit. 
!!$
!!$
!!$      ! Initialize the output variables:
!!$      u = 0.
!!$      v = 0. 
!!$
!!$      ! Open the file. 
!!$      ! NF90_NOWRITE tells netCDF we want READ-ONLY access to the file.
!!$      !WRITE(*,*) 'loading file: ', FILE_NAME
!!$      CALL check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
!!$      
!!$
!!$      ! Get the variable id of u_geo_wind and v_geo_wind
!!$      CALL check( nf90_inq_varid(ncid, "uwnd", nc_u_id) )
!!$      CALL check( nf90_inq_varid(ncid, "vwnd", nc_v_id) )
!!$      CALL check( nf90_inq_varid(ncid, "time", time_id) )
!!$
!!$      ! Compute the time index corresponding to the given date
!!$      ! 1 : 1979-01-01:00
!!$      ! 2 : 1979-01-01:06
!!$      ! 3 : 1979-01-01:12
!!$      ! etc...
!!$
!!$      CALL check( nf90_get_var(ncid, time_id, time0) )
!!$      n_units = datetime2num(dt, 'hours', since)
!!$
!!$      time_index = INT((n_units-time0)/steps)+1
!!$
!!$      ! Check that the index is correct, ie that the time corresponds to the datetime given.
!!$      CALL check( nf90_get_var(ncid, time_id, timei, start=(/time_index/)))
!!$      IF (n_units /= timei) STOP
!!$
!!$
!!$      start = (/ 1,1,time_index /)
!!$      count = (/ 522, 442, 1 /)
!!$
!!$      CALL check(nf90_get_var(ncid, nc_u_id, u(0:521, 0:441), start=start, count=count))
!!$      CALL check(nf90_get_var(ncid, nc_v_id, v(0:521, 0:441), start=start, count=count))
!!$      CALL check( nf90_close(ncid) )
!!$
!!$    END SUBROUTINE nc_load_NARR_wind
!!$      
!!$
!!$
!!$
!!$
!!$    SUBROUTINE nc_create_output_file(file_name, id, overwrite)
!!$      ! Create a netCDF file with  to store modeling results. 
!!$      !
!!$      ! This routine gets an ID number for the output file and calls one
!!$      ! subroutine for each variable that is to be stored. 
!!$
!!$      USE netcdf
!!$      USE params, ONLY: nx, ny
!!$      IMPLICIT NONE
!!$     
!!$      TYPE( netcdf_id_type), intent (inout) :: id
!!$      CHARACTER(len=*), intent(in) :: file_name
!!$      LOGICAL, INTENT(in), OPTIONAL :: overwrite
!!$      
!!$      LOGICAL :: o_overwrite
!!$
!!$      ! This is the name of the data file we will create.
!!$      ! CHARACTER (len = *), PARAMETER:: FILE_NAME = "output/current_run.nc"
!!$     
!!$      ! We are writing 2D data, a 12 x 6 grid.
!!$      INTEGER, PARAMETER :: NDIMS = 3
!!$
!!$     
!!$      ! When we create netCDF files, variables and dimensions, we get back
!!$      ! an ID for each one.
!!$      INTEGER :: ncid, dimids(NDIMS)
!!$     
!!$      ! Always check the return code of every netCDF function call. In
!!$      ! this example program, wrapping netCDF calls with "call check()"
!!$      ! makes sure that any return which is not equal to nf90_noerr (0)
!!$      ! will print a netCDF error message and exit.
!!$      
!!$      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
!!$      ! overwrite this file, if it already exists.
!!$
!!$      o_overwrite = .FALSE.
!!$      IF (PRESENT(overwrite)) o_overwrite = overwrite
!!$
!!$      IF (o_overwrite) THEN
!!$         CALL check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
!!$      ELSE
!!$         CALL check( nf90_create(FILE_NAME, NF90_NOCLOBBER, ncid) )
!!$      END IF
!!$
!!$      id%file = ncid
!!$
!!$     
!!$      ! There are 3 locations on the grid: 
!!$      !   * the tracer point (center of the cell) (nx, ny)
!!$      !   * the C-grid u component (left center of the cell (nx+1, ny)
!!$      !   * the C-grid v component (bottom center of the cell (nx, ny+1)
!!$      !
!!$      ! If we include a boundary around the domain, the dimensions increase by 2.
!!$
!!$      ! To avoid defining different dimensions for each one of those, let's 
!!$      ! use the same dimensions (x, y, time) with some conventions:
!!$      !   * The grid is nx+2 by ny+2
!!$      !   * Latitudes and longitudes are given at the corners of the grid. 
!!$      !   * Scalars A and h are stored on (1:nx, 1:ny)
!!$      !   * Scalars To, Ta are stored on (0:nx+1, 0:ny+1)
!!$      !   * U_comp are stored on (0:nx, 1:ny)
!!$      !   * V_comp are stored on (1:nx, 0:ny)
!!$      !   * Empty locations are replaced by masked values = -99999.0
!!$
!!$       ! Define the dimensions for variables. NetCDF will hand back an ID for each.
!!$       CALL check( nf90_def_dim(ncid, "x", NX+2, id%dim_x) )
!!$       CALL check( nf90_def_dim(ncid, "y", NY+2, id%dim_y) )
!!$       CALL check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, id%dim_time ))
!!$       
!!$
!!$
!!$       CALL check( nf90_def_var(ncid, 'time', NF90_REAL, id%dim_time, id%time) )
!!$       CALL check( nf90_put_att(ncid, id%time, 'units', & 
!!$            TRIM(id%units) // ' since ' //datetime_str(id%since)) )
!!$       CALL check( nf90_put_att(ncid, id%time, 'step', id%step) )
!!$
!!$       ! The dimids array is used to pass the IDs of the dimensions of
!!$       ! the variables. Note that in fortran arrays are stored in
!!$       ! column-major format.
!!$       dimids =  (/ id%dim_x, id%dim_y, id%dim_time /)
!!$
!!$
!!$       ! Define the scalar variables.
!!$       call check( nf90_def_var(ncid, "h", NF90_REAL, dimids, id%h) )
!!$       CALL check( nf90_put_att(ncid, id%h, 'units', 'meter') )
!!$       CALL check (nf90_put_att(ncid, id%h, 'long_name', 'ice thickness') ) 
!!$
!!$       call check( nf90_def_var(ncid, "A", NF90_REAL, dimids, id%A) )
!!$       CALL check( nf90_put_att(ncid, id%A, 'units', '%') )
!!$       CALL check (nf90_put_att(ncid, id%A, 'long_name', 'ice covered area') ) 
!!$       
!!$
!!$       ! Define the vector variables
!!$       CALL check( nf90_def_var(ncid, 'u', NF90_REAL, dimids, id%u) )
!!$       CALL check( nf90_put_att(ncid, id%u, 'units', 'm/s') )
!!$       CALL check (nf90_put_att(ncid, id%u, 'long_name', 'x-component of the ice velocity') ) 
!!$
!!$       CALL check( nf90_def_var(ncid, 'v', NF90_REAL, dimids, id%v) )
!!$       CALL check( nf90_put_att(ncid, id%v, 'units', 'm/s') )
!!$       CALL check (nf90_put_att(ncid, id%v, 'long_name', 'y-component of the ice velocity') ) 
!!$
!!$       ! Define thermodynamical variables
!!$       CALL check( nf90_def_var(ncid, 'Ta', NF90_REAL, dimids, id%Ta) )
!!$       CALL check( nf90_def_var(ncid, 'To', NF90_REAL, dimids, id%To) )
!!$       CALL check( nf90_def_var(ncid, 'Ti', NF90_REAL, dimids, id%Ti) )
!!$
!!$       CALL check( nf90_def_var(ncid, 'Qsh_io', NF90_REAL, dimids, id%Qsh_io) )
!!$       CALL check( nf90_def_var(ncid, 'Qoa', NF90_REAL, dimids, id%Qoa) )
!!$
!!$
!!$       CALL nc_store_parameters(ncid)
!!$
!!$       ! End define mode. This tells netCDF we are done defining metadata.
!!$       call check( nf90_enddef(ncid) )
!!$
!!$
!!$     END SUBROUTINE nc_create_output_file
!!$
!!$     SUBROUTINE nc_save(id, dt)
!!$       ! Save the ice thickess to file. 
!!$       ! The netCDF ids are defined in id, and the time in datetime type dt.
!!$       !
!!$       ! TODO: apply the mask
!!$
!!$       USE netcdf
!!$       USE dynvariables, ONLY: h, A, uice, vice
!!$       USE thermovariables, ONLY: Ta, Ti, To, Qsh_io, Qoa
!!$       USE params, ONLY : nx, ny
!!$       IMPLICIT NONE
!!$
!!$       TYPE(netcdf_id_type), INTENT(inout) :: id
!!$       TYPE(datetime_type), INTENT(in) :: dt
!!$
!!$       INTEGER, PARAMETER :: ndims=3  ! Number of dimensions of the variables (y, x, time)
!!$       INTEGER, DIMENSION(ndims) :: start, count  ! Lower and upper bounds for indices.
!!$       REAL :: n_units, index_units, time0
!!$       INTEGER :: time_index, time_shape
!!$       
!!$       n_units = datetime2num(dt, id%units, id%since)
!!$
!!$       CALL check( nf90_inquire_dimension(id%file, id%dim_time, len=time_shape) )
!!$       IF (time_shape == 0) THEN
!!$          time_index = 1
!!$       ELSE
!!$          CALL check( nf90_get_var(id%file, id%time, time0) )
!!$          time_index = INT((n_units-time0)/id%step)+1
!!$       END IF
!!$
!!$       start = (/ 1,1,time_index /)
!!$       count = (/ nx+2, ny+2, 1 /)
!!$       
!!$       CALL check( nf90_put_var(id%file, id%time, n_units, start=(/time_index/)) ) 
!!$       CALL check( nf90_put_var(id%file, id%h, h, start=start, count=count) )
!!$       CALL check( nf90_put_var(id%file, id%A, A, start=start, count=count) )
!!$
!!$       CALL check (nf90_put_var(id%file, id%u, uice(0:nx+1, 0:ny+1), start=start, count=count) )
!!$       CALL check (nf90_put_var(id%file, id%v, vice(0:nx+1, 0:ny+1), start=start, count=count) )
!!$
!!$       CALL check( nf90_put_var(id%file, id%Ti, Ti, start=start, count=count) )
!!$       CALL check( nf90_put_var(id%file, id%Ta, Ta, start=start, count=count) )
!!$       CALL check( nf90_put_var(id%file, id%To, To, start=start, count=count) )
!!$       CALL check( nf90_put_var(id%file, id%Qsh_io, Qsh_io, start=start, count=count) )
!!$       CALL check( nf90_put_var(id%file, id%Qoa, Qoa, start=start, count=count) )
!!$
!!$       CALL check ( nf90_sync(id%file) )
!!$
!!$     END SUBROUTINE nc_save
!!$
!!$
!!$     SUBROUTINE nc_restart(file_name, dt, u, v, h, A, Ti, Ta, To, Qsh_io, Qoa)
!!$       ! Load a netCDF output file and return the values of the variables needed for 
!!$       ! a model restart. 
!!$       ! 
!!$       ! Input
!!$       ! -----
!!$       ! file_name : string
!!$       !   Path of the netCDF to use as a restart. 
!!$       ! dt : datetime_type
!!$       !   Date and time at which the fields are desired. 
!!$       !
!!$       ! Output
!!$       ! ------
!!$       ! u : 2D array
!!$       !   The x-component of the ice velocity [m/s].
!!$       ! v : 2D array
!!$       !   The y-component of the ice velocity [m/s].
!!$       ! h : 2D array
!!$       !   The ice thickness [m].
!!$       ! ...
!!$       USE params, ONLY: nx, ny
!!$
!!$       CHARACTER(len=*), intent(in) :: file_name
!!$       TYPE(datetime_type), INTENT(in) :: dt
!!$       REAL(KIND=RP), DIMENSION(:,:), INTENT(inout) :: u, v, h, A, Ti, Ta, To, Qsh_io, Qoa
!!$
!!$       TYPE(netcdf_id_type) :: id
!!$
!!$       INTEGER, PARAMETER :: ndims=3  ! Number of dimensions of the variables (y, x, time)
!!$       INTEGER, DIMENSION(ndims) :: start, count  ! Lower and upper bounds for indices.
!!$       INTEGER :: time_index
!!$            
!!$       id = nc_load_ids(file_name)
!!$
!!$
!!$       time_index = FLOOR(datetime2num(dt, id%units, id%since)/id%step)+1
!!$
!!$       start = (/ 1,1,time_index /)
!!$       count = (/ nx+2, ny+2, 1 /)
!!$       WRITE(LOG_OUT,*) '---------------------------------'
!!$       WRITE(LOG_OUT,*) 'Reading fields from restart file.'
!!$       WRITE(LOG_OUT,*) '---------------------------------'
!!$       CALL check( nf90_get_var(id%file, id%h, h, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%A, A, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%u, u, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%v, v, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%Ti, Ti, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%Ta, Ta, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%To, To, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%Qsh_io, Qsh_io, start=start, count=count) )
!!$       CALL check( nf90_get_var(id%file, id%Qoa, Qoa, start=start, count=count) )
!!$       WRITE(LOG_OUT,*) 'Done.'
!!$
!!$     END SUBROUTINE nc_restart
!!$
!!$
!!$     FUNCTION nc_last_entry(file_name) 
!!$       ! Return the datetime of the last entry in the file. 
!!$       ! This assumes that the time is in units of hours.
!!$       !
!!$       ! Input
!!$       ! -----
!!$       ! file_name : str
!!$       !   Name of the netCDF file containing the results.
!!$       !  
!!$       ! Output
!!$       ! ------
!!$       ! nc_last_entry : datetime_type
!!$       !   Datetime of the last_entry in the file.
!!$       CHARACTER(len=*), INTENT(in):: file_name
!!$       TYPE(datetime_type):: nc_last_entry, since
!!$       TYPE(netcdf_id_type) :: id
!!$       INTEGER:: time_index(1)
!!$       REAL(KIND=RP):: time
!!$       CHARACTER(len=50):: units
!!$       CHARACTER(len=20):: datestr
!!$       
!!$       id = nc_load_ids(file_name)
!!$       CALL check( nf90_inquire_dimension(id%file, id%dim_time, len=time_index(1)) )
!!$       CALL check( nf90_get_var(id%file, id%time, time, start=time_index) )
!!$       CALL check( nf90_get_att(id%file, id%time, 'units', units) )
!!$       units = ADJUSTR(units)
!!$       datestr = units(31:50)
!!$       since = str2dt(datestr)
!!$       nc_last_entry = since + delta_init(hours=int(time))
!!$
!!$     END FUNCTION nc_last_entry
!!$       
!!$     FUNCTION nc_load_ids(file_name, mode) RESULT (id)
!!$       ! Load the netcdf file containing the results from a previous run and retrieve
!!$       ! the id of the file, dimensions are variables. 
!!$       ! 
!!$       ! Input
!!$       ! -----
!!$       ! file_name : string
!!$       !   Path of the netCDF file to load. 
!!$       ! mode : integer 
!!$       !   NF90_NOWRITE, NF90_WRITE, NF90_SHARE. 
!!$       ! Output
!!$       ! ------
!!$       ! id : netcdf_id_type
!!$       !   ID of the file, dimensions and variables.
!!$
!!$       IMPLICIT NONE
!!$
!!$       CHARACTER(len = *), INTENT(in) :: FILE_NAME 
!!$       TYPE(netcdf_id_type) :: id
!!$       INTEGER, INTENT(in), optional :: mode
!!$       INTEGER :: omode
!!$       character(len=50) :: units
!!$
!!$       omode = NF90_SHARE
!!$       IF (PRESENT(mode)) omode = mode
!!$
!!$       CALL check (nf90_open(file_name, omode, id%file) )
!!$
!!$       CALL check( nf90_inq_dimid(id%file, 'time', id%dim_time) )
!!$
!!$
!!$       CALL check( nf90_inq_dimid(id%file, 'x', id%dim_x) )
!!$       CALL check( nf90_inq_dimid(id%file, 'y', id%dim_y) )
!!$       
!!$       CALL check( nf90_inq_varid(id%file, 'time', id%time) ) 
!!$       ! Get since, unit and step.
!!$       CALL check( nf90_get_att(id%file, id%time, 'units', units) )
!!$       units = ADJUSTR(units)
!!$       id%since = str2dt(units(31:50))
!!$       id%units = TRIM(ADJUSTL(units(:INDEX(units, 'since')-1)))
!!$
!!$       CALL check( nf90_get_att(id%file, id%time, 'step', id%step) )
!!$
!!$       CALL check( nf90_inq_varid(id%file, 'h', id%h) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'A', id%A) ) 
!!$       
!!$       CALL check( nf90_inq_varid(id%file, 'u', id%u) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'v', id%v) ) 
!!$
!!$       CALL check( nf90_inq_varid(id%file, 'Ti', id%Ti) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'Ta', id%Ta) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'To', id%To) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'Qsh_io', id%Qsh_io) ) 
!!$       CALL check( nf90_inq_varid(id%file, 'Qoa', id%Qoa) ) 
!!$       
!!$     END FUNCTION nc_load_ids
!!$
!!$
!!$     SUBROUTINE nc_store_parameters(ncid)
!!$       ! Set simulation parameters as attributes of a netCDF file.
!!$       !
!!$       ! :param ncid: ID of the netCDF file where those attributes will be stored. 
!!$
!!$       USE params
!!$       USE constants
!!$       USE dyndim
!!$       USE options
!!$       use thermodim
!!$
!!$       IMPLICIT NONE
!!$
!!$       INTEGER, INTENT(in) :: ncid
!!$
!!$       ! Params
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'ksuperloop', ksuperloop))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'ksor', ksor))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'klsor', klsor))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'wsor', wsor))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'wlsor', wlsor))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'beta', beta))
!!$
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'tthermo', tthermo))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'tdyn', tdyn))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'tdyn_h', tdyn_h))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'tadv', tadv))
!!$
!!$       
!!$       ! Constants
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'relhum', relhum))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'theta_a', theta_a))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'theta_w', theta_w))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Tif', Tif))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Tof', Tof))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'albedoo', albedoo))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'albedoi', albedoi))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'emisice', emisice))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'emisatml', emisatml))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'AMR', AMR))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'hmin', hmin))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Cpwater', Cpwater))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'rhowater', rhowater))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'rhoice', rhoice))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Lfusion', Lfusion))
!!$
!!$       ! Dyndim
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Cdw', Cdw))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Cda', Cda))
!!$
!!$       ! Options
!!$!       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Dynamic', DYNAMIC))
!!$!       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Thermodyn', Thermodyn))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'AirTemp', AirTemp))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'OcnTemp', OcnTemp))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Current', Current))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'BndyCond', BndyCond))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'DragLaw', DragLaw))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Rheology', Rheology))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'linearization', linearization))
!!$       
!!$       ! Thermodim 
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Klat_ia', Klat_ia))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Klat_oa', Klat_oa))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Ksens_ai', Ksens_ai))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Ksens_ao', Ksens_ao))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Ksens_io', Ksens_io))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kice', Kice))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kocn', Kocn))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kadvo', Kadvo))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kemis_i', Kemis_i))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kemis_al', Kemis_al))
!!$       CALL check( nf90_put_att(ncid, NF90_GLOBAL, 'Kemis_o', Kemis_o))
!!$      
!!$     END SUBROUTINE nc_store_parameters
!!$
!!$

!!$
!!$
!!$    SUBROUTINE ocn_Tclim
!!$
!!$!************************************************************************
!!$!     Subroutine load_data: load forcing data (wind stress, air temperature,
!!$!       water temperature and ocean current.
!!$!
!!$!     Revision History
!!$!     ----------------
!!$!
!!$!     Ver             Date (dd-mm-yy)        Author
!!$!
!!$!     V01             14-05-97               L.-B. Tremblay
!!$!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!!$!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!!$!
!!$!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!!$!     -------   Montreal, Quebec, Canada!!$!     Email   :  bruno.tremblay@mcgill.ca
!!$!
!!$!************************************************************************
!!$
!!$
!!$      USE PARAMS
!!$      USE OPTIONS
!!$      USE THERMOVARIABLES
!!$      USE THERMOFORCING
!!$      USE MASKS
!!$
!!$      implicit none
!!$
!!$
!!$      character(LEN=60) fname1
!!$
!!$      integer i, j,  kmo, mo
!!$
!!$      if ( Thermodyn ) then
!!$
!!$!------------------------------------------------------------------------
!!$!     load ocean temperature created by Tocn_clim_gen.m
!!$!------------------------------------------------------------------------
!!$
!!$      do kmo = 0, 13
!!$
!!$         mo = kmo
!!$
!!$         if ( kmo .eq. 0 ) then
!!$            mo = 12
!!$         elseif ( kmo .eq. 13 ) then
!!$            mo = 1
!!$         endif
!!$
!!$         fname1 = ''
!!$         
!!$         if ( OcnTemp .eq. 'MonthlyClim' .or.                     &
!!$                   OcnTemp .eq. 'calculated' ) then
!!$
!!$         fname1  = 'forcing/Tocn/monthlyclim/Tocn'
!!$   
!!$            write (fname1,'(a29,i2.2)') fname1, mo
!!$
!!$            open(unit = 30, file = fname1, status = 'unknown')
!!$         
!!$            do j = 0 , ny+1
!!$               read(30,*) ( To_clim(i,j,kmo), i = 0, nx+1 )
!!$            enddo
!!$         
!!$            close(30)
!!$
!!$         endif
!!$
!!$            
!!$         do j = 0, ny+1
!!$            do i = 0, nx+1
!!$
!!$               To_clim(i,j,kmo) = ( To_clim(i,j,kmo) + 273.15e0 ) &
!!$                                 * maskC(i,j)
!!$
!!$!     user specified ocean temperature  (non-dim)
!!$!     WARNING specify To not To_clim
!!$
!!$               if ( OcnTemp .eq. 'specified' ) then
!!$                  To_clim(i,j,kmo) = ( -1.8e0 + 273.15e0 ) * maskC(i,j)
!!$               endif
!!$                  
!!$            enddo
!!$         enddo
!!$         
!!$      enddo
!!$
!!$      endif
!!$
!!$      return
!!$    END SUBROUTINE ocn_Tclim
!!$      
!!$
!!$
!!$
!!$
!!$    SUBROUTINE load_wind (date) 
!!$!------------------------------------------------------------------------
!!$!     load wind data for a given date and time
!!$!------------------------------------------------------------------------
!!$
!!$
!!$      USE PARAMS
!!$      USE DYNFORCING
!!$      USE DYNVARIABLES
!!$      USE DYNDIM
!!$      USE OPTIONS
!!$      USE CONSTANTS
!!$      USE MASKS
!!$      USE datetime, ONLY: datetime_type, datetime_str
!!$
!!$
!!$      implicit none
!!$
!!$
!!$      TYPE(datetime_type), INTENT(in) :: date
!!$
!!$      character(LEN=6) date1, date2, datetemp
!!$      CHARACTER(LEN=60) fname1, fname2, var, path
!!$
!!$      integer i, j, year
!!$      integer time, time1, time2
!!$      real wm
!!$
!!$      if ( Wind .eq. '6hours' ) then
!!$         
!!$         path   = 'forcing/wind/6hours/'
!!$
!!$         IF ( tdyn_h .EQ. 6 ) THEN
!!$            
!!$            var    = 'uwgeo'                   
!!$
!!$            fname1 = ''
!!$            
!!$            write  ( fname1, '(a20,i4.4,"/",a5,a6,"_",i2.2,".dat")' ) &
!!$                     path,date%year,var, datetime_str(date), date%hour
!!$            print *, 'opening/reading ', fname1
!!$            open   ( unit = 20, file = fname1, status = 'old' )
!!$            
!!$
!!$
!!$            read(20,*)                                                &
!!$                 ( ( uair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
!!$
!!$            close(20)
!!$
!!$            var    = 'vwgeo'                   
!!$
!!$            fname2 = ''
!!$            
!!$            write  ( fname2, '(a20,i4.4,"/",a5,a6,"_",i2.2,".dat")' ) &
!!$                     path,date%year,var, datetime_str(date),date%hour
!!$
!!$            print *, 'opening/reading ', fname2
!!$            OPEN   ( unit = 21, file = fname2, status = 'old' )
!!$          
!!$          
!!$
!!$            read(21,*)                                            &
!!$                 ( ( vair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
!!$
!!$            close(21)
!!$         ENDIF
!!$
!!$
!!$      elseif ( Wind .eq. '60yrs_clim' ) then 
!!$
!!$      fname1 = 'forcing/wind/60yrs_clim/uwgeo_clim.dat'
!!$
!!$         open ( unit = 20, file = fname1, status = 'old' )
!!$            
!!$         print *, 'opening/reading ', fname1
!!$
!!$         read(20,*)                                               &
!!$              ( ( uair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
!!$         
!!$         close(20)
!!$
!!$      fname2 = 'forcing/wind/60yrs_clim/vwgeo_clim.dat'
!!$            
!!$         open ( unit = 21, file = fname2, status = 'old' )
!!$          
!!$         print *, 'opening/reading ', fname2
!!$
!!$         read(21,*)                                               &
!!$              ( ( vair(i,j), i = 1, nx+1 ), j = 1, ny+1 )
!!$
!!$         close(21)
!!$
!!$      elseif ( Wind .eq. 'specified' ) then
!!$
!!$         wm = 10.0e0
!!$
!!$         do i = 1, nx+1
!!$            do j = 1, ny+1
!!$
!!$               uair(i,j) = 10e0 
!!$               vair(i,j) =  0e0
!!$
!!$!               call random_number(rdnumb)
!!$!               uair(i,j) = wm * (rdnumb - 0.5e0)
!!$!               call random_number(rdnumb)
!!$!               vair(i,j) = wm * (rdnumb - 0.5e0)
!!$
!!$            enddo
!!$         enddo
!!$
!!$      endif
!!$    END SUBROUTINE load_wind
!!$
!!$    SUBROUTINE print_parameters(unit)
!!$      ! Print parameters in unit.
!!$      USE OPTIONS
!!$      USE RHEO
!!$      USE params
!!$
!!$      INTEGER, intent(in) :: unit
!!$      CHARACTER(len=30) :: rfmt, ifmt, sfmt, lfmt
!!$
!!$      rfmt = "(A20, ': ', es8.2)"
!!$      ifmt = "(A20, ': ', i0)"
!!$      sfmt = "(A20, ': ', A)"
!!$      lfmt = "(A20, ': ', l1)"
!!$
!!$      WRITE(unit,*) 
!!$      WRITE(unit,*) 'Parameters'
!!$      WRITE(unit,*) '----------'
!!$      WRITE(unit,sfmt) 'Rheology', Rheology
!!$      WRITE(unit,rfmt) 'Pstar', Pstar
!!$      WRITE(unit,rfmt) 'zetamin', zetamin
!!$      WRITE(unit,rfmt) 'ellipticity', 1./SQRT(ell_2)
!!$      WRITE(unit,sfmt) 'linearization', linearization
!!$      WRITE(unit,lfmt) 'Dynamic', Dynamic
!!$      WRITE(unit,lfmt) 'Thermodyn', Thermodyn
!!$      WRITE(unit,sfmt) 'Current', Current
!!$      WRITE(unit,sfmt) 'Wind', Wind
!!$      WRITE(unit,sfmt) 'AirTemp', AirTemp
!!$      WRITE(unit,sfmt) 'OcnTemp', OcnTemp
!!$      WRITE(unit,ifmt) 'Superloops', ksuperloop 
!!$      WRITE(unit,ifmt) 'Precond it.', klsor
!!$      WRITE(unit,rfmt) 'Precond relax.', wlsor
!!$      WRITE(unit,*) 
!!$      WRITE(unit,ifmt) 'thermo time step', int(tthermo/3600.)
!!$      WRITE(unit,ifmt) 'dynamic time step', tdyn_h
!!$      WRITE(unit,ifmt) 'advection time step', int(tadv/3600.)
!!$      WRITE(unit,*) 
!!$
!!$    END SUBROUTINE print_parameters
!!$
!!$

    SUBROUTINE check(status)
      INTEGER, INTENT ( in) :: status

      
      IF(status /= nf90_noerr) THEN
         PRINT *, TRIM(nf90_strerror(status))
         STOP 2
      END IF
    END SUBROUTINE check


  END MODULE IO
