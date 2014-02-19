!*************************************************************************
!     subroutine shortwave:
!       Calculates the daily average short wave solar energy input as a function of 
!       latitude, time and cloud cover after Zillman (1972). 
!
!
!       Q0 = S cos(Z)^2 / (cos(Z) + 2.7) e x 10^-5 + 1.085 cos(Z) + 0.1)
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


      subroutine shortwave (date)
        USE datetime, ONLY: datetime_type, datetime_str, datetime_str_6
        use datetime, only: julianday
      implicit none

      include 'parameter.h'

      include 'CB_const.h'
      include 'CB_ThermoForcing.h'
      include 'CB_ThermoVariables.h'
      
      TYPE(datetime_type), INTENT(in) :: date
      integer i, j, l, Jday, npts

      character(LEN=6) datestr
       
      double precision pi, deg2rad, tilt, decli, sindecli, cosdecli
      double precision factor, dt, sinsin, coscos, cosZ
      double precision albedoa, coalbedo

      pi       = 4d0 * atan( 1d0 )
      deg2rad  = pi / 180d0
      albedoa  = 0.26d0         ! % of reflected SW radiation (atm)
      coalbedo = (1d0 - albedoa)


      npts     = 24                ! # points for the avg SW rad
      tilt     = 23.44 * deg2rad   ! Earth axis tilt [rad]

      datestr = datetime_str_6(date)
      
      Jday = int(julianday(date))

!------------------------------------------------------------------------
!    compute the declination angle and cosine of the zenith angle
!------------------------------------------------------------------------

      decli    = tilt * cos((172d0 - Jday) * pi / 180d0)

      sindecli = sin(decli)
      cosdecli = cos(decli)

      factor = 1d0 / (npts*1d0)

      do j = 1, ny
         do i = 1, nx
            Qsw(i,j) = 0d0
         enddo
      enddo
      
      
      do l = 1, npts
         
         dt = 24d0 / (npts * 1d0)
         
         do j = 1, ny
            do i = 1, nx
               
!------------------------------------------------------------------------
!     cos (lat - decli) = sin(lat) sin(decli) + cos(lat) cos(decli)
!------------------------------------------------------------------------

               sinsin = sinlat(i,j) * sindecli
               coscos = coslat(i,j) * cosdecli
               
               cosZ   = max (sinsin + coscos * cos((12-l*dt) * pi/12d0), &
                         0d0)

               Qsw(i,j) = Qsw(i,j) + S0 * coalbedo* factor * cosZ ** 2 / &
                         ((cosZ + 2.7) * Pvap(i,j) * 1d-05 + &
                         1.085d0 * cosZ + 0.1)
              
            enddo
         enddo
         
      enddo
      
      return
      end
      

