!*************************************************************************
!     subroutine thermodynamic: 
!       computes the the thermodynamic ice growth 
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             29-07-97               L.-B. Tremblay
!     V2.0            16-10-06               L.-B. Tremblay & JF Lemieux
!     V3.0            30-01-08               JF Lemieux & L.-B. Tremblay
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!************************************************************************


      subroutine thermodynamic (date)
        USE datetime, ONLY: datetime_type, datetime_str, datetime_str_6
      implicit none

      include 'parameter.h'
      include 'CB_Dyndim.h'
      include 'CB_Thermodim.h'
      include 'CB_DynVariables.h'
      include 'CB_ThermoVariables.h'
      include 'CB_mask.h'
      include 'CB_const.h'

      TYPE(datetime_type), INTENT(in) :: date
      integer day, month, i, j
      
      character(LEN=6) datestr

      double precision tiny, Qtot, Source_h, Source_A

      tiny = 1d-07   !jfl to be consistent with the nodim model

      datestr = datetime_str_6(date)
      day = date%day
      month = date%month

      
      
      ! The interpolation of the ocean temperature in bc_get has daily accuracy, 
      ! but this only affects boundary conditions. 
      call bc_get (date)

      
      call atmosphere (date)

      call oceanTadv ! advection of ocean T
     
      call oceanTthermo

      call HeatFluxes
      

      do i = 1, nx
         do j = 1, ny
                     

!------------------------------------------------------------------------
!     Thermodynamic source term of ice thickness
!------------------------------------------------------------------------

! the concentration is taken into account in Heatfluxes.f

! Bruno, the if on source_h looks useless. Can it be removed ?

            Qtot = Qia(i,j) - Qsh_io(i,j) + Qoa_f(i,j) 

            if ( Qtot .gt. 0d0 .and.                                 &
                      (month .eq. 12 .or. month .eq. 1 .or.          &
                       month .eq. 2  .or. month .eq. 3 .or.          &
                       month .eq. 4  .or. month .eq. 5 .or.          &
                       month .eq. 6) ) then 


            Source_h = (Qia(i,j) + Qoa_f(i,j) - Qsh_io(i,j) )/(rhoice * Lfusion)


            else
               
               Source_h = ( Qia(i,j) + Qoa_f(i,j) - Qsh_io(i,j) )/   &
                          (rhoice * Lfusion)
            endif

!------------------------------------------------------------------------
!     Thermodynamic source term of ice concentration
!------------------------------------------------------------------------

            if ( Qtot .gt. 0d0 .and.                                 &
                      (month .eq. 12 .or. month .eq. 1 .or.          &
                       month .eq. 2  .or. month .eq. 3 .or.          &
                       month .eq. 4  .or. month .eq. 5 .or.          &
                       month .eq. 6) ) then 


               Source_A =  Qoa_f(i,j) / (hmin * rhoice * Lfusion)

            else

               Source_A = Qoa_f(i,j) / (hmin * rhoice * Lfusion)

            endif


            if ( Source_h .lt. 0d0 )                                 &
               Source_A = Source_A + tracer(i,j,2) * Source_h  /     &
                          ( 2d0 * max( tracer(i,j,1), tiny ) )


!------------------------------------------------------------------------
!     Update ice thickness and ice concentration ( after SA(h_t-1) )
!------------------------------------------------------------------------


            tracer(i,j,1)  = tracer(i,j,1) + Deltat * Source_h 

            tracer(i,j,1)  = max( tracer(i,j,1), 0d0 ) * maskC(i,j)

            tracer(i,j,2)  = tracer(i,j,2) + Deltat * Source_A

            tracer(i,j,2)  = max( min( tracer(i,j,2), 1d0 ), 0d0 ) * &
                             maskC(i,j)
            
         enddo
      enddo


      return
      end



