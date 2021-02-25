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


      subroutine thermo_source_terms (date, htp, Atp)
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

      double precision, intent(in) :: htp(0:nx+1,0:ny+1), Atp(0:nx+1,0:ny+1)
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

      call HeatFluxes (htp, Atp)
      
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
               Source_A = Source_A + Atp(i,j) * Source_h  /     &
                          ( 2d0 * max( htp(i,j), tiny ) )


!------------------------------------------------------------------------
!     Update ice thickness and ice concentration ( after SA(h_t-1) )
!------------------------------------------------------------------------

            if ( maskC(i,j) == 1 ) then

               Sh(i,j)=Source_h
               SA(i,j)=Source_A
               
            else

               Sh(i,j)=0d0
               SA(i,j)=0d0

            endif
            
         enddo
      enddo

      return
    end subroutine thermo_source_terms

    subroutine  dh_dA_thermo (htp, Atp)

      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_ThermoVariables.h'
      include 'CB_options.h'

      integer i, j

      double precision, intent(inout) :: htp(0:nx+1,0:ny+1), Atp(0:nx+1,0:ny+1)

      do i = 1, nx
         do j = 1, ny

            if ( maskC(i,j) == 1 ) then

               if ( adv_scheme == 'upwind' .or. adv_scheme == 'upwindRK2' ) then

                  htp(i,j) = htp(i,j) + Deltat * Sh(i,j)
                  Atp(i,j) = Atp(i,j) + Deltat * SA(i,j)

               elseif ( adv_scheme == 'semilag' ) then

                  print *, 'NOT CODED YET'
                  stop

               endif

               htp(i,j) = max(htp(i,j), 0d0)
               Atp(i,j) = max(Atp(i,j), 0d0)
               Atp(i,j) = min(Atp(i,j), 1d0)
               
            else
               
               Sh(i,j)=0d0
               SA(i,j)=0d0

            endif

         enddo
      enddo

    end subroutine dh_dA_thermo
