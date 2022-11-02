!************************************************************************
!     Subroutine atmosphere: 
!       Calculates the heat transfers between the surface and the atmosphere
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

      subroutine HeatFluxes (htp, Atp)

      use ICE_ALBEDO
      implicit none
      include 'parameter.h'
      include 'CB_ThermoForcing.h'
      include 'CB_ThermoVariables.h'
      include 'CB_Thermodim.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'

      integer i, j

      double precision, intent(in) :: htp(0:nx+2,0:ny+2), Atp(0:nx+2,0:ny+2)

      double precision tiny, tiny2, absorpatm1, absorpatm2, absorpatm
      double precision qsi, qso, qa, Bi, Bi1, Ai1
      double precision hmax, delT, relvelocity
      double precision esat
      external esat

      tiny       = 1d-8
      tiny2      = 2.73d-6    !jfl to be consistent with the nodim model

      absorpatm1 = 0.25d0     ! atm absorptivity to SW radiation (A = 1)
      absorpatm2 = 0.50d0     ! atm absorptivity to SW radiation (A = 0)
      
!------------------------------------------------------------------------
!     Calculate the ice-atmosphere (Qia = Qia), ocean-atmosphere (Qoa) 
!       ocean-ice heat fluxes, and ice surface temperature (Ti).
!     Heat fluxes are positive upward.
!     To at the previous time step is used to calculate 
!       ice/ocean-atmosphere heat fluxes.
!     The ice base temperature is equal to Tof.
!------------------------------------------------------------------------


      do i = 1, nx
         do j = 1, ny

            if ( maskC(i,j) .eq. 0 ) goto 10

            absorpatm= -(absorpatm2 - absorpatm1) * Atp(i,j) + absorpatm2

            qsi      = esat (Ti(i,j), 1d0) 
            qso      = esat (To(i,j), 0d0)
            qa       = esat (Ta(i,j), Atp(i,j)) * relhum

!------------------------------------------------------------------------
!     Coefficients of the linearized total heat flux between the ice and
!     the atmosphere ( Qia = Ai Ti + Bi )
!------------------------------------------------------------------------

            Bi =  -Ksens_ai * speeda(i,j) * Ta(i,j)                     &
                       - Kemis_al * Ta(i,j) ** 4                        &
                       + max(Klat_ia * speeda(i,j)                      &
                            * ( qsi - qa ), 0d0 )                       &
                       - Qsw(i,j) * (1d0 - absorpatm) *                 &
                       (1d0 - albedo_manabe(Ti(i,j), htp(i,j)))
            
! jfl: Ti is the old value and we are solving for delta T         

            Bi1 = Bi + Ksens_ai * speeda(i,j) * Ti(i,j)                 &
                     + Kemis_i * Ti(i,j) ** 4

            Ai1 = Ksens_ai * speeda(i,j) + 4d0 * Kemis_i * Ti(i,j) ** 3 
 
!------------------------------------------------------------------------
!     From the equation Qia(Ti,Tib,To) = Qia(Ti,Ta), calculate Ti 
!
!     if Ti > Tif 
!        Ti = Tif
!      -(Qia - Qcond), melting at the top surface
!       -Qcond,        melting at the bottom surface
!        Qia,          total melting on both surfaces
!     if Ti < Tif
!        Qia,          freezing of ice at the bottom surface
!        
!------------------------------------------------------------------------


            hmax = max(htp(i,j) / max(Atp(i,j), tiny), tiny)

            delT=(-Bi1 - Kice * (Ti(i,j)-Tof) / hmax )/(Kice/hmax + Ai1)

            Ti(i,j) = min(Ti(i,j)+delT, Tif)

! net upward flux over ice

            Qia(i,j) = ( Ksens_ai * speeda(i,j) * Ti(i,j)               &
                         + Kemis_i * Ti(i,j) ** 4 + Bi ) * Atp(i,j)

!------------------------------------------------------------------------
!     Ocean-atmosphere heat flux, leads to ice formation if To = Tof(S)
!------------------------------------------------------------------------

! net upward flux over water

            Qoa(i,j) = ( Kemis_o * To(i,j) ** 4                         &
                          - Kemis_al * Ta(i,j) ** 4                     &
                          - Qsw(i,j) * (1d0 - albedoo) *                &
                            (1d0 - absorpatm)                           &
                          + max(Klat_oa * speeda(i,j)                   &
                                 * ( qso - qa ), 0d0)                   &
                          - Ksens_ao * speeda(i,j)                      &
                                 * ( Ta(i,j) - To(i,j) ) )              &
                          * (1d0 - Atp(i,j))

            Qoa_f(i,j) = 0d0

            if ( abs( To(i,j) - Tof) .lt. tiny2 .and.                   &
                      Qoa(i,j) .gt. 0d0) then
               Qoa_f(i,j) = Qoa(i,j)
               Qoa(i,j) = 0d0

            endif
            
!     Now form ice whenever there is heat loss from the leads even though 
!     the ocean temperature... when the heat lost is negative (i.e. from 
!     atm to ocn) then the heat flux should warm up the ocean... CHANGE IT BACK! 
!     IT ONLY MADE SENSE TO HAVE THAT WHEN WE HAD FOCN... NOW I KILLED THAT
!     BECAUSE IT WAS TOO AWKWARD... AND SO I NEED TO BRING THIS THING BACK AS 
!     WELL, IF NOT, QOA IN OPEN WATERS FORMS ICE, THE WATER GETS WARM BUT 
!     THE QSH_IO WILL INLY PREVENT THE H TO GROW SLOWER, NOT THE AICE!! AND SO 
!     WHEN RUNNING WITH OVERTURN YOU END UP WITH THIN ICE, BUT THE AICE IS HUGE
!     IN THE NORTHERN NORTH ATLANTIC.!

!            if ( Qoa(i,j) .gt. 0d0) then
!               Qoa_f(i,j) = Qoa(i,j)
!               Qoa(i,j) = 0d0

!            endif


!------------------------------------------------------------------------
!     Ocean-ice sensible heat flux
!------------------------------------------------------------------------

! ask Bruno...where do we calculate speediw? WE NEED TO CALCULATE JFL

! change 0.005...0.0005

            relvelocity  = max(speediw(i,j), 0.005d0)

            Qsh_io(i,j)  = max (Ksens_io * relvelocity *                &
                      (To(i,j) - Tof), 0d0) * Atp(i,j)

!------------------------------------------------------------------------
!     Define various heat fluxes for plotting purposes
! jfl modify these they are not good since it is with dims!!!
!------------------------------------------------------------------------
!!$
!!$            Qsws(i,j) = S0 * (1d0 - albedoo) * (1d0 - absorpatm) *      &
!!$                        Qsw(i,j) * (1d0 - A(i,j)) +                     &
!!$                        S0 * (1d0 - albedo_manabe(Ti(i,j), h(i,j))) *    &
!!$                        (1d0 - absorpatm) * Qsw(i,j) * A(i,j)
!!$
!!$            Qlwup(i,j)= emisocn * To(i,j) ** 4 * (1d0 - A(i,j)) +       &
!!$                        emisice * Ti(i,j) ** 4 * A(i,j)
!!$            Qlwdo(i,j)= emisatml * Ta(i,j) ** 4
!!$
!!$            Qsh(i,j)  =-Ksens_ao * speeda(i,j) *                        &
!!$                        ( Ta(i,j) - To(i,j) )  * (1d0 - A(i,j)) -       &
!!$                        Ksens_ai * speeda(i,j) *                        &
!!$                        ( Ta(i,j) - Ti(i,j) ) * A(i,j)
!!$            Qlh(i,j)  = max(-Klat_oa * speeda(i,j) * ( qa - qso ), 0d0) &
!!$                        * (1d0 - A(i,j)) +                              &
!!$                        max(-Klat_ia * speeda(i,j) * ( qa - qsi ), 0d0) &
!!$                        * A(i,j)

 10         continue


         enddo
      enddo


      return
    end subroutine HeatFluxes

!************************************************************************
!     function Tofreeze:
!       Freezing point temperature as a function of salinity.
!       Assumed linear between (T,S) = (0,0) and (-1.8C, 35ppt)
!************************************************************************


!      function Tofreeze (S)

!      include 'parameter.f'

!      Tofreeze = -1.8d0 * S / 35 + 273.15d0
!      Tofreeze = Tofreeze / Temp

!      return
!      end

