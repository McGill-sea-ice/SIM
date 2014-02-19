!*************************************************************************
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

      subroutine oceanTthermo 

      implicit none

      include 'parameter.h'
      include 'CB_Thermodim.h'
      include 'CB_ThermoVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j

      double precision relvelocity, Qnet



      do j = 1, ny
         do i = 1, nx
            
!------------------------------------------------------------------------
!     Net heat gain by the mixed layer
!------------------------------------------------------------------------
            
! ask Bruno...speediw is 0.0...therefore 0.005 is always chosen...should we change this value?

            relvelocity   = max( speediw(i,j), 0.005d0 )

            Qnet          = -Qoa(i,j)  - Qsh_io(i,j) 

!------------------------------------------------------------------------
!     New ocean temperature due to thermo
!------------------------------------------------------------------------


            if (OcnTemp .eq. 'calculated') then

               To(i,j) = To(i,j) + Deltat * Qnet / Kadvo
               To(i,j) = max( To(i,j), Tof ) * maskC(i,j) 

            endif
            

         enddo
      enddo
     
      return
    end subroutine oceanTthermo







