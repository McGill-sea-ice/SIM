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

      subroutine oceanTadv ()

      implicit none

      include 'parameter.h'
      include 'CB_Thermodim.h'
      include 'CB_ThermoForcing.h'
      include 'CB_ThermoVariables.h'
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'


      integer i, j, ii, jj

      double precision Qnet

      double precision fluxxtot(nx+1,ny), fluxytot(nx,ny+1),             &
                advfluxx(nx+1,ny), advfluxy(nx,ny+1),                    &
                diffluxx(nx+1,ny), diffluxy(nx,ny+1)

!------------------------------------------------------------------------
!     Diffusion and advection heat fluxes
!------------------------------------------------------------------------

      do j = 1, ny
         do i = 1, nx+1

           
            diffluxx(i,j) = -Kocn * ( To(i,j) - To(i-1,j) ) / Deltax
            diffluxx(i,j) = diffluxx(i,j)                                &
                         * max( maskC(i,j) + maskC(i-1,j) - 1, 0 )


            ii       = i - max ( 0, int(sign(1d0, uwatnd(i,j)))) 

            advfluxx(i,j) = Kadvo * Hocn * To(ii,j) * uwatnd(i,j) 

            fluxxtot(i,j) = diffluxx(i,j) + advfluxx(i,j)


         enddo
      enddo

      do j = 1, ny+1
         do i = 1, nx

            diffluxy(i,j) = -Kocn * ( To(i,j) - To(i,j-1) ) / Deltax
            diffluxy(i,j) = diffluxy(i,j)                                &
                         * max(maskC(i,j) + maskC(i,j-1) - 1, 0)
      
            jj       = j - max ( 0, int( sign( 1d0, vwatnd(i,j) ) ) ) 

            advfluxy(i,j) = Kadvo * Hocn * To(i,jj)  * vwatnd(i,j) 

            fluxytot(i,j) = diffluxy(i,j) + advfluxy(i,j)
            
            
         enddo
      enddo

!------------------------------------------------------------------------
!     New ocean temperature due to advection and diffusion
!------------------------------------------------------------------------

      do j = 1, ny
         do i = 1, nx
            
            Qadvdiff(i,j) = ( fluxxtot(i+1,j) - fluxxtot(i,j) ) / Deltax &
                          + ( fluxytot(i,j+1) - fluxytot(i,j) ) / Deltax

            Qnet          = - Qadvdiff(i,j)

            if (OcnTemp .eq. 'calculated') then

               To(i,j) = To(i,j) + Deltat * Qnet / (Kadvo * Hocn)
!               To(i,j) = max( To(i,j), Tof ) * maskC(i,j) 

            endif
            

         enddo
      enddo
     
      return
      end







