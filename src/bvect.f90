! calculates the complete b vector for forcing including air stress + sstilt 
! (R1, R2) with pressure gradient (bu_ind, bv_ind) and terms related to the 
! water drag

      subroutine bvect (utp,vtp,rhs)

      use solver_choice
      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'

      integer i, j

      double precision &
                       uwavg(0:nx+2,0:ny+2), & ! uw evaluated at the v-loc
                       vwavg(0:nx+2,0:ny+2)    ! vw evaluated at the u-loc

      double precision speed1p, speed2p, rhs(nvar) 
      double precision uavg, vavg
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)


      do j = 1, ny
         do i = 1, nx+1

            if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

               vwavg(i,j) = ( vwatnd(i,j)   + vwatnd(i,j+1) &
                            + vwatnd(i-1,j) + vwatnd(i-1,j+1) ) / 4d0
            endif

         enddo
      enddo

      do j = 1, ny+1
         do i = 1, nx

            if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then
               
               uwavg(i,j) = ( uwatnd(i,j)   + uwatnd(i+1,j) &
                            + uwatnd(i,j-1) + uwatnd(i+1,j-1) ) / 4d0
            endif

         enddo
      enddo


       do j = 1, ny
         do i = 1, nx+1

            if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

            vavg = ( vtp(i,j)   + vtp(i,j+1) &
                   + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

            speed1p = sqrt( ( utp(i,j) - uwatnd(i,j) ) ** 2 &
                              + ( vavg - vwavg(i,j)  ) ** 2 )


            CdwC1(i,j) = max(Cdw * speed1p, 1d-10)

            endif

         enddo
      enddo


       do j = 1, ny+1
         do i = 1, nx

            if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then

            uavg = ( utp(i,j)   + utp(i+1,j) &
                   + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

            speed2p = sqrt( ( uavg - uwavg(i,j)  ) ** 2 &
                              + ( vtp(i,j) - vwatnd(i,j) ) ** 2 )


            CdwC2(i,j) = max(Cdw * speed2p, 1d-10)

            endif

         enddo
      enddo


      do j = 1, ny
         do i = 1, nx+1

            if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

            if (solver .le. 2) then ! Picard or JFNK
               
               bu(i,j) = bu_ind(i,j) - ( P(i,j) - P(i-1,j) ) / Deltax + & ! P is the replacement pressure 
                         CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                         vwavg(i,j)  * sintheta_w   )

            elseif (solver .eq. 3) then ! EVP solver
 
               bu(i,j) = R1(i,j) + &
                         CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                         vwavg(i,j)  * sintheta_w   )

            endif

            else

               bu(i,j) = 0d0

            endif

         enddo
      enddo


      do j = 1, ny+1
         do i = 1, nx
            
            if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then

            if (solver .le. 2) then ! Picard or JFNK

               bv(i,j) = bv_ind(i,j) - ( P(i,j) - P(i,j-1) ) / Deltax + & ! P is the replacement pressure
                         CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                         uwavg(i,j)  * sintheta_w   )

            elseif (solver .eq. 3) then ! EVP solver
               
               bv(i,j) = R2(i,j) + &
                         CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                         uwavg(i,j)  * sintheta_w   )
            
            endif

            else

               bv(i,j) = 0d0
            
            endif

         enddo
      enddo

!-----------------------------------------------------------------------------
!   Set bu and bv to 0.0 at the 'appropriate' open bc
!-----------------------------------------------------------------------------

      do j = 1, ny+1

         bu(1,j)    = 0.0d0   ! Bering Strait
         bu(nx+1,j) = 0.0d0   !

         bv(nx+1,j)    = 0.0d0 

      enddo

      do i = 1, nx+1

         bv(i,1)    = 0.0d0   ! North Atlantic
         bv(i,ny+1) = 0.0d0   ! no open bc in current configuration
            
         bu(i,ny+1)    = 0.0d0 

      enddo

      call transformer (bu,bv,rhs,1)

      return
    end subroutine bvect


!-----------------------------------------------------------------------------
! calculates R1 + component of pressure gradient + tendency term un1 part                                                                         
    subroutine bvect_ind

      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_const.h'
      include 'CB_DynForcing.h'
      include 'CB_options.h'

      integer i, j
      double precision hvert

      do j = 1, ny
         do i = 1, nx+1

            hvert = ( h(i,j) + h(i-1,j) ) / 2d0

            bu_ind(i,j) = R1(i,j) - &
                          rhof * hvert * ( vwater(i,j) + vwater(i,j+1) ) / 2d0

            if ( BDF .eq. 0 ) then
               bu_ind(i,j) = bu_ind(i,j) + (rhoice * hvert * un1(i,j)) / Deltat
            elseif ( BDF .eq. 1 ) then
               bu_ind(i,j) = bu_ind(i,j) + &
                    (rhoice * hvert * ( 2d0*un1(i,j) - 0.5d0*un2(i,j) ) ) / Deltat
            endif

         enddo
      enddo

      do j = 1, ny+1
         do i = 1, nx

            hvert = ( h(i,j) + h(i,j-1) ) / 2d0

            bv_ind(i,j) = R2(i,j) + &
                          rhof * hvert * ( uwater(i,j) + uwater(i+1,j) ) / 2d0

            if ( BDF .eq. 0 ) then
               bv_ind(i,j) = bv_ind(i,j) + (rhoice * hvert * vn1(i,j)) / Deltat
            elseif ( BDF .eq. 1 ) then
               bv_ind(i,j) = bv_ind(i,j) + &
                    (rhoice * hvert * ( 2d0*vn1(i,j) - 0.5d0*vn2(i,j) ) ) / Deltat
            endif

         enddo
      enddo

      return
    end subroutine bvect_ind
