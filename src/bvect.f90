! calculates the complete b vector for forcing including air stress + sstilt 
! (R1, R2) with pressure gradient (bu_ind, bv_ind) and terms related to the 
! water drag

! stress at BC option: to simplify things, bvect_ind is not modified. All the 
! changes are made in bvect. 

      subroutine bvect (utp,vtp,rhs)

      use solver_choice
      use basal_param

      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'
      include 'CB_bathymetry.h'
      include 'CB_options.h'
      include 'CB_stressBC.h'

      integer i, j

      double precision &
                       uwavg(0:nx+2,0:ny+2), & ! uw evaluated at the v-loc
                       vwavg(0:nx+2,0:ny+2)    ! vw evaluated at the u-loc

      double precision speed1p, speed2p, rhs(nvar) 
      double precision uavg, vavg
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision A_at_u, Ae, bathy_at_u, CBfactor
      double precision h_at_u, hc, minA, h1,h2, A1,A2, alpha
      double precision v1,v2,va,vb,vc,vd,v_at_u, u_at_v

      minA=0.01d0
      alpha = 1d06

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


       do j = 1, ny ! u comp
         do i = 1, nx+1

            if ( maskB(i,j) + maskB(i,j+1) .gt. 0 ) then

               vavg = ( vtp(i,j)   + vtp(i,j+1) &
                    + vtp(i-1,j) + vtp(i-1,j+1) ) / 4d0

               if (stressBC .and. i .eq. 1) then

                vavg = ( vtp(i,j) + vtp(i,j+1) ) / 2d0

               endif

               if (stressBC .and. i .eq. nx+1) then

                vavg = ( vtp(i-1,j) + vtp(i-1,j+1) ) / 2d0

               endif

               speed1p = sqrt( ( utp(i,j) - uwatnd(i,j) ) ** 2 &
                       + ( vavg - vwavg(i,j)  ) ** 2 )
            
!               CdwC1(i,j) = max(Cdw * speed1p, 1d-10)
               CdwC1(i,j) = Cdw * speed1p
            
	      if (BasalStress) then
	      
               bathy_at_u = min(bathy(i-1,j), bathy(i,j))

               h1=h(i-1,j)
               h2=h(i,j)
               A1=A(i-1,j)
               A2=A(i,j)

               h_at_u = (h1+h2)/2d0 + (h1/2d0)*tanh(alpha*(h1-h2)) + & 
                                      (h2/2d0)*tanh(alpha*(h2-h1))

               A_at_u = (A1+A2)/2d0 + (A1/2d0)*tanh(alpha*(A1-A2)) + &
                                      (A2/2d0)*tanh(alpha*(A2-A1))

               if ( A_at_u .gt. minA ) then 
                  hc = ( A_at_u * bathy_at_u ) / k1
               else
                  hc = 10000d0
               endif
                                    
               Cbasal1(i,j)=0d0

               if (h_at_u .gt. hc) then

                  va=abs(vtp(i-1,j+1))
                  vb=abs(vtp(i,j+1))
                  vc=abs(vtp(i-1,j))
                  vd=abs(vtp(i,j))

                  v1 = (va+vb)/2d0 - (va/2d0)*tanh(alpha*(va-vb)) - &
                       (vb/2d0)*tanh(alpha*(vb-va))

                  v2 = (vc+vd)/2d0 - (vc/2d0)*tanh(alpha*(vc-vd)) - &
                       (vd/2d0)*tanh(alpha*(vd-vc))
                  
                  v_at_u = (v1+v2)/2d0 - (v1/2d0)*tanh(alpha*(v1-v2)) - &
                       (v2/2d0)*tanh(alpha*(v2-v1))

                  speed1p=sqrt( utp(i,j)**2 + v_at_u**2 )

                  Cbfactor=k2/(speed1p+umin)
!                     Cbasal1(i,j) = Cbfactor * (h_at_u -hc) * dexp(-CC * (1d0 - A_at_u))
                  Cbasal1(i,j)=Cbfactor * (h_at_u - hc) * dexp(-CC * (1d0 - A_at_u))
               else
                  Cbasal1(i,j)=0d0
               endif
               
               endif
               
            endif
            
         enddo
      enddo

      if (stressBC .and. clipping) bu(nx+1,1)=-2.76d-09

       do j = 1, ny+1 ! v comp
         do i = 1, nx

            if ( maskB(i,j) + maskB(i+1,j) .gt. 0 ) then

               uavg = ( utp(i,j)   + utp(i+1,j) &
                    + utp(i,j-1) + utp(i+1,j-1) ) / 4d0

               if (stressBC .and. j .eq. 1) then

                uavg = ( utp(i,j) + utp(i+1,j) ) / 2d0

               endif

               if (stressBC .and. j .eq. ny+1) then

                uavg = ( utp(i,j-1) + utp(i+1,j-1) ) / 2d0

               endif

               speed2p = sqrt( ( uavg - uwavg(i,j)  ) ** 2 &
                       + ( vtp(i,j) - vwatnd(i,j) ) ** 2 )

!               CdwC2(i,j) = max(Cdw * speed2p, 1d-10)
               CdwC2(i,j) = Cdw * speed2p

              if (BasalStress) then
               
               bathy_at_u = min(bathy(i,j), bathy(i,j-1)) ! in fact at v
 
               h1=h(i,j-1)
               h2=h(i,j)
               A1=A(i,j-1)
               A2=A(i,j)

               h_at_u = (h1+h2)/2d0 + (h1/2d0)*tanh(alpha*(h1-h2)) + &
                                      (h2/2d0)*tanh(alpha*(h2-h1))

               A_at_u = (A1+A2)/2d0 + (A1/2d0)*tanh(alpha*(A1-A2)) + &
                                      (A2/2d0)*tanh(alpha*(A2-A1))

               if ( A_at_u .gt. minA ) then
                  hc = ( A_at_u * bathy_at_u ) / k1
               else
                  hc = 10000d0
               endif

               Cbasal2(i,j)=0d0

               if (h_at_u .gt. hc) then

                  va=abs(utp(i,j))
                  vb=abs(utp(i,j-1))
                  vc=abs(utp(i+1,j))
                  vd=abs(utp(i+1,j-1))

                  v1 = (va+vb)/2d0 - (va/2d0)*tanh(alpha*(va-vb)) - &
                       (vb/2d0)*tanh(alpha*(vb-va))

                  v2 = (vc+vd)/2d0 - (vc/2d0)*tanh(alpha*(vc-vd)) - &
                       (vd/2d0)*tanh(alpha*(vd-vc))

                  u_at_v = (v1+v2)/2d0 - (v1/2d0)*tanh(alpha*(v1-v2)) - &
                       (v2/2d0)*tanh(alpha*(v2-v1))

                  speed2p=sqrt( u_at_v**2 + vtp(i,j)**2 )

                  Cbfactor=k2/(speed2p+umin)
                  !                  Cbasal2(i,j) = Cbfactor * (h_at_u -hc) * dexp(-CC * (1d0 - A_at_u))
                  Cbasal2(i,j)=Cbfactor * (h_at_u - hc) * dexp(-CC * (1d0 - A_at_u))
               else
                  Cbasal2(i,j)=0d0
               endif
               
               endif
               
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

             if (stressBC) then
             
              if (i .eq. 1) then ! West

               bu(i,j)= bu_ind(i,j) - ( P(i,j) / Deltaxh ) + & 
                        CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                        vwavg(i,j)  * sintheta_w   ) - sigmaW(j)/ Deltaxh + &
                        ( tauW(j+1)-tauW(j) ) / Deltax
              endif

              if (i .eq. nx+1) then ! East

              bu(i,j) = bu_ind(i,j) + ( P(i-1,j) / Deltaxh ) + & 
                        CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                        vwavg(i,j)  * sintheta_w   ) + sigmaE(j)/ Deltaxh + &
                        ( tauE(j+1)-tauE(j) ) / Deltax
            
              endif

              if (j .eq. 1 .and. i .ne. 1 .and. i .ne. nx+1) then ! South

              bu(i,j) = bu_ind(i,j) - ( P(i,j) - P(i-1,j) ) / Deltax + & 
                        CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                        vwavg(i,j)  * sintheta_w   ) - tauS(i) / Deltax   

              endif

              if (j .eq. ny .and. i .ne. 1 .and. i .ne. nx+1) then ! North

              bu(i,j) = bu_ind(i,j) - ( P(i,j) - P(i-1,j) ) / Deltax + & 
                        CdwC1(i,j) * ( uwatnd(i,j) * costheta_w - &
                        vwavg(i,j)  * sintheta_w   ) + tauN(i) / Deltax   

              endif

             endif

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

             if (stressBC) then             

              if (j .eq. 1) then ! South
       
               bv(i,j) = bv_ind(i,j) - ( P(i,j) ) / Deltaxh + & 
                        CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                        uwavg(i,j)  * sintheta_w   ) - sigmaS(i)/ Deltaxh + &
                        ( tauS(i+1)-tauS(i) ) / Deltax

              endif

             if (j .eq. ny+1) then ! North

              bv(i,j) = bv_ind(i,j) + ( P(i,j-1) ) / Deltaxh + & 
                        CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                        uwavg(i,j)  * sintheta_w   ) + sigmaN(i)/ Deltaxh + &
                        ( tauN(i+1)-tauN(i) ) / Deltax

             endif

             if (i .eq. 1 .and. j .ne. 1 .and. j .ne. ny+1) then ! West

              bv(i,j) = bv_ind(i,j) - ( P(i,j) - P(i,j-1) ) / Deltax + & 
                        CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                        uwavg(i,j)  * sintheta_w   ) - tauW(j) / Deltax

             endif

             if (i .eq. nx .and. j .ne. 1 .and. j .ne. ny+1) then ! East

              bv(i,j) = bv_ind(i,j) - ( P(i,j) - P(i,j-1) ) / Deltax + & 
                        CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                        uwavg(i,j)  * sintheta_w   ) + tauE(j) / Deltax

             endif

            elseif (solver .eq. 3) then ! EVP solver
               
               bv(i,j) = R2(i,j) + &
                         CdwC2(i,j) * ( vwatnd(i,j) * costheta_w + &
                         uwavg(i,j)  * sintheta_w   )
            
             endif
            endif

            else

               bv(i,j) = 0d0
            
            endif

         enddo
      enddo

!-----------------------------------------------------------------------------
!   Set bu and bv to 0.0 at the 'appropriate' open bc
!-----------------------------------------------------------------------------
     if (.not. stressBC) then
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
     endif

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
