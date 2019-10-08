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

      double precision, intent(in) :: htp(0:nx+2,0:ny+2), Atp(0:nx+2,0:ny+2)
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
      include 'CB_semilag.h'

      integer i, j, caseSL, iloc, jloc
      integer isw, jsw, inw, jnw, ine, jne, ise, jse !SL SouthWest=sw, nw, ne, se corners

      double precision, intent(inout) :: htp(0:nx+2,0:ny+2), Atp(0:nx+2,0:ny+2)
      double precision                :: xdn1, ydn1, rhsthermoh, rhsthermoA
      double precision                :: fsw, fnw, fne, fse
      double precision                :: fxsw, fxnw, fxne, fxse, fysw, fynw, fyne, fyse
      double precision                :: fxysw, fxynw, fxyne, fxyse
      double precision                :: fx, fy, fxy, cubic_interp         ! functions

      if ( adv_scheme == 'upwind' .or. adv_scheme == 'upwindRK2' ) then

         do i = 1, nx
            do j = 1, ny

               if ( maskC(i,j) == 1 ) then
                  htp(i,j) = htp(i,j) + Deltat * Sh(i,j)
                  Atp(i,j) = Atp(i,j) + Deltat * SA(i,j)
                  htp(i,j) = max(htp(i,j), 0d0)
                  Atp(i,j) = max(Atp(i,j), 0d0)
                  Atp(i,j) = min(Atp(i,j), 1d0)
               endif

            enddo
         enddo

      elseif ( adv_scheme == 'semilag' ) then

         do i = 1, nx
            do j = 1, ny   

               caseSL=3
               if (maskC(i,j) .eq. 1) then
                  caseSL=1 ! SL                                                                                           
                  if (i .lt. 3 .or. i .gt. nx-2 .or. j .lt. 3 .or. j .gt. ny-2) then
                     caseSL=2 ! upwind                                                                                    
                  else
                     do jloc=j-2, j+2
                        do iloc=i-2, i+2
                           if (maskC(iloc,jloc) .eq. 0) caseSL=2 ! upwind                                                 
                        enddo
                     enddo
                  endif
               endif

               if (caseSL == 1) then ! SL advection

! identify coordinates of 4 corners. These corners are 4 tracer points.                                                     

! Tnw--------Tne                                                                                                            
!  |         |                                                                                                              
!  |         |                                                                                                              
!  |         |                                                                                                              
! Tsw--------Tse                                                                                                            

! xd and yd are distances in the interval [0,1]. They are calc from the sw corner                                           
! xdn1, ydn1 is the position of the particle at n-1 (with respect to the sw corner)                                         

                  if (alpmx(i,j) .ge. 0) then  ! particle coming from the West (u .ge. 0)                                      
                     xdn1 = 1d0 - alpmx(i,j) / Deltax
                     if (alpmy(i,j) .ge. 0) then ! particle coming from the South (v .ge. 0)                                   
                        isw=i-1
                        jsw=j-1
                        inw=i-1
                        jnw=j
                        ine=i
                        jne=j
                        ise=i
                        jse=j-1
                        ydn1 = 1d0 - alpmy(i,j) / Deltax
                     else                     ! particle coming from the North (v .lt. 0)      
                        isw=i-1
                        jsw=j
                        inw=i-1
                        jnw=j+1
                        ine=i
                        jne=j+1
                        ise=i
                        jse=j
                        ydn1 = -1d0*alpmy(i,j) / Deltax
                     endif
                  else                      ! particle coming from the East (u .lt. 0)                                      
                     xdn1 = -1d0*alpmx(i,j) / Deltax
                     if (alpmy(i,j) .ge. 0) then ! particle coming from the South (v .ge. 0)                                   
                        isw=i
                        jsw=j-1
                        inw=i
                        jnw=j
                        ine=i+1
                        jne=j
                        ise=i+1
                        jse=j-1
                        ydn1 = 1d0 - alpmy(i,j) / Deltax
                     else                     ! particle coming from the North (v .lt. 0)                                   
                        isw=i
                        jsw=j
                        inw=i
                        jnw=j+1
                        ine=i+1
                        jne=j+1
                        ise=i+1
                        jse=j
                        ydn1 = -1d0*alpmy(i,j) / Deltax
                     endif
                  endif

! find rhsthermoh using cubic interpolation
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation                                

                  fsw=Sh(isw,jsw)
                  fnw=Sh(inw,jnw)
                  fne=Sh(ine,jne)
                  fse=Sh(ise,jse)
                  fxsw=fx(fse, Sh(isw-1,jsw), 2d0)
                  fxnw=fx(fne, Sh(inw-1,jnw), 2d0)
                  fxne=fx(Sh(ine+1,jne), fnw, 2d0)
                  fxse=fx(Sh(ise+1,jse), fsw, 2d0)
                  fysw=fy(fnw, Sh(isw,jsw-1), 2d0)
                  fynw=fy(Sh(inw,jnw+1), fsw, 2d0)
                  fyne=fy(Sh(ine,jne+1), fse, 2d0)
                  fyse=fy(fne, Sh(ise,jse-1), 2d0)
                  fxysw=fxy(fne, Sh(inw-1,jnw), &
                            Sh(ise,jse-1),Sh(isw-1,jsw-1),4d0)
                  fxynw=fxy(Sh(ine,jne+1), Sh(inw-1,jnw+1), &
                            fse,Sh(isw-1,jsw), 4d0)
                  fxyne=fxy(Sh(ine+1,jne+1),Sh(inw,jnw+1), &
                            Sh(ise+1,jse), fsw, 4d0)
                  fxyse=fxy(Sh(ine+1,jne), fnw, &
                            Sh(ise+1,jse-1),Sh(isw,jsw-1),4d0)

                  rhsthermoh=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                           fxsw, fxnw, fxne, fxse,  &
                                           fysw, fynw, fyne, fyse,  &
                                           fxysw,fxynw,fxyne,fxyse, &
                                           xdn1, ydn1 )

! find rhsthermoA using cubic interpolation
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation 

                  fsw=SA(isw,jsw)
                  fnw=SA(inw,jnw)
                  fne=SA(ine,jne)
                  fse=SA(ise,jse)
                  fxsw=fx(fse, SA(isw-1,jsw), 2d0)
                  fxnw=fx(fne, SA(inw-1,jnw), 2d0)
                  fxne=fx(SA(ine+1,jne), fnw, 2d0)
                  fxse=fx(SA(ise+1,jse), fsw, 2d0)
                  fysw=fy(fnw, SA(isw,jsw-1), 2d0)
                  fynw=fy(SA(inw,jnw+1), fsw, 2d0)
                  fyne=fy(SA(ine,jne+1), fse, 2d0)
                  fyse=fy(fne, SA(ise,jse-1), 2d0)
                  fxysw=fxy(fne, SA(inw-1,jnw), &
                            SA(ise,jse-1),SA(isw-1,jsw-1),4d0)
                  fxynw=fxy(SA(ine,jne+1), SA(inw-1,jnw+1), &
                            fse,SA(isw-1,jsw), 4d0)
                  fxyne=fxy(SA(ine+1,jne+1),SA(inw,jnw+1), &
                            SA(ise+1,jse), fsw, 4d0)
                  fxyse=fxy(SA(ine+1,jne), fnw, &
                            SA(ise+1,jse-1),SA(isw,jsw-1),4d0)

                  rhsthermoA=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                           fxsw, fxnw, fxne, fxse,  &
                                           fysw, fynw, fyne, fyse,  &
                                           fxysw,fxynw,fxyne,fxyse, &
                                           xdn1, ydn1 )
!------------------------------------------------------------------------                                                   
! find output h and A at time level n                                                                                       
!------------------------------------------------------------------------ 
                  
                  htp(i,j) = htp(i,j) + 2d0*Deltat*rhsthermoh
                  Atp(i,j) = Atp(i,j) + 2d0*Deltat*rhsthermoA

               elseif (caseSL == 2) then ! upwind for special cases 
                  
                  htp(i,j) = htp(i,j) + Deltat * Sh(i,j)
                  Atp(i,j) = Atp(i,j) + Deltat * SA(i,j)

               endif

               htp(i,j) = max(htp(i,j), 0d0)
               Atp(i,j) = max(Atp(i,j), 0d0)
               Atp(i,j) = min(Atp(i,j), 1d0)

            enddo
         enddo
         
      endif

    end subroutine dh_dA_thermo
