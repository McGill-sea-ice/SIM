!***********************************************************************
!     subroutine damage: 
!       computes the level of damage inside the ice
!       should be computed after 
!
!     Revision History
!     ----------------
!
!     Ver             Date (dd-mm-yy)        Author
!
!     V01             19-04-22               A. Savard
!
!     Address : Dept. of Atmospheric and Oceanic Sciences, McGill University
!     -------   Montreal, Quebec, Canada
!     Email   :  bruno.tremblay@mcgill.ca
!
!***********************************************************************


      subroutine dam_source_terms (damtp)
      
      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_DamageVariables.h'
      include 'CB_mask.h'
      include 'CB_const.h'

      integer i, j
      double precision, intent(in) :: damtp(0:nx+1,0:ny+1)
      double precision Source_dam

      !call ViscousCoefficient(uice,vice) ! update damSS verify this

      do i = 1, nx
         do j = 1, ny

!-----------------------------------------------------------------------
!     Damage source term verify this (td and th)
!-----------------------------------------------------------------------

            Source_dam = (damSS(i,j) - damtp(i,j)) / td - damtp(i,j) / th

!-----------------------------------------------------------------------
!     Update ice damage
!-----------------------------------------------------------------------

            if ( maskC(i,j) == 1 ) then

               Sdam(i,j)=Source_dam
               
            else

               Sdam(i,j)=0d0

            endif
            
         enddo
      enddo

      return
    end subroutine dam_source_terms

    subroutine  ddam_damage (damtp)

      implicit none

      include 'parameter.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DamageVariables.h'
      include 'CB_options.h'
      include 'CB_semilag.h'

      integer i, j, caseSL, iloc, jloc
      integer isw, jsw, inw, jnw, ine, jne, ise, jse !SL SouthWest=sw, nw, ne, se corners

      double precision, intent(inout) :: damtp(0:nx+1,0:ny+1)
      double precision                :: xdn1, ydn1, rhsdam
      double precision                :: fsw, fnw, fne, fse
      double precision                :: fxsw, fxnw, fxne, fxse, fysw, fynw, fyne, fyse
      double precision                :: fxysw, fxynw, fxyne, fxyse
      double precision                :: fx, fy, fxy, cubic_interp         ! functions

      if ( adv_scheme == 'upwind' .or. adv_scheme == 'upwindRK2' ) then

         do i = 1, nx
            do j = 1, ny

               if ( maskC(i,j) == 1 ) then
                  damtp(i,j) = damtp(i,j) + Deltat * Sdam(i,j)
                  damtp(i,j) = max(damtp(i,j), 0d0)
                  damtp(i,j) = min(damtp(i,j), 1d0)
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

! find rhsdam using cubic interpolation
! PREPARATION: set 4 corners values and compute derivatives required for cubic interpolation                                

                  fsw=Sdam(isw,jsw)
                  fnw=Sdam(inw,jnw)
                  fne=Sdam(ine,jne)
                  fse=Sdam(ise,jse)
                  fxsw=fx(fse, Sdam(isw-1,jsw), 2d0)
                  fxnw=fx(fne, Sdam(inw-1,jnw), 2d0)
                  fxne=fx(Sdam(ine+1,jne), fnw, 2d0)
                  fxse=fx(Sdam(ise+1,jse), fsw, 2d0)
                  fysw=fy(fnw, Sdam(isw,jsw-1), 2d0)
                  fynw=fy(Sdam(inw,jnw+1), fsw, 2d0)
                  fyne=fy(Sdam(ine,jne+1), fse, 2d0)
                  fyse=fy(fne, Sdam(ise,jse-1), 2d0)
                  fxysw=fxy(fne, Sdam(inw-1,jnw), &
                            Sdam(ise,jse-1),Sdam(isw-1,jsw-1),4d0)
                  fxynw=fxy(Sdam(ine,jne+1), Sdam(inw-1,jnw+1), &
                            fse,Sdam(isw-1,jsw), 4d0)
                  fxyne=fxy(Sdam(ine+1,jne+1),Sdam(inw,jnw+1), &
                            Sdam(ise+1,jse), fsw, 4d0)
                  fxyse=fxy(Sdam(ine+1,jne), fnw, &
                            Sdam(ise+1,jse-1),Sdam(isw,jsw-1),4d0)

                  rhsdam=cubic_interp( fsw,  fnw,  fne,  fse,   &
                                           fxsw, fxnw, fxne, fxse,  &
                                           fysw, fynw, fyne, fyse,  &
                                           fxysw,fxynw,fxyne,fxyse, &
                                           xdn1, ydn1 )

!------------------------------------------------------------------------                                                   
! find output h and A at time level n                                                                                       
!------------------------------------------------------------------------ 
                  
                  damtp(i,j) = damtp(i,j) + 2d0*Deltat*rhsdam

               elseif (caseSL == 2) then ! upwind for special cases 
                  
                  damtp(i,j) = damtp(i,j) + Deltat * Sdam(i,j)

               endif

               damtp(i,j) = max(damtp(i,j), 0d0)
               damtp(i,j) = min(damtp(i,j), 1d0)

            enddo
         enddo
         
      endif

    end subroutine ddam_damage