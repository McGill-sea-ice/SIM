!************************************************************************                                  
!     Subroutine ViscousCoefficient: computes the non-linear viscous                                         
!     coefficient at the centers (zetaC, etaC) and at the nodes (etaB) of the                                
!     grid. All the strain rate calculations are second order accurate.                                      
!     The ice strength (Pp) is collocated with zetaC, etaC.
!                                                                                                           
!     We use the method described in Bouillon et al. 2013 (the EVP revisited):    
!                                                                                                            
!     1) ep12 are calculated at the grid nodes (ep12b).                                                      
!     2) ep12c**2 = ( ep12b**2 + ep12b**2 + ep12b**2 + ep12b**2 / 4).                                        
!     3) zetaC = Pp / 2f(ep11,ep22,ep12c**2)                                                                 
!     4) etaB = ( etaC + etaC + etaC + etaC ) / 4                                                            
!                                                                                                            
!     note : etaC = e^-2 * zetaC                                                                             
!                                                                                                            
!     We also apply the replacement closure of Kreyscher et al 2000 (eq.7).                                   
!                                                                                                            
!     JF Lemieux, 15 Febrary 2013                                                                            
!                                                                                                            
!************************************************************************     


      subroutine ViscousCoefficient(utp,vtp)
        use ellipse
        implicit none

      include 'parameter.h'
      include 'CB_Dyndim.h'
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_options.h'

      integer i, j, rheo, summaskC

      double precision dudx, dvdy, dudy, dvdx, deno, denomin, denoT
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)
      double precision ep12(0:nx+2,0:ny+2), meanep12sq

      denomin = 2d-09 ! Hibler, 1979

      rheo = Rheology ! define local variable to speed up the code

!------------------------------------------------------------------------
!     free slip boundary condition:
!       d(v_tangential)/d(normal) = 0 & v_normal = 0, at close boundary
!------------------------------------------------------------------------


      if ( BndyCond .eq. 'freeslip' ) then

               print *, 'WARNING: freeslip option not fully tested yet'

!------------------------------------------------------------------------
!     strain-rates calculation
!     no-slip boundary condition :  v_normal = 0 and v_tangential = 0 
!------------------------------------------------------------------------

      elseif ( BndyCond .eq. 'noslip' ) then


         do i = 0, nx+1
            do j = 0, ny+1

               etaC(i,j)  = 0d0
               zetaC(i,j) = 0d0 
               etaB(i,j)  = 0d0
	       ep12(i,j)  = 0d0
 
            enddo
         enddo


!------------------------------------------------------------------------
!     ep12 calculation at the grid node (see p.2-118 PDF notebook)
!------------------------------------------------------------------------

         
         ! for sig12B (defined at the node)
         do j = 1, ny+1 
            do i = 1, nx+1

               summaskC = maskC(i-1,j) + maskC(i,j) + & 
                          maskC(i,j-1) + maskC(i-1,j-1)

               if (summaskC .ge. 2) then

                  if (summaskC .eq. 4) then
! oo
! oo normal
                     dudy = ( utp(i,j) - utp(i,j-1) ) / Deltax !case 1 
                     dvdx = ( vtp(i,j) - vtp(i-1,j) ) / Deltax
                     
                  elseif (summaskC .eq. 3) then 

                     if (maskC(i-1,j) .eq. 0) then !case 2
! xo
! oo    
                        dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                        dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax

                     elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           
                        dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                        dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                     elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
                        dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                        dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax

                     elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      
                        dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                        dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                     else

                        print *, 'wowowo1'
                        stop

                     endif !summaskC .eq. 3

                  elseif (summaskC .eq. 2) then !case 7

                     if (maskC(i-1,j) .eq. 0 .and. &
                         maskC(i-1,j-1) .eq. 0) then
! xo
! xo
                        dudy = 0d0
                        dvdx = (3d0*vtp(i,j) - vtp(i+1,j)/3d0) / Deltax
                        
                     elseif(maskC(i,j) .eq. 0 .and. & !case 6
                            maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           
                        dudy = 0d0
                        dvdx = (-3d0*vtp(i-1,j) + vtp(i-2,j)/3d0) / Deltax
                        
                     elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                            maskC(i,j) .eq. 0) then
! xx
! oo
                        dudy = (-3d0*utp(i,j-1) + utp(i,j-2)/3d0) / Deltax
                        dvdx = 0d0

                     elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                            maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                        dudy = (3d0*utp(i,j) - utp(i,j+1)/3d0) / Deltax
                        dvdx = 0d0

                     elseif(maskC(i-1,j) .eq. 0 .and. & !case 15
                            maskC(i,j-1) .eq. 0) then
! xo                                                                    
! ox
			
                     elseif(maskC(i,j) .eq. 0 .and. & !case 16
                            maskC(i-1,j-1) .eq. 0) then
! ox                      	                    
! xo                    ! don't do anything (could be improved)                                                     
                     else

                     print *, 'wowowo2'
                     stop

                     endif  !summaskC .eq. 2

                  else
                     print *, 'wowowo3'
                     stop

                  endif !if (summaskC .eq. 4) then...

		  ep12(i,j) = dudy + dvdx
                  
	      endif !summaskC .ge. 2 

            enddo
         enddo


!------------------------------------------------------------------------
!     Shear and bulk viscosity calculation at the grid center
!------------------------------------------------------------------------   

         do i = 1, nx
            do j = 1, ny

               
               if ( maskC(i,j) .eq. 1 ) then
                  
                  dudx = ( utp(i+1,j) - utp(i,j) ) / Deltax
                  dvdy = ( vtp(i,j+1) - vtp(i,j) ) / Deltax
                  
                  if ( rheo .eq. 1 ) then ! ellipse, jfl p.892

		     meanep12sq = ( (ep12(i,j))**2 + (ep12(i,j+1))**2+ & 
				    (ep12(i+1,j+1))**2 + (ep12(i+1,j))**2 )/4d0

                     deno = sqrt(( dudx **2 + dvdy **2 )*(1.0d0 + ell_2) &
                                 + 2.0d0 * dudx * dvdy * (1.0d0 - ell_2) &
                                 + ell_2 * meanep12sq )

                     if ( regularization .eq. 'tanh' ) then

                        deno = max( deno, 1d-20 )
                     
                        zetaC(i,j) =  ( Pp(i,j)/denomin ) &
                                     *( tanh(denomin*(1/deno)))

                     elseif ( regularization .eq. 'Kreyscher' ) then
                        
                        denoT = deno + denomin

                        zetaC(i,j) = Pp(i,j)/denoT 

                     else
                        
                        print *, 'WRONG REGULARIZATION'
                        stop
                        
                     endif

                     P(i,j) = zetaC(i,j)*deno ! replacement pressure 

                     etaC(i,j)  = zetaC(i,j) * ell_2

                        
                  elseif ( rheo .eq. 2 ) then ! triangle, jfl p.1124

                     stop
                   
                  endif

                  
               endif
               
               
            enddo
         enddo

!------------------------------------------------------------------------
!     Set zetaC and etaC to 0.0 at the open boundaries (see p.32-33 EC-2)
!------------------------------------------------------------------------

         do i = 1, nx+1

            if (maskC(i,0) .eq. 1) then
               etaC(i,1)  = 0d0
               zetaC(i,1) = 0d0
            endif
            
            if (maskC(i,ny+1) .eq. 1) then
               etaC(i,ny)  = 0d0
               zetaC(i,ny) = 0d0
            endif
            
         enddo
         
         do j = 1, ny+1
            
            if (maskC(0,j) .eq. 1) then   
               etaC(1,j)  = 0d0
               zetaC(1,j) = 0d0
            endif
            
            if (maskC(nx+1,j) .eq. 1) then  
               etaC(nx,j)  = 0d0
               zetaC(nx,j) = 0d0
            endif           

         enddo

!------------------------------------------------------------------------
!     Shear viscosity calculation at the grid node (see p.2-118 PDF notebook)
!------------------------------------------------------------------------

         
         ! for sig12B (defined at the node)
         do j = 1, ny+1 
            do i = 1, nx+1

               summaskC = maskC(i-1,j) + maskC(i,j) + & 
                          maskC(i,j-1) + maskC(i-1,j-1)

               if (summaskC .ge. 2) then

                  if (summaskC .eq. 4) then
! oo
! oo normal
                     
                     etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 4d0
                     
                     
                  elseif (summaskC .eq. 3) then 

                     if (maskC(i-1,j) .eq. 0) then !case 2
! xo
! oo    
                        
		      etaB(i,j) = ( etaC(i,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 3d0

                     elseif (maskC(i,j) .eq. 0) then !case 3
! ox
! oo           

		      etaB(i,j) = ( etaC(i-1,j) + etaC(i,j-1) + etaC(i-1,j-1) ) / 3d0

                     elseif (maskC(i,j-1) .eq. 0) then !case 5
! oo                                                          
! ox
    
		      etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i-1,j-1) ) / 3d0

                     elseif (maskC(i-1,j-1) .eq. 0) then !case 4
! oo                                                            
! xo      

		      etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) + etaC(i,j-1) ) / 3d0
		      
                     else

                        print *, 'wowowo1'
                        stop

                     endif !summaskC .eq. 3

                  elseif (summaskC .eq. 2) then !case 7

                     if (maskC(i-1,j) .eq. 0 .and. &
                         maskC(i-1,j-1) .eq. 0) then
! xo
! xo

		      etaB(i,j) = ( etaC(i,j) + etaC(i,j-1) ) / 2d0

                     elseif(maskC(i,j) .eq. 0 .and. & !case 6
                            maskC(i,j-1) .eq. 0) then
! ox
! ox                                                                           

                        etaB(i,j) = ( etaC(i-1,j) + etaC(i-1,j-1) ) / 2d0

                     elseif(maskC(i-1,j) .eq. 0 .and. & !case 8           
                            maskC(i,j) .eq. 0) then
! xx
! oo

		      etaB(i,j) = ( etaC(i-1,j-1) + etaC(i,j-1) ) / 2d0

                     elseif(maskC(i,j-1) .eq. 0 .and. & !case 9
                            maskC(i-1,j-1) .eq. 0) then
! oo                                                                     
! xx                                          
                        
		      etaB(i,j) = ( etaC(i-1,j) + etaC(i,j) ) / 2d0  
  
                     elseif(maskC(i-1,j) .eq. 0 .and. & !case 15
                            maskC(i,j-1) .eq. 0) then
! xo                                                                    
! ox                    ! don't do anything (could be improved)                                        
                        
			
                     elseif(maskC(i,j) .eq. 0 .and. & !case 16
                            maskC(i-1,j-1) .eq. 0) then
! ox                                                                      
! xo                    ! don't do anything (could be improved)                                                     
			
                     else

                     print *, 'wowowo2'
                     stop

                     endif  !summaskC .eq. 2

                  else
                     print *, 'wowowo3'
                     stop

                  endif !if (summaskC .eq. 4) then...
                  
	      endif !summaskC .ge. 2 

            enddo
         enddo
         
!------------------------------------------------------------------------
!     Set etaB to 0.0 at the open boundaries (see p.32-33 EC-2)
!------------------------------------------------------------------------

	 etaB(1,1)=0d0
	 etaB(1,ny+1)=0d0
         do j = 2, ny 

	    summaskC = maskC(0,j) + maskC(1,j) + & 
                       maskC(1,j-1) + maskC(0,j-1)

            if (summaskC .ge. 2) then

               if (summaskC .ge. 3) then
		  etaB(1,j) = 0d0 
	       elseif (summaskC .eq. 2) then
		  if (maskC(0,j) .eq. 0 .and. maskC(1,j) .eq. 0) then
		    etaB(1,j) = 0d0
		  elseif (maskC(0,j-1) .eq. 0 .and. maskC(1,j-1) .eq. 0) then
		    etaB(1,j) = 0d0
		  elseif (maskC(0,j) .eq. 0 .and. maskC(0,j-1) .eq. 0) then
		    ! don't do anything...it can be non zero
		  else
		    print *, 'wowowo' ! ocean on the left is not possible
		  endif	
		
	       endif	
						    					    
	    endif
	    
	 enddo
	 
	 etaB(nx+1,1)=0d0
	 etaB(nx+1,ny+1)=0d0
         do j = 2, ny 

	    summaskC = maskC(nx,j) + maskC(nx+1,j) + & 
                       maskC(nx+1,j-1) + maskC(nx,j-1)

            if (summaskC .ge. 2) then

               if (summaskC .ge. 3) then
		  etaB(nx+1,j) = 0d0 
	       elseif (summaskC .eq. 2) then
		  if (maskC(nx,j) .eq. 0 .and. maskC(nx+1,j) .eq. 0) then
		    etaB(nx+1,j) = 0d0
		  elseif (maskC(nx,j-1) .eq. 0 .and. maskC(nx+1,j-1) .eq. 0) then
		    etaB(nx+1,j) = 0d0
		  elseif (maskC(nx+1,j) .eq. 0 .and. maskC(nx+1,j-1) .eq. 0) then
		    ! don't do anything...it can be non zero
		  else
		    print *, 'wowowo' ! ocean on the right is not possible
		  endif	
		
	       endif	
						    					    
	    endif
	    
	 enddo
	 
         do i = 2, nx 

	    summaskC = maskC(i-1,1) + maskC(i,1) + & 
                       maskC(i-1,0) + maskC(i,0)

            if (summaskC .ge. 2) then

               if (summaskC .ge. 3) then
		  etaB(i,1) = 0d0 
	       elseif (summaskC .eq. 2) then
		  if (maskC(i-1,1) .eq. 0 .and. maskC(i-1,0) .eq. 0) then
		    etaB(i,1) = 0d0
		  elseif (maskC(i,1) .eq. 0 .and. maskC(i,0) .eq. 0) then
		    etaB(i,1) = 0d0
		  elseif (maskC(i-1,0) .eq. 0 .and. maskC(i,0) .eq. 0) then 
		    ! don't do anything...it can be non zero
		  else
		    print *, 'wowowo' ! ocean below is not possible
		  endif	
		
	       endif	
						    					    
	    endif
	    
	 enddo
	 
	 do i = 2, nx 

	    summaskC = maskC(i-1,ny+1) + maskC(i,ny+1) + & 
                       maskC(i-1,ny) + maskC(i,ny)

            if (summaskC .ge. 2) then

               if (summaskC .ge. 3) then
		  etaB(i,ny+1) = 0d0 
	       elseif (summaskC .eq. 2) then
		  if (maskC(i-1,ny+1) .eq. 0 .and. maskC(i-1,ny) .eq. 0) then
		    etaB(i,ny+1) = 0d0
		  elseif (maskC(i,ny+1) .eq. 0 .and. maskC(i,ny) .eq. 0) then
		    etaB(i,ny+1) = 0d0
		  elseif (maskC(i-1,ny+1) .eq. 0 .and. maskC(i,ny+1) .eq. 0) then
		    ! don't do anything...it can be non zero
		  else
		    print *, 'wowowo' ! ocean above is not possible
		  endif	
		
	       endif	
						    					    
	    endif
	    
	 enddo

      endif
      
      return
    end subroutine ViscousCoefficient



