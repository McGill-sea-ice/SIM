
!****************************************************************************
!     subroutine UVsolveNR_B: 
!       computes the free-drift ice velocities on the B-grid keeping the
!       non-linear water drag term, using Newton-Raphson method. 
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


      subroutine UVsolveNR_B

      implicit none

      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_Dyndim.h'
      include 'CB_const.h'
      include 'CB_mask.h'
      include 'CB_DynForcing.h'

      double precision tol, AvgSpeed, Cdw_l, deno, error
      double precision A11, A12, A21, A22, DetA, G1, G2, du, dv

      double precision                            &
                speed (0:nx+2,0:ny+2),            & ! ice speed on B-grid
                hnode (nx+1,ny+1)            ! ice thickness on B-grid


      integer i, j, k, nmax, uflag, vflag
      integer ierroru, jerroru, ierrorv, jerrorv
    
      tol    = 1d-07       ! could be increased to 1d-04
      nmax   = 20          ! maximum number of iteration
      k      = 0           ! counter for # of iteration

!------------------------------------------------------------------------
!     Initialize velocity field 
!------------------------------------------------------------------------

      do i = 0, nx+2
         do j = 0, ny+2

            uice(i,j) = 0d0
            vice(i,j) = 0d0

         enddo
      enddo


!------------------------------------------------------------------------
!     Calculate the relative linear free drift solution (u_r = u_i - uw^g) 
!     on the B-grid to use as initial guess for the non-linear free drift 
!     solution.
!     .        [A] {u_i - u_w} = R
!
!     note that the du/dt term is not included in the calculation here. 
!------------------------------------------------------------------------

      do i = 1, nx+1
         do j = 1, ny+1

            AvgSpeed  = 0.1d0           ! mean speed [m/s]

            Cdw_l     = Cdw * AvgSpeed

            deno = max( ( maskC(i,j)   + maskC(i-1,j)                  &
                        + maskC(i,j-1) + maskC(i-1,j-1) ) * 1d0,       &
                          1d-04 )
            hnode(i,j) = ( h(i,j) + h(i-1,j) + h(i,j-1) + h(i-1,j-1) ) &
                         / deno


            A11  =  Cdw_l * costheta_w
            A21  =  Cdw_l * sintheta_w + rhof * hnode(i,j) 
            A12  = -A21
            A22  =  A11

            DetA = A11 * A22 - A12 * A21 

!------------------------------------------------------------------------
!     free drift ice velocities relative to geostrophic currents (u_r)
!------------------------------------------------------------------------

            uice(i,j) = (  A22 * R1n(i,j) - A12 * R2n(i,j) ) / DetA
            vice(i,j) = ( -A21 * R1n(i,j) + A11 * R2n(i,j) ) / DetA

         enddo
      enddo




!------------------------------------------------------------------------
!     Solve for the relative ice velocities on the B-grid using 
!     Newton-Raphson method:
!
!     .                     2
!     Gi(U + dU) = Gi(U) + Sum (dGi/dUj) dUj + ... = 0
!     .                    j=1
!     
!     {dU} = [J]^-1 {-G}, where Jij = dGi/dUj  and U = u_i^fd - u_w^g.
!     Ref. Numerical recipes, pp381.
!------------------------------------------------------------------------



 5    continue

      error   = 0d0   !   \  
      ierroru = 0     !   |   Initialization of variables:
      jerroru = 0     !    >   error between successive iterations &  
      ierrorv = 0     !   |   location (i,j) where max error occurs
      jerrorv = 0     !   /



      do i = 1, nx+1
         do j = 1, ny+1


!------------------------------------------------------------------------
!     Compute the inverse Jacobian J-1 = [A] and function G on the B-grid: 
!------------------------------------------------------------------------


            speed(i,j) = max(                                         &
                              sqrt( uice(i,j) ** 2 + vice(i,j) ** 2 ) &
                            , 1d-31  & !jfl to be consist with nodim model
                            )


            A11 =  Cdw / speed(i,j) *                                  &
                     (                                                 &
                       ( 2d0 * vice(i,j) ** 2 + uice(i,j) ** 2 ) *     &
                           costheta_w                                  &
                       + vice(i,j) * uice(i,j) * sintheta_w         )

            A22 =  Cdw / speed(i,j) *                                  &
                     (                                                 &
                       ( 2d0 * uice(i,j) ** 2 + vice(i,j) ** 2 ) *     &
                           costheta_w                                  &
                       - uice(i,j) * vice(i,j) * sintheta_w         )

            A21 = -Cdw / speed(i,j) *                                  &
                     (                                                 &
                       ( 2d0 * uice(i,j) ** 2 + vice(i,j) ** 2 ) *     &
                           sintheta_w                                  &
                       + uice(i,j) * vice(i,j) * costheta_w         )  &
                 - rhof * hnode(i,j) 

            A12 =  Cdw / speed(i,j) *                                  &
                     (                                                 &
                       ( 2d0 * vice(i,j) ** 2 + uice(i,j) ** 2 ) *     &
                           sintheta_w -                                &
                         uice(i,j) * vice(i,j) * costheta_w         )  &
                 + rhof * hnode(i,j)



            G1 = R1n(i,j)                                              &
                      - Cdw * speed(i,j) * ( uice(i,j) * costheta_w    &
                                           - vice(i,j) * sintheta_w )  &
                      + rhof *  hnode(i,j) * vice(i,j)


            G2 = R2n(i,j)                                              &
                      - Cdw * speed(i,j) * ( vice(i,j) * costheta_w    &
                                           + uice(i,j) * sintheta_w )  &
                      - rhof * hnode(i,j) * uice(i,j)


!------------------------------------------------------------------------
!     Calculate the relative ice velocities: u_i = u_i + du
!     Note: Freedrift equations are algebraic, they do not support BC's.
!     .     If BC are artificially applied during the solution, the 
!     .     quadratic convergence is lost. 
!     .     Boundary condition on the velocities are imposed at the end.
!------------------------------------------------------------------------

            detA = A11 * A22 - A12 * A21



            du   = ( A11 * G1 + A12 * G2 ) / detA  
            dv   = ( A21 * G1 + A22 * G2 ) / detA  


            uice(i,j) = uice(i,j) + du
            vice(i,j) = vice(i,j) + dv


!------------------------------------------------------------------------
!     Maximum error between succesive iterations
!------------------------------------------------------------------------

            if (abs(du) .gt. error) then
               error = abs(du)
               ierroru = i
               jerroru = j
               uflag   = 1
               vflag   = 0
            endif
            if (abs(dv) .gt. error) then
               error = abs(dv)
               ierrorv = i
               jerrorv = j
               uflag   = 0
               vflag   = 1
            endif


         enddo
      enddo



      k = k + 1   


!------------------------------------------------------------------------
!     If the maximum # of iteration is reached, print the error and exit
!------------------------------------------------------------------------

      if (k .eq. nmax) then

         print *, 'WARNING: Max # iteration reached in UVsolveNR_B:'
         if ( uflag .eq. 1)                                            &
            print *, 'max error on u = ', error,error**2,              &
                     ' at (i,j) = ', ierroru, jerroru
         if ( vflag .eq. 1)                                            &
            print *, 'max error on v = ', error,error**2,              &
                     ' at (i,j) =', ierrorv, jerrorv

         goto 10 

      endif 


      if ( error .gt. tol ) goto 5  
      
 10   continue

!      print *, '# of iteration in UVsolveNR_B', k


!------------------------------------------------------------------------
!     Interpolate the ice velocities (u_i) on the C-grid, and apply
!     the boundary condition u_n = 0.
!------------------------------------------------------------------------


      do i = 1, nx+1
         do j = 1, ny

            uice(i,j) = ( uice(i,j) + uice(i,j+1) ) / 2d0

            uice(i,j) = ( uice(i,j) + uwatnd(i,j) )                    &
                          * min ( maskB(i,j) + maskB(i,j+1), 1)

         enddo
         uice(i,ny+1) = 0d0              ! outside domain, u set to 0
      enddo

      do j = 1, ny+1
         do i = 1, nx

            vice(i,j) = ( vice(i,j) + vice(i+1,j) ) / 2d0

            vice(i,j) = ( vice(i,j) + vwatnd(i,j) )                    &
                          * min ( maskB(i,j) + maskB(i+1,j), 1)

         enddo
         vice(nx+1,j) = 0d0              ! outside domain, v set to 0
      enddo


!------------------------------------------------------------------------
!     Apply open boundary conditions ( du/dn = 0 )
!------------------------------------------------------------------------

      do j = 1, ny

         if ( maskB(1,j) + maskB(1,j+1) .ge. 1 )                       &
                   uice(1,j)    = ( 4d0 * uice(2,j)   - uice(3,j)    ) &
                                  / 3d0
         if ( maskB(nx+1,j) + maskB(nx+1,j+1) .ge. 1 )                 &
                   uice(nx+1,j) = ( 4d0 * uice(nx,j)  - uice(nx-1,j) ) &
                                  / 3d0

      enddo

      do i = 1, nx

         if ( maskB(i,1) + maskB(i+1,1) .ge. 1 )                      &
                   vice(i,1)    = ( 4d0 * vice(i,2)  - vice(i,3)    ) &
                                  / 3d0
         if ( maskB(i,ny+1) + maskB(i+1,ny+1) .ge. 1 )                &
                   vice(i,ny+1) = ( 4d0 * vice(i,ny) - vice(i,ny-1) ) &
                                  / 3d0
      enddo

      return
      end
      

