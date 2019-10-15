
      subroutine PrepFGMRES_NK(k,res,resk_1,tot_its, x, Fu)
      
      use numerical_VP
      
      implicit none

      include 'parameter.h' 
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_Dyndim.h'

      integer  icode, iter, maxits, iout, tot_its
      integer, intent(in) :: k

      double precision, intent(in) :: res
      double precision, intent(inout) :: Fu(nvar)
      double precision  eps, resk_1, gamma, epsilon
      double precision  sol(nvar)
      double precision  vv(nvar,img1), wk(nvar,img), Funeg(nvar)
      double precision  wk1(nvar), wk2(nvar),x(nvar)

!------------------------------------------------------------------------
!     This routine solves J(u)du = -F(u) where u = u^k, du = du^k using the
!     Jacobian free Newton Krylov method. The Krylov method is the precon-
!     ditioned FGMRES method. A linesearch globalization method can be used
!     Some usefull references are:
!
!     Lemieux et al., J.of.Comput.Physics, 2012.
!     Lemieux et al., J.of.Comput.Physics, 2010.
!     Knoll and Keyes, J.of.Comput.Physics, 2004.
!     Eisenstat and Walker, SIAM J.Optimization, 1994.
!     Eisenstat and Walker, SIAM J.Sci.Comput., 1996.
!
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!     Making of the RHS vector : -F(u)
!------------------------------------------------------------------------

      CdwC1f = CdwC1  ! the matrix M(u^k) for the precond should stay the same 
      CdwC2f = CdwC2  ! when solving J(u^k)du = -F(u^k). Because Jacfreevec 
      etaCf  = etaC   ! needs F(u^k +eps*v), this changes the visc coeff and
      zetaCf = zetaC  ! the CdwC. This is why we fixe them here for the precond.
      etaBf  = etaB   ! precond therefore uses CdwC1f...etaBf where f = fixed
      Pf     = P

      Funeg = Fu
      Fu = -1d0 * Fu ! mult by -1 because we solve Jdu = -F(u)
      
!------------------------------------------------------------------------
!     Initial guess vector: because we solve for du and du is just a 
!     correction, we set the initial guess to zero
!------------------------------------------------------------------------

      sol = 0d0

!------------------------------------------------------------------------
!     Choose the forcing term (gamma)
!------------------------------------------------------------------------

      call forcing_term(k,res,resk_1,gamma)
      resk_1 = res
      
!------------------------------------------------------------------------
!      Begining of FGMRES method    
!------------------------------------------------------------------------

      print *, 'L2norm is', k, res, gamma

      eps = gamma * res ! setting the tolerance for fgmres

      maxits = 50  ! should be divisible by img (15)   
      iout   = 0    ! set higher than 0 to have res(ite)

      icode = 0
10    CONTINUE

      call fgmres (nvar,img,Fu,sol,iter,vv,wk,wk1,wk2, &
           eps,maxits,iout,icode,tot_its)

      IF ( icode == 1 ) THEN
         CALL precondLSOR (wk1, wk2)
!         CALL precondSOR (wk1, wk2)
!         CALL precondJacobi (wk1, wk2)  
!          CALL precond_EVPC (wk1, wk2)
         GOTO 10
      ELSEIF ( icode >= 2 ) THEN
         epsilon = 1d-07
         call JacfreeVec (wk1, wk2, Funeg, epsilon) ! approximates Jv
         GOTO 10
      ENDIF

!------------------------------------------------------------------------
!      End of FGMRES method    
!------------------------------------------------------------------------

      if (tot_its .eq. maxits) then
         print *,'WARNING: FGMRES has not converged'
         print*, 'Please check the precond relaxation param (wlsor or wsor).'
!         stop
      endif

! icode = 0 means that fgmres has finished and sol contains the app. solution

!------------------------------------------------------------------------
!      Find kth iterate
!------------------------------------------------------------------------

      if ( k .lt. klinesearch ) then
         x = x + sol ! u^k = u^k-1 + du  
      else
         call linesearch(sol, x, res, Fu) !u^k=u^k-1 +a*du (a=0.125,0.25,0.5 or 1)
      endif

      call transformer(uice,vice,x,0)

      return
    end subroutine PrepFGMRES_NK
      

    subroutine forcing_term(k,res,resk_1,gamma)
      use numerical_VP
      implicit none

      integer, intent(in) :: k
      
      double precision, intent(in) :: res, resk_1
      double precision :: gamma_ini, phi_e, alp_e
      double precision, intent(out) :: gamma
      
      gamma_ini = 0.99d0
      phi_e     = 1d0
      alp_e     = 1.5d0 !      alp_e = (1d0 + 5d0**0.5d0)/2d0 !2d0                                                                                  

      if (k .eq. 1) then

         gamma = gamma_ini

      elseif (k .gt. 1 .and. res .gt. res_t) then

         gamma = gamma_ini

      else
         gamma = phi_e * (res/resk_1)**alp_e ! Eisenstat, 1996, equ.(2.6)      
         gamma = min(gamma_ini,gamma)
         gamma = max(0.1d0,gamma)
      endif

    end subroutine forcing_term

    subroutine locate_max_residual(k,Fu_tp)

! Finds the location of the maximum residual (at the u-point)     

      implicit none
      include 'parameter.h'

      integer, intent(in) :: k
      integer :: i,j,imaxu,jmaxu,imaxv,jmaxv,n

      double precision, intent(in) :: Fu_tp(nvar)
      double precision :: Fmaxu, Fmaxv

      Fmaxu=0d0
      Fmaxv=0d0
      imaxu=0
      jmaxu=0
      imaxv=0
      jmaxv=0

      do j = 1, ny
         do i = 1, nx+1

            n = i+(j-1)*(nx+1)
            if ( abs(Fu_tp(n)) .gt. Fmaxu) then
               imaxu=i
               jmaxu=j
               Fmaxu=abs(Fu_tp(n))
            endif

         enddo
      enddo

      do j = 1, ny+1
         do i = 1, nx

            n = i+(j-1)*nx+(nx+1)*ny
            if ( abs(Fu_tp(n)) .gt. Fmaxv) then
               imaxv=i
               jmaxv=j
               Fmaxv=abs(Fu_tp(n))
            endif
         enddo
      enddo
         
      print *, 'k, Fmaxu, imaxu, jmaxu: ', k,Fmaxu,imaxu,jmaxu
      print *, 'k, Fmaxv, imaxv, jmaxv: ', k,Fmaxv,imaxv,jmaxv
      
    end subroutine locate_max_residual

    subroutine linesearch(sol,x, res, Fu)

! Apply the linesearch method     
      use datetime, only: datetime_type
      implicit none
      include 'parameter.h'
      include 'CB_options.h'
      include 'CB_DynVariables.h'

      integer :: l

      type(datetime_type) :: dummy

      double precision, intent(in) :: sol(nvar)
      double precision, intent(inout) :: x(nvar), Fu(nvar)
      double precision, intent(in) :: res
      double precision :: aLS, resnew, xini(nvar), rhs(nvar)
      double precision utp(0:nx+2,0:ny+2), vtp(0:nx+2,0:ny+2)

      xini = x

      do l = 1, 4

<<<<<<< HEAD
         aLS = 1d0/(2d0**(1d0*(l-1)))
         x = xini + aLS*sol
         call transformer (utp,vtp,x,0)
         if ( IMEX .gt. 0 ) then ! IMEX method 1 or 2   
            call advection ( un1, vn1, utp, vtp, hn2, An2, hn1, An1, h, A )
            if (Rheology .eq. 3) then ! calculate the damage factor
               dam  = dam1
               damB = damB1
               call stress_strain_MEB(uice,vice,dummy,0,0)
            else
               call Ice_strength()
            endif
            call bvect_ind
         endif
         call ViscousCoefficient(utp,vtp)
         call bvect (utp,vtp,rhs)
         call Funk (x,rhs,Fu)
         resnew = sqrt(DOT_PRODUCT(Fu,Fu))

         if (resnew .lt. res) exit
         print *, 'LINESEARCH', l, res, resnew

      enddo

    end subroutine linesearch
