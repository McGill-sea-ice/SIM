
      subroutine PrepFGMRES(k,res,tot_its,sol,rhs)

      implicit none

      include 'parameter.h' 
      include 'CB_DynVariables.h'
      include 'CB_DynForcing.h'
      include 'CB_Dyndim.h'

      INTEGER  icode, iter, maxits, iout, k, tot_its

      double precision, intent(in) :: res,rhs(nvar)
      double precision  eps
      double precision  sol(nvar), Flin(nvar)
      double precision  vv(nvar,img1), wk(nvar,img)
      double precision  wk1(nvar), wk2(nvar)

!------------------------------------------------------------------------
!     Save Cdw and viscous coeffs for the preconditioner
!------------------------------------------------------------------------

      CdwC1f = CdwC1 ! this is done in order to use the same preconditioners
      CdwC2f = CdwC2 ! as used in the JFNK solver. The precond used the   
      etaCf  = etaC  ! f (fixed) version of these variables. For the standard
      zetaCf = zetaC ! solver they dont change during a OL ite but they change 
      etaBf  = etaB  ! during a Newton ite (because of the calc of F(u+eps*v)
      Pf     = P

!------------------------------------------------------------------------
!      Begining of FGMRES method    
!------------------------------------------------------------------------
      
      print *, 'L2norm is', k, res

      if (k .le. 5) then
         eps = res / 5d0
      else
         eps = res / 25d0
      endif

      maxits = 50  ! should be divisible by img (15)   
      iout   = 0    ! set higher than 0 to have res(ite)

      icode = 0
 10   CONTINUE

      call fgmres (nvar,img,rhs,sol,iter,vv,wk,wk1,wk2, &
                   eps,maxits,iout,icode,tot_its)

      IF ( icode == 1 ) THEN
         CALL precondLSOR ( wk1, wk2)
!         CALL precondSOR ( wk1, wk2)  
!         CALL precondJacobi ( wk1, wk2)          
!         CALL precond_EVPC (wk1, wk2)     
!        wk2=wk1
         GOTO 10
      ELSEIF ( icode >= 2 ) THEN
         CALL MATVEC ( wk1, wk2) 
         GOTO 10
      ENDIF

!------------------------------------------------------------------------
!      End of FGMRES method    
!------------------------------------------------------------------------

      if (tot_its .eq. maxits) then
         print *,'WARNING: FGMRES has not converged'
         print*, 'Please check the preconditioner relaxation parameter (wlsor or wsor).'
!         stop
      endif

! icode = 0 means that fgmres has finished and sol contains the app. solution

!------------------------------------------------------------------------
!     From vector to matrix: puts sol in uice and vice
!------------------------------------------------------------------------

      call transformer(uice,vice,sol,0)

      return
    end subroutine PrepFGMRES
      




