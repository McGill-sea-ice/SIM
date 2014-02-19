MODULE UTILS

  IMPLICIT NONE

  INTEGER, PARAMETER :: RP = SELECTED_REAL_KIND(12)
  REAL(KIND=RP), PARAMETER :: PI = 3.14159265
  
  


  INTERFACE average
     MODULE PROCEDURE average_2d
  END INTERFACE

CONTAINS
  

  ELEMENTAL FUNCTION radian(angle)
    ! Return the angle [degree] in radians.
    REAL(KIND=RP), INTENT(in) :: angle
    REAL(KIND=RP) :: radian
    radian = MOD(angle/180.*pi, 2*pi)
  END FUNCTION radian


  PURE FUNCTION linspace(a, b, n)
    ! Return a vector of evenly spaced numbers.
    !
    ! Input
    ! -----
    ! a : real
    !   Starting point.
    ! b : real
    !   Ending point (inclusive).
    ! n : integer
    !   Number of points.

    REAL(KIND=RP), INTENT(in) :: a, b
    INTEGER, INTENT(in) :: n
    REAL(KIND=RP) :: linspace(n)
    REAL(KIND=RP) :: dx
    INTEGER :: t(n), i

    dx = (b-a)/(n-1.)
    t = (/ (i, i=0,n) /)
    linspace = a + dx*t

  END FUNCTION linspace
    

  PURE  FUNCTION arange(a, dx, n)
    ! Return a vector of evenly spaced numbers. 
    ! 
    ! Input 
    ! -----
    ! a : real
    !   Starting point,
    ! dx : real
    !   Distance between each point.
    ! n : integer
    !   Number of points.
    ! 
    REAL(KIND=RP), INTENT(in) :: a, dx
    INTEGER, INTENT(in) :: n
    REAL(KIND=RP) :: arange(n)
    INTEGER :: t(n), i

    t = (/ (i, i=0,n) /)
    arange = a + dx * t
  END FUNCTION arange
    


  PURE SUBROUTINE meshgrid(x, y, mx, my)
    ! Return 2D arrays constructed by stacking
    ! x and y.
    
    REAL(KIND=RP), DIMENSION(:), INTENT(in):: x,y
    REAL(KIND=RP), DIMENSION(SIZE(x), SIZE(y)), INTENT(out) :: mx, my

    mx = SPREAD(x, 2, SIZE(y))
    my = TRANSPOSE(SPREAD(y, 2, SIZE(x)))

  END SUBROUTINE meshgrid


  FUNCTION vec_radian(angle)
    ! Only for benchmarking purposes.
    REAL(KIND=RP), INTENT(in) :: angle(:)
    REAL(KIND=RP) :: vec_radian(SIZE(angle))
    
    vec_radian = MOD(angle/180.*pi, 2*pi)
  END FUNCTION vec_radian


  ELEMENTAL  FUNCTION divided_difference(x0,x1,f0,f1) RESULT (dd)
    ! Return the divided difference:
    !
    !      f1-f0
    ! dd = ----- .
    !      x1-x0
    !
    REAL(KIND=RP), INTENT(in) :: x0,x1,f0,f1
    REAL(KIND=RP) :: dd
    dd = (f1-f0)/(x1-x0)
  END FUNCTION divided_difference
  

   ELEMENTAL FUNCTION finite_difference_2nd(x0, x1, x2, f0, f1, f2, x) RESULT (d)
    ! Compute the derivative of a function f at position x, 
    ! assuming fi = f(xi).
    !
    ! This is a second order scheme, so three known values are
    ! needed.
    !
    ! Input
    ! -----
    ! x0,x1,x2: REAL
    !   Positions.
    ! f0, f2, f2: REAL
    !   Value of function evaluated at x0,x1,x2.
    !
    ! Output
    ! ------
    ! d : REAL
    !   Function derivative f'(x) evaluated at x. 
    
    REAL(KIND=RP), INTENT(in):: x0,x1,x2,f0,f1,f2,x
    REAL(KIND=RP):: d, c0, c1, c2
    
    c0 = divided_difference(x0, x1, f0, f1)
    c1 = divided_difference(x1, x2, f1, f2)
    c2 = divided_difference(x0, x2, c0, c1)

    d = c0 + (2.*x - x1 - x0) * c2

  END FUNCTION finite_difference_2nd


  ELEMENTAL FUNCTION mid_central_difference(f1, f2, delta) RESULT (der)
    ! Return the derivative at the mid point: (f2-f1)/delta.
    !
    ! That is, if f1 and f2 are defined at x1 and x2, return the derivative
    ! of f in 
    !
    !  -------f1----X-----f2---------
    !          |  delta   |

    REAL(KIND=RP), INTENT(in) :: f1, f2, delta
    REAL(KIND=RP) :: der

    der = (f2-f1)/delta
  END FUNCTION mid_central_difference



  ELEMENTAL FUNCTION central_difference(f0, f1, f2, delta1, delta2) RESULT (der)
    ! Return the derivative at the central point (1).
    !
    ! ------f0-------f1-------f2------
    !       | delta1 | delta2 |

    REAL(KIND=RP), INTENT(in) :: f0, f1, f2, delta1, delta2
    REAL(KIND=RP) :: der


    der = ( (f1-f0)*delta2/delta1 + (f2-f1)*delta1/delta2 )/(delta1+delta2)

  END FUNCTION central_difference


  PURE FUNCTION mid_derivative(f, delta, mask) RESULT (der)
    ! Compute the first derivative of function f assuming points are
    ! equally spaced. The derivative is evaluated using at the mid point
    ! using (f(i+1) - f(i)) / delta.
    !
    ! The derivative is set to 0.0 where the mask is False. 
    !
    ! Input
    ! -----
    ! f : real array (n+1)
    !   Function evaluated at n+1 equally spaced locations.
    ! delta : real
    !   Distance between each location.
    ! mask : logical array (n)
    !   Whether or not to evaluate the derivative at that point. 
    !
    ! Output 
    ! ------
    ! der : real array (n)
    !   First derivative. Set to 0. if mask is False.

    REAL(KIND=RP), DIMENSION(:), INTENT(in) :: f
    LOGICAL, DIMENSION(SIZE(f)-1), INTENT(in) :: mask
    REAL(KIND=RP), INTENT(in) :: delta
    REAL(KIND=RP), DIMENSION(SIZE(mask)) :: der
    INTEGER :: n

    n = SIZE(mask)

    WHERE (mask)
       der = mid_central_difference(f(1:n), f(2:n+1), delta)
    ELSEWHERE
       der = 0.0
    END WHERE

  END FUNCTION mid_derivative


  PURE FUNCTION no_slip_derivative(f, delta, mask) RESULT (der)
    ! Compute the first derivative of function f assuming points are
    ! equally spaced. The derivative is evaluated at the same position 
    ! as f using central differences for interior points and forward 
    ! and backward differences for the left- and rightmost points 
    ! respectively.  
    !
    ! Masked values are handled by assuming a no splip condition.   
    ! That is, if for some location `i` the mask is False, then the value of f(i)
    ! is ignored and the derivative at `i+1` is computed assuming f(i+1/2) = 0.
    ! 
    ! 
    !   0    i    1   i+1   2    i+2   3
    ! --|----+----|----+----|-----+----|-----
    !        M    0. f(i+1) |  f(i+2)  |  
    !
    !
    ! Input
    ! -----
    ! f : real array(n)
    !   Value of a FUNCTION evaluated at regular intervals.
    ! delta : real
    !   Length of the interval. 
    ! mask : logical array (n)
    !   Indicate whether or not the value of f is to be taken into account. 
    !   When computing ice velocity derivative, the mask should be true for 
    !   ocean cells and false for land cells. 
    ! 
    ! Output
    ! ------
    ! der : real array(n)
    !   The first derivative evaluated at the same location as f. 
    
    REAL(KIND=RP), DIMENSION(:), INTENT(in) :: f
    LOGICAL, DIMENSION(SIZE(f)), INTENT(in) :: mask
    REAL(KIND=RP), INTENT(in) :: delta
    REAL(KIND=RP), DIMENSION(SIZE(mask)):: der

    INTEGER :: n
    REAL(KIND=RP), DIMENSION(SIZE(f)) :: fmasked
    REAL(KIND=RP), DIMENSION(SIZE(f)-2) :: delta1, delta2

    n = SIZE(mask)

    ! Set masked values to 0. This enforces the no-slip condition. 
    WHERE(mask)
       fmasked = f
    ELSEWHERE
       fmasked = 0.0_RP
    END WHERE

    ! Leftmost element: Forward difference
    IF (mask(1) .AND. mask(2)) THEN
       IF (mask(3)) THEN
          der(1) =  finite_difference_2nd(0.0_RP, 1.0_RP, 2.0_RP, fmasked(1), fmasked(2), fmasked(3), 0.0_RP)
       ELSE
          der(1) =  finite_difference_2nd(0.0_RP, 1.0_RP, 1.5_RP, fmasked(1), fmasked(2), fmasked(3), 0.0_RP)
       ENDIF
    ELSE
       der(1) = 0.0_RP
    END IF



    ! Middle points, central difference
    WHERE(mask(:n-2))
       delta1 = 1.0_RP       ! Ocean in the left grid cell.
    ELSEWHERE
       delta1 = 0.5_RP       ! Land in the left grid cell. 
    END WHERE

    WHERE (mask(3:))
       delta2 = 1.0_RP      ! Ocean in the right grid cell
    ELSEWHERE
       delta2 = 0.5_RP      ! Land in the right grid cell
    END WHERE

    WHERE (mask(2:n-1))
       der(2:n-1) = central_difference(fmasked(1:n-2), fmasked(2:n-1), fmasked(3:n), delta1, delta2)
    ELSEWHERE
       der(2:n-1) = 0.0_RP
    END WHERE

    ! Rightmost element: Backward difference
    IF (mask(n) .AND. mask(n-1)) THEN
       IF (mask(n-2)) THEN
          der(n) = finite_difference_2nd(-2.0_RP, -1.0_RP, 0.0_RP, fmasked(n-2), fmasked(n-1), fmasked(n), 0.0_RP)
       ELSE
          der(n) = finite_difference_2nd(-1.5_RP, -1.0_RP, 0.0_RP, fmasked(n-2), fmasked(n-1), fmasked(n), 0.0_RP)
       ENDIF
    ELSE
       der(n) = 0.0_RP
    END IF

    der = der/delta


  END FUNCTION no_slip_derivative


  SUBROUTINE C_grid_no_slip_gradient(u, v, deltax, deltay, mask, dudx, dvdy, dudy, dvdx)
    ! Return the derivatives of u and v along x and y. 
    ! Assume u and v are located on a staggered C grid. 
    !
    !   +---v---+
    !   |       |
    !   u   n   u
    !   |       |
    !   +---v---+   -> x
    !
    !  The derivatives are returned at the node point. 
    !
    ! Input
    ! -----
    ! u : real array (nx+1, ny)
    !   X-component of some function, eg. the ice velocity.
    ! v : real array (nx, ny+1)
    !   Y-component of some function. 
    ! deltax : real
    !   The grid spacing in the x direction. 
    ! deltay : real
    !   The grid spacing in the y direction. 
    ! mask : logical array (nx, ny)
    !   Where values must be computed, eg. ocean cells. 
    ! 
    ! Output
    ! ------
    ! dudx : real array (nx,ny)
    !   The partial derivative du/dx evaluated at the node point. 
    ! dvdy : real array (nx, ny)
    !   The partial derivative dv/dy.
    ! dudy : real array (nx, ny)
    !   The partial derivative du/dy.
    ! dvdx : real array (nx, ny)
    !   The partial derivative dv/dx.
    

    REAL(KIND=RP), DIMENSION(:,:), INTENT(in) :: u,v
    REAL(KIND=RP), INTENT(in) :: deltax, deltay
    LOGICAL, DIMENSION(SIZE(v,1),SIZE(u,2)), INTENT(in) :: mask
    REAL(KIND=RP), DIMENSION(SIZE(v,1), SIZE(u,2)), INTENT(out) :: dudx, dvdy, dudy, dvdx

    INTEGER :: i, j, nx, ny
    
    nx = SIZE(mask, 1)
    ny = SIZE(mask, 2)
 

    FORALL (j=1:ny)
       ! du/dx
       dudx(1:nx, j) = mid_derivative(u(:,j), deltax, mask(:, j))

       ! dv/dx
       dvdx(1:nx,j) = no_slip_derivative((v(:,j) + v(:,j))/2., deltax, mask(:,j))
    END FORALL


    FORALL (i=1:nx)
       ! dv/dy
       dvdy(i, 1:ny) = mid_derivative(v(i,:), deltay, mask(i,:))
      
       !du/dy
       dudy(i, 1:ny) = no_slip_derivative((u(i, :) + u(i+1, :))/2., deltay, mask(i, :))
    END FORALL
     
  END SUBROUTINE C_grid_no_slip_gradient


  SUBROUTINE C_grid_gradient(s, deltax, deltay, mask, dsdx, dsdy)
    ! Compute the gradient of a scalar quantity defined at the C-grid node
    ! on the C-grid staggered vector locations. 
    !
    ! Input
    ! -----
    ! s : 2D real array
    !   Scalar input. 
    ! deltax : real
    !   Grid width along the x axis. 
    ! deltay : real
    !   Grid width along the y axis.
    ! mask : 2D logical array
    !   Where the scalar is defined. 
    !
    ! Output
    ! ------
    ! dsdx : 2D real array
    !   Partial derivative ds/dx.
    !  dsdy : 2D real array 
    !   Partial derivative ds/dy.

    REAL(KIND=RP), DIMENSION(:,:), INTENT(in) :: s
    REAL(KIND=RP), INTENT(in) :: deltax, deltay
    LOGICAL, DIMENSION(:,:), INTENT(in) :: mask
    REAL(KIND=RP), DIMENSION(SIZE(s,1)-1, SIZE(s,2)), INTENT(out) :: dsdx
    REAL(KIND=RP), DIMENSION(SIZE(s,1), SIZE(s,2)-1), INTENT(out) :: dsdy
    
    INTEGER :: i,j, nx, ny

    nx = SIZE(s, 1)
    ny = SIZE(s, 2)

    FORALL (i=1:nx)
       dsdy(i, :) = mid_derivative(s(i, :), deltay, mask(i, 1:ny-1) .AND. mask(i,2:ny))
    END FORALL
       
    FORALL (j=1:ny)
       dsdx(:,j) = mid_derivative(s(:,j), deltax, mask(1:nx+1, j) .AND. mask(2:nx, j))
    END FORALL
    

  END SUBROUTINE C_grid_gradient



  FUNCTION average_2d(x, mask)
    ! Return the average of 2D array x.
    ! Only the values of x where mask is .True. are counted.
    REAL(KIND=RP), INTENT(in) :: x(:,:)
    LOGICAL, INTENT(in), OPTIONAL :: mask(:,:)
    REAL(KIND=RP) :: average_2d, tmp
    
    IF (PRESENT(mask)) THEN
       tmp = SUM(x, mask=mask)
       average_2d = tmp/COUNT(mask)
    ELSE
       average_2d = SUM(x)/SIZE(x)
    ENDIF
  END FUNCTION average_2d


  ELEMENTAL FUNCTION coriolis_parameter(lat)
    ! Return the Coriolis parameter.
    ! 
    ! .. math::
    !     f = 2\omega\sin\phi
    ! 
    ! where `\omega=7.29211538E-5` rad/s is Earth's angular velocity 
    ! and `\phi` the latitude. 
    !
    ! Input
    ! -----
    ! lat : real
    !   Latitude [degrees]
    !
    ! Output
    ! ------
    ! coriolis_parameter : real
    !   2 \omega \sin \phi
    
    REAL(KIND=RP), INTENT(in) :: lat
    REAL(KIND=RP) :: coriolis_parameter

    coriolis_parameter = 2.0 * 7.29211538E-5 * SIN(radian(lat))

  END FUNCTION coriolis_parameter



  PURE  FUNCTION determinant_2x2(A)
    ! Return the determinant of a 2 by 2 matrix:
    !
    ! ..  math::
    !     \det(A) = A_{11}A_{22} - A_{21}A_{12}

    REAL(KIND=RP), DIMENSION(2,2), INTENT(in) :: A 
    REAL(KIND=RP) :: determinant_2x2

    determinant_2x2 = A(1,1) * A(2,2) - A(2,1) * A(1,2)

  END FUNCTION determinant_2x2

  
  PURE FUNCTION find(x,a)
    ! Return the smallest index `i` for which x(i) == a.
    !
    ! INPUT
    ! -----
    ! x : 1D integer array
    !   The array to search
    ! a : integer
    !   The value to match
    !
    ! OUTPUT
    ! ------
    ! where : integer
    !   The smallest index such that x(i) == a. If a is not found in x, return lbound(x)-1.

    INTEGER, INTENT(in) :: x(:), a
    INTEGER :: find, i

    DO i=LBOUND(x,1), UBOUND(x,1)
       IF (x(i) == a) THEN
          find = i
          RETURN
       END IF
    END DO

    ! Case was not found in x
    find = LBOUND(x,1) - 1

  END FUNCTION find




  ELEMENTAL FUNCTION interpolate(x1, x2, pos) RESULT (out)
    ! Interpolate between two values.
    !
    ! INPUT
    ! -----
    ! x1 : real
    !   First reference value.
    ! x2 : real
    !   Second reference value.
    ! pos : real
    !   Relative distance from x1 to x1 [0.,1.]
    !
    ! OUTPUT
    ! ------
    ! out : real
    !   The value at pos interpolated from x1 and x2:  :math:`x_1*(1-pos) + x_2*pos`

    REAL(KIND=RP), INTENT(in) :: x1, x2, pos
    REAL(KIND=RP) :: out
    out = x1*(1.-pos) + x2*pos
    
  END FUNCTION interpolate


  PURE  FUNCTION inverse_2x2(A)
    ! Return the inverse of a 2 by 2 matrix:
    !
    ! .. math::
    !    A^{-1} = \frac{1}{\det A} \left[ \begin{matrix}
    !               A_{22}  & -A_{12} \\
    !              -A_{21} & A_{11}
    !              \end{matrix} \right].

    REAL(KIND=RP), DIMENSION(2,2), INTENT(in) :: A
    REAL(KIND=RP), DIMENSION(2,2) :: inverse_2x2
    REAL(KIND=RP) :: idet

    idet = 1./determinant_2x2(A)

    inverse_2x2(1,1) =  idet * A(2,2)
    inverse_2x2(1,2) = -idet * A(1,2)
    inverse_2x2(2,1) = -idet * A(2,1)
    inverse_2x2(2,2) =  idet * A(1,1)
    
  END FUNCTION inverse_2x2


  SUBROUTINE B_to_C_grid(uin, vin, uout, vout)
    ! Convert B-grid vectors to C-grid vectors by a simple averaging.
    !
    !
    !
    !   +---v(i,j)--+
    !   |           |
    !   |           |
    !  u(i,j)      u(i+1,j)
    !   |           | 
    !   |           |
    !  u,v--v(i,j+1)+   -> x
    !
    !
    !
    ! Input
    ! -----
    ! uin : 2D real array (nx. ny)
    !   Input x-vector at the B-grid node.
    ! vin : 2D real array (nx, ny)
    !   Input y-vector at the B-grid node.
    !
    ! Output
    ! ------
    ! uout : 2D real array (nx, ny-1)
    !  Output x-vector at the C-grid x location.
    ! vout : 2D real array (nx-1, ny)
    !  Output y-vector at the C-grid y location.
    
    REAL(KIND=RP), DIMENSION(:, :), INTENT(in) :: uin
    REAL(KIND=RP), DIMENSION(SIZE(uin,1), SIZE(uin,2)), INTENT(in) :: vin
    REAL(KIND=RP), INTENT(out) :: uout(SIZE(uin,1), SIZE(uin,2)-1), vout(SIZE(vin,1)-1, SIZE(vin,2))
    INTEGER :: nx, ny

    nx = SIZE(uin, 1)
    ny = SIZE(uin, 2)
    
    uout = (uin(:, 1:ny-1) + uin(:, 2:ny))/2.
    vout = (vin(1:nx-1, :) + vin(2:nx, :))/2.

  END SUBROUTINE B_to_C_grid


  SUBROUTINE C_to_B_grid(uin, vin, uout, vout)
    ! Convert C-grid vectors to B-grid vectors by a simple averaging.
    ! WATCH OUT. THIS MUST BE TESTED THOROUGHLY BEFORE BEING USED.
    ! NOT READY FOR CONSUMPTION.
    !
    !
    !
    !   +---v(i,j)--+
    !   |           |
    !   |           |
    !  u(i,j)      u(i+1,j)
    !   |           | 
    !   |           |
    !   +-v(i,j+1)--B   -> x
    !
    !
    !
    ! Input
    ! -----
    ! uin : 2D real array (nx+1. ny)
    !   Input x-vector at the C-grid staggered locations.
    ! vin : 2D real array (nx, ny+1)
    !   Input y-vector at the C-grid staggered locations.
    !
    ! Output
    ! ------
    ! uout : 2D real array (nx-1, ny-1)
    !  Output x-vector at the B-grid node.
    ! vout : 2D real array (nx-1, ny-1)
    !  Output y-vector at the B-grid node.
    
    REAL(KIND=RP), DIMENSION(:, :), INTENT(in) :: uin
    REAL(KIND=RP), DIMENSION(SIZE(uin,1)-1, SIZE(uin,2)+1), INTENT(in) :: vin
    REAL(KIND=RP), DIMENSION(SIZE(vin,1)-1, SIZE(uin,2)-1), INTENT(inout) :: uout, vout
    INTEGER :: nx, ny

    nx = SIZE(vin, 1)
    ny = SIZE(uin, 2)
    
    uout = (uin(2:nx, 1:ny-1) + uin(2:nx, 2:ny))/2._RP ! Average along y.
    vout = (vin(1:nx-1, 2:ny) + vin(2:nx, 2:ny))/2._RP ! Average along x.


  END SUBROUTINE C_to_B_grid
    





!!$  SUBROUTINE vec2mat (xu,xv,x)
!!$    ! transformer case 0
!!$    
!!$    USE PARAMS
!!$    IMPLICIT NONE
!!$    
!!$    REAL(KIND=RP), DIMENSION(0:nx+2,0:ny+2), INTENT(OUT) :: xu,xv
!!$    REAL(KIND=RP), INTENT(IN) :: x(nvar)
!!$    INTEGER i, j
!!$    
!!$    DO j = 1, ny
!!$       DO i = 1, nx+1
!!$          
!!$          xu(i,j) =  x(i+(j-1)*(nx+1))
!!$          
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    DO j = 1, ny+1
!!$       DO i = 1, nx
!!$          
!!$          xv(i,j) = x(i+(j-1)*nx+(nx+1)*ny) 
!!$          
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    RETURN
!!$  END SUBROUTINE vec2mat
!!$  
!!$    
!!$  SUBROUTINE mat2vec(xu, xv, out)
!!$    ! transformer case 1
!!$    USE PARAMS
!!$    IMPLICIT NONE
!!$    
!!$    REAL(KIND=RP), DIMENSION(0:nx+2,0:ny+2), INTENT(IN) :: xu,xv
!!$    REAL(KIND=RP), INTENT(INOUT) :: out(nvar)
!!$    INTEGER i, j, k
!!$
!!$    DO j = 1, ny
!!$       DO i = 1, nx+1
!!$          k = i+(j-1)*(nx+1)   
!!$          out(k) = xu(i,j)  
!!$          
!!$       ENDDO
!!$    ENDDO
!!$      
!!$    DO j = 1, ny+1
!!$       DO i = 1, nx
!!$          k = i+(j-1)*nx+(nx+1)*ny  
!!$          out(k) = xv(i,j)  
!!$       ENDDO
!!$    ENDDO
!!$    
!!$  END SUBROUTINE mat2vec
!!$    

  SUBROUTINE mypack(m1, m2, vec)
    ! Pack matrices m1 and m2 in vector vec. 

    REAL(KIND=RP), DIMENSION(:,:), INTENT(in) :: m1, m2
    REAL(KIND=RP), DIMENSION(SIZE(m1) + SIZE(m2)), INTENT(out) :: vec

    vec(1:SIZE(m1)) = PACK(m1, .TRUE.)
    vec(SIZE(m1)+1:) = PACK(m2, .TRUE.)

  END SUBROUTINE mypack




END MODULE UTILS
  
