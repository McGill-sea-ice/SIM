! forms the velocity field used to linearize the system of equations at ite k
! A(ul^k)u^k = b(ul^k)

      subroutine uv_linearization (uk2, vk2, ul, vl)

      implicit none
      
      include 'parameter.h'
      include 'CB_DynVariables.h'
      include 'CB_options.h'

      double precision, intent(in) :: uk2(0:nx+2,0:ny+2), vk2(0:nx+2,0:ny+2)
      double precision, intent(inout) :: ul(0:nx+2,0:ny+2), vl(0:nx+2,0:ny+2)

      if (linearization .eq. 'Zhang') then

         ul = ( uk2 + uice ) / 2d0 ! these are vectors, uice = u^k-1
         vl = ( vk2 + vice ) / 2d0

      elseif (linearization .eq. 'Tremblay') then
         
         ul = ( ul + uice ) / 2d0  ! these are vectors, uice = u^k-1
         vl = ( vl + vice ) / 2d0

      endif

      return
    end subroutine uv_linearization







