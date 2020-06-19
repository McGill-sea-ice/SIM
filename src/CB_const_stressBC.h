!========================================================================
!     Common block const_stressBC: program constants
!========================================================================

      integer iseed

      double precision sig11bc, sig22bc, sig12bc, hlevel, bpfactor, ecboat


      common/const_stressBC/             &
                sig11bc,        & ! normal stress BC (constant)
	        sig22bc,        & ! normal stress BC (constant)
		sig12bc,        & ! shear  stress BC (constant)
		hlevel,         & ! level ice thickness (stress BC)
                bpfactor,       & ! factor for boat ice strength
                ecboat,         & ! e aspect ratio for boat contour
		iseed             ! used by random nb generator
