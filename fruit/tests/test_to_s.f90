PROGRAM test_to_s

USE fruit_util, only : to_s


LOGICAL :: lscalar, lshort1d(5), llong1d(31) 
INTEGER :: iscalar, ishort1d(5), ilong1d(11), ishort2d(4,4), ilong2d(20,20)
INTEGER :: imixed1(4, 11), imixed2(11,4)
REAL :: rscalar, rshort1d(5), rlong1d(11), rshort2d(5,5), rlong2d(20,20)
DOUBLE PRECISION :: dscalar, dshort1d(5), dlong1d(11), dshort2d(5,5), dlong2d(20,20)
COMPLEX :: cscalar, cshort1d(5), clong1d(11), cshort2d(5,5), clong2d(20,20)


lscalar = .TRUE.
lshort1d = .TRUE.
llong1d = .FALSE.

iscalar = 5
ishort1d = (/ (i-2, i=1,5) /)
ilong1d =   (/ (i-2, i=1,11) /)
ishort2d = RESHAPE((/ (i-2, i=1,16) /), SHAPE(ishort2d))
ilong2d = RESHAPE((/ (i-2, i=1,400) /), SHAPE(ilong2d))
imixed1 = RESHAPE((/ (i-2, i=1,44) /), SHAPE(imixed1))
imixed2 = RESHAPE((/ (i-2, i=1,44) /), SHAPE(imixed2))

rscalar = -5.2
rshort1d = (/ (i-2.45, i=1,5) /)
rlong1d = 2.2
rshort2d = RESHAPE((/ (i-2.45, i=1,25) /), SHAPE(rshort2d))
rlong2d = 4.2

dscalar = 5.2
dshort1d = 1.2
dlong1d = 2.2
dshort2d = 3.2
dlong2d = 4.2

cscalar = 5.2
cshort1d = 1.2
clong1d = 2.2
cshort2d = 3.2
clong2d = 4.2


WRITE(*,*) 'logical scalar: ', trim(to_s(lscalar))
WRITE(*,*) 'short 1D logical: ', TRIM(to_s(lshort1d))
WRITE(*,*) 'long 1D logical: ', TRIM(to_s(llong1d))

WRITE(*,*) 'scalar integer: ', trim(to_s (iscalar))
WRITE(*,*) 'scalar real: ', trim(to_s (rscalar))
WRITE(*,*) 'scalar dp: ', trim(to_s (dscalar))
WRITE(*,*) 'scalar complex: ', trim(to_s (cscalar))

WRITE(*,*) 'short 1D int: ', trim(to_s (ishort1d))
WRITE(*,*) 'short 1D real: ', trim(to_s (rshort1d))


WRITE(*,*) 'long 1D int: ', trim(to_s (ilong1d))
WRITE(*,*) 'long 1D real: ', trim(to_s (rlong1d))


WRITE(*,*) 'short 2D int: ', trim(to_s (ishort2d))
WRITE(*,*) 'short 2D real: ', trim(to_s (rshort2d))

WRITE(*,*) 'long 2D int: ', trim(to_s (ilong2d))
WRITE(*,*) 'long 2D real: ', trim(to_s (rlong2d))

WRITE(*,*) 'mixed1 2D int: ', trim(to_s (imixed1)) 
WRITE(*,*) 'mixed2 2D int: ', trim(to_s (imixed2)) 


END PROGRAM test_to_s
