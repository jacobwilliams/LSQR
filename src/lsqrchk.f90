!***************************************************************************************************
!>
!  These routines are for testing algorithms [[LSQR]] and [[CRAIG]]
!  (Paige and Saunders, ACM TOMS, 1982).
!
!  * [[acheck]] tests if products of the form y := y + Ax and x := x + A'y
!    seem to refer to the same A.
!  * [[xcheck]] tests if a solution from [[LSQR]] or [[CRAIG]] seems to be
!    correct.
!
!  [[acheck]] and [[xcheck]] may be helpful to all users of [[LSQR]] and [[CRAIG]].
!
!### History
!  * 10 Feb 1992: acheck revised and xcheck implemented.
!  * 27 May 1993: acheck and xcheck kept separate from test problems.
!
!### Author
!  * Michael Saunders, Dept of Operations Research, Stanford University.

   module lsqr_check

   use lsqr_kinds
   use lsqpblas_module

   implicit none

   contains
!***************************************************************************************************

!***************************************************************************************************
!>
!  Checks the two modes of aprod for LSQR and CRAIG.
!
!  acheck may be called to test the user-written subroutine
!  aprod required by LSQR and CRAIG.  For some m x n matrix A,
!  aprod with mode = 1 and 2 supplies LSQR and CRAIG with products
!  of the form
!     y := y + Ax  and  x := x + A'y
!  respectively, where A' means A(transpose).
!  acheck tries to verify that A and A' refer to the same matrix.
!
!### Method
!  We cook up some "unlikely" vectors x and y of unit length
!  and test if  y'(y + Ax)  =  x'(x + A'y).
!
!### History
! * 04 Sep 1991  Initial design and code.
!   Michael Saunders, Dept of Operations Research,
!   Stanford University
! * 10 Feb 1992  aprod and eps added as parameters.
!   tol defined via power.

subroutine acheck( m, n, nout, aprod, eps, &
                   v, w, x, y, inform )

implicit none

external :: aprod
integer,intent(in)   :: m       !! No. of rows of A.
integer,intent(in)   :: n       !! No. of columns of A.
integer,intent(in)   :: nout    !! A file number for printed output.
integer,intent(out)  :: inform  !! Error indicator.
                                !! inform = 0 if aprod seems to be
                                !! consistent.
                                !! inform = 1 otherwise.
real(wp),intent(in)  :: eps     !! The machine precision.
real(wp)             :: v(n)
real(wp)             :: w(m)
real(wp)             :: x(n)
real(wp)             :: y(m)

real(wp),parameter :: power = 0.5_wp   !! eps**power is used as the tolerance for judging
                                       !! whether `y'(y + Ax)  =  x'(x + A'y)`
                                       !! to sufficient accuracy.
                                       !! power should be in the range (0.25, 0.9) say.
                                       !! For example, power = 0.75 means that we are happy
                                       !! if three quarters of the available digits agree.
                                       !! power = 0.5 seems a reasonable requirement
                                       !! (asking for half the digits to agree).

integer :: i, j
real(wp) :: alfa, beta, t, test1, test2, test3, tol

tol = eps**power
if (nout /= 0) write(nout, '(A)') &
   'Enter acheck. Test of aprod for LSQR and CRAIG'

! ==================================================================
! Cook up some "unlikely" vectors x and y of unit length.
! ==================================================================
t = one
do j = 1, n
   t    = t + one
   x(j) = sqrt( t )
end do

t = one
do i = 1, m
   t    = t + one
   y(i) = one / sqrt( t )
end do

alfa = dnrm2 ( n, x, 1 )
beta = dnrm2 ( m, y, 1 )
call dscal ( n, (one / alfa), x, 1 )
call dscal ( m, (one / beta), y, 1 )

! ==================================================================
! Test if y'(y + Ax) = x'(x + A'y).
! ==================================================================

! First set w = y + Ax, v = x + A'y.

call dcopy ( m, y, 1, w, 1 )
call dcopy ( n, x, 1, v, 1 )
call aprod ( 1, m, n, x, w )
call aprod ( 2, m, n, v, y )

! Now set alfa = y'w, beta = x'v.

alfa  = ddot  ( m, y, 1, w, 1 )
beta  = ddot  ( n, x, 1, v, 1 )
test1 = abs( alfa - beta )
test2 = one  +  abs( alfa )  +  abs( beta )
test3 = test1 / test2

! See if alfa and beta are essentially the same.

if ( test3 <= tol ) then
   inform = 0
   if (nout /= 0) write(nout, '(A,1X,E10.1)') &
      'aprod seems OK. Relative error =', test3
else
   inform = 1
   if (nout /= 0) write(nout, '(A,1X,E10.1)') &
      'aprod seems incorrect. Relative error =', test3
end if

end subroutine acheck
!***************************************************************************************************

!***************************************************************************************************
!>
!  Tests if x solves a certain least-squares problem.
!
!  xcheck computes residuals and norms associated with the
!  vector x and the least-squares problem solved by LSQR or CRAIG.
!  It determines whether x seems to be a solution to any of three
!  possible systems:
!
!  1.  Ax = b
!  2.  min norm(Ax - b)
!  3.  min norm(Ax - b)^2 + damp^2 * norm(x)^2
!
!### History
! * 07 Feb 1992  Initial design and code.
!   Michael Saunders, Dept of Operations Research,
!   Stanford University.

subroutine xcheck( m, n, nout, aprod, anorm, damp, eps, &
                   b, u, v, w, x, &
                   inform, test1, test2, test3 )

implicit none

external :: aprod !! The subroutine defining A.
                  !! See LSQR or CRAIG.
integer,intent(in)   :: m  !! The number of rows in A.
integer,intent(in)   :: n  !! The number of columns in A.
integer,intent(in)   :: nout  !! A file number for printed output.
                              !! If nout = 0, nothing is printed.
integer,intent(out)  :: inform   !! inform = 0 if b = 0 and x = 0.
                                 !! inform = 1, 2 or 3 if x seems to
                                 !! solve systems 1 2 or 3 above.
real(wp),intent(in)  :: anorm    !! An estimate of norm(A) or
                                 !! norm( A, delta*I ) if delta > 0.
                                 !! Normally this will be available
                                 !! from LSQR or CRAIG.
real(wp),intent(in)  :: damp     !! Possibly defines a damped problem.
real(wp),intent(in)  :: eps      !! Machine precision.
real(wp),intent(out) :: test1    !! These are dimensionless quantities
                                 !! that should be "small" if x does
                                 !! seem to solve one of the systems.
                                 !! "small" means less than
                                 !! tol = eps**power, where power is
                                 !! defined as a parameter below.
real(wp),intent(out) :: test2    !! These are dimensionless quantities
                                 !! that should be "small" if x does
                                 !! seem to solve one of the systems.
                                 !! "small" means less than
                                 !! tol = eps**power, where power is
                                 !! defined as a parameter below.
real(wp),intent(out) :: test3    !! These are dimensionless quantities
                                 !! that should be "small" if x does
                                 !! seem to solve one of the systems.
                                 !! "small" means less than
                                 !! tol = eps**power, where power is
                                 !! defined as a parameter below.
real(wp),intent(in)  :: b(m)     !! The right-hand side of Ax = b etc.
real(wp),intent(out) :: u(m)     !! On exit, u = r (where r = b - Ax).
real(wp),intent(out) :: v(n)     !! On exit, v = A'r.
real(wp),intent(out) :: w(n)     !! On exit, w = A'r - damp^2 x.
real(wp),intent(in)  :: x(n)     !! The given estimate of a solution.

real(wp),parameter :: power = 0.5_wp

integer :: j
real(wp) :: bnorm, dampsq, rho1, rho2, sigma1, sigma2, &
            tol, snorm, xnorm, xsnorm

dampsq = damp**2
tol    = eps**power

! Compute u = b - Ax via u = -b + Ax, u = -u.
! This is usual residual vector r.

call dcopy ( m, b, 1, u, 1 )
call dscal ( m, (-one), u, 1 )
call aprod ( 1, m, n, x, u )
call dscal ( m, (-one), u, 1 )

! Compute v = A'u via v = 0, v = v + A'u.

do j = 1, n
   v(j)  = zero
end do
call aprod ( 2, m, n, v, u )

! Compute w = A'u - damp**2 * x.
! This will be close to zero in all cases
! if x is close to a solution.

call dcopy ( n, v, 1, w, 1 )
if (damp /= zero) then
   do j = 1, n
      w(j)  = w(j)  -  dampsq * x(j)
   end do
end if

! Compute the norms associated with b, x, u, v, w.

bnorm  = dnrm2 ( m, b, 1 )
xnorm  = dnrm2 ( n, x, 1 )
rho1   = dnrm2 ( m, u, 1 )
sigma1 = dnrm2 ( n, v, 1 )
if (nout /= 0) then
   write(nout, '(A)') 'Enter xcheck. Does x solve Ax = b, etc?'
   write(nout, '(A,E10.3)')      ' damp        =', damp
   write(nout, '(A,E10.3)')      ' norm(x)     =', xnorm
   write(nout, '(A,E15.8,A)')    ' norm(r)     =', rho1,   ' = rho1'
   write(nout, '(A,E10.3,5X,A)') ' norm(A''r)  =', sigma1, ' = sigma1'
end if

if (damp == zero) then
   rho2   = rho1
   sigma2 = sigma1
else
   rho2   = sqrt( rho1**2  +  dampsq * xnorm**2 )
   sigma2 = dnrm2 ( n, w, 1 )
   snorm  = rho1 / damp
   xsnorm = rho2 / damp
   if (nout /= 0) then
      write(nout, '(A,E10.3)')      ' norm(s)         =', snorm
      write(nout, '(A,E10.3)')      ' norm(x,s)       =', xsnorm
      write(nout, '(A,E15.8,A)')    ' norm(rbar)      =', rho2,    ' = rho2'
      write(nout, '(A,E10.3,5X,A)') ' norm(Abar''rbar) =', sigma2, ' = sigma2'
   end if
end if

! See if x seems to solve Ax = b or min norm(Ax - b)
! or the damped least-squares system.

if (bnorm == zero  .and.  xnorm == zero) then
   inform = 0
   test1  = zero
   test2  = zero
   test3  = zero
else
   inform = 4
   test1  = rho1 / (bnorm  +  anorm * xnorm)
   test2  = zero
   if (rho1 > zero) test2  = sigma1 / (anorm * rho1)
   test3  = test2
   if (rho2 > zero) test3  = sigma2 / (anorm * rho2)
   if (test3 <= tol) inform = 3
   if (test2 <= tol) inform = 2
   if (test1 <= tol) inform = 1
end if

if (nout /= 0) then
   write(nout, '(A,I2)')      ' inform          =', inform
   write(nout, '(A,E10.3)')   ' tol             =', tol
   write(nout, '(A,E10.3,A)') ' test1           =', test1, ' (Ax = b)'
   write(nout, '(A,E10.3,A)') ' test2           =', test2, ' (least-squares)'
   write(nout, '(A,E10.3,A)') ' test3           =', test3, ' (damped least-squares)'
end if

end subroutine xcheck
!***************************************************************************************************

!***************************************************************************************************
   end module lsqr_check
!***************************************************************************************************
