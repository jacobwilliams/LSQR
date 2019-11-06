*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File lsqrchk  fortran   (double precision)
*
*     acheck   xcheck
*
*     These routines are for testing algorithms LSQR and CRAIG
*     (Paige and Saunders, ACM TOMS, 1982).
*
*     acheck  tests if products of the form y := y + Ax and x := x + A'y
*             seem to refer to the same A.
*     xcheck  tests if a solution from LSQR or CRAIG seems to be
*             correct.
*     acheck and xcheck may be helpful to all users of LSQR and CRAIG.
*
*     10 Feb 1992:  acheck revised and xcheck implemented.
*     27 May 1993:  acheck and xcheck kept separate from test problems.
*
*     Michael Saunders, Dept of Operations Research, Stanford University.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine acheck( m, n, nout, aprod, eps,
     $                   leniw, lenrw, iw, rw,
     $                   v, w, x, y,
     $                   inform )

      implicit           none
      external           aprod
      integer            m, n, nout, inform, leniw, lenrw
      integer            iw(leniw)
      double precision   eps
      double precision   v(n), w(m), x(n), y(m), rw(lenrw)

* One-liner: acheck checks the two modes of aprod for LSQR and CRAIG.
*
* Purpose:   acheck may be called to test the user-written subroutine
*     aprod required by LSQR and CRAIG.  For some m x n matrix A,
*     aprod with mode = 1 and 2 supplies LSQR and CRAIG with products
*     of the form
*        y := y + Ax  and  x := x + A'y
*     respectively, where A' means A(transpose).
*     acheck tries to verify that A and A' refer to the same matrix.
*
* Method:    We cook up some "unlikely" vectors x and y of unit length
*     and test if  y'(y + Ax)  =  x'(x + A'y).
*
* Arguments:
*     Arg       Dim     Type    I/O/S Description
*
*     m                 Integer I     No. of rows of A.
*     n                 Integer I     No. of columns of A.
*     nout              Integer I     A file number  for printed output.
*     aprod             External      The routine to be tested.
*                                     See LSQR or CRAIG.
*     eps               Double  I     The machine precision.
*     leniw             Integer I     These four parameters are passed
*     lenrw             Integer I     to aprod but not otherwise used.
*     iw        leniw   Integer I     LSQR and CRAIG pass them to aprod
*     rw        lenrw   Double  I     in the same way.
*     v         n       Double      S
*     w         m       Double      S
*     x         n       Double      S
*     y         m       Double      S
*     inform            Integer   O   Error indicator.
*                                     inform = 0 if aprod seems to be
*                                     consistent.
*                                     inform = 1 otherwise.
*
* Error Handling:  See inform above.
*
* Parameter Constants:
*     Param   Type   Description
*     one     double 1.0d+0
*     power   double eps**power is used as the tolerance for judging
*                    whether    y'(y + Ax)  =  x'(x + A'y)
*                    to sufficient accuracy.
*                    power should be in the range (0.25, 0.9) say.
*                    For example, power = 0.75 means that we are happy
*                    if three quarters of the available digits agree.
*                    power = 0.5 seems a reasonable requirement
*                    (asking for half the digits to agree).
*                    
*
* Files Used:
*     nout       O   Screen diagnostics
* 
* Procedures:
*     dcopy    BLAS
*     ddot     BLAS
*     dnrm2    BLAS
*     dscal    BLAS
*
* Environment:
*     FORTRAN 77 with implicit none.
*
* History:
*     04 Sep 1991  Initial design and code.
*                  Michael Saunders, Dept of Operations Research,
*                  Stanford University
*     10 Feb 1992  aprod and eps added as parameters.
*                  tol defined via power.
C-----------------------------------------------------------------------

*     Procedures.

      external           dcopy, ddot, dnrm2, dscal
      double precision                 ddot, dnrm2

*     Local constants.

      double precision   one,          power
      parameter         (one = 1.0d+0, power = 0.5d+0)

*     Local variables.

      integer            i, j
      double precision   alfa, beta, t, test1, test2, test3, tol


*     Execution.

      tol    = eps**power
      if (nout .gt. 0) write(nout, 1000)

*     ==================================================================
*     Cook up some "unlikely" vectors x and y of unit length.
*     ==================================================================
      t     = one
      do 20 j = 1, n
         t    = t + one
         x(j) = sqrt( t )
   20 continue
      
      t     = one
      do 30 i = 1, m
         t    = t + one
         y(i) = one / sqrt( t )
   30 continue
      
      alfa   =   dnrm2 ( n, x, 1 )
      beta   =   dnrm2 ( m, y, 1 )
      call dscal ( n, (one / alfa), x, 1 )
      call dscal ( m, (one / beta), y, 1 )
      
*     ==================================================================
*     Test if       y'(y + Ax)  =  x'(x + A'y).
*     ==================================================================

*     First set    w = y + Ax,    v = x + A'y.

      call dcopy ( m, y, 1, w, 1 )
      call dcopy ( n, x, 1, v, 1 )
      call aprod ( 1, m, n, x, w, leniw, lenrw, iw, rw )
      call aprod ( 2, m, n, v, y, leniw, lenrw, iw, rw )
      
*     Now set      alfa = y'w,    beta = x'v.
      
      alfa    =   ddot  ( m, y, 1, w, 1 )
      beta    =   ddot  ( n, x, 1, v, 1 )
      test1   =   abs( alfa - beta )            
      test2   =   one  +  abs( alfa )  +  abs( beta )
      test3   =   test1 / test 2

*     See if alfa and beta are essentially the same.

      if ( test3 .le. tol ) then
         inform = 0
         if (nout .gt. 0) write(nout, 1010) test3
      else
         inform = 1
         if (nout .gt. 0) write(nout, 1020) test3
      end if

      return

 1000 format(//' Enter acheck.     Test of aprod for LSQR and CRAIG')
 1010 format(  ' aprod seems OK.   Relative error =', 1p, e10.1)
 1020 format(  ' aprod seems incorrect.   Relative error =', 1p, e10.1)
              
*     end of acheck
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine xcheck( m, n, nout, aprod, anorm, damp, eps, 
     $                   leniw, lenrw, iw, rw,
     $                   b, u, v, w, x,
     $                   inform, test1, test2, test3 )

      implicit           none
*     implicit           double precision (a-h,o-z)
      external           aprod
      integer            m, n, nout, inform, leniw, lenrw
      integer            iw(leniw)
      double precision   anorm, damp, eps, test1, test2, test3
      double precision   rw(lenrw), b(m), u(m), v(n), w(n), x(n)

*     ------------------------------------------------------------------
* One-liner: xcheck tests if x solves a certain least-squares problem.
*
* Purpose:   xcheck computes residuals and norms associated with the
*     vector x and the least-squares problem solved by LSQR or CRAIG.
*     It determines whether x seems to be a solution to any of three
*     possible systems:  1.  Ax = b
*                        2.  min norm(Ax - b)
*                        3.  min norm(Ax - b)^2 + damp^2 * norm(x)^2.
*
* Arguments:
*     Arg       Dim     Type    I/O/S Description
*
*     m                 Integer I     The number of rows in A.
*     n                 Integer I     The number of columns in A.
*     nout              Integer I     A file number for printed output.
*                                     If nout = 0, nothing is printed.
*     aprod             External      The subroutine defining A.
*                                     See LSQR or CRAIG.
*     anorm             Double  I     An estimate of norm(A) or
*                                     norm( A, delta*I ) if delta > 0.
*                                     Normally this will be available
*                                     from LSQR or CRAIG.  
*     damp              Double  I     Possibly defines a damped problem.
*     eps               Double  I     Machine precision.
*     leniw             Integer I     These four parameters are passed
*     lenrw             Integer I     to aprod but not otherwise used.
*     iw        leniw   Integer I     LSQR and CRAIG pass them to aprod
*     rw        lenrw   Double  I     in the same way.
*     b         m       Double  I     The right-hand side of Ax = b etc.
*     u         m       Double    O   On exit, u = r (where r = b - Ax).
*     v         n       Double    O   On exit, v = A'r.
*     w         n       Double    O   On exit, w = A'r - damp^2 x.
*     x         n       Double  I     The given estimate of a solution.
*     inform            Integer   O   inform = 0 if b = 0 and x = 0.
*                                     inform = 1, 2 or 3 if x seems to
*                                     solve systems 1 2 or 3 above.
*     test1             Double    O   These are dimensionless quantities
*     test2             Double    O   that should be "small" if x does
*     test3             Double    O   seem to solve one of the systems.
*                                     "small" means less than
*                                     tol = eps**power, where power is
*                                     defined as a parameter below.
*
* Files Used:
*     nout       O   Screen diagnostics
* 
* Procedures:
*     dcopy    BLAS
*     dnrm2    BLAS
*     dscal    BLAS
*
* Environment:
*     FORTRAN 77 with implicit none.
*
* History:
*     07 Feb 1992  Initial design and code.
*                  Michael Saunders, Dept of Operations Research,
*                  Stanford University.
*     ------------------------------------------------------------------

*     Procedures.

      external           dcopy, dnrm2, dscal
      double precision          dnrm2

*     Local constants.

      double precision   zero,          one,          power
      parameter         (zero = 0.0d+0, one = 1.0d+0, power = 0.5d+0)

*     Local variables.

      integer            j
      double precision   bnorm, dampsq, rho1, rho2, sigma1, sigma2, 
     $                   tol, snorm, xnorm, xsnorm


*     Execution.

      dampsq = damp**2
      tol    = eps**power

*     Compute  u = b - Ax   via   u = -b + Ax,  u = -u.
*     This is usual residual vector r.

      call dcopy ( m, b, 1, u, 1 )
      call dscal ( m, (-one), u, 1 )
      call aprod ( 1, m, n, x, u, leniw, lenrw, iw, rw )
      call dscal ( m, (-one), u, 1 )

*     Compute  v = A'u   via   v = 0,  v = v + A'u. 

      do 100 j = 1, n
         v(j)  = zero
  100 continue
      call aprod ( 2, m, n, v, u, leniw, lenrw, iw, rw )

*     Compute  w = A'u  -  damp**2 * x.
*     This will be close to zero in all cases
*     if  x  is close to a solution.

      call dcopy ( n, v, 1, w, 1 )
      if (damp .ne. zero) then
         do 200 j = 1, n
            w(j)  = w(j)  -  dampsq * x(j)
  200    continue
      end if

*     Compute the norms associated with  b, x, u, v, w.

      bnorm  = dnrm2 ( m, b, 1 )
      xnorm  = dnrm2 ( n, x, 1 )
      rho1   = dnrm2 ( m, u, 1 )
      sigma1 = dnrm2 ( n, v, 1 )
      if (nout .gt. 0) write(nout, 2200) damp, xnorm, rho1, sigma1

      if (damp .eq. zero) then
         rho2   = rho1
         sigma2 = sigma1
      else
         rho2   = sqrt( rho1**2  +  dampsq * xnorm**2 )
         sigma2 = dnrm2 ( n, w, 1 )
         snorm  = rho1 / damp
         xsnorm = rho2 / damp
         if (nout .gt. 0) write(nout, 2300) snorm, xsnorm, rho2, sigma2
      end if

*     See if x seems to solve Ax = b  or  min norm(Ax - b)
*                    or the damped least-squares system.

      if (bnorm .eq. zero  .and.  xnorm .eq. zero) then
         inform = 0
         test1  = zero
         test2  = zero
         test3  = zero
      else
         inform = 4
         test1  = rho1 / (bnorm  +  anorm * xnorm)
         test2  = zero
         if (rho1 .gt. zero) test2  = sigma1 / (anorm * rho1)
         test3  = test2
         if (rho2 .gt. zero) test3  = sigma2 / (anorm * rho2)

         if (test3 .le. tol) inform = 3
         if (test2 .le. tol) inform = 2
         if (test1 .le. tol) inform = 1
      end if

      if (nout .gt. 0)
     $   write(nout, 3000) inform, tol, test1, test2, test3
      return

 2200 format(1p
     $ // ' Enter xcheck.     Does x solve Ax = b, etc?'
     $ /  '    damp            =', e10.3
     $ /  '    norm(x)         =', e10.3
     $ /  '    norm(r)         =', e15.8, ' = rho1'
     $ /  '    norm(A''r)       =',e10.3, '      = sigma1')
 2300 format(1p
     $ /  '    norm(s)         =', e10.3
     $ /  '    norm(x,s)       =', e10.3
     $ /  '    norm(rbar)      =', e15.8, ' = rho2'
     $ /  '    norm(Abar''rbar) =',e10.3, '      = sigma2')
 3000 format(1p
     $ /  '    inform          =', i2
     $ /  '    tol             =', e10.3
     $ /  '    test1           =', e10.3, ' (Ax = b)'
     $ /  '    test2           =', e10.3, ' (least-squares)'
     $ /  '    test3           =', e10.3, ' (damped least-squares)')

*     end of xcheck
      end
