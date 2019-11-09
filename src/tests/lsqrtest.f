*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File lsqrtst  fortran  (double precision)
*
*     aprod    aprod1   aprod2   hprod    lstp     test     MAIN
*
*     These routines define a class of least-squares test problems
*     for testing algorithms LSQR and CRAIG
*     (Paige and Saunders, ACM TOMS, 1982).
*
*     1982---1991:  Various versions implemented.
*     06 Feb 1992:  Test-problem generator lstp generalized to allow
*                   any m and n.  lstp is now the same as the generator
*                   for LSQR and CRAIG.
*     30 Nov 1993:  Modified lstp.
*                   For a while, damp = 0 implied r = damp*s = 0.
*                   This was a result of generating x and s.
*                   Reverted to generating x and r as in LSQR paper.
*
*     Michael Saunders, Dept of Operations Research, Stanford University.
*------------------------------------------------------------------------

      subroutine aprod ( mode, m, n, x, y, leniw, lenrw, iw, rw )

      implicit           double precision (a-h,o-z)
      integer            mode, m, n, leniw, lenrw
      integer            iw(leniw)
      double precision   x(n), y(m), rw(lenrw)

*     ------------------------------------------------------------------
*     This is the matrix-vector product routine required by subroutines
*     LSQR and CRAIG for a test matrix of the form  A = HY*D*HZ.
*     The quantities defining D, HY, HZ are in the work array rw,
*     followed by a work array w.  These are passed to aprod1 and aprod2
*     in order to make the code readable.
*     ------------------------------------------------------------------

      integer            locd, lochy, lochz, locw, maxmn, minmn

      maxmn  = max( m, n )
      minmn  = min( m, n )
      locd   = 1
      lochy  = locd  + minmn
      lochz  = lochy + m
      locw   = lochz + n

      if (mode .eq. 1) then
         call aprod1( m, n, maxmn, minmn, x, y,
     $                rw(locd), rw(lochy), rw(lochz), rw(locw) )
      else
         call aprod2( m, n, maxmn, minmn, x, y,
     $                rw(locd), rw(lochy), rw(lochz), rw(locw) )
      end if

*     end of aprod
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine aprod1( m, n, maxmn, minmn, x, y, d, hy, hz, w )

      implicit           double precision (a-h,o-z)
      integer            m, n, maxmn, minmn
      double precision   x(n), y(m), d(minmn), hy(m), hz(n), w(maxmn)

*     ------------------------------------------------------------------
*     aprod1  computes  y = y + A*x  for subroutine aprod,
*     where A is a test matrix of the form  A = HY*D*HZ,
*     and the latter matrices HY, D, HZ are represented by
*     input vectors with the same name.
*     ------------------------------------------------------------------

      integer            i
      double precision   zero
      parameter        ( zero = 0.0d+0 )

      call hprod ( n, hz, x, w )

      do 100 i = 1, minmn
         w(i)  = d(i) * w(i)
  100 continue

      do 200 i = n + 1, m
         w(i)  = zero
  200 continue

      call hprod ( m, hy, w, w )

      do 600 i = 1, m
         y(i)  = y(i) + w(i)
  600 continue

*     end of aprod1
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine aprod2( m, n, maxmn, minmn, x, y, d, hy, hz, w )

      implicit           double precision (a-h,o-z)
      integer            m, n, maxmn, minmn
      double precision   x(n), y(m), d(minmn), hy(m), hz(n), w(maxmn)

*     ------------------------------------------------------------------
*     aprod2  computes  x = x + A(t)*y  for subroutine aprod,
*     where  A  is a test matrix of the form  A = HY*D*HZ,
*     and the latter matrices  HY, D, HZ  are represented by
*     input vectors with the same name.
*     ------------------------------------------------------------------

      integer            i
      double precision   zero
      parameter        ( zero = 0.0d+0 )

      call hprod ( m, hy, y, w )

      do 100 i = 1, minmn
         w(i)  = d(i)*w(i)
  100 continue

      do 200 i = m + 1, n
         w(i)  = zero
  200 continue

      call hprod ( n, hz, w, w )

      do 600 i = 1, n
         x(i)  = x(i) + w(i)
  600 continue

*     end of aprod2
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine hprod ( n, hz, x, y )

      implicit           double precision (a-h,o-z)
      integer            n
      double precision   hz(n), x(n), y(n)

*     ------------------------------------------------------------------
*     hprod  applies a Householder transformation stored in  hz
*     to get  y = ( I - 2*hz*hz(transpose) ) * x.
*     ------------------------------------------------------------------

      integer            i
      double precision   s

      s      = 0.0
      do 100 i = 1, n
         s     = hz(i) * x(i)  +  s
  100 continue

      s      = s + s
      do 200 i = 1, n
         y(i)  = x(i)  -  s * hz(i)
  200 continue

*     end of hprod
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lstp  ( m, n, maxmn, minmn, nduplc, npower, damp, x,
     $                   b, d, hy, hz, w, acond, rnorm )

      implicit           double precision (a-h,o-z)
      integer            m, n, maxmn, minmn, nduplc, npower
      double precision   damp, acond, rnorm
      double precision   b(m), x(n), d(minmn), hy(m), hz(n), w(maxmn)

*     ------------------------------------------------------------------
*     lstp  generate a sparse least-squares test problem of the form
*                (   A    )*x = ( b ) 
*                ( damp*I )     ( 0 )
*     for solution by LSQR, or a sparse underdetermined system
*                   Ax + damp*s = b
*     for solution by CRAIG.  The matrix A is m by n and is
*     constructed in the form  A = HY*D*HZ,  where D is an m by n
*     diagonal matrix, and HY and HZ are Householder transformations.
*
*     m and n may contain any positive values.
*     The first 8 parameters are input to lstp.  The last 8 are output.
*     If m .ge. n  or  damp = 0, the true solution is x as given.
*     Otherwise, x is modified to contain the true solution.
*
*     Functions and subroutines
*
*     testprob           aprod1, hprod
*     blas               dnrm2 , dscal
*     ------------------------------------------------------------------

*     Intrinsics and local variables

      intrinsic          cos,  min, sin, sqrt
      integer            i, j
      double precision   dnrm2
      double precision   alfa, beta, dampsq, fourpi, t
      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

*     ------------------------------------------------------------------
*     Make two vectors of norm 1.0 for the Householder transformations.
*     fourpi  need not be exact.
*     ------------------------------------------------------------------
      minmn  = min( m, n )
      dampsq = damp**2
      fourpi = 4.0 * 3.141592
      alfa   = fourpi / m
      beta   = fourpi / n

      do 100 i = 1, m
         hy(i) = sin( i * alfa )
  100 continue

      do 200 i = 1, n
         hz(i) = cos( i * beta )
  200 continue

      alfa   = dnrm2 ( m, hy, 1 )
      beta   = dnrm2 ( n, hz, 1 )
      call dscal ( m, (- one / alfa), hy, 1 )
      call dscal ( n, (- one / beta), hz, 1 )
*
*     ------------------------------------------------------------------
*     Set the diagonal matrix  D.  These are the singular values of  A.
*     ------------------------------------------------------------------
      do 300 i = 1, minmn
         j     = (i - 1 + nduplc) / nduplc
         t     =  j * nduplc
         t     =  t / minmn
         d(i)  =  t**npower
  300 continue

      acond  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
      acond  = sqrt( acond )

*     ------------------------------------------------------------------
*     Set the true solution   x.
*     It must be of the form  x = Z ( w )  for some  w.
*                                   ( 0 )
*     ------------------------------------------------------------------
      call hprod ( n, hz, x, w )

      do 420 i = m + 1, n
         w(i)  = zero
  420 continue

      call hprod ( n, hz, w, x )

*     Solve D r1bar = dampsq x1bar 
*     where r1bar and x1bar are both in w.

      do 440 i = 1, minmn
         w(i)  = dampsq * w(i) / d(i)
  440 continue

*     Set r2bar to be anything.  (It is empty if m .le. n)
*     Then form r = Y rbar (again in w).

      do 450 i = minmn+1, m
         w(i)  = one
  450 continue    
      
      call Hprod ( m, hy, w, w )

*     Compute the rhs    b = r  +  Ax.

      rnorm  = dnrm2 ( m, w, 1 )
      call dcopy ( m, w, 1, b, 1 )
      call aprod1( m, n, maxmn, minmn, x, b, d, hy, hz, w )

*     end of lstp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine test  ( m, n, nduplc, npower, damp )

      implicit           double precision (a-h,o-z)
      integer            m, n, nduplc, npower
      double precision   damp

*     ------------------------------------------------------------------
*     This is an example driver routine for running LSQR.
*     It generates a test problem, solves it, and examines the results.
*     Note that subroutine aprod must be declared external
*     if it is used only in the call to LSQR (and acheck).
*
*     1982---1991:  Various versions implemented.
*     04 Sep 1991:  "wantse" added to argument list of LSQR,
*                   making standard errors optional.
*     10 Feb 1992:  Revised for use with lsqrchk fortran.
*
*     Michael Saunders, Dept of Operations Research, Stanford University.
*------------------------------------------------------------------------

      external           acheck, aprod, xcheck
      integer            inform, istop, itnlim, j, nout
      logical            wantse
      double precision   dnrm2

      parameter        ( maxm = 2000,  maxn = 2000, mxmn = 2000 )
      double precision   b(maxm),  u(maxm),
     $                   v(maxn),  w(mxmn), x(maxn),
     $                   se(maxn), xtrue(maxn), y(mxmn)
      double precision   atol, btol, conlim,
     $                   anorm, acond, rnorm, arnorm,
     $                   enorm, etol, xnorm

      parameter        ( leniw = 1,  lenrw = 10000 )
      integer            iw(leniw)
      double precision   rw(lenrw)
      integer            locd, lochy, lochz, locw, ltotal

      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      character*34       line
      data               line
     $                 /'----------------------------------'/


      if (m .gt. maxm  .or.  n .gt. maxn) go to 800

*     Set the output unit and the machine precision.

      nout   = 6
      eps    = 2.22d-16

*     Set the desired solution xtrue.
*     For least-squares problems, this is it.
*     For underdetermined systems, lstp may alter it.

      do 100 j = 1, n
*        xtrue(j) = one
         xtrue(j) = j * 0.1d+0
  100 continue

*     Generate the specified test problem.
*     The workspace array  iw  is not needed in this application.
*     The workspace array  rw  is used for the following vectors:
*        d(minmn), hy(m), hz(n), w(maxmn).
*     The vectors  d, hy, hz  will define the test matrix A.
*     w is needed for workspace in aprod1 and aprod2.

      maxmn  = max( m, n )
      minmn  = min( m, n )
      locd   = 1
      lochy  = locd  + minmn
      lochz  = lochy + m
      locw   = lochz + n
      ltotal = locw  + maxmn - 1
      if (ltotal .gt. lenrw) go to 900

      call lstp  ( m, n, maxmn, minmn, nduplc, npower, damp, xtrue,
     $             b, rw(locd), rw(lochy), rw(lochz), rw(locw),
     $             acond, rnorm )

      write(nout, 1000) line, line,
     $                  m, n, nduplc, npower, damp, acond, rnorm,
     $                  line, line

*     Check that aprod generates y + Ax and x + A'y consistently.

      call acheck( m, n, nout, aprod, eps,
     $             leniw, lenrw, iw, rw,
     $             v, w, x, y,
     $             inform )

      if (inform .gt. 0) then
         write(nout, '(a)') ' Check eps and power in subroutine acheck'
         stop
      end if

*     Solve the problem defined by aprod, damp and b.
*     Copy the rhs vector b into u  (LSQR will overwrite u)
*     and set the other input parameters for LSQR.
*     We ask for standard errors only if they are well-defined.

      call dcopy ( m, b, 1, u, 1 )
*---  wantse = m .gt. n  .or.  damp .gt. zero
      wantse = .false.
      atol   = eps**0.99
      btol   = atol
      conlim = 1000.0 * acond
      itnlim = 4*(m + n + 50)

      call LSQR  ( m, n, aprod, damp, wantse,
     $             leniw, lenrw, iw, rw,
     $             u, v, w, x, se,
     $             atol, btol, conlim, itnlim, nout,
     $             istop, itn, anorm, acond, rnorm, arnorm, xnorm )

*     Examine the results.

*-      if (damp .eq. zero) then
*-         write(nout, 2000)       xnorm, rnorm, arnorm
*-      else
*-         write(nout, 2100) damp, xnorm, rnorm, arnorm
*-      end if

      call xcheck( m, n, nout, aprod, anorm, damp, eps,
     $             leniw, lenrw, iw, rw,
     $             b, u, v, w, x,
     $             inform, test1, test2, test3 )

*     Print the solution and standard error estimates from  LSQR.

      nprint = min( m, n, 8 )
      write(nout, 2500)    (j, x(j) , j = 1, nprint)
      if ( wantse ) then
         write(nout, 2600) (j, se(j), j = 1, nprint)
      end if

*     Print a clue about whether the solution looks OK.

      do 500 j = 1, n
         w(j)  = x(j) - xtrue(j)
  500 continue
      wnorm    = dnrm2 ( n, w    , 1 )
      xnorm    = dnrm2 ( n, xtrue, 1 )
      enorm    = wnorm / (one + xnorm)
      etol     = 0.001
      if (enorm .le. etol) then
         write(nout, 3000) enorm
      else
         write(nout, 3100) enorm
      end if
      return

*     m or n too large.

  800 write(nout, 8000)
      return

*     Not enough workspace.

  900 write(nout, 9000) ltotal
      return

 1000 format(1p
     $ // 1x, 2a
     $ /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )'
     $ // ' Condition no. =', e12.4,  '     Residual function =', e17.9
     $ /  1x, 2a)
 2000 format(1p
     $ // ' We are solving    min norm(Ax - b)    with no damping.'
     $ // ' Estimates from LSQR:'
     $ /  '    norm(x)         =', e10.3, ' = xnorm'
     $ /  '    norm(r)         =', e10.3, ' = rnorm'
     $ /  '    norm(A''r)       =', e10.3, ' = arnorm')
 2100 format(1p
     $ // ' We are solving    min norm(Ax - b)    with damp =', e10.3
     $ /  '                           (damp*x)'
     $ // ' Estimates from LSQR:'
     $ /  '    norm(x)         =', e10.3, ' = xnorm'
     $ /  '    norm(rbar)      =', e10.3, ' = rnorm'
     $ /  '    norm(Abar''rbar) =', e10.3, ' = arnorm')
 2500 format(//' Solution  x:' / 4(i6, g14.6))
 2600 format(/ ' Standard errors  se:' / 4(i6, g14.6))
 3000 format(1p / ' LSQR  appears to be successful.',
     $        '     Relative error in  x  =', e10.2)
 3100 format(1p / ' LSQR  appears to have failed.  ',
     $        '     Relative error in  x  =', e10.2)
 8000 format(/ ' XXX  m or n is too large.')
 9000 format(/ ' XXX  Insufficient workspace.',
     $        '  The length of  rw  should be at least', i6)

*     end of test
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*     MAIN PROGRAM

      program            main
      implicit           double precision (a-h,o-z)

      open( 6, file='LSQR.LIS', status='NEW' )

      nbar   = 1000
      nduplc = 40

      m      = 2*nbar
      n      = nbar
      do  ndamp = 2, 7
         npower = ndamp
         damp   = 10.0d+0**(-ndamp-6)
         call test  ( m, n, nduplc, npower, damp )
      end do

      m      = nbar
      n      = nbar
      do  ndamp = 2, 7
         npower = ndamp
         damp   = 10.0d+0**(-ndamp-6)
         call test  ( m, n, nduplc, npower, damp )
      end do

      m      = nbar
      n      = 2*nbar
      do  ndamp = 2, 7
         npower = ndamp
         damp   = 10.0d+0**(-ndamp-6)
         call test  ( m, n, nduplc, npower, damp )
      end do

*     end of main program for testing LSQR
      end
