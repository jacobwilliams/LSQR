*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File lsqr   fortran   (double precision)
*
*     LSQR     d2norm
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine LSQR  ( m, n, aprod, damp, wantse,
     $                   leniw, lenrw, iw, rw,
     $                   u, v, w, x, se,
     $                   atol, btol, conlim, itnlim, nout,
     $                   istop, itn, anorm, acond, rnorm, arnorm, xnorm)

      implicit           double precision (a-h,o-z)
      external           aprod
      logical            wantse
      integer            m, n, leniw, lenrw, itnlim, nout, istop, itn
      integer            iw(leniw)
      double precision   rw(lenrw), u(m), v(n), w(n), x(n), se(*),
     $                   atol, btol, conlim, damp,
     $                   anorm, acond, rnorm, arnorm, xnorm
*     ------------------------------------------------------------------
*
*     LSQR  finds a solution x to the following problems:
*
*     1. Unsymmetric equations --    solve  A*x = b
*
*     2. Linear least squares  --    solve  A*x = b
*                                    in the least-squares sense
*
*     3. Damped least squares  --    solve  (   A    )*x = ( b )
*                                           ( damp*I )     ( 0 )
*                                    in the least-squares sense
*
*     where A is a matrix with m rows and n columns, b is an
*     m-vector, and damp is a scalar.  (All quantities are real.)
*     The matrix A is intended to be large and sparse.  It is accessed
*     by means of subroutine calls of the form
*
*                call aprod ( mode, m, n, x, y, leniw, lenrw, iw, rw )
*
*     which must perform the following functions:
*
*                If mode = 1, compute  y = y + A*x.
*                If mode = 2, compute  x = x + A(transpose)*y.
*
*     The vectors x and y are input parameters in both cases.
*     If  mode = 1,  y should be altered without changing x.
*     If  mode = 2,  x should be altered without changing y.
*     The parameters leniw, lenrw, iw, rw may be used for workspace
*     as described below.
*
*     The rhs vector b is input via u, and subsequently overwritten.
*
*
*     Note:  LSQR uses an iterative method to approximate the solution.
*     The number of iterations required to reach a certain accuracy
*     depends strongly on the scaling of the problem.  Poor scaling of
*     the rows or columns of A should therefore be avoided where
*     possible.
*
*     For example, in problem 1 the solution is unaltered by
*     row-scaling.  If a row of A is very small or large compared to
*     the other rows of A, the corresponding row of ( A  b ) should be
*     scaled up or down.
*
*     In problems 1 and 2, the solution x is easily recovered
*     following column-scaling.  Unless better information is known,
*     the nonzero columns of A should be scaled so that they all have
*     the same Euclidean norm (e.g., 1.0).
*
*     In problem 3, there is no freedom to re-scale if damp is
*     nonzero.  However, the value of damp should be assigned only
*     after attention has been paid to the scaling of A.
*
*     The parameter damp is intended to help regularize
*     ill-conditioned systems, by preventing the true solution from
*     being very large.  Another aid to regularization is provided by
*     the parameter acond, which may be used to terminate iterations
*     before the computed solution becomes very large.
*
*     Note that x is not an input parameter.
*     If some initial estimate x0 is known and if damp = 0,
*     one could proceed as follows:
*
*       1. Compute a residual vector     r0 = b - A*x0.
*       2. Use LSQR to solve the system  A*dx = r0.
*       3. Add the correction dx to obtain a final solution x = x0 + dx.
*
*     This requires that x0 be available before and after the call
*     to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
*     to solve A*x = b and k2 iterations to solve A*dx = r0.
*     If x0 is "good", norm(r0) will be smaller than norm(b).
*     If the same stopping tolerances atol and btol are used for each
*     system, k1 and k2 will be similar, but the final solution x0 + dx
*     should be more accurate.  The only way to reduce the total work
*     is to use a larger stopping tolerance for the second system.
*     If some value btol is suitable for A*x = b, the larger value
*     btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
*
*     Preconditioning is another way to reduce the number of iterations.
*     If it is possible to solve a related system M*x = b efficiently,
*     where M approximates A in some helpful way
*     (e.g. M - A has low rank or its elements are small relative to
*     those of A), LSQR may converge more rapidly on the system
*           A*M(inverse)*z = b,
*     after which x can be recovered by solving M*x = z.
*
*     NOTE: If A is symmetric, LSQR should not be used!
*     Alternatives are the symmetric conjugate-gradient method (cg)
*     and/or SYMMLQ.
*     SYMMLQ is an implementation of symmetric cg that applies to
*     any symmetric A and will converge more rapidly than LSQR.
*     If A is positive definite, there are other implementations of
*     symmetric cg that require slightly less work per iteration
*     than SYMMLQ (but will take the same number of iterations).
*
*
*     Notation
*     --------
*
*     The following quantities are used in discussing the subroutine
*     parameters:
*
*     Abar   =  (   A    ),          bbar  =  ( b )
*               ( damp*I )                    ( 0 )
*
*     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
*
*     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
*            =  norm( rbar )
*
*     relpr  =  the relative precision of floating-point arithmetic
*               on the machine being used.  On most machines,
*               relpr is about 1.0e-7 and 1.0d-16 in single and double
*               precision respectively.
*
*     LSQR  minimizes the function rnorm with respect to x.
*
*
*     Parameters
*     ----------
*
*     m       input      m, the number of rows in A.
*
*     n       input      n, the number of columns in A.
*
*     aprod   external   See above.
*
*     damp    input      The damping parameter for problem 3 above.
*                        (damp should be 0.0 for problems 1 and 2.)
*                        If the system A*x = b is incompatible, values
*                        of damp in the range 0 to sqrt(relpr)*norm(A)
*                        will probably have a negligible effect.
*                        Larger values of damp will tend to decrease
*                        the norm of x and reduce the number of 
*                        iterations required by LSQR.
*
*                        The work per iteration and the storage needed
*                        by LSQR are the same for all values of damp.
*
*     wantse  input      A logical variable to say if the array se(*)
*                        of standard error estimates should be computed.
*                        If m .gt. n  or  damp .gt. 0,  the system is
*                        overdetermined and the standard errors may be
*                        useful.  (See the first LSQR reference.)
*                        Otherwise (m .le. n  and  damp = 0) they do not
*                        mean much.  Some time and storage can be saved
*                        by setting wantse = .false. and using any
*                        convenient array for se(*), which won't be
*                        touched.
*
*     leniw   input      The length of the workspace array iw.
*     lenrw   input      The length of the workspace array rw.
*     iw      workspace  An integer array of length leniw.
*     rw      workspace  A real array of length lenrw.
*
*             Note:  LSQR  does not explicitly use the previous four
*             parameters, but passes them to subroutine aprod for
*             possible use as workspace.  If aprod does not need
*             iw or rw, the values leniw = 1 or lenrw = 1 should
*             be used, and the actual parameters corresponding to
*             iw or rw  may be any convenient array of suitable type.
*
*     u(m)    input      The rhs vector b.  Beware that u is
*                        over-written by LSQR.
*
*     v(n)    workspace
*
*     w(n)    workspace
*
*     x(n)    output     Returns the computed solution x.
*
*     se(*)   output     If wantse is true, the dimension of se must be
*             (maybe)    n or more.  se(*) then returns standard error
*                        estimates for the components of x.
*                        For each i, se(i) is set to the value
*                           rnorm * sqrt( sigma(i,i) / t ),
*                        where sigma(i,i) is an estimate of the i-th
*                        diagonal of the inverse of Abar(transpose)*Abar
*                        and  t = 1      if  m .le. n,
*                             t = m - n  if  m .gt. n  and  damp = 0,
*                             t = m      if  damp .ne. 0.
*
*                        If wantse is false, se(*) will not be touched.
*                        The actual parameter can be any suitable array
*                        of any length.
*
*     atol    input      An estimate of the relative error in the data
*                        defining the matrix A.  For example,
*                        if A is accurate to about 6 digits, set
*                        atol = 1.0e-6 .
*
*     btol    input      An estimate of the relative error in the data
*                        defining the rhs vector b.  For example,
*                        if b is accurate to about 6 digits, set
*                        btol = 1.0e-6 .
*
*     conlim  input      An upper limit on cond(Abar), the apparent
*                        condition number of the matrix Abar.
*                        Iterations will be terminated if a computed
*                        estimate of cond(Abar) exceeds conlim.
*                        This is intended to prevent certain small or
*                        zero singular values of A or Abar from
*                        coming into effect and causing unwanted growth
*                        in the computed solution.
*
*                        conlim and damp may be used separately or
*                        together to regularize ill-conditioned systems.
*
*                        Normally, conlim should be in the range
*                        1000 to 1/relpr.
*                        Suggested value:
*                        conlim = 1/(100*relpr)  for compatible systems,
*                        conlim = 1/(10*sqrt(relpr)) for least squares.
*
*             Note:  If the user is not concerned about the parameters
*             atol, btol and conlim, any or all of them may be set
*             to zero.  The effect will be the same as the values
*             relpr, relpr and 1/relpr respectively.
*
*     itnlim  input      An upper limit on the number of iterations.
*                        Suggested value:
*                        itnlim = n/2   for well-conditioned systems
*                                       with clustered singular values,
*                        itnlim = 4*n   otherwise.
*
*     nout    input      File number for printed output.  If positive,
*                        a summary will be printed on file nout.
*
*     istop   output     An integer giving the reason for termination:
*
*                0       x = 0  is the exact solution.
*                        No iterations were performed.
*
*                1       The equations A*x = b are probably
*                        compatible.  Norm(A*x - b) is sufficiently
*                        small, given the values of atol and btol.
*
*                2       damp is zero.  The system A*x = b is probably
*                        not compatible.  A least-squares solution has
*                        been obtained that is sufficiently accurate,
*                        given the value of atol.
*
*                3       damp is nonzero.  A damped least-squares
*                        solution has been obtained that is sufficiently
*                        accurate, given the value of atol.
*
*                4       An estimate of cond(Abar) has exceeded
*                        conlim.  The system A*x = b appears to be
*                        ill-conditioned.  Otherwise, there could be an
*                        error in subroutine aprod.
*
*                5       The iteration limit itnlim was reached.
*
*     itn     output     The number of iterations performed.
*
*     anorm   output     An estimate of the Frobenius norm of  Abar.
*                        This is the square-root of the sum of squares
*                        of the elements of Abar.
*                        If damp is small and if the columns of A
*                        have all been scaled to have length 1.0,
*                        anorm should increase to roughly sqrt(n).
*                        A radically different value for anorm may
*                        indicate an error in subroutine aprod (there
*                        may be an inconsistency between modes 1 and 2).
*
*     acond   output     An estimate of cond(Abar), the condition
*                        number of Abar.  A very high value of acond
*                        may again indicate an error in aprod.
*
*     rnorm   output     An estimate of the final value of norm(rbar),
*                        the function being minimized (see notation
*                        above).  This will be small if A*x = b has
*                        a solution.
*
*     arnorm  output     An estimate of the final value of
*                        norm( Abar(transpose)*rbar ), the norm of
*                        the residual for the usual normal equations.
*                        This should be small in all cases.  (arnorm
*                        will often be smaller than the true value
*                        computed from the output vector x.)
*
*     xnorm   output     An estimate of the norm of the final
*                        solution vector x.
*
*
*     Subroutines and functions used              
*     ------------------------------
*
*     USER               aprod
*     LSQR               d2norm
*     BLAS               dcopy, dnrm2, dscal (see Lawson et al. below)
*
*
*     Precision
*     ---------
*
*     The number of iterations required by LSQR will usually decrease
*     if the computation is performed in higher precision.
*     At least 15-digit arithmetic should normally be used.
*     To convert LSQR and D2NORM between single and double precision,
*     change
*                        double precision
*                        dcopy, dnrm2, dscal
*     to the appropriate FORTRAN and BLAS equivalents.
*     Also change 'd+' or 'e+' in the parameter statement.
*
*
*     References
*     ----------
*
*     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
*          linear equations and sparse least squares,
*          ACM Transactions on Mathematical Software 8, 1 (March 1982),
*          pp. 43-71.
*
*     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
*          linear equations and least-squares problems,
*          ACM Transactions on Mathematical Software 8, 2 (June 1982),
*          pp. 195-209.
*
*     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
*          Basic linear algebra subprograms for Fortran usage,
*          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
*          pp. 308-323 and 324-325.
*     ------------------------------------------------------------------
*
*
*     LSQR development:
*     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
*     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
*     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
*                     if ( (one + dabs(t)) .le. one ) GO TO 200
*                  from loop 200.  The test was an attempt to reduce
*                  underflows, but caused w(i) not to be updated.
*     17 Mar 1989: First F77 version.
*     04 May 1989: Bug (David Gay, AT&T).  When the second beta is zero,
*                  rnorm = 0 and
*                  test2 = arnorm / (anorm * rnorm) overflows.
*                  Fixed by testing for rnorm = 0.
*     05 May 1989: Sent to "misc" in netlib.
*     14 Mar 1990: Bug (John Tomlin via IBM OSL testing).
*                  Setting rhbar2 = rhobar**2 + dampsq can give zero
*                  if rhobar underflows and damp = 0.
*                  Fixed by testing for damp = 0 specially.
*     15 Mar 1990: Converted to lower case.
*     21 Mar 1990: d2norm introduced to avoid overflow in numerous
*                  items like  c = sqrt( a**2 + b**2 ).
*     04 Sep 1991: wantse added as an argument to LSQR, to make
*                  standard errors optional.  This saves storage and
*                  time when se(*) is not wanted.
*     13 Feb 1992: istop now returns a value in [1,5], not [1,7].
*                  1, 2 or 3 means that x solves one of the problems
*                  Ax = b,  min norm(Ax - b)  or  damped least squares.
*                  4 means the limit on cond(A) was reached.
*                  5 means the limit on iterations was reached.
*     07 Dec 1994: Keep track of dxmax = max_k norm( phi_k * d_k ).
*                  So far, this is just printed at the end.
*                  A large value (relative to norm(x)) indicates
*                  significant cancellation in forming
*                  x  =  D*f  =  sum( phi_k * d_k ).
*                  A large column of D need NOT be serious if the
*                  corresponding phi_k is small.
*     27 Dec 1994: Include estimate of alfa_opt in iteration log.
*                  alfa_opt is the optimal scale factor for the
*                  residual in the "augmented system", as described by
*                  A. Bjorck (1992),
*                  Pivoting and stability in the augmented system method,
*                  in D. F. Griffiths and G. A. Watson (eds.),
*                  "Numerical Analysis 1991",
*                  Proceedings of the 14th Dundee Conference,
*                  Pitman Research Notes in Mathematics 260,
*                  Longman Scientific and Technical, Harlow, Essex, 1992.
*                  
*
*     Michael A. Saunders                  mike@sol-michael.stanford.edu
*     Dept of Operations Research          na.Msaunders@na-net.ornl.gov
*     Stanford University
*     Stanford, CA 94305-4022              (415) 723-1875
*-----------------------------------------------------------------------

      external           d2norm, dnrm2, dcopy, dscal
      double precision   d2norm, dnrm2

*     Local variables

      logical            damped, extra
      integer            i, maxdx, nchar, nconv, nstop
      double precision   alpha, beta, bnorm,
     $                   cs, cs1, cs2, ctol,
     $                   delta, dknorm, dnorm, dxk, dxmax,
     $                   gamma, gambar, phi, phibar, psi,
     $                   res2, rho, rhobar, rhbar1,
     $                   rhs, rtol, sn, sn1, sn2,
     $                   t, tau, temp, test1, test2, test3,
     $                   theta, t1, t2, t3, xnorm1, z, zbar

      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )

      character*14       enter, exit              
      character*53       msg(0:5)

      data               enter /' Enter LSQR.  '/
      data               exit  /' Exit  LSQR.  '/
      data               msg
     $ / 'The exact solution is  x = 0',
     $   'A solution to Ax = b was found, given atol, btol',
     $   'A least-squares solution was found, given atol',
     $   'A damped least-squares solution was found, given atol',
     $   'Cond(Abar) seems to be too large, given conlim',
     $   'The iteration limit was reached' /
*-----------------------------------------------------------------------


*     Initialize.

      if (nout .gt. 0) then
         write(nout, 1000) enter, m, n, damp, wantse,
     $                     atol, conlim, btol, itnlim
      end if

      damped =   damp .gt. zero
      extra  =  .true.     ! true for extra printing below.
      itn    =   0
      istop  =   0
      nstop  =   0
      maxdx  =   0
      ctol   =   zero
      if (conlim .gt. zero) ctol = one / conlim
      anorm  =   zero
      acond  =   zero
      dnorm  =   zero
      dxmax  =   zero
      res2   =   zero
      psi    =   zero
      xnorm  =   zero
      xnorm1 =   zero
      cs2    = - one
      sn2    =   zero
      z      =   zero

*     ------------------------------------------------------------------
*     Set up the first vectors u and v for the bidiagonalization.
*     These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
*     ------------------------------------------------------------------
      do 10  i = 1, n
         v(i)  =  zero
         x(i)  =  zero
   10 continue

      if ( wantse ) then
         do 20  i = 1, n
            se(i) =  zero
   20    continue
      end if

      alpha  =   zero
      beta   =   dnrm2 ( m, u, 1 )

      if (beta .gt. zero) then
         call dscal ( m, (one / beta), u, 1 )
         call aprod ( 2, m, n, v, u, leniw, lenrw, iw, rw )
         alpha  =   dnrm2 ( n, v, 1 )
      end if

      if (alpha .gt. zero) then
         call dscal ( n, (one / alpha), v, 1 )
         call dcopy ( n, v, 1, w, 1 )
      end if

      arnorm =   alpha * beta
      if (arnorm .eq. zero) go to 800

      rhobar =   alpha
      phibar =   beta
      bnorm  =   beta
      rnorm  =   beta

      if (nout   .gt.  0  ) then
         if ( damped ) then
             write(nout, 1300)
         else
             write(nout, 1200)
         end if
         test1  = one
         test2  = alpha / beta

         if ( extra ) then
            write(nout, 1400)
         end if
         write(nout, 1500) itn, x(1), rnorm, test1, test2
         write(nout, 1600)
      end if


*     ==================================================================
*     Main iteration loop.
*     ==================================================================
  100 itn    = itn + 1

*     ------------------------------------------------------------------
*     Perform the next step of the bidiagonalization to obtain the
*     next  beta, u, alpha, v.  These satisfy the relations
*                beta*u  =  A*v  -  alpha*u,
*               alpha*v  =  A(transpose)*u  -  beta*v.
*     ------------------------------------------------------------------
      call dscal ( m, (- alpha), u, 1 )
      call aprod ( 1, m, n, v, u, leniw, lenrw, iw, rw )
      beta   =   dnrm2 ( m, u, 1 )

*     Accumulate  anorm = || Bk ||
*                       =  sqrt( sum of  alpha**2 + beta**2 + damp**2 ).

      temp   =   d2norm( alpha, beta )
      temp   =   d2norm( temp , damp )
      anorm  =   d2norm( anorm, temp )

      if (beta .gt. zero) then
         call dscal ( m, (one / beta), u, 1 )
         call dscal ( n, (- beta), v, 1 )
         call aprod ( 2, m, n, v, u, leniw, lenrw, iw, rw )
         alpha  =   dnrm2 ( n, v, 1 )
         if (alpha .gt. zero) then
            call dscal ( n, (one / alpha), v, 1 )
         end if
      end if

*     ------------------------------------------------------------------
*     Use a plane rotation to eliminate the damping parameter.
*     This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
*     ------------------------------------------------------------------
      rhbar1 = rhobar
      if ( damped ) then
         rhbar1 = d2norm( rhobar, damp )
         cs1    = rhobar / rhbar1
         sn1    = damp   / rhbar1
         psi    = sn1 * phibar
         phibar = cs1 * phibar
      end if                   

*     ------------------------------------------------------------------
*     Use a plane rotation to eliminate the subdiagonal element (beta)
*     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
*     ------------------------------------------------------------------
      rho    =   d2norm( rhbar1, beta )
      cs     =   rhbar1 / rho
      sn     =   beta   / rho
      theta  =   sn * alpha
      rhobar = - cs * alpha
      phi    =   cs * phibar
      phibar =   sn * phibar
      tau    =   sn * phi

*     ------------------------------------------------------------------
*     Update  x, w  and (perhaps) the standard error estimates.
*     ------------------------------------------------------------------
      t1     =   phi   / rho
      t2     = - theta / rho
      t3     =   one   / rho
      dknorm =   zero

      if ( wantse ) then
         do 200  i =  1, n
            t      =  w(i)
            x(i)   =  t1*t  +  x(i)
            w(i)   =  t2*t  +  v(i)
            t      = (t3*t)**2
            se(i)  =  t     +  se(i)
            dknorm =  t     +  dknorm
  200    continue
      else
         do 220  i =  1, n
            t      =  w(i)
            x(i)   =  t1*t  +  x(i)
            w(i)   =  t2*t  +  v(i)
            dknorm = (t3*t)**2  +  dknorm
  220    continue
      end if

*     ------------------------------------------------------------------
*     Monitor the norm of d_k, the update to x.
*     dknorm = norm( d_k )
*     dnorm  = norm( D_k ),        where   D_k = (d_1, d_2, ..., d_k )
*     dxk    = norm( phi_k d_k ),  where new x = x_k + phi_k d_k.
*     ------------------------------------------------------------------
      dknorm = sqrt( dknorm )
      dnorm  = d2norm( dnorm, dknorm )
      dxk    = abs( phi * dknorm )
      if (dxmax .lt. dxk ) then
          dxmax   =  dxk
          maxdx   =  itn
      end if

*     ------------------------------------------------------------------
*     Use a plane rotation on the right to eliminate the
*     super-diagonal element (theta) of the upper-bidiagonal matrix.
*     Then use the result to estimate  norm(x).
*     ------------------------------------------------------------------
      delta  =   sn2 * rho
      gambar = - cs2 * rho
      rhs    =   phi    - delta * z
      zbar   =   rhs    / gambar
      xnorm  =   d2norm( xnorm1, zbar  )
      gamma  =   d2norm( gambar, theta )
      cs2    =   gambar / gamma
      sn2    =   theta  / gamma
      z      =   rhs    / gamma
      xnorm1 =   d2norm( xnorm1, z     )

*     ------------------------------------------------------------------
*     Test for convergence.
*     First, estimate the norm and condition of the matrix  Abar,
*     and the norms of  rbar  and  Abar(transpose)*rbar.
*     ------------------------------------------------------------------
      acond  =   anorm * dnorm
      res2   =   d2norm( res2 , psi    )
      rnorm  =   d2norm( res2 , phibar )
      arnorm =   alpha * abs( tau )

*     Now use these norms to estimate certain other quantities,
*     some of which will be small near a solution.

      alfopt =   sqrt( rnorm / (dnorm * xnorm) )
      test1  =   rnorm /  bnorm
      test2  =   zero
      if (rnorm .gt. zero) test2 = arnorm / (anorm * rnorm)
      test3  =   one   /  acond
      t1     =   test1 / (one  +  anorm * xnorm / bnorm)
      rtol   =   btol  +  atol *  anorm * xnorm / bnorm

*     The following tests guard against extremely small values of
*     atol, btol  or  ctol.  (The user may have set any or all of
*     the parameters  atol, btol, conlim  to zero.)
*     The effect is equivalent to the normal tests using
*     atol = relpr,  btol = relpr,  conlim = 1/relpr.

      t3     =   one + test3
      t2     =   one + test2
      t1     =   one + t1
      if (itn .ge. itnlim) istop = 5
      if (t3  .le. one   ) istop = 4
      if (t2  .le. one   ) istop = 2
      if (t1  .le. one   ) istop = 1

*     Allow for tolerances set by the user.

      if (test3 .le. ctol) istop = 4
      if (test2 .le. atol) istop = 2
      if (test1 .le. rtol) istop = 1

*     ------------------------------------------------------------------
*     See if it is time to print something.
*     ------------------------------------------------------------------
      if (nout  .le.  0       ) go to 600
      if (n     .le. 40       ) go to 400
      if (itn   .le. 10       ) go to 400
      if (itn   .ge. itnlim-10) go to 400
      if (mod(itn,10) .eq. 0  ) go to 400
      if (test3 .le.  2.0*ctol) go to 400
      if (test2 .le. 10.0*atol) go to 400
      if (test1 .le. 10.0*rtol) go to 400
      if (istop .ne.  0       ) go to 400
      go to 600

*     Print a line for this iteration.
*     "extra" is for experimental purposes.

  400 if ( extra ) then
         write(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
     $                   , phi, dknorm, dxk, alfopt
      else
         write(nout, 1500) itn, x(1), rnorm, test1, test2, anorm, acond
      end if
      if (mod(itn,10) .eq. 0) write(nout, 1600)

*     ------------------------------------------------------------------
*     Stop if appropriate.
*     The convergence criteria are required to be met on  nconv
*     consecutive iterations, where  nconv  is set below.
*     Suggested value:  nconv = 1, 2  or  3.
*     ------------------------------------------------------------------
  600 if (istop .eq. 0) then
         nstop  = 0
      else
         nconv  = 1
         nstop  = nstop + 1
         if (nstop .lt. nconv  .and.  itn .lt. itnlim) istop = 0
      end if
      if (istop .eq. 0) go to 100

*     ==================================================================
*     End of iteration loop.
*     ==================================================================

*     Finish off the standard error estimates.

      if ( wantse ) then
         t    =   one
         if (m .gt. n)  t = m - n
         if ( damped )  t = m
         t    =   rnorm / sqrt( t )
      
         do 700  i = 1, n
            se(i)  = t * sqrt( se(i) )
  700    continue
      end if

*     Decide if istop = 2 or 3.
*     Print the stopping condition.

  800 if (damped  .and.  istop .eq. 2) istop = 3
      if (nout .gt. 0) then
         write(nout, 2000) exit, istop, itn,
     $                     exit, anorm, acond,
     $                     exit, bnorm, xnorm,
     $                     exit, rnorm, arnorm
         write(nout, 2100) exit, dxmax, maxdx,
     $                     exit, dxmax/(xnorm + 1.0d-20)
         write(nout, 3000) exit, msg(istop)
      end if

  900 return

*     ------------------------------------------------------------------
 1000 format(// 1p, a, '     Least-squares solution of  Ax = b'
     $    / ' The matrix  A  has', i7, ' rows   and', i7, ' columns'
     $    / ' damp   =', e22.14, 3x,        'wantse =', l10
     $    / ' atol   =', e10.2, 15x,        'conlim =', e10.2
     $    / ' btol   =', e10.2, 15x,        'itnlim =', i10)
 1200 format(// '   Itn       x(1)           Function',
     $   '     Compatible   LS        Norm A    Cond A')
 1300 format(// '   Itn       x(1)           Function',
     $   '     Compatible   LS     Norm Abar Cond Abar')
 1400 format(80x, '    phi    dknorm   dxk  alfa_opt')
 1500 format(1p, i6, 2e17.9, 4e10.2, e9.1, 3e8.1)
 1600 format(1x)
 2000 format(/ 1p, a, 5x, 'istop  =', i2,   15x, 'itn    =', i8
     $       /     a, 5x, 'anorm  =', e12.5, 5x, 'acond  =', e12.5
     $       /     a, 5x, 'bnorm  =', e12.5, 5x, 'xnorm  =', e12.5
     $       /     a, 5x, 'rnorm  =', e12.5, 5x, 'arnorm =', e12.5)
 2100 format(  1p, a, 5x, 'max dx =', e8.1 , ' occurred at itn ', i8,
     $       /     a, 5x, '       =', e8.1 , '*xnorm' )
 3000 format( a, 5x, a )

*     End of LSQR
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      function          d2norm( a, b )

      double precision  d2norm, a, b

*     ------------------------------------------------------------------
*     d2norm  returns  sqrt( a**2 + b**2 )  with precautions
*     to avoid overflow.
*
*     21 Mar 1990: First version.
*     ------------------------------------------------------------------

      intrinsic         abs, sqrt
      double precision  scale
      double precision  zero
      parameter       ( zero = 0.0d+0 )

      scale  = abs( a ) + abs( b )
      if (scale .eq. zero) then
         d2norm = zero
      else
         d2norm = scale * sqrt( (a/scale)**2   +  (b/scale)**2 )
      end if

*     end of d2norm
      end
