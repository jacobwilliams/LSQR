!***************************************************************************************************
!>
!  Test module for [[lsqr]].
!
!  These routines define a class of least-squares test problems
!  for testing algorithms LSQR and CRAIG
!  (Paige and Saunders, ACM TOMS, 1982).
!
!### Author
! *  Michael Saunders, Dept of Operations Research, Stanford University.
!
!### History
!  * 1982---1991:  Various versions implemented.
!  * 06 Feb 1992:  Test-problem generator lstp generalized to allow
!    any m and n.  lstp is now the same as the generator
!    for LSQR and CRAIG.
!  * 30 Nov 1993:  Modified lstp.
!    For a while, damp = 0 implied r = damp*s = 0.
!    This was a result of generating x and s.
!    Reverted to generating x and r as in LSQR paper.
!  * 12 Nov 2019 : Jacob Williams : significant refactoring into modern Fortran.

    module lsqrtest_module

    use lsqr_kinds
    use lsqpblas_module, only: dnrm2,dcopy,dscal
    use lsqr_module,     only: lsqr_solver

    implicit none

    private

    integer,parameter :: lenrw = 10000

    type,extends(lsqr_solver) :: test_solver
        integer :: nout = -1  !! output unit for printing
        real(wp),dimension(lenrw) :: rw  !! workspace array
    contains
        procedure :: aprod => aprod_test_solver
        procedure :: test
        procedure :: aprod1
        procedure :: aprod2
        procedure :: lstp
    end type

    public :: lsqr_test

    contains
!***************************************************************************************************

!***************************************************************************************************
!>
!  Unit test.

    subroutine lsqr_test()

    integer :: iunit,nbar,nduplc,n,m,ndamp,npower
    real(wp) :: damp
    type(test_solver) :: solver

    open( newunit=iunit, file='LSQR.LIS', status='REPLACE' )

    solver%nout = iunit

    nbar   = 1000
    nduplc = 40

    m = 2*nbar
    n = nbar
    do ndamp = 2, 7
        npower = ndamp
        damp   = 10.0_wp**(-ndamp-6)
        call solver%test( m, n, nduplc, npower, damp )
    end do

    m = nbar
    n = nbar
    do ndamp = 2, 7
        npower = ndamp
        damp   = 10.0_wp**(-ndamp-6)
        call solver%test( m, n, nduplc, npower, damp )
    end do

    m = nbar
    n = 2*nbar
    do ndamp = 2, 7
        npower = ndamp
        damp   = 10.0_wp**(-ndamp-6)
        call solver%test( m, n, nduplc, npower, damp )
    end do

    close(iunit)

    end subroutine lsqr_test
!***************************************************************************************************

!***************************************************************************************************
!>
!  This is an example driver routine for running LSQR.
!  It generates a test problem, solves it, and examines the results.
!  Note that subroutine aprod must be declared external
!  if it is used only in the call to LSQR (and acheck).
!
!### History
!  * 1982---1991:  Various versions implemented.
!  * 04 Sep 1991:  "wantse" added to argument list of LSQR,
!    making standard errors optional.
!  * 10 Feb 1992:  Revised for use with lsqrchk fortran.
!  * 31 Mar 2005: changed atol = eps**0.666667 to eps*0.99
!    to increase accuracy of the solution.  LSQR appears to be
!    successful on all 18 test problems except 5 and 6
!    (which are over-determined and too ill-conditioned to
!    permit any correct digits).
!    The output from an Intel Xeon system with g77 is in LSQR.LIS.
!    The two "appears to have failed" messages are no cause for alarm.

!  Michael Saunders, Dept of Operations Research, Stanford University.

    subroutine test  ( me, m, n, nduplc, npower, damp )

    class(test_solver),intent(inout) :: me
    integer :: m, n, nduplc, npower
    real(wp) :: damp

    integer,parameter :: maxm = 2000
    integer,parameter :: maxn = 2000
    integer,parameter :: mxmn = 2000
    real(wp),parameter :: eps = epsilon(1.0_wp) !! machine precision
    character(len=34),parameter :: line = '----------------------------------'

    integer :: inform, istop, itnlim, j, itn, maxmn, minmn, nprint
    integer :: locd, lochy, lochz, locw, ltotal
    logical :: wantse
    real(wp) :: b(maxm),  u(maxm), &
                v(maxn),  w(mxmn), x(maxn), &
                se(maxn), xtrue(maxn), y(mxmn)
    real(wp) :: atol, btol, conlim, &
                anorm, acond, rnorm, arnorm, &
                enorm, etol, xnorm, test1, test2, test3, wnorm

    if (m > maxm  .or.  n > maxn) then
        ! m or n too large.
        write(me%nout, 8000)
        return
    end if

    ! Set the desired solution xtrue.
    ! For least-squares problems, this is it.
    ! For underdetermined systems, lstp may alter it.

    do j = 1, n
        ! xtrue(j) = one
        xtrue(j) = j * 0.1_wp
    end do

    ! Generate the specified test problem.
    ! The workspace array  rw  is used for the following vectors:
    !    d(minmn), hy(m), hz(n), w(maxmn).
    ! The vectors  d, hy, hz  will define the test matrix A.
    ! w is needed for workspace in aprod1 and aprod2.

    maxmn  = max( m, n )
    minmn  = min( m, n )
    locd   = 1
    lochy  = locd  + minmn
    lochz  = lochy + m
    locw   = lochz + n
    ltotal = locw  + maxmn - 1
    if (ltotal > lenrw) then
        ! Not enough workspace.
        write(me%nout, 9000) ltotal
    end if

    call me%lstp  ( m, n, maxmn, minmn, nduplc, npower, damp, xtrue, &
                b, me%rw(locd), me%rw(lochy), me%rw(lochz), me%rw(locw), &
                acond, rnorm )

    write(me%nout, 1000) line, line, &
                         m, n, nduplc, npower, damp, acond, rnorm, &
                         line, line

    ! Check that aprod generates y + Ax and x + A'y consistently.
    call me%acheck( m, n, me%nout, eps, v, w, x, y, inform )

    if (inform > 0) then
        write(me%nout, '(a)') ' Check eps and power in subroutine acheck'
        stop
    end if

    ! Solve the problem defined by aprod, damp and b.
    ! Copy the rhs vector b into u  (LSQR will overwrite u)
    ! and set the other input parameters for LSQR.
    ! We ask for standard errors only if they are well-defined.

    call dcopy ( m, b, 1, u, 1 )
    !---  wantse = m > n  .or.  damp > zero
    wantse = .false.
    atol   = eps**0.99_wp
    btol   = atol
    conlim = 1000.0_wp * acond
    itnlim = 4*(m + n + 50)

    call me%lsqr(  m, n, damp, wantse, &
                   u, v, w, x, se, &
                   atol, btol, conlim, itnlim, me%nout, &
                   istop, itn, anorm, acond, rnorm, arnorm, xnorm )

    ! Examine the results.

    !-  if (damp == zero) then
    !-     write(nout, 2000)       xnorm, rnorm, arnorm
    !-  else
    !-     write(nout, 2100) damp, xnorm, rnorm, arnorm
    !-  end if

    call me%xcheck( m, n, me%nout, anorm, damp, eps, &
                    b, u, v, w, x, &
                    inform, test1, test2, test3 )

    ! Print the solution and standard error estimates from  LSQR.

    nprint = min( m, n, 8 )
    write(me%nout, 2500)    (j, x(j) , j = 1, nprint)
    if ( wantse ) then
        write(me%nout, 2600) (j, se(j), j = 1, nprint)
    end if

    ! Print a clue about whether the solution looks OK.

    do j = 1, n
        w(j)  = x(j) - xtrue(j)
    end do
    wnorm    = dnrm2 ( n, w    , 1 )
    xnorm    = dnrm2 ( n, xtrue, 1 )
    enorm    = wnorm / (one + xnorm)
    etol     = 0.001_wp
    if (enorm <= etol) then
        write(me%nout, 3000) enorm
    else
        write(me%nout, 3100) enorm
    end if
    return

    1000 format(1p &
    // 1x, 2a &
    /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )' &
    // ' Condition no. =', e12.4,  '     Residual function =', e17.9 &
    /  1x, 2a)
    2000 format(1p &
    // ' We are solving    min norm(Ax - b)    with no damping.' &
    // ' Estimates from LSQR:' &
    /  '    norm(x)         =', e10.3, ' = xnorm' &
    /  '    norm(r)         =', e10.3, ' = rnorm' &
    /  '    norm(A''r)       =', e10.3, ' = arnorm')
    2100 format(1p &
    // ' We are solving    min norm(Ax - b)    with damp =', e10.3 &
    /  '                           (damp*x)' &
    // ' Estimates from LSQR:' &
    /  '    norm(x)         =', e10.3, ' = xnorm' &
    /  '    norm(rbar)      =', e10.3, ' = rnorm' &
    /  '    norm(Abar''rbar) =', e10.3, ' = arnorm')
    2500 format(//' Solution  x:' / 4(i6, g14.6))
    2600 format(/ ' Standard errors  se:' / 4(i6, g14.6))
    3000 format(1p / ' LSQR  appears to be successful.', &
            '     Relative error in  x  =', e10.2)
    3100 format(1p / ' LSQR  appears to have failed.  ', &
            '     Relative error in  x  =', e10.2)
    8000 format(/ ' XXX  m or n is too large.')
    9000 format(/ ' XXX  Insufficient workspace.', &
            '  The length of  rw  should be at least', i6)

    end subroutine test
!***************************************************************************************************

!***************************************************************************************************
!>
!   This is the matrix-vector product routine required by subroutines
!   LSQR and CRAIG for a test matrix of the form  A = HY*D*HZ.
!   The quantities defining D, HY, HZ are in the work array rw,
!   followed by a work array w.  These are passed to aprod1 and aprod2
!   in order to make the code readable.

    subroutine aprod_test_solver (me, mode, m, n, x, y)

    class(test_solver),intent(inout) :: me
    integer,intent(in) :: mode
    integer,intent(in) :: m
    integer,intent(in) :: n
    real(wp),dimension(:),intent(inout) :: x   !! dimension n
    real(wp),dimension(:),intent(inout) :: y   !! dimension m

    integer :: locd, lochy, lochz, locw, maxmn, minmn

    maxmn  = max( m, n )
    minmn  = min( m, n )
    locd   = 1
    lochy  = locd  + minmn
    lochz  = lochy + m
    locw   = lochz + n

    if (mode == 1) then
        call me%aprod1( m, n, maxmn, minmn, x, y, &
                        me%rw(locd), me%rw(lochy), me%rw(lochz), me%rw(locw) )
    else
        call me%aprod2( m, n, maxmn, minmn, x, y, &
                        me%rw(locd), me%rw(lochy), me%rw(lochz), me%rw(locw) )
    end if

    end subroutine aprod_test_solver
!***************************************************************************************************

!***************************************************************************************************
!>
!  aprod1  computes  y = y + A*x  for subroutine aprod,
!  where A is a test matrix of the form  A = HY*D*HZ,
!  and the latter matrices HY, D, HZ are represented by
!  input vectors with the same name.

    subroutine aprod1( me, m, n, maxmn, minmn, x, y, d, hy, hz, w )

    class(test_solver),intent(inout) :: me
    integer :: m, n, maxmn, minmn
    real(wp) :: x(n), y(m), d(minmn), hy(m), hz(n), w(maxmn)

    integer :: i

    call hprod ( n, hz, x, w )

    do i = 1, minmn
        w(i)  = d(i) * w(i)
    end do

    do i = n + 1, m
        w(i)  = zero
    end do

    call hprod ( m, hy, w, w )

    do i = 1, m
        y(i)  = y(i) + w(i)
    end do

    end subroutine aprod1
!***************************************************************************************************

!***************************************************************************************************
!>
!  aprod2  computes  x = x + A(t)*y  for subroutine aprod,
!  where  A  is a test matrix of the form  A = HY*D*HZ,
!  and the latter matrices  HY, D, HZ  are represented by
!  input vectors with the same name.

    subroutine aprod2( me, m, n, maxmn, minmn, x, y, d, hy, hz, w )

    class(test_solver),intent(inout) :: me
    integer :: m, n, maxmn, minmn
    real(wp) :: x(n), y(m), d(minmn), hy(m), hz(n), w(maxmn)

    integer :: i

    call hprod ( m, hy, y, w )

    do i = 1, minmn
        w(i)  = d(i)*w(i)
    end do

    do i = m + 1, n
        w(i)  = zero
    end do

    call hprod ( n, hz, w, w )

    do i = 1, n
        x(i)  = x(i) + w(i)
    end do

    end subroutine aprod2
!***************************************************************************************************

!***************************************************************************************************
!>
!  hprod  applies a Householder transformation stored in  hz
!  to get  y = ( I - 2*hz*hz(transpose) ) * x.

    subroutine hprod ( n, hz, x, y )

    integer :: n
    real(wp) :: hz(n), x(n), y(n)

    integer :: i
    real(wp) :: s

    s = zero
    do i = 1, n
        s = hz(i) * x(i)  +  s
    end do

    s = s + s
    do i = 1, n
        y(i)  = x(i)  -  s * hz(i)
    end do

    end subroutine hprod
!***************************************************************************************************

!***************************************************************************************************
!>
!  lstp  generate a sparse least-squares test problem of the form
!             (   A    )*x = ( b )
!             ( damp*I )     ( 0 )
!  for solution by LSQR, or a sparse underdetermined system
!                Ax + damp*s = b
!  for solution by CRAIG.  The matrix A is m by n and is
!  constructed in the form  A = HY*D*HZ,  where D is an m by n
!  diagonal matrix, and HY and HZ are Householder transformations.
!
!  m and n may contain any positive values.
!  The first 8 parameters are input to lstp.  The last 8 are output.
!  If m >= n  or  damp = 0, the true solution is x as given.
!  Otherwise, x is modified to contain the true solution.

    subroutine lstp  ( me, m, n, maxmn, minmn, nduplc, npower, damp, x, &
                       b, d, hy, hz, w, acond, rnorm )

    class(test_solver),intent(inout) :: me
    integer :: m, n, maxmn, minmn, nduplc, npower
    real(wp) :: damp, acond, rnorm
    real(wp) :: b(m), x(n), d(minmn), hy(m), hz(n), w(maxmn)

    integer :: i, j
    real(wp) alfa, beta, dampsq, t

    real(wp),parameter :: fourpi = 4.0_wp * acos(-1.0_wp)   !4.0*3.141592

    ! ------------------------------------------------------------------
    ! Make two vectors of norm 1.0 for the Householder transformations.
    ! fourpi  need not be exact.
    ! ------------------------------------------------------------------
    minmn  = min( m, n )
    dampsq = damp**2
    alfa   = fourpi / m
    beta   = fourpi / n

    do i = 1, m
        hy(i) = sin( i * alfa )
    end do

    do i = 1, n
        hz(i) = cos( i * beta )
    end do

    alfa   = dnrm2 ( m, hy, 1 )
    beta   = dnrm2 ( n, hz, 1 )
    call dscal ( m, (- one / alfa), hy, 1 )
    call dscal ( n, (- one / beta), hz, 1 )

    ! ------------------------------------------------------------------
    ! Set the diagonal matrix  D.  These are the singular values of  A.
    ! ------------------------------------------------------------------
    do i = 1, minmn
        j     = (i - 1 + nduplc) / nduplc
        t     =  j * nduplc
        t     =  t / minmn
        d(i)  =  t**npower
    end do

    acond  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
    acond  = sqrt( acond )

    ! ------------------------------------------------------------------
    ! Set the true solution   x.
    ! It must be of the form  x = Z ( w )  for some  w.
    !                               ( 0 )
    ! ------------------------------------------------------------------
    call hprod ( n, hz, x, w )

    do i = m + 1, n
        w(i)  = zero
    end do

    call hprod ( n, hz, w, x )

    ! Solve D r1bar = dampsq x1bar
    ! where r1bar and x1bar are both in w.

    do i = 1, minmn
        w(i)  = dampsq * w(i) / d(i)
    end do

    ! Set r2bar to be anything.  (It is empty if m <= n)
    ! Then form r = Y rbar (again in w).

    do i = minmn+1, m
        w(i)  = one
    end do

    call Hprod ( m, hy, w, w )

    ! Compute the rhs    b = r  +  Ax.

    rnorm  = dnrm2 ( m, w, 1 )
    call dcopy ( m, w, 1, b, 1 )
    call me%aprod1( m, n, maxmn, minmn, x, b, d, hy, hz, w )

    end subroutine lstp
!***************************************************************************************************

!***************************************************************************************************
    end module lsqrtest_module
!***************************************************************************************************