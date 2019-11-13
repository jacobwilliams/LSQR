!***************************************************************************************************
!>
!  this file contains BLAS routines required by subroutines [[lsqr]] and [[acheck]].
!
!### History
!  * Jacob Williams : 11/8/2019 : using modernized versions of these routines.

    module lsqpblas_module

    use lsqr_kinds, only: wp, zero, one

    implicit none

    private

    public :: dcopy,ddot,dnrm2,dscal

    contains
!***************************************************************************************************

!***************************************************************************************************
    subroutine dcopy(n,dx,incx,dy,incy)

    integer incx,incy,n
    real(wp) dx(*),dy(*)

    integer i,ix,iy,m,mp1

    if (n<=0) return
    if (incx==1 .and. incy==1) then
        ! code for both increments equal to 1
        ! clean-up loop
        m = mod(n,7)
        if (m/=0) then
            do i = 1,m
                dy(i) = dx(i)
            end do
            if (n<7) return
        end if
        mp1 = m + 1
        do i = mp1,n,7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
        end do
    else
        ! code for unequal increments or equal increments
        ! not equal to 1
        ix = 1
        iy = 1
        if (incx<0) ix = (-n+1)*incx + 1
        if (incy<0) iy = (-n+1)*incy + 1
        do i = 1,n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
        end do
    end if

    end subroutine dcopy
!***************************************************************************************************

!***************************************************************************************************
    real(wp) function ddot(n,dx,incx,dy,incy)

    integer incx,incy,n
    real(wp) dx(*),dy(*)

    real(wp) dtemp
    integer i,ix,iy,m,mp1

    ddot = zero
    dtemp = zero
    if (n<=0) return
    if (incx==1 .and. incy==1) then
        ! code for both increments equal to 1
        ! clean-up loop
        m = mod(n,5)
        if (m/=0) then
            do i = 1,m
                dtemp = dtemp + dx(i)*dy(i)
            end do
            if (n<5) then
                ddot=dtemp
                return
            end if
        end if
        mp1 = m + 1
        do i = mp1,n,5
            dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
        end do
    else
        ! code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if (incx<0) ix = (-n+1)*incx + 1
        if (incy<0) iy = (-n+1)*incy + 1
        do i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
        end do
    end if
    ddot = dtemp

    end function ddot
!***************************************************************************************************

!***************************************************************************************************
    real(wp) function dnrm2(n,x,incx)

    integer incx,n
    real(wp) x(*)

    real(wp) absxi,norm,scale,ssq
    integer ix

    if (n<1 .or. incx<1) then
        norm = zero
    else if (n==1) then
        norm = abs(x(1))
    else
        scale = zero
        ssq = one

        ! the following loop is equivalent to this call to the lapack
        ! auxiliary routine:
        ! call dlassq( n, x, incx, scale, ssq )

        do ix = 1,1 + (n-1)*incx,incx
            if (x(ix)/=zero) then
                absxi = abs(x(ix))
                if (scale<absxi) then
                    ssq = one + ssq* (scale/absxi)**2
                    scale = absxi
                else
                    ssq = ssq + (absxi/scale)**2
                end if
            end if
        end do
        norm = scale*sqrt(ssq)
    end if

    dnrm2 = norm

    end function dnrm2
!***************************************************************************************************

!***************************************************************************************************
    subroutine dscal(n,da,dx,incx)

    real(wp) da
    integer incx,n
    real(wp) dx(*)

    integer i,m,mp1,nincx

    if (n<=0 .or. incx<=0) return
    if (incx==1) then
        ! code for increment equal to 1
        ! clean-up loop
        m = mod(n,5)
        if (m/=0) then
            do i = 1,m
                dx(i) = da*dx(i)
            end do
            if (n<5) return
        end if
        mp1 = m + 1
        do i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
        end do
    else
        ! code for increment not equal to 1
        nincx = n*incx
        do i = 1,nincx,incx
            dx(i) = da*dx(i)
        end do
    end if

    end subroutine dscal
!***************************************************************************************************

!***************************************************************************************************
    end module lsqpblas_module
!***************************************************************************************************