!***************************************************************************************************
!>
!  Main program for EZ test.

    program main

    use lsqr_kinds
    use lsqr_module,     only: lsqr_solver_ez
    use iso_fortran_env, only: output_unit

    implicit none

    ! define a 3x3 dense system to solve:
    integer,parameter :: m = 3 !! number of rows in `A` matrix
    integer,parameter :: n = 3 !! number of columns in `A` matrix
    real(wp),dimension(m),parameter :: b = real([1,2,3],wp)
    integer,dimension(m*n),parameter :: icol = [1,1,1,2,2,2,3,3,3]
    integer,dimension(m*n),parameter :: irow = [1,2,3,1,2,3,1,2,3]
    real(wp),dimension(m*n),parameter :: a = real([1,4,7,2,5,88,3,66,9],wp)
    real(wp),parameter :: damp = zero
    real(wp),dimension(m,1),parameter :: b_vec = b
    real(wp),dimension(m,n),parameter :: a_mat = reshape(a,[3,3])

    type(lsqr_solver_ez) :: solver
    real(wp),dimension(n) :: x
    real(wp),dimension(1) :: se !! not used
    real(wp),dimension(n,1) :: x_vec

    call solver%initialize(m,n,a,irow,icol,&
                            itnlim = 100,&
                            nout   = output_unit)
    call solver%solve(b,damp,x)

    x_vec(1:3,1) = x

    write(*,*) ''
    write(*,'(1P,A,*(E16.6))') 'x       = ', x
    write(*,'(1P,A,*(E16.6))') 'A*x     = ', matmul(a_mat, x_vec)
    write(*,'(1P,A,*(E16.6))') 'A*x - b = ', matmul(a_mat, x_vec) - b_vec

    end program main
!***************************************************************************************************