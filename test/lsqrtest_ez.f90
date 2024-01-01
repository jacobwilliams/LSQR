!***************************************************************************************************
!>
!  Main program for EZ test.

    program main

    use lsqr_kinds
    use lsqr_module,     only: lsqr_solver_ez
    use iso_fortran_env, only: output_unit

    implicit none

    call test_1()
    call test_2()

    contains 

    subroutine test_1()

        ! define a 3x3 dense system to solve:
        integer,parameter :: m = 3 !! number of rows in `A` matrix
        integer,parameter :: n = 3 !! number of columns in `A` matrix
        real(wp),dimension(m),parameter :: b = real([1,2,3],wp)
        integer,dimension(m*n),parameter :: icol = [1,1,1,2,2,2,3,3,3]
        integer,dimension(m*n),parameter :: irow = [1,2,3,1,2,3,1,2,3]
        real(wp),dimension(m*n),parameter :: a = real([1,4,7,2,5,88,3,66,9],wp)
        real(wp),parameter :: damp = zero
        real(wp),dimension(m,1),parameter :: b_vec = reshape(b,[m,1])
        real(wp),dimension(m,n),parameter :: a_mat = reshape(a,[3,3])

        type(lsqr_solver_ez) :: solver
        real(wp),dimension(n) :: x
        real(wp),dimension(n,1) :: x_vec
        integer :: istop

        call solver%initialize(m,n,a,irow,icol,&
                                itnlim = 100,&
                                nout   = output_unit)
        call solver%solve(b,damp,x,istop)

        x_vec(1:3,1) = x

        write(*,*) ''
        write(*,*) 'istop = ', istop
        write(*,*) ''
        write(*,'(1P,A,*(E16.6))') 'x       = ', x
        write(*,'(1P,A,*(E16.6))') 'A*x     = ', matmul(a_mat, x_vec)
        write(*,'(1P,A,*(E16.6))') 'A*x - b = ', matmul(a_mat, x_vec) - b_vec

        if (any(abs(matmul(a_mat, x_vec) - b_vec) > 1.0e-12)) error stop 'TEST FAILED'

    end subroutine test_1

    subroutine test_2()

        ! another test case (with n>m). compare to scipy.

        !!```
        !! >>> a
        !! array([[  4. ,   5. ,  66. ,   0.1],
        !!        [  1. ,  -3. ,   8. ,  -9. ],
        !!        [ 11. ,   3. , -87. ,   2. ]])
        !! >>> b
        !! array([1, 2, 3])
        !! >>> scipy.sparse.linalg.lsqr(a, b)
        !! (array([ 0.26437473,  0.04901579, -0.00426183, -0.21297414]), 1, 3, 5.5785963493386424e-12, 5.5785963493386424e-12, 110.70234866523838, 15.316189089999897, 6.119548932941366e-10, 0.343034538979173, array([0., 0., 0., 0.]))
        !! >>> 
        !!```

        integer,parameter :: m = 3 !! number of rows in `A` matrix
        integer,parameter :: n = 4 !! number of columns in `A` matrix
        real(wp),dimension(m),parameter :: b = real([1,2,3],wp)
        integer,dimension(m*n),parameter :: icol = [1,1,1,2,2,2,3,3,3,4,4,4]
        integer,dimension(m*n),parameter :: irow = [1,2,3,1,2,3,1,2,3,1,2,3]
        real(wp),dimension(m*n),parameter :: a = [4.1_wp,1.1_wp,11.1_wp,&
                                                  5.1_wp,-3.1_wp,3.1_wp,&
                                                  66.1_wp,8.1_wp,-87.1_wp,&
                                                  0.1_wp,-9.1_wp,2.1_wp]
        real(wp),parameter :: damp = zero
        real(wp),dimension(m,1),parameter :: b_vec = reshape(b,[m,1])
        real(wp),dimension(m,n),parameter :: a_mat = reshape(a,[3,4])

        type(lsqr_solver_ez) :: solver
        real(wp),dimension(n) :: x
        real(wp),dimension(n,1) :: x_vec
        integer :: istop

        call solver%initialize(m,n,a,irow,icol,&
                                itnlim = 100,&
                                nout   = output_unit)
        call solver%solve(b,damp,x,istop)

        x_vec(1:4,1) = x

        write(*,*) ''
        write(*,*) 'istop = ', istop
        write(*,*) ''
        write(*,'(1P,A,*(E16.6))') 'x       = ', x
        write(*,'(1P,A,*(E16.6))') 'A*x     = ', matmul(a_mat, x_vec)
        write(*,'(1P,A,*(E16.6))') 'A*x - b = ', matmul(a_mat, x_vec) - b_vec

        if (any(abs(matmul(a_mat, x_vec) - b_vec) > 1.0e-12)) error stop 'TEST FAILED'

    end subroutine test_2


    end program main
!***************************************************************************************************