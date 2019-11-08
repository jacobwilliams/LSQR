!***************************************************************************************************
!>
!  Module for LSQR kinds and parameters

    module lsqr_kinds

    use iso_fortran_env, only: wp => real64  ! double precision kinds

    implicit none

    public

    ! parameters:
    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp

!***************************************************************************************************
    end module lsqr_kinds
!***************************************************************************************************