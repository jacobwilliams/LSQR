!***************************************************************************************************
!>
!  Module for LSQR kinds and parameters
!
!### History
!  * Jacob Williams : 8 Nov 2019 : created module

    module lsqr_kinds

    use iso_fortran_env 

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: wp = real32   !! default real kind [4 bytes]
#elif REAL64
    integer,parameter,public :: wp = real64   !! default real kind [8 bytes]
#elif REAL128
    integer,parameter,public :: wp = real128  !! default real kind [16 bytes]
#else
    integer,parameter,public :: wp = real64   !! default real kind [8 bytes]
#endif

    ! parameters:
    real(wp),parameter,public :: zero = 0.0_wp
    real(wp),parameter,public :: one  = 1.0_wp

!***************************************************************************************************
    end module lsqr_kinds
!***************************************************************************************************