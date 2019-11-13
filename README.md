### Brief Description

A modern Fortran edition of LSQR, a conjugate-gradient type method for solving sparse linear equations and sparse least-squares problems.

The original Fortran 77 version of this algorithm can be found here: https://web.stanford.edu/group/SOL/software/lsqr/

The updated version has been significantly refactored.

### Compiling

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`LSQR.fobis`) is provided that can build the library and examples. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f LSQR.fobis -mode tests-gnu`
  * To build all the examples using ifort: `FoBiS.py build -f LSQR.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f LSQR.fobis -mode static-gnu`
  * To build a static library using ifort: `FoBiS.py build -f LSQR.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`

  To generate the documentation using [ford](https://github.com/cmacmackin/ford), run: ```FoBis.py rule --execute makedoc -f LSQR.fobis```