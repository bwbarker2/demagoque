# Process this file with autoconf to produce a configure script.

AC_CONFIG_AUX_DIR([build-aux])
AC_INIT([demagoque], m4_esyscmd_s([echo 0.01r`bzr revno`]), [barker@nscl.msu.edu],[],[[https://launchpad.net/demagoque]])
AC_CONFIG_SRCDIR([src/dmtdhf.f90])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_FC
AC_PROG_F77
AC_LANG([Fortran])
AC_FC_SRCEXT([f90])
AC_CHECK_LIB([fftw3],[dfftw_execute],[],[AC_MSG_ERROR([library 'fftw3' not found])])
AC_CHECK_LIB([lapack],[zgeev],[],[AC_MSG_ERROR([library 'lapack' not found])])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
