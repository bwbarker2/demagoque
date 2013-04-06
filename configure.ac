# Process this file with autoconf to produce a configure script.

# we actually need this version to get fortran module support
AC_PREREQ([2.69])

AC_CONFIG_AUX_DIR([build-aux])
AC_INIT([demagoque], m4_esyscmd_s([echo 0.01r`bzr revno`]), [barker@nscl.msu.edu],[],[[https://launchpad.net/demagoque]])

# aims to remove default -O2 option from FCFLAGS while still allowing command-line override
if test -n "${FCFLAGS}"; then
  user_set_fcflags=yes
fi

AC_CONFIG_SRCDIR([src/dmtdhf.f90])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PROG_AWK
AC_PROG_FC
AC_PROG_F77
AC_PROG_SED

# finds compiler flag to include Fortran 90 module info from another directory,
# stores in FC_MODINC. Used in tests/ to find mod files in src/
# \TODO add error message if compiler doesn't have module directory options
AC_FC_MODULE_FLAG

# actually removes the default -O2 option from FCFLAGS
if test X$user_set_fcflags != Xyes; then
 FCFLAGS=`echo ${FCFLAGS} | ${SED} s/\ -O2//`
 FFLAGS=`echo ${FFLAGS} | ${SED} s/\ -O2//`
fi

if test X$FC == Xgfortran; then
 FCFLAGS="$FCFLAGS -Wall -Wno-surprising"
 echo 'checking for gfortran compiler options... '
 filetest=`date +%N`
 touch ${filetest}.f90
 echo -n '   fcheck=all... '
 if $FC -c -fcheck=all ${filetest}.f90 &> /dev/null; then
  HAVE_FCHECK='yes'
  FCFLAGS="$FCFLAGS -fcheck=all"
  FFLAGS="$FFLAGS -fcheck=all"
  echo 'yes'
 else
  echo 'no'
  echo -n '   fbounds-check...'
  if $FC -c -fbounds-check ${filetest}.f90 &> /dev/null; then
   FCFLAGS="$FCFLAGS -fbounds-check"
   FFLAGS="$FFLAGS -fbounds-check"
   echo 'yes'
  else
   echo 'no'
  fi
 fi
 rm ${filetest}*
fi
 
AC_LANG([Fortran])
AC_FC_SRCEXT([f90])
AC_CHECK_LIB([fftw3],[dfftw_execute],[],[AC_MSG_ERROR([library 'fftw3' not found])])
AC_CHECK_LIB([lapack],[zgeev],[],[AC_MSG_ERROR([library 'lapack' not found])])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
