srcdir = $(top_srcdir)/src

check_PROGRAMS=bmath_dZeroBrent_run

TESTS=$(check_PROGRAMS)
#XFAIL_TESTS=bmath_dZeroBrent

bmath_dZeroBrent_run_SOURCES = $(srcdir)/prec_def.f90 ${srcdir}/bmath.f90 $(srcdir)/bexception.f90 bmath_dZeroBrent_functions.F90 bmath_dZeroBrent_run.F90

mostlyclean-local:
	rm -f *.mod

