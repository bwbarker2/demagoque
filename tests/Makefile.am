srcdir = $(top_srcdir)/src

AM_FCFLAGS=$(FC_MODINC). $(FC_MODINC)../src

check_PROGRAMS=bmath_dZeroBrent_run class_WfSquareWell_run \
               test_bmath_LDiffRichardson

TESTS=$(check_PROGRAMS)
#XFAIL_TESTS=bmath_dZeroBrent

bmath_dZeroBrent_run_SOURCES = $(srcdir)/prec_def.f90 \
                               $(srcdir)/specfun.f90 \
                               ${srcdir}/bmath.f90 \
                               $(srcdir)/bexception.f90 \
                               $(srcdir)/iso_varying_string.f90 \
                               $(srcdir)/phys_cons.f90 \
                               bmath_dZeroBrent_functions.F90 \
                               bmath_dZeroBrent_run.F90

class_WfSquareWell_run_SOURCES = $(srcdir)/prec_def.f90 \
                                 $(srcdir)/class_Wavefunction.F90 \
                                 $(srcdir)/class_WfSquareWell.F90 \
                                 $(srcdir)/specfun.f90 \
                                 $(srcdir)/bmath.f90 \
                                 $(srcdir)/bexception.f90 \
                                 $(srcdir)/iso_varying_string.f90 \
                                 $(srcdir)/phys_cons.f90 \
                                 class_WfSquareWell_run.F90

test_bmath_LDiffRichardson_SOURCES = $(srcdir)/prec_def.f90 \
                                     $(srcdir)/bexception.f90 \
                                     $(srcdir)/specfun.f90 \
                                     $(srcdir)/bmath.f90 \
                                     $(srcdir)/iso_varying_string.f90 \
                                     $(srcdir)/phys_cons.f90 \
                                     test_bmath_LDiffRichardson.F90

mostlyclean-local:
	rm -f *.mod

