srcdir = $(top_srcdir)/src

check_PROGRAMS=bmath_dZeroBrent_run class_SquareWellState_run

TESTS=$(check_PROGRAMS)
#XFAIL_TESTS=bmath_dZeroBrent

bmath_dZeroBrent_run_SOURCES = $(srcdir)/prec_def.f90 \
                               ${srcdir}/bmath.f90 \
                               $(srcdir)/bexception.f90 \
                               $(srcdir)/phys_cons.f90 \
                               bmath_dZeroBrent_functions.F90 \
                               bmath_dZeroBrent_run.F90

class_SquareWellState_run_SOURCES = $(srcdir)/prec_def.f90 \
                                    $(srcdir)/class_SquareWellState.F90 \
                                    $(srcdir)/bmath.f90 \
                                    $(srcdir)/bexception.f90 \
                                    $(srcdir)/phys_cons.f90 \
                                    class_SquareWellState_run.F90

mostlyclean-local:
	rm -f *.mod

