AM_FCFLAGS = -ggdb -O3
AM_FCFLAGS_f90 = -g -O3
AM_FFLAGS = -ggdb -O3

bin_PROGRAMS = demagoque
demagoque_SOURCES = iso_varying_string.f95 \
          bexception.f90 bmath.f90 \
          bstring.f90 prec_def.f90 \
          cons_laws.f90 input_parameters.f90 phys_cons.f90 lib_lapack.f90 \
          lib_fftw.f90 mesh.f90 time.f90 formatting.f90 dmtdhf.f90 ener.f90 \
          initial.f90 integra.f90 \
          interp.f \
          outAnalHarmonic.f90 output.f90 \
          renormalizeDM.f90 skyrme_params.f90 \
          time_evol.f90 wfnho.f90

#demagoque_CPPFLAGS = -fbounds-check -Wall -ggdb -O3
demagoque_FCFLAGS = -ggdb -O3

mostlyclean-local:
	rm -f *.mod

#SUFFIXES = .mod
#
#.o.mod:
#	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) $<
#
#.obj.mod:
#	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) `$(CYGPATH_W) '$<'`
#
