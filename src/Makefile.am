AM_CPPFLAGS = -fbounds-check -Wall -ggdb -O3

bin_PROGRAMS = demagoque
demagoque_SOURCES = bmath.f90 bstring.f90 cons_laws.f90 dmtdhf.f90 ener.f90 \
          initial.f90 input_parameters.f90 integra.f90 \
          interp.f lib_fftw.f90 lib_lapack.f90 mesh.f90 \
          outAnalHarmonic.f90 output.f90 phys_cons.f90 \
          prec_def.f90 renormalizeDM.f90 skyrme_params.f90 \
          time.f90 time_evol.f90 wfnho.f90

demagoque_CPPFLAGS = -fbounds-check -Wall -ggdb -O3
demagoque_FCFLAGS = -fbounds-check -Wall -ggdb -O3

#SUFFIXES = .mod
#
#.o.mod:
#	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) $<
#
#.obj.mod:
#	$(FCCOMPILE) -c -o $@ $(FCFLAGS_f90) `$(CYGPATH_W) '$<'`
#
