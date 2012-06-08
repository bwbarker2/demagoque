AM_FCFLAGS = -ggdb -O3
AM_FCFLAGS_f90 = $(AM_FCFLAGS)
AM_FFLAGS = -ggdb -O3

bin_PROGRAMS = demagoque
demagoque_SOURCES = iso_varying_string.f90 prec_def.f90 phys_cons.f90 \
          bexception.f90 bmath.f90 \
          bstring.f90 \
          cons_laws.f90 class_Wavefunction.F90 class_SuperWavefunction.F90 \
          class_WfSquareWell.F90 class_WfKronigPenney.F90 \
          input_parameters.f90 lib_lapack.f90 \
          lib_fftw.f90 mesh.f90 time.f90 formatting.f90 dmtdhf.f90 ener.f90 \
          initial.f90 integra.f90 \
          interp.f \
          outAnalHarmonic.f90 output.f90 \
          renormalizeDM.f90 skyrme_params.f90 \
          time_evol.f90 wfnho.f90

#demagoque_CPPFLAGS = -fbounds-check -Wall -ggdb -O3
#demagoque_FCFLAGS = -ggdb -std=f2008 -warn -check all -stand f03 -standard-semantics -O3 #-fdefault-real-8 -fdefault-double-8
demagoque_FCFLAGS = -ggdb -std=f2008 -fcheck=all -Wall -O3 -Wno-surprising

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
