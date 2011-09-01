#!/bin/bash

gfortran cons_laws.f90 formatting.f90 input_parameters.f90 mesh.f90 bmath.f90 \
         bstring.f90 time.f90 dmtdhf.f90 \
         ener.f90 initial.f90 integra.f90 \
         interp.f lib_fftw.f90 lib_lapack.f90 \
         outAnalHarmonic.f90 output.f90 phys_cons.f90 \
         prec_def.f90 renormalizeDM.f90 \
         skyrme_params.f90 time_evol.f90 wfnho.f90 \
         -lfftw3_threads -lpthread -lfftw3 -llapack -lm \
         -fbounds-check -Wall -ggdb \
         -O3 \
         -o time.x
