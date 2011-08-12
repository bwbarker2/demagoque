#!/bin/bash

gfortran mesh.f90 bmath.f90 bstring.f90 common.f90 dmtdhf.f90 \
         ener.f90 initial.f90 input_parameters.f90 integra.f90 \
         interp.f lib_fftw.f90 lib_lapack.f90 \
         outAnalHarmonic.f90 output.f90 renormalizeDM.f90 \
         time_evol.f90 wfnho.f90 \
         -lfftw3_threads -lpthread -lfftw3 -llapack -lm \
         -fbounds-check -Wall -ggdb \
         -O3 \
         -o time.x
