#!/bin/bash

gfortran common.f90 dmtdhf.f90 ener.f90 initial.f90 integra.f90 \
         interp.f lib_fftw.f90 lib_lapack.f90 mesh.f90 \
         outAnalHarmonic.f90 output.f90 renormalizeDM.f90 \
         time_evol.f90 wfnho.f90 \
         -lfftw3_threads -lpthread -lfftw3 -llapack -lm \
         -Wall -ggdb \
         -O0 \
         -o time.x
