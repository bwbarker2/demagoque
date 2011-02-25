#!/bin/bash

gfortran common.f90 dmtdhf.f90 ener.f90 fft_fftw.f90 initial.f90 integra.f90 \
         interp.f lib_lapack.f90 mesh.f90 \
         outAnalHarmonic.f90 output.f90 renormalizeDM.f90 \
         time_evol.f90 wfnho.f90 \
         -lfftw3 -llapack -lm \
         -Wall -ggdb \
         -o time.x
