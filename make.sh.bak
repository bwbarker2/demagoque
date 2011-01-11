#!/bin/bash
# fake make

COMP=gfortran
#COMP=ifort

FFT=fftw
#FFT=pack

modules="common.f90"

main="bmath.f90 dmtdhf.f90 common.f90 ener.f90 initial.f90 integra.f90 output.f90 time_evol.f90 wfnho.f90 fft_fftw.f90 interp.f mesh.f90"

exe=time.x

ARCH=`uname -m`

#options="-lm ../../NAG/fll3a21dgl/lib/libnag_nag.a"
#options="-Llib -lfftpack5"
#options="-warn"
options="-Wall -ggdb -lfftw3"

$COMP -c $modules

$COMP $main $options -o $exe
