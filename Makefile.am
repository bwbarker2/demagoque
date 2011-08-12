AM_CPPFLAGS = -fbounds-check -Wall -ggdb -O3

bin_PROGRAMS = demagoque
demagoque_SOURCES = bmath.f90 bstring.f90 common.f90 dmtdhf.c90 ener.f90 \
                    initial.f90 input_parameters.f90 integra.f90 \
                    interp.f90 lib_fftw.f90 lib_lapack.f90 mesh.f90 \
                    outAnalHarmonic.f90 output.f90 renormalizeDM.f90 \
                    time_evol.f90 wfnho.f90

demagoque_CPPFLAGS = -lfftw3_threads -lpthread -lfftw3 -llapack -lm \
                     -fbounds-check -Wall -ggdb -O3

