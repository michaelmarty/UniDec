#!/usr/bin/env bash
gcc UniDec.c -o ../bin/unideclinux -I/usr/local/include -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lm -lfftw3 -lhdf5 -lhdf5_hl -fopenmp -O3 -w
