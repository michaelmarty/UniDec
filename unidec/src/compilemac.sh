#!/usr/bin/env bash
/opt/local/usr/local/bin/gcc UniDec.c -fopenmp -lm -lfftw3 -lhdf5 -lhdf5_hl -O3 -o ../bin/unidecmac -std=c99 -arch x86_64

