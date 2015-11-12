#!/usr/bin/env bash
gcc -Wall -c pythonfunctions.c -O3 -fopenmp -fPIC
gcc -shared -o libmypfunc.so pythonfunctions.o -O3 -fopenmp
cp libmypfunc.so ../../unidec_bin/libmypfunc.so
