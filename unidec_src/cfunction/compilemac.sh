#!/usr/bin/env bash
/opt/local/usr/local/bin/gcc -Wall -c pythonfunctions.c -O3 -fopenmp
/opt/local/usr/local/bin/gcc -shared -o libmypfunc.dylib pythonfunctions.o -O3 -fopenmp
cp libmypfunc.dylib ..\..\unidec_bin\libmypfunc.dylib
