gcc -Wall -c pythonfunctions.c -O3 -fopenmp
gcc -shared -o libmypfunc32.dll pythonfunctions.o -O3 -fopenmp
copy libmypfunc32.dll ..\..\unidec_bin\libmypfunc32.dll
