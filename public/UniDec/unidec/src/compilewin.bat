gcc "-IC:\\MinGW\\fftw" "-IC:\\Python\\UniDec" -O3 -fopenmp -Wall -c -o UniDec.o "..\\UniDec.c" 
gcc "-LC:\\MinGW\\fftw" -fopenmp -o UniDec.exe UniDec.o -lm -lfftw3-3.dll
cp UniDec.exe "..\\" 