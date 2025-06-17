###MacOS Makefile for IsoGen and IsoDec
###Important notes for making this script run:
###Run this command in Powershell to Install WSL, it will be used to compile the .c files into the corresponding .so files

#wsl --install

###Run these commands to obtain the gcc compiler and make utility

#sudo apt update
#sudo apt install build-essential

#inside CLion. Go to Settings->ToolChains -> Add
#Build Tool: /usr/bin/make
#C Compiler: /usr/bin/gcc
#C++ Compiler /usr/bin/g++

#Even if this flags as not correctly finding WSL, it will still auto recognize the created Makefile. Allowing
#you to create the .so for isogen and isodec

# #navigate to the current directory in Ubuntu terminal:

#cd /mnt/c/Python/UniDecDev/unidec/IsoDec/src_cmake

#Then run the files to create the .so

#make clean

#make



CC = gcc
CFLAGS = -fPIC -O3
LDFLAGS = -lfftw3 -lm
FFTW_DIR = ./fftw
FFTW_SOURCE = /usr/lib/x86_64-linux-gnu/libfftw3.so

all: libisogen.so libisodeclib.so copy_fftw

libisogen.so: isogendep.c isogenpep.c isogenrna.c
	$(CC) $(CFLAGS) -shared -o libisogen.so \
	isogendep.c isogenpep.c isogenrna.c -I$(FFTW_DIR) $(LDFLAGS)

libisodeclib.so: isodeclib.c
	$(CC) $(CFLAGS) -shared -o libisodeclib.so \
	isodeclib.c -I. -L. -lisogen $(LDFLAGS)

copy_fftw:
	cp $(FFTW_SOURCE) ./libfftw3.so

clean:
	rm -f *.so *.o libfftw3.so
