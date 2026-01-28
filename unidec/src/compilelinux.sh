#gcc UniDec.c -o ../bin/unideclinux -I/usr/local/include -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -lm -lfftw3 -lhdf5 -lhdf5_hl -fopenmp -O3 -w
# Build with CMake
# cmake -S . -B cmake-build-release-wsl -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -G "Unix Makefiles" --fresh
cmake -S . -B cmake-build-release-wsl -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=../../bin -G "Unix Makefiles" --fresh
# Install
cmake --build cmake-build-release-wsl --target all -- -j$(nproc)
