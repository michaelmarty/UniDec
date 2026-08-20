FROM python:3.13.9-bookworm

# C dependencies
RUN apt-get update && apt-get install -y libfftw3-dev hdf5-tools libhdf5-serial-dev libhdf5-dev cmake --no-install-recommends && rm -rf /var/lib/apt/lists/*

# Copy the source tree
COPY . /opt/UniDec

# Compile C
WORKDIR /opt/UniDec/unidec/src/
RUN /opt/UniDec/unidec/src/compilelinux.sh

WORKDIR /opt/UniDec/
ENV PATH="$PATH:/opt/UniDec/unidec/bin"

# Install UniDec and its external IsoDec, IsoGen, and importer dependencies.
RUN python -m pip install --no-cache-dir .

# Test everything
ENV TESTFILE="/opt/UniDec/unidec/bin/TestSpectra/test_1.txt"
ENV TESTFILE2="/opt/UniDec/unidec/bin/TestSpectra/test_2.txt"

RUN python -m unidec $TESTFILE
RUN python -m isodec $TESTFILE2 --centroided


