FROM python:3.11.9-bookworm

# C dependencies
RUN apt-get update && apt-get install -y libfftw3-dev hdf5-tools libhdf5-dev cmake --no-install-recommends && rm -rf /var/lib/apt/lists/*

# To Run GUI
RUN pip install setuptools attrdict attrdict3 numpy pymzml networkx natsort h5py matplotlib scipy pyteomics mpld3 pandas plotly numba
# RUN pip install wxpython

# Copy in stuff
RUN mkdir /opt/UniDec
COPY . /opt/UniDec

# Compile C
WORKDIR /opt/UniDec/unidec/src/
RUN /opt/UniDec/unidec/src/compilelinux.sh

# Compile IsoDec C
WORKDIR /opt/UniDec/unidec/IsoDec/src_cmake/
RUN /opt/UniDec/unidec/IsoDec/src_cmake/compilelinux.sh

# Install Python
WORKDIR /opt/UniDec/
RUN python setupdocker.py install

# Add the UniDec bin folder to the path
ENV PATH="$PATH:/opt/UniDec/unidec/bin"

# Test everything
ENV TESTFILE="/opt/UniDec/unidec/bin/TestSpectra/test_1.txt"
ENV TESTFILE2="/opt/UniDec/unidec/bin/TestSpectra/test_2.txt"
