FROM python:3.13.9-bookworm

# C dependencies
RUN apt-get update && apt-get install -y libfftw3-dev hdf5-tools libhdf5-serial-dev libhdf5-dev cmake --no-install-recommends && rm -rf /var/lib/apt/lists/*

# Copy in stuff
RUN mkdir /opt/UniDec
COPY . /opt/UniDec

# Compile C
WORKDIR /opt/UniDec/unidec/src/
RUN /opt/UniDec/unidec/src/compilelinux.sh

# Compile IsoDec C
WORKDIR /opt/UniDec/unidec/IsoDec/src_cmake/
RUN /opt/UniDec/unidec/IsoDec/src_cmake/compilelinux.sh

# Test that it works
WORKDIR /opt/UniDec/
# Add the UniDec bin folder to the path
ENV PATH="$PATH:/opt/UniDec/unidec/bin"
#RUN unideclinux

# To Run GUI
RUN pip install setuptools attrdict attrdict3 numpy pymzml h5py networkx natsort matplotlib scipy pyteomics mpld3 pandas plotly numba lxml lxml_html_clean
# RUN pip install wxpython

# Install Python
RUN python setupdocker.py install

# Test everything
ENV TESTFILE="/opt/UniDec/unidec/bin/TestSpectra/test_1.txt"
ENV TESTFILE2="/opt/UniDec/unidec/bin/TestSpectra/test_2.txt"

RUN python -m unidec $TESTFILE
RUN python -m unidec.IsoDec $TESTFILE2 -precentroided


