FROM python:3.9-bullseye

# C dependencies
RUN apt-get update && apt-get install -y libgtk-3-dev libwxgtk3.0-gtk3-dev libfftw3-dev hdf5-tools libhdf5-dev --no-install-recommends && rm -rf /var/lib/apt/lists/*

# To Run GUI
RUN pip install setuptools attrdict attrdict3 numpy pymzml networkx natsort h5py matplotlib scipy pyteomics mpld3 pandas
# RUN pip install wxpython

# Copy in stuff
RUN mkdir /opt/UniDec3
COPY . /opt/UniDec3

# Compile C
WORKDIR /opt/UniDec3/unidec_src/UniDec/
RUN /opt/UniDec3/unidec_src/unidec/compilelinux.sh

# Install Python
WORKDIR /opt/UniDec3/
RUN python setupdocker.py install

ENV TESTFILE /opt/UniDec3/unidec_bin/TestSpectra/test_1.txt
