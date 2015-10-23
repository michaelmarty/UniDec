from setuptools import setup

__author__ = 'Michael.Marty'

setup(name="UniDec",
      version="0.3",
      description="Universal Deconvolution of Mass and Ion Mobility Spectra",
      author="Michael Marty",
      author_email="michael.marty@chem.ox.ac.uk",
      url="http://unidec.chem.ox.ac.uk/",
      packages=['../UniDec'],
      install_requires=["numpy", "matplotlib", "scipy", "natsort", "twython", "pymzml", "networkx"]
      )
