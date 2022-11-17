from setuptools import setup, find_packages

setup(
    name="UniDec",
    version="5.2.0",
    author="Michael Marty",
    author_email="mtmarty@arizona.edu",
    keywords="UniDec Mass Spectrometry Deconvolution Ion Mobility",
    description='Universal Deconvolution of Electrospray Mass Spectrometry Data',
    url="https://github.com/michaelmarty/UniDec",

    install_requires=["numpy>=1.16", "scipy>=1.2", "matplotlib>=3.1", "wxpython>=4.0", "mpld3", 'pandas',
                      "pymzml", "natsort", "networkx", "h5py", "pypubsub", "multiplierz", "massql", "pymsfilereader",
                      "pyimzml", "pythonnet", "pyteomics", "lxml", "molmass", "rdkit"],
    python_requires='>=3.6',

    packages=find_packages(exclude=["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*"]),
    include_package_data=True,
    package_data={"": ["logo.ico", "readme.md", "UniDecLogoMR.png", "LICENSE"]},
    exclude_package_data={"": ["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*"]},
    download_url='https://github.com/michaelmarty/UniDec/archive/refs/tags/v.5.1.1.tar.gz',
    classifiers=[],
)
