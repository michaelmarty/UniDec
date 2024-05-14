from setuptools import setup, find_packages

setup(
    setup_requires=["setuptools-git"],
    name="unidec",
    author="Michael Marty",
    author_email="mtmarty@arizona.edu",
    keywords="UniDec Mass Spectrometry Deconvolution Ion Mobility",
    description='Universal Deconvolution of Electrospray Mass Spectrometry Data',
    url="https://github.com/michaelmarty/UniDec",
    python_requires='>=3.7',

    install_requires=["numpy>=1.16", "scipy>=1.2", "matplotlib>=3.1", "lxml", 'pandas',
                      "pymzml", "natsort", "networkx", "h5py", "pyimzml", "pyteomics", "lxml_html_clean"],

    packages=find_packages(exclude=["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*"]),
    include_package_data=True,
    package_data={"": ["readme.md", "LICENSE"]},
    exclude_package_data={"": ["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*"]}
)
