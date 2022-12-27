from setuptools import setup, find_packages

setup(
    setup_requires=["setuptools-git"],
    name="unidec",
    #version="5.2.1.1",
    author="Michael Marty",
    author_email="mtmarty@arizona.edu",
    keywords="unidec Mass Spectrometry Deconvolution Ion Mobility",
    description='Universal Deconvolution of Electrospray Mass Spectrometry Data',
    url="https://github.com/michaelmarty/UniDec",
    python_requires='>=3.7',
    packages=find_packages(#"unidec",
                            #include=["Launcher", "modules", "metaunidec", "iFAMS", "LipiDec", ""],
                           exclude=["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*", "test_GUI.py"]
        ),
    include_package_data=True,
    package_data={"": ["logo.ico", "readme.md", "UniDecLogoMR.png", "LICENSE"]},
    exclude_package_data={"": ["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*", ".gitignore"],
                          },
    download_url='https://github.com/michaelmarty/UniDec/archive/refs/tags/v.5.2.1.tar.gz',
    classifiers=[],
)
