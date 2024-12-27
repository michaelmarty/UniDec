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
    packages=find_packages(
                           exclude=["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*", "test_GUI.py"]
        ),
    include_package_data=True,
    package_data={"": [ "readme.md", "LICENSE"]},
    exclude_package_data={"": ["Scripts", "Scripts.*", "*.Scripts", "*.Scripts.*", ".gitignore", "unidec_doc"],
                          },

    classifiers=[],
)
