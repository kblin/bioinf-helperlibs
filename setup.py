import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "helperlibs",
    version = "0.1.0",
    author = "Kai Blin",
    author_email = "kai.blin@biotech.uni-tuebingen.de",
    description = ("A collection of bioinformatics-related helper functions"),
    license = "GPL",
    keywords = "bioinformatics",
    url = "https://github.com/kblin/bioinf-helperlibs/wiki",
    packages=['helperlibs'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3",
    ],
)
