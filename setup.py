import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = 'undefined'
for line in open(os.path.join('helperlibs', '__init__.py'), 'r'):
    if line.startswith('version'):
        exec(line.strip())

setup(
    name = "helperlibs",
    version = version,
    author = "Kai Blin",
    author_email = "kai.blin@biotech.uni-tuebingen.de",
    description = ("A collection of bioinformatics-related helper functions"),
    license = "GPL",
    keywords = "bioinformatics",
    url = "https://github.com/kblin/bioinf-helperlibs/wiki",
    packages=['helperlibs', 'helperlibs.bio', 'helperlibs.wrappers',
              'helperlibs.tests', 'helperlibs.tests.bio',
              'helperlibs.tests.wrappers'],
    install_requires=['BioPython>=1.62'],
    tests_require=['unittest2','minimock','nose'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)
