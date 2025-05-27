import os
import sys
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


version = 'undefined'
for line in open(os.path.join('helperlibs', '__init__.py'), 'r'):
    if line.startswith('version'):
        exec(line.strip())

short_description = "A collection of bioinformatics-related helper functions"
long_description = open('README.md').read()

tests_require = [
    'flake8',
    'minimock',
    'pytest',
    'pytest-cover',
]


biopython = "Biopython>=1.76"
if sys.version_info.major == 3 and sys.version_info.minor == 5:
    biopython += ",<1.77"


setup(
    name="helperlibs",
    version=version,
    author="Kai Blin",
    author_email="kblin@biosustain.dtu.dk",
    description=short_description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="GPL",
    keywords="bioinformatics",
    url="https://github.com/kblin/bioinf-helperlibs/wiki",
    packages=['helperlibs', 'helperlibs.bio', 'helperlibs.wrappers',
              'helperlibs.tests', 'helperlibs.tests.bio',
              'helperlibs.tests.wrappers'],
    install_requires=[biopython],
    tests_require=tests_require,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    extras_require={
        'testing': tests_require,
    }
)
