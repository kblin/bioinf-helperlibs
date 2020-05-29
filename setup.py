import os
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
    'minimock',
    'pytest',
    'pytest-cover'
]


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
    install_requires=['BioPython>=1.62,<1.77'],
    tests_require=tests_require,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
    extras_require={
        'testing': tests_require,
    }
)
