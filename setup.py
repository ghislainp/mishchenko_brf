#!/usr/bin/env python

"""The setup script."""
import os

import setuptools
from setuptools import find_packages
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


requirements = ["numpy", "scipy"]

test_requirements = ['pytest>=3', ]


FORTRAN_SOURCE_DIR = 'src'                     # fortran source files

# Define the fortran extensions to make
ext_modules = [
    Extension(name="mishchenko_brf.lib.refl", sources=[os.path.join(FORTRAN_SOURCE_DIR, 'refl.f')]),
    Extension(name="mishchenko_brf.lib.spher", sources=[os.path.join(FORTRAN_SOURCE_DIR, 'spher.f')])
]


setup(
    author="Ghislain Picard",
    author_email='ghislain.picard@univ-grenoble-alpes.fr',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: Free For Educational Use',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Compute BRDF using Fortran Michael Mishchenko's radiative transfer solver",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='mishchenko_brf',
    name='mishchenko_brf',
    packages=find_packages(include=['mishchenko_brf', 'mishchenko_brf.*']),
    ext_modules=ext_modules,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ghislainp/mishchenko_brf',
    version='0.1.0',
    zip_safe=False,
)
