#!/usr/bin/env python3

# cuttlery - 'Codon Usage Table Tools-lery.'
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
#
# This file is part of cuttlery.
#
# cuttlery is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cuttlery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pauvre.  If not, see <http://www.gnu.org/licenses/>.

# Tutorials on how to setup python package here:
#   - http://python-packaging.readthedocs.io/en/latest/testing.html
#   - https://jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/
#   - https://packaging.python.org/tutorials/distributing-packages/

# Use `pip install -e .` from the git repo root dir to make an editable installation

import os
from setuptools import setup, find_packages

version_py = os.path.join(os.path.dirname(__file__), 'cuttlery', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','').strip()
print(version)


setup(name='cuttlery',
      requires=['python (>3.0)'],
      version=version,
      description='Codon Usage Table Tools-lery.',
      long_description="""
          'cuttlery' is a package for determining the codon usage table of input
          sequences. Instead of calculating an absolute codon usage frequency
          for all observed genes as other packages do, cuttlery handles codon
          usage frequency as a distribution for each codon. The cuttlery package
          also contains tools for plotting codon usage information, calculating
          if an unknown ORF is likely to belong to coding or noncoding sequence,
          calculating nucleotide diversity of protein alignments and plotting
          the results, as well as looking for synonymous and nonsynonymous
          mutation heterogeneity in protein alignments.
          https://github.com/conchoecia/cuttlery
      """,

      url='https://github.com/conchoecia/cuttlery',
      author='Darrin Schultz',
      author_email='dts@ucsc.edu',
      classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Operating System :: POSIX :: Linux',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Intended Audience :: Science/Research'
          ],
      license='GPLv3',
      provides=['cuttlery'],
      packages=find_packages(),
      install_requires=[
          "matplotlib >= 2.0.2",
          "biopython >= 1.68",
          "pandas >= 0.20.1",
          "numpy >= 1.12.1"
      ],
      entry_points={
            'console_scripts': ['cuttlery=cuttlery.cuttlery_main:main'],
        },
      zip_safe=False,
      include_package_data=True)
