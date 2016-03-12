#!/usr/bin/env python

'''
:author: Iyad Kandalaft <iyad.kandalaft@canada.ca>
:organization: Agriculture and Agri-Foods Canada
:group: Microbial Biodiversity Bioinformatics
:contact: mbb@agr.gc.ca 
:license: LGPL v3
'''

from setuptools import setup

setup(name='pybold',
      version='0.1',
      description='Barcode of Life Public API consumers in Python 2.7',
      author='Iyad Kandalaft',
      author_email='ik@iyadk.com',
      license='LGPL v3',
      url='http://packages.python.org/pybold',
      packages=['pybold', 'tests'],
      classifiers=[
                   "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
                   "Programming Language :: Python :: 2.7",
                   "Topic :: Scientific/Engineering :: Bio-Informatics",
                   "Development Status :: 4 - Beta"
                   ],
     )