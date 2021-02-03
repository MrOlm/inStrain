#!/usr/bin/env python

from setuptools import setup, find_packages

from inStrain._version import __version__

setup(name='inStrain',
      version=__version__,
      description='Calculation of strain-level metrics',
      url='https://github.com/MrOlm/inStrain',
      author='Matt Olm and Alex Crits-Christoph',
      author_email='mattolm@berkeley.edu',
      license='MIT',
      package_data={'inStrain': ['helper_files/NullModel.txt']},
      #packages=['inStrain'],
      packages=find_packages(),
      scripts=['bin/inStrain'],
      python_requires='>=3.4.0',
      install_requires=[
          'numpy',
          'pandas>=0.25,!=1.1.3',
          'seaborn',
          'matplotlib',
          'biopython<=1.74',
          'scikit-learn',
          'pytest',
          'tqdm',
          'pysam>=0.15', # This sets a requirement for python 3.7 for now, but so be it. pysam v0.9 (which works on python 3.8) has a broken iterator (has no stop)
          'networkx',
          'h5py',
          'psutil',
          'lmfit',
          'numba'
      ],
      zip_safe=False)
