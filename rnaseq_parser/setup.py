#!/usr/bin/env python

from setuptools import setup

setup(name='RNAseqParser',
      version='0.1',
      description='Add data from RNAseq pipeline to molgenis database',
      author='Niek de Klein',
      author_email='niekdeklein@gmail.com',
      url='https://github.com/molgenis/ngs-utils/tree/master/rnaseq_parser',
      license="?",
      packages=['RNAseqParser'],
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3'
      ],
    keywords='genotyping RNAseq ENA Molgenis',
    install_requires=['pyyaml',
                      'requests',
                      'PBKDF2',
                      'crypto',
                      'configparser'],
    test_suite='nose.collector',
    tests_require=['nose'],
     )
