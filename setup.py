#!/usr/bin/env python3

"""
MultiQC_PGx is a plugin for MultiQC developed for the LUMC PGx
project. It adds additional statistics related to the phasing of the target
regions, which are not suitable to add as full MultiQC modules.
"""

from setuptools import setup, find_packages

version = '0.0.1'

setup(
    name = 'multiqc_PGx',
    version = version,
    author = 'Redmar van den Berg',
    description = 'MultiQC plugin for the LUMC PGx project',
    long_description = __doc__,
    keywords = 'bioinformatics',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [ 'multiqc' ],
    entry_points = {
        'multiqc.cli_options.v1': [
            'add_header = multiqc_pgx.cli:add_header'
        ],
        'multiqc.modules.v1': [
            'target_phasing = multiqc_pgx.modules.target_phasing:MultiqcModule'
        ]
    }
)
